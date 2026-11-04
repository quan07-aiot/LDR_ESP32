/**
 * BÀI TẬP LỚN: HỆ THỐNG THU THẬP & PHÂN TÍCH PHỔ TÍN HIỆU CẢM BIẾN LDR (ESP32-S3)
 * 1. Thu thập mẫu LDR qua Web Serial
 * 2. Hai bộ lọc LTI rời rạc: IIR bậc một và FIR trung bình trượt
 * 3. FFT Radix-2, so sánh phổ thô / phổ đã lọc
 * 4. Biểu đồ, phiên ghi 30 giây và xuất báo cáo
 */
"use strict";

const BAUD_RATE = 115200;
const MAX_SAMPLES = 6000;
const MAX_FFT_SAMPLES = 4096;
const samples = []; // Lưu trữ { time, raw, filtered }

let port = null, reader = null, reading = false;
let textBuffer = "", prevDeviceTime = null, elapsedTime = 0;
let filterState = { iir: null, window: [] };
let renderPending = false;

// Trạng thái phiên ghi 30s
const RECORD_MS = 30000;
let isRecording = false, recStart = 0, recSamples = [];
let isFrozen = false, lastFftResult = null, lastRawSpectrum = null, lastTableTime = 0;

const $ = id => document.getElementById(id);

// =====================================================================
// BỘ LỌC LTI RỜI RẠC: IIR BẬC MỘT HOẶC FIR TRUNG BÌNH TRƯỢT
// =====================================================================
function filterSample(raw, state = filterState) {
  const mode = $("filterMode").value;
  const alpha = parseFloat($("filterAlpha").value) || 0.2;
  const winSize = parseInt($("windowSize").value) || 5;

  if (mode === "rc_iir") {
    // Hệ thức sai phân IIR: y[n] = αx[n] + (1 - α)y[n-1].
    state.iir = state.iir === null ? raw : alpha * raw + (1 - alpha) * state.iir;
    return state.iir;
  }

  // FIR: đáp ứng xung hữu hạn, trung bình N mẫu đầu vào gần nhất.
  state.window.push(raw);
  if (state.window.length > winSize) state.window.shift();
  return state.window.reduce((sum, value) => sum + value, 0) / state.window.length;
}

function recomputeFilter(data) {
  const state = { iir: null, window: [] };
  data.forEach(sample => { sample.filtered = filterSample(sample.raw, state); });
  return state;
}

function refreshFilterSettings() {
  const mode = $("filterMode").value;
  $("alphaControl").classList.toggle("hidden", mode !== "rc_iir");
  $("windowControl").classList.toggle("hidden", mode !== "moving_avg");

  if (mode === "rc_iir") {
    const alpha = parseFloat($("filterAlpha").value);
    $("filterTitle").textContent = "Thông thấp bậc 1 · IIR";
    $("filterEquation").textContent = `y[n] = ${alpha.toFixed(2)}x[n] + ${(1 - alpha).toFixed(2)}y[n−1]`;
    $("filterDescription").textContent = "Đầu ra dùng mẫu hiện tại và đầu ra trước đó. α nhỏ làm suy giảm mạnh hơn các thành phần tần số cao.";
  } else {
    const count = parseInt($("windowSize").value, 10);
    $("filterTitle").textContent = "Trung bình trượt · FIR";
    $("filterEquation").textContent = `y[n] = (x[n] + … + x[n−${count - 1}]) / ${count}`;
    $("filterDescription").textContent = "Lấy trung bình N mẫu gần nhất, không dùng đầu ra trước. N lớn làm đường tín hiệu mượt hơn nhưng phản ứng chậm hơn.";
  }

  filterState = recomputeFilter(samples);
  if (recSamples.length) recomputeFilter(recSamples);
  if (samples.length) $("latestFiltered").textContent = Math.round(samples[samples.length - 1].filtered);
  scheduleRender();
}

// =====================================================================
// KHỐI 2: THUẬT TOÁN BIẾN ĐỔI FOURIER NHANH FFT (COOLEY-TUKEY RADIX-2)
// =====================================================================
function fftCooleyTukey(real, imag) {
  const n = real.length;
  // 1. Bit-reversal permutation (Đảo bit thứ tự mẫu)
  for (let i = 0, j = 0; i < n - 1; i++) {
    if (i < j) {
      [real[i], real[j]] = [real[j], real[i]];
      [imag[i], imag[j]] = [imag[j], imag[i]];
    }
    let k = n >> 1;
    while (k <= j) { j -= k; k >>= 1; }
    j += k;
  }
  // 2. Cooley-Tukey butterfly computations (Phép tính cánh bướm)
  for (let len = 2; len <= n; len <<= 1) {
    const half = len >> 1, ang = -2 * Math.PI / len;
    const wStepR = Math.cos(ang), wStepI = Math.sin(ang);
    for (let i = 0; i < n; i += len) {
      let wR = 1, wI = 0;
      for (let k = 0; k < half; k++) {
        const pEven = i + k, pOdd = pEven + half;
        const tr = wR * real[pOdd] - wI * imag[pOdd];
        const ti = wR * imag[pOdd] + wI * real[pOdd];
        real[pOdd] = real[pEven] - tr;
        imag[pOdd] = imag[pEven] - ti;
        real[pEven] += tr;
        imag[pEven] += ti;
        const nxtR = wR * wStepR - wI * wStepI;
        wI = wR * wStepI + wI * wStepR;
        wR = nxtR;
      }
    }
  }
}

// Phân tích phổ: Khử DC + Cửa sổ Hann + FFT + Dò đỉnh
function computeSpectrum(data, sourceKey) {
  const count = data.length;
  if (count < 2) return null;
  const duration = (data[count - 1].time - data[0].time) / 1000;
  if (duration <= 0 || !Number.isFinite(duration)) return null;
  const fs = (count - 1) / duration;

  // 1. Khử thành phần DC Offset (Zero-mean) & áp dụng cửa sổ Hann
  const mean = data.reduce((s, d) => s + d[sourceKey], 0) / count;
  const real = new Float32Array(count), imag = new Float32Array(count);
  let hannSum = 0;
  for (let i = 0; i < count; i++) {
    const h = 0.5 * (1 - Math.cos(2 * Math.PI * i / (count - 1)));
    real[i] = (data[i][sourceKey] - mean) * h;
    hannSum += h;
  }

  // 2. Biến đổi Fourier
  fftCooleyTukey(real, imag);

  // 3. Tính biên độ đơn biên (Single-sided Magnitude)
  const scale = hannSum > 0 ? 2 / hannSum : 2 / count;
  const spectrum = [];
  let totalEnergy = 0;

  for (let bin = 0; bin <= count / 2; bin++) {
    const freq = bin * fs / count;
    let amp = Math.hypot(real[bin], imag[bin]) * scale;
    if (bin === 0 || bin === count / 2) amp *= 0.5;
    spectrum.push({ bin, freq, amp });
    if (bin > 0) totalEnergy += amp * amp;
  }

  // 4. Thuật toán tìm đỉnh cực đại (Peak Detection)
  const peaks = [];
  for (let b = 1; b < count / 2; b++) {
    if (spectrum[b].amp > 0.5 && spectrum[b].amp >= spectrum[b - 1].amp && spectrum[b].amp >= spectrum[b + 1].amp) {
      const p = spectrum[b];
      peaks.push({
        freq: p.freq, amp: p.amp, period: p.freq > 0 ? 1 / p.freq : 0,
        powerDb: 20 * Math.log10(p.amp + 1e-6),
        energyPercent: totalEnergy > 0 ? (p.amp * p.amp / totalEnergy) * 100 : 0
      });
    }
  }
  peaks.sort((a, b) => b.amp - a.amp);

  // Gán nhãn phân loại đỉnh
  if (peaks.length > 0) {
    const f0 = peaks[0].freq;
    peaks[0].type = "dominant"; peaks[0].label = "Tín hiệu chủ đạo (Dominant)";
    for (let i = 1; i < peaks.length; i++) {
      const ratio = peaks[i].freq / f0;
      if (Math.abs(peaks[i].freq - 50) < 1.5) { peaks[i].type = "noise"; peaks[i].label = "Thành phần gần 50 Hz"; }
      else if (Math.abs(ratio - 2) < 0.15) { peaks[i].type = "harmonic"; peaks[i].label = "Hài bậc 2 (2× f₀)"; }
      else if (Math.abs(ratio - 3) < 0.15) { peaks[i].type = "harmonic"; peaks[i].label = "Hài bậc 3 (3× f₀)"; }
      else if (peaks[i].energyPercent > 5) { peaks[i].type = "secondary"; peaks[i].label = "Đỉnh phụ rõ rệt"; }
      else { peaks[i].type = "noise"; peaks[i].label = "Nhiễu dao động nhỏ"; }
    }
  }

  return { fs, duration, count, spectrum, peaks, dominantPeak: peaks[0] || null };
}

// =====================================================================
// KHỐI 3: VẼ ĐỒ THỊ ĐỒNG BỘ CHUẨN MATLAB (CANVAS HELPER)
// =====================================================================
function prepareCanvas(canvas) {
  const w = Math.max(1, canvas.clientWidth), h = Math.max(1, canvas.clientHeight);
  const ratio = window.devicePixelRatio || 1;
  canvas.width = Math.round(w * ratio);
  canvas.height = Math.round(h * ratio);
  const ctx = canvas.getContext("2d");
  ctx.setTransform(ratio, 0, 0, ratio, 0, 0);
  ctx.fillStyle = "#08131f";
  ctx.fillRect(0, 0, w, h);
  return { ctx, w, h };
}

function drawAxes(canvas, { xMax, xStep, yMin, yMax, yStep, xLabel, yLabel, left = 66, right = 20, top = 18, bottom = 40 }) {
  const { ctx, w, h } = prepareCanvas(canvas);
  const pw = Math.max(1, w - left - right), ph = Math.max(1, h - top - bottom);
  const getX = x => left + (x / xMax) * pw;
  const getY = y => top + (1 - (y - yMin) / (yMax - yMin)) * ph;

  ctx.font = "10px 'DM Mono', Consolas, monospace";
  ctx.lineWidth = 1;
  ctx.strokeStyle = "rgba(133, 164, 188, .12)";
  ctx.fillStyle = "#6e8498";
  for (let y = yMin, index = 0; y <= yMax + yStep * 0.15; y += yStep, index++) {
    const yp = getY(y);
    ctx.beginPath(); ctx.moveTo(left, yp); ctx.lineTo(left + pw, yp); ctx.stroke();
    ctx.textAlign = "right";
    ctx.fillText(Number.isInteger(y) ? String(y) : y.toFixed(1), left - 9, yp + 3);
  }
  for (let x = 0; x <= xMax + xStep * 0.15; x += xStep) {
    const xp = getX(x);
    ctx.beginPath(); ctx.moveTo(xp, top); ctx.lineTo(xp, top + ph); ctx.stroke();
    ctx.textAlign = "center";
    ctx.fillText(Number.isInteger(x) ? String(x) : x.toFixed(1), xp, top + ph + 17);
  }
  ctx.strokeStyle = "rgba(164, 190, 215, .28)";
  ctx.strokeRect(left, top, pw, ph);
  ctx.fillStyle = "#8094a7";
  ctx.textAlign = "center";
  ctx.fillText(xLabel, left + pw / 2, h - 8);
  ctx.save(); ctx.translate(16, top + ph / 2); ctx.rotate(-Math.PI / 2);
  ctx.fillText(yLabel, 0, 0); ctx.restore();

  return { ctx, w, h, left, top, pw, ph, getX, getY, clipBox: () => { ctx.beginPath(); ctx.rect(left, top, pw, ph); ctx.clip(); } };
}

function drawEmptyChart(canvas, message) {
  const { ctx, w, h } = prepareCanvas(canvas);
  ctx.fillStyle = "#6e8498";
  ctx.textAlign = "center";
  ctx.font = "12px 'Manrope', sans-serif";
  ctx.fillText(message, w / 2, h / 2);
}

// 1. Vẽ đồ thị miền thời gian
function drawTimeChart() {
  const data = (isFrozen && recSamples.length > 1) ? recSamples : samples;
  if (data.length < 2) return drawEmptyChart($("timeChart"), "Kết nối ESP32 để xem tín hiệu theo thời gian");

  const winSec = parseInt($("timeWindowSelect").value) || 30;
  const tEnd = data[data.length - 1].time, tStart = tEnd - winSec * 1000;
  const pts = data.filter(d => d.time >= tStart);
  if (pts.length < 2) return drawEmptyChart($("timeChart"), "Đang chờ thêm mẫu...");

  let yMin = Infinity, yMax = -Infinity;
  pts.forEach(d => { yMin = Math.min(yMin, d.raw, d.filtered); yMax = Math.max(yMax, d.raw, d.filtered); });
  const margin = Math.max(8, (yMax - yMin) * 0.16);
  yMin = Math.max(0, Math.floor((yMin - margin) / 10) * 10);
  yMax = Math.min(4095, Math.ceil((yMax + margin) / 10) * 10);
  if (yMax <= yMin) yMax = yMin + 10;
  const yStep = Math.max(1, Math.ceil((yMax - yMin) / 5 / 10) * 10);

  const p = drawAxes($("timeChart"), {
    xMax: winSec, xStep: winSec <= 10 ? 2 : 5, yMin, yMax, yStep,
    xLabel: "Thời gian (s)", yLabel: "Giá trị ADC"
  });
  p.ctx.save(); p.clipBox();
  const toX = t => (t - tStart) / 1000;

  p.ctx.beginPath(); p.ctx.strokeStyle = "#ffb45d"; p.ctx.lineWidth = 1.05; p.ctx.globalAlpha = .72;
  pts.forEach((d, i) => i === 0 ? p.ctx.moveTo(p.getX(toX(d.time)), p.getY(d.raw)) : p.ctx.lineTo(p.getX(toX(d.time)), p.getY(d.raw)));
  p.ctx.stroke();

  p.ctx.beginPath(); p.ctx.strokeStyle = "#40e0c0"; p.ctx.lineWidth = 1.8;
  p.ctx.shadowColor = "rgba(64, 224, 192, .48)"; p.ctx.shadowBlur = 7;
  pts.forEach((d, i) => i === 0 ? p.ctx.moveTo(p.getX(toX(d.time)), p.getY(d.filtered)) : p.ctx.lineTo(p.getX(toX(d.time)), p.getY(d.filtered)));
  p.ctx.stroke(); p.ctx.shadowBlur = 0;
  p.ctx.restore();
}

// So sánh phổ FFT của cùng một đoạn tín hiệu trước và sau lọc.
function drawFrequencyChart() {
  const data = (isFrozen && recSamples.length >= 16) ? recSamples : samples;
  const requestedCount = Math.min(data.length, isFrozen ? MAX_FFT_SAMPLES : MAX_SAMPLES);
  const count = 2 ** Math.floor(Math.log2(requestedCount));
  if (count < 16) {
    $("sampleCount").textContent = `${data.length}`;
    $("spectrumMeta").textContent = `Đang chờ đủ mẫu · cần ít nhất 16 mẫu, hiện có ${data.length}`;
    return drawEmptyChart($("frequencyChart"), "Đang thu mẫu để dựng phổ FFT...");
  }

  const windowData = data.slice(-count);
  const rawRes = computeSpectrum(windowData, "raw");
  const filteredRes = computeSpectrum(windowData, "filtered");
  if (!rawRes || !filteredRes) return drawEmptyChart($("frequencyChart"), "Không thể tính phổ từ đoạn tín hiệu này");
  lastFftResult = filteredRes;
  lastRawSpectrum = rawRes;

  $("sampleRate").innerHTML = `${filteredRes.fs.toFixed(1)} <small>Hz</small>`;
  $("sampleCount").textContent = count.toLocaleString("vi-VN");
  $("resolutionLabel").textContent = `Δf = ${(filteredRes.fs / count).toFixed(3)} Hz/bin`;
  const segmentLabel = isFrozen
    ? `${count.toLocaleString("vi-VN")} mẫu FFT cuối / ${data.length.toLocaleString("vi-VN")} mẫu phiên`
    : `${count.toLocaleString("vi-VN")} mẫu FFT`;
  $("spectrumMeta").textContent = `${segmentLabel} · Fs = ${filteredRes.fs.toFixed(1)} Hz · Cửa sổ Hann`;

  if (filteredRes.dominantPeak) {
    const peak = filteredRes.dominantPeak;
    $("peakFreq").innerHTML = `${peak.freq.toFixed(2)} <small>Hz</small>`;
    $("peakBadge").textContent = `${peak.freq.toFixed(2)} Hz · biên độ ${peak.amp.toFixed(2)} ADC`;
  } else {
    $("peakFreq").innerHTML = "— <small>Hz</small>";
    $("peakBadge").textContent = "Chưa phát hiện đỉnh rõ";
  }

  if (!isFrozen && performance.now() - lastTableTime > 1000) {
    lastTableTime = performance.now();
    renderSpectrumTable(filteredRes, false);
  }

  const maxF = parseFloat($("freqSpanSelect").value) || 100;
  let maxA = 0;
  rawRes.spectrum.forEach(s => { if (s.freq <= maxF) maxA = Math.max(maxA, s.amp); });
  filteredRes.spectrum.forEach(s => { if (s.freq <= maxF) maxA = Math.max(maxA, s.amp); });
  maxA = Math.max(1, Math.ceil(maxA * 1.18 / 5) * 5);
  const yStep = Math.max(1, Math.ceil(maxA / 5 / 5) * 5);
  const p = drawAxes($("frequencyChart"), {
    xMax: maxF, xStep: maxF <= 20 ? 2 : maxF <= 50 ? 5 : 10,
    yMin: 0, yMax: maxA, yStep, xLabel: "Tần số (Hz)", yLabel: "Biên độ ADC"
  });

  p.ctx.save(); p.clipBox();
  const drawSpectrum = (result, color, glow, fill) => {
    const points = result.spectrum.filter(point => point.freq <= maxF);
    if (points.length < 2) return;
    p.ctx.beginPath();
    p.ctx.moveTo(p.getX(points[0].freq), p.getY(0));
    points.forEach(point => p.ctx.lineTo(p.getX(point.freq), p.getY(point.amp)));
    p.ctx.lineTo(p.getX(points[points.length - 1].freq), p.getY(0));
    p.ctx.closePath();
    if (fill) {
      const gradient = p.ctx.createLinearGradient(0, p.top, 0, p.top + p.ph);
      gradient.addColorStop(0, "rgba(64,224,192,.22)");
      gradient.addColorStop(1, "rgba(64,224,192,0)");
      p.ctx.fillStyle = gradient;
      p.ctx.fill();
    }
    p.ctx.beginPath();
    points.forEach((point, index) => index === 0 ? p.ctx.moveTo(p.getX(point.freq), p.getY(point.amp)) : p.ctx.lineTo(p.getX(point.freq), p.getY(point.amp)));
    p.ctx.strokeStyle = color;
    p.ctx.lineWidth = fill ? 1.8 : 1.15;
    p.ctx.globalAlpha = fill ? 1 : .78;
    p.ctx.shadowColor = glow;
    p.ctx.shadowBlur = fill ? 9 : 2;
    p.ctx.stroke();
    p.ctx.shadowBlur = 0;
    p.ctx.globalAlpha = 1;
  };

  // Mốc lưới điện giúp đối chiếu hiện tượng nhấp nháy ở đèn chiếu sáng.
  if (maxF >= 50) {
    const x50 = p.getX(50);
    p.ctx.beginPath(); p.ctx.setLineDash([4, 5]); p.ctx.strokeStyle = "rgba(187,143,255,.75)"; p.ctx.lineWidth = 1;
    p.ctx.moveTo(x50, p.top); p.ctx.lineTo(x50, p.top + p.ph); p.ctx.stroke(); p.ctx.setLineDash([]);
    p.ctx.fillStyle = "#c7a8ff"; p.ctx.font = "10px 'DM Mono', Consolas, monospace"; p.ctx.textAlign = "center";
    p.ctx.fillText("50 Hz", x50, p.top + 12);
  }

  drawSpectrum(rawRes, "#ffb45d", "rgba(255,180,93,.25)", false);
  drawSpectrum(filteredRes, "#40e0c0", "rgba(64,224,192,.55)", true);

  const peak = filteredRes.dominantPeak;
  if (peak && peak.freq <= maxF) {
    const px = p.getX(peak.freq), py = p.getY(peak.amp);
    p.ctx.beginPath(); p.ctx.arc(px, py, 4, 0, Math.PI * 2);
    p.ctx.fillStyle = "#d6fff4"; p.ctx.shadowColor = "#40e0c0"; p.ctx.shadowBlur = 12; p.ctx.fill(); p.ctx.shadowBlur = 0;
    p.ctx.fillStyle = "#b9f6e8"; p.ctx.font = "10px 'DM Mono', Consolas, monospace";
    p.ctx.textAlign = px > p.left + p.pw * .72 ? "right" : "left";
    p.ctx.fillText(`${peak.freq.toFixed(2)} Hz`, px + (p.ctx.textAlign === "right" ? -8 : 8), Math.max(p.top + 16, py - 10));
  }
  p.ctx.restore();
}

function scheduleRender() {
  if (renderPending) return;
  renderPending = true;
  requestAnimationFrame(() => {
    renderPending = false;
    drawTimeChart();
    drawFrequencyChart();
  });
}

// =====================================================================
// KHỐI 4: BẢNG PHỔ, GHI 30S & XUẤT BÁO CÁO (CSV/PNG)
// =====================================================================
function renderSpectrumTable(res, isDone) {
  const tbody = $("spectrumTableBody");
  $("tableModeBadge").textContent = isDone ? "30 GIÂY" : "LIVE";
  if (res) {
    $("spectrumTableMeta").textContent = isDone
      ? `FFT ${res.count} mẫu cuối của phiên ${recSamples.length.toLocaleString("vi-VN")} mẫu · Fs = ${res.fs.toFixed(1)} Hz · Δf = ${(res.fs / res.count).toFixed(3)} Hz/bin`
      : `Chế độ thời gian thực: ${res.count} mẫu · Fs = ${res.fs.toFixed(1)} Hz`;
  }
  if (!res || !res.peaks.length) {
    tbody.innerHTML = `<tr><td colspan="7" class="empty-table">Chưa phát hiện đỉnh phổ rõ trong đoạn mẫu này.</td></tr>`;
    return;
  }

  tbody.innerHTML = res.peaks.slice(0, 8).map((p, i) => `
    <tr class="${i === 0 ? 'top-peak' : ''}">
      <td><strong class="peak-rank">${String(i + 1).padStart(2, "0")}</strong></td>
      <td><strong>${p.freq.toFixed(2)} Hz</strong></td>
      <td>${p.period.toFixed(2)} s</td>
      <td>${p.amp.toFixed(1)}</td>
      <td>${p.powerDb.toFixed(1)} dB</td>
      <td>${p.energyPercent.toFixed(1)} %</td>
      <td><span class="badge-tag tag-${p.type}">${p.label}</span></td>
    </tr>`).join("");
}

function start30sRecording() {
  if (!reading || !port) return alert("Vui lòng kết nối ESP32 trước khi ghi 30s!");
  isRecording = true; isFrozen = false; recStart = performance.now(); recSamples = [];
  $("btnRecord30s").textContent = "Đang ghi · bấm để dừng";
  $("btnRecord30s").className = "btn-record recording";
  $("btnResetRecord").classList.remove("hidden");
  $("btnExportCsv").disabled = true;
}

function finish30sRecording() {
  isRecording = false; isFrozen = true;
  $("btnRecord30s").textContent = "Ghi phiên mới";
  $("btnRecord30s").className = "btn-record";
  $("btnResetRecord").textContent = "Quay lại live";
  $("btnExportCsv").disabled = true;
  $("recordProgressBar").style.width = "100%";
  $("recordTimerLabel").textContent = `Đã ghi xong · ${recSamples.length.toLocaleString("vi-VN")} mẫu`;

  const count = 2 ** Math.floor(Math.log2(Math.min(recSamples.length, MAX_FFT_SAMPLES)));
  if (count >= 16) {
    const windowData = recSamples.slice(-count);
    const rawRes = computeSpectrum(windowData, "raw");
    const filteredRes = computeSpectrum(windowData, "filtered");
    if (rawRes && filteredRes) {
      lastRawSpectrum = rawRes;
      lastFftResult = filteredRes;
      renderSpectrumTable(filteredRes, true);
      $("btnExportCsv").disabled = false;
    }
  }
  scheduleRender();
}

function resetRecording() {
  isRecording = false; isFrozen = false; recSamples = [];
  $("btnRecord30s").textContent = "Bắt đầu ghi";
  $("btnRecord30s").className = "btn-record";
  $("btnResetRecord").classList.add("hidden");
  $("btnExportCsv").disabled = true;
  $("recordProgressBar").style.width = "0%";
  $("recordTimerLabel").textContent = "Sẵn sàng ghi · 00:00 / 00:30";
  $("recordSampleCount").textContent = "0 mẫu";
  scheduleRender();
}

function exportCsv() {
  if (!lastFftResult) return alert("Chưa có dữ liệu phổ để xuất file CSV!");
  const res = lastFftResult;
  let csv = "\uFEFF--- THONG TIN PHIEN PHAN TICH PHO LDR ---\r\n";
  csv += `Ngay do,${new Date().toLocaleString("vi-VN")}\r\n`;
  csv += `Fs (Hz),${res.fs.toFixed(2)},So mau,${res.count},Do phan giai delta_f,${(res.fs / res.count).toFixed(4)}\r\n\r\n`;
  csv += "Hang,Tan so (Hz),Chu ky (s),Bien do sau loc,Cong suat (dB),Ty le nang luong (%),Phan loai\r\n";
  res.peaks.forEach((p, i) => {
    csv += `#${i + 1},${p.freq.toFixed(3)},${p.period.toFixed(3)},${p.amp.toFixed(2)},${p.powerDb.toFixed(2)},${p.energyPercent.toFixed(2)}%,${p.label}\r\n`;
  });
  csv += "\r\nBin,Tan so (Hz),Bien do tho,Bien do sau loc\r\n";
  res.spectrum.forEach((s, i) => { csv += `${s.bin},${s.freq.toFixed(3)},${(lastRawSpectrum?.spectrum[i]?.amp ?? 0).toFixed(3)},${s.amp.toFixed(3)}\r\n`; });

  const a = document.createElement("a");
  a.href = URL.createObjectURL(new Blob([csv], { type: "text/csv;charset=utf-8;" }));
  a.download = `Bang_Pho_LDR_30s_${Date.now()}.csv`;
  a.click();
}

function exportPng(canvasId, defaultName) {
  const a = document.createElement("a");
  a.href = $(canvasId).toDataURL("image/png");
  a.download = defaultName;
  a.click();
}

// =====================================================================
// KHỐI 5: GIAO TIẾP WEB SERIAL & SỰ KIỆN
// =====================================================================
function receiveLine(line) {
  line = line.trim();
  if (!line) return;

  let raw = NaN, devTime = null;
  try {
    const d = JSON.parse(line);
    if (d && typeof d === "object") {
      raw = Number(d.raw ?? d.ldr);
      if (d.timeMs !== undefined) devTime = Number(d.timeMs);
    }
  } catch (_) {}

  if (!Number.isFinite(raw)) {
    raw = Number(line.replace(/^LDR\s*[:=]\s*/i, ""));
  }
  if (!Number.isFinite(raw)) return;

  let time;
  if (devTime !== null) {
    if (prevDeviceTime === null) prevDeviceTime = devTime;
    else {
      const diff = (devTime - prevDeviceTime + 4294967296) % 4294967296;
      if (diff > 2147483648) { samples.length = 0; elapsedTime = 0; }
      else elapsedTime += Math.max(1, diff);
      prevDeviceTime = devTime;
    }
    time = elapsedTime;
  } else {
    time = performance.now();
  }

  const s = { time, raw, filtered: filterSample(raw) };
  samples.push(s);
  if (samples.length > MAX_SAMPLES) samples.shift();

  if (isRecording) {
    recSamples.push({ ...s });
    const ms = performance.now() - recStart;
    $("recordProgressBar").style.width = Math.min(100, (ms / RECORD_MS) * 100).toFixed(1) + "%";
    $("recordTimerLabel").textContent = `Đang ghi · còn ${(Math.max(0, RECORD_MS - ms) / 1000).toFixed(1)} giây`;
    $("recordSampleCount").textContent = `${recSamples.length.toLocaleString("vi-VN")} mẫu`;
    if (ms >= RECORD_MS) finish30sRecording();
  }

  $("latestRaw").textContent = Math.round(raw);
  $("latestFiltered").textContent = Math.round(s.filtered);
  scheduleRender();
}

async function connectSerial() {
  if (!("serial" in navigator)) return alert("Vui lòng dùng Chrome hoặc Edge trên máy tính để hỗ trợ Web Serial!");
  $("connect").disabled = true;
  try {
    port = await navigator.serial.requestPort();
    await port.open({ baudRate: BAUD_RATE });
    samples.length = 0;
    recSamples = [];
    filterState = { iir: null, window: [] };
    prevDeviceTime = null;
    elapsedTime = 0;
    isRecording = false;
    isFrozen = false;
    lastFftResult = null;
    lastRawSpectrum = null;
    $("btnExportCsv").disabled = true;
    $("btnRecord30s").textContent = "Bắt đầu ghi";
    $("btnRecord30s").className = "btn-record";
    $("btnResetRecord").classList.add("hidden");
    $("sampleRate").innerHTML = "— <small>Hz</small>";
    $("sampleCount").textContent = "—";
    $("latestRaw").textContent = "—";
    $("latestFiltered").textContent = "—";
    $("peakFreq").innerHTML = "— <small>Hz</small>";
    $("peakBadge").textContent = "Chưa có đỉnh";
    $("resolutionLabel").textContent = "Δf = — Hz/bin";
    $("spectrumTableMeta").textContent = "Kết nối ESP32 để bắt đầu phân tích trực tiếp.";
    $("spectrumTableBody").innerHTML = `<tr><td colspan="7" class="empty-table">Đang chờ dữ liệu phổ.</td></tr>`;
    $("recordProgressBar").style.width = "0%";
    $("recordTimerLabel").textContent = "Sẵn sàng ghi · 00:00 / 00:30";
    $("recordSampleCount").textContent = "0 mẫu";
    reading = true;
    $("connect").textContent = "Ngắt kết nối";
    $("status").textContent = "Đã kết nối · 115200 baud";
    $("status").className = "status-pill connected";

    const decoder = new TextDecoder();
    reader = port.readable.getReader();
    while (reading) {
      const { value, done } = await reader.read();
      if (done) break;
      textBuffer += decoder.decode(value, { stream: true });
      const lines = textBuffer.split("\n");
      textBuffer = lines.pop();
      lines.forEach(receiveLine);
    }
  } catch (err) {
    $("status").textContent = err.message || "Không kết nối được";
    $("status").className = "status-pill";
  } finally {
    try { reader?.releaseLock(); } catch (_) {}
    $("connect").disabled = false;
  }
}

async function disconnectSerial() {
  reading = false;
  if (isRecording) resetRecording();
  try { await reader?.cancel(); } catch (_) {}
  try { await port?.close(); } catch (_) {}
  port = null;
  $("connect").textContent = "Kết nối ESP32";
  $("status").textContent = "Đã ngắt kết nối";
  $("status").className = "status-pill";
}

// Đăng ký sự kiện giao diện
$("connect").addEventListener("click", () => port ? disconnectSerial() : connectSerial());
$("btnRecord30s").addEventListener("click", () => isRecording ? finish30sRecording() : start30sRecording());
$("btnResetRecord").addEventListener("click", resetRecording);
$("btnExportCsv").addEventListener("click", exportCsv);
$("btnExportTimePng")?.addEventListener("click", () => exportPng("timeChart", "Tin_hieu_thoi_gian_ADC_LDR.png"));
$("btnExportPlotPng")?.addEventListener("click", () => exportPng("frequencyChart", "Pho_tan_so_so_sanh_Raw_Filtered.png"));

$("filterAlpha").addEventListener("input", e => { $("alphaVal").textContent = parseFloat(e.target.value).toFixed(2); refreshFilterSettings(); });
$("windowSize").addEventListener("input", e => { $("windowVal").textContent = e.target.value; refreshFilterSettings(); });
$("filterMode").addEventListener("change", refreshFilterSettings);
$("timeWindowSelect").addEventListener("change", scheduleRender);
$("freqSpanSelect")?.addEventListener("change", scheduleRender);
window.addEventListener("resize", scheduleRender);

refreshFilterSettings();
scheduleRender();

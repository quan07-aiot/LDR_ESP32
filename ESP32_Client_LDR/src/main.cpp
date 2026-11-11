#include <Arduino.h>

// ======================== LDR SENSOR ========================
// LDR analog input. Chọn GPIO tương thích ADC trên bo mạch đang sử dụng.
const int LDR_SENSOR_PIN = 5;

// Lấy mẫu cách đều 5 ms: Fs = 200 Hz, giới hạn Nyquist là 100 Hz.
const uint32_t SAMPLE_INTERVAL_US = 5000;
uint32_t lastSampleMicros = 0;

// ======================== LED CHASING ========================
// 3 LED nối vào GPIO 16, 17, 18 – chạy vòng lặp khi trời tối.
const int LED_PINS[] = {16, 17, 18};
const int NUM_LEDS   = 3;

// Ngưỡng phân biệt sáng/tối (ADC 12-bit: 0-4095).
// Giá trị LDR thấp = nhiều ánh sáng, giá trị cao = tối.
// Điều chỉnh ngưỡng này phù hợp với mạch thực tế của bạn.
const uint16_t DARK_THRESHOLD = 2000;

// Tốc độ chuyển LED (ms) – mỗi LED sáng trong khoảng thời gian này.
const uint32_t LED_CHASE_INTERVAL_MS = 300;

uint32_t lastLedChangeMs = 0;   // Thời điểm đổi LED lần cuối
int      currentLedIndex = 0;   // LED đang sáng (0, 1, 2)
bool     isDark          = false; // Trạng thái tối/sáng hiện tại

// Tắt tất cả LED
void allLedsOff() {
  for (int i = 0; i < NUM_LEDS; i++) {
    digitalWrite(LED_PINS[i], LOW);
  }
}

void setup() {
  Serial.begin(115200);

  // ---------- LDR ----------
  analogReadResolution(12); // ADC 12-bit: 0-4095.
#if defined(ADC_11db)
  analogSetAttenuation(ADC_11db);
#elif defined(ADC_12db)
  analogSetAttenuation(ADC_12db);
#endif

  pinMode(LDR_SENSOR_PIN, INPUT);

  // ---------- LEDs ----------
  for (int i = 0; i < NUM_LEDS; i++) {
    pinMode(LED_PINS[i], OUTPUT);
    digitalWrite(LED_PINS[i], LOW);
  }

  delay(100);

  lastSampleMicros = micros();
  lastLedChangeMs  = millis();

  Serial.println("{\"status\":\"ESP32 LDR + LED Chase Ready\",\"fs\":200}");
}

void loop() {
  const uint32_t currentMicros = micros();
  const uint32_t currentMs     = millis();

  // ======================== LDR Sampling (giữ nguyên) ========================
  if (currentMicros - lastSampleMicros >= SAMPLE_INTERVAL_US) {
    // Giữ nhịp lấy mẫu đều; chỉ gửi giá trị ADC gốc để web áp dụng bộ lọc.
    lastSampleMicros += SAMPLE_INTERVAL_US;
    const uint16_t rawSample = analogRead(LDR_SENSOR_PIN);
    const uint32_t sampleTimeMs = millis();

    // Cập nhật trạng thái sáng/tối dựa trên giá trị LDR
    isDark = (rawSample >= DARK_THRESHOLD);

    Serial.printf("{\"raw\":%u,\"timeMs\":%lu,\"dark\":%s}\n",
                  static_cast<unsigned int>(rawSample),
                  static_cast<unsigned long>(sampleTimeMs),
                  isDark ? "true" : "false");
  }

  // ======================== LED Chase ========================
  if (isDark) {
    // Trời tối → chạy vòng lặp LED: 16 → 17 → 18 → 16 ...
    if (currentMs - lastLedChangeMs >= LED_CHASE_INTERVAL_MS) {
      lastLedChangeMs = currentMs;

      // Tắt tất cả LED trước
      allLedsOff();

      // Bật LED hiện tại
      digitalWrite(LED_PINS[currentLedIndex], HIGH);

      // Chuyển sang LED tiếp theo (vòng lặp)
      currentLedIndex = (currentLedIndex + 1) % NUM_LEDS;
    }
  } else {
    // Trời sáng → tắt hết LED, reset vị trí
    allLedsOff();
    currentLedIndex = 0;
    lastLedChangeMs = currentMs; // Reset timer để khi tối lại bắt đầu ngay
  }
}

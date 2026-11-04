#include <Arduino.h>

// LDR analog input. Chọn GPIO tương thích ADC trên bo mạch đang sử dụng.
const int LDR_SENSOR_PIN = 5;

// Lấy mẫu cách đều 5 ms: Fs = 200 Hz, giới hạn Nyquist là 100 Hz.
const uint32_t SAMPLE_INTERVAL_US = 5000;
uint32_t lastSampleMicros = 0;

void setup() {
  Serial.begin(115200);

  analogReadResolution(12); // ADC 12-bit: 0-4095.
#if defined(ADC_11db)
  analogSetAttenuation(ADC_11db);
#elif defined(ADC_12db)
  analogSetAttenuation(ADC_12db);
#endif

  pinMode(LDR_SENSOR_PIN, INPUT);
  delay(100);

  lastSampleMicros = micros();
  Serial.println("{\"status\":\"ESP32 LDR Ready\",\"fs\":200}");
}

void loop() {
  const uint32_t currentMicros = micros();
  if (currentMicros - lastSampleMicros < SAMPLE_INTERVAL_US) return;

  // Giữ nhịp lấy mẫu đều; chỉ gửi giá trị ADC gốc để web áp dụng bộ lọc.
  lastSampleMicros += SAMPLE_INTERVAL_US;
  const uint16_t rawSample = analogRead(LDR_SENSOR_PIN);
  const uint32_t sampleTimeMs = millis();

  Serial.printf("{\"raw\":%u,\"timeMs\":%lu}\n",
                static_cast<unsigned int>(rawSample),
                static_cast<unsigned long>(sampleTimeMs));
}

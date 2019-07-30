#include <TinyGPS++.h>
#include <SoftwareSerial.h>

#define GPSbaud 9600

TinyGPSPlus GPS;
SoftwareSerial GPSbaudrate(7, 8);

int stp1 = 6;
int dir1 = 5;
int stp2 = 4;
int dir2 = 3;
int ena = 2;

float Az;
float Alt;

String bar;

void setup() {
  pinMode(stp1, OUTPUT);
  pinMode(dir1, OUTPUT);
  pinMode(stp2, OUTPUT);
  pinMode(dir2, OUTPUT);
  pinMode(ena, OUTPUT);
  Serial.begin(9600);
  GPSbaudrate.begin(GPSbaud);
  //delay(30000);

}
void loop() {
  while(true) 
  {
    smartDelay (1000);   
    //Serial.println(GPS.location.lat());
    if (GPS.location.lat() != 0 & GPS.location.lng() != 0) 
    {
      Serial.print (GPS.location.lat(), 6);
      Serial.print (" ");
      Serial.println (GPS.location.lng(), 6);
      bar = Serial.readStringUntil('\n');
      if(bar == "foo"){
       // smartDelay (1000);
        break;
      }
      break;
    }
    break;   
  }
  while(true){
  if(Serial.available() > 0){
   String input = Serial.readStringUntil('\n');
   Az = getValue(input,' ', 0).toFloat(); // ento azimut
   Alt = getValue(input,' ',1).toFloat(); // ento visota
  }
}
 


  analogWrite (ena, 0); // работа двгателей, вставить метод (Тима)
  for (int i = 0; i < 100; i++) //вставить переменную для первого двигателя (6 пин)
  {
    digitalWrite(stp1, HIGH);
    delay (10);
    digitalWrite(stp1, LOW);
    delay (10);
  }
  delay (2000);
  for (int i = 0; i < 200; i++) //вставить переменную для вторго двигателя (4 пин)
  {
    digitalWrite(stp2, HIGH);
    delay (10);
    digitalWrite(stp2, LOW);
    delay (10);
  }
  while (true) {}
}

String getValue(String data, char separator, int index)
{
    int found = 0;
    int strIndex[] = { 0, -1 };
    int maxIndex = data.length() - 1;

    for (int i = 0; i <= maxIndex && found <= index; i++) {
        if (data.charAt(i) == separator || i == maxIndex) {
            found++;
            strIndex[0] = strIndex[1] + 1;
            strIndex[1] = (i == maxIndex) ? i+1 : i;
        }
    }
    return found > index ? data.substring(strIndex[0], strIndex[1]) : "";
}

static void smartDelay(unsigned long ms)
{
  unsigned long start = millis();
  do {
    while (GPSbaudrate.available())
      GPS.encode(GPSbaudrate.read());
  }
  while (millis() - start < ms);
}

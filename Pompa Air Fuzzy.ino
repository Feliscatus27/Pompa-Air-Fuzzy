#include<math.h>

float A, B;
int sel_;

float a1, b1a, b1b, c1;

float L1, L2, L3, L4, L5, L6, L7;
float M1, M2, M3, M4, M5, M6, M7;

float fudeket[4] = {0, 0, 3, 6};
float fulumayan[3] = {3, 7, 11};
float fujauh[4] = {7, 11, 15, 15};

float fylambat[3] = {0, 0, 120}; 
float fysedang[3] = {25, 127.5, 235};
float fycepat[3] = {155, 255, 255};

float inputf, outputf;

//Membership function untuk dekat
float FiN() {
  if (inputf < fudeket[2]) {
    return 1;
  }
  else if (inputf >= fudeket[2] && inputf <= fudeket[3]) {
    return (fudeket[3] - inputf) / (fudeket[3] - fudeket[2]);
  }
  else if (inputf > fudeket[3]) {
    return 0;
  }
}

//Membership function untuk Lumayan
float FiZ() {

  if (inputf < fulumayan[0]) {
    return 0;
  }
  else if (inputf >= fulumayan[0] && inputf <= fulumayan[1]) {
    return (inputf - fulumayan[0]) / (fulumayan[1] - fulumayan[0]);
  }
  else if (inputf >= fulumayan[1] && inputf <= fulumayan[2]) {
    return (fulumayan[2] - inputf) / (fulumayan[2] - fulumayan[1]);
  }
  else if (inputf > fulumayan[2]) {
    return 0;
  }

}

//Membership function untuk jauh
float FiP() {
  if (inputf < fujauh[0]) {
    return 0;
  }
  else if (inputf >= fujauh[0] && inputf <= fujauh[1]) {
    return (inputf - fujauh[0]) / (fujauh[1] - fujauh[0]);
  }
  else if (inputf > fujauh[1]) {
    return 1;
  }
}

//Membership Function untuk output lambat
float FoN() {
  if (outputf < fylambat[1]) {
    return 1;
  }
  else if (outputf >= fylambat[1] && outputf <= fylambat[2]) {
    return (fylambat[2] - outputf) / (fylambat[2] - fylambat[1]);
  }
  else if (outputf > fylambat[2]) {
    return 0;
  }
}

//Membership Function untuk output sedang
float FoZ() {

  if (outputf < fysedang[0]) {
    return 0;
  }
  else if (outputf >= fysedang[0] && outputf <= fysedang[1]) {
    return (outputf - fysedang[0]) / (fysedang[1] - fysedang[0]);
  }
  else if (outputf >= fysedang[1] && outputf <= fysedang[2]) {
    return (fysedang[2] - outputf) / (fysedang[2] - fysedang[1]);
  }
  else if (outputf > fysedang[2]) {
    return 0;
  }

}

//Membership Function untuk output cepat
float FoP() {
  if (outputf < fycepat[0]) {
    return 0;
  }
  else if (outputf >= fycepat[0] && outputf <= fycepat[1]) {
    return (outputf - fycepat[0]) / (fycepat[1] - fycepat[0]);
  }
  else if (outputf > fycepat[1]) {
    return 1;
  }
}

void implikasi () {
  //sesuai dengan rule
  // if deket then lambat                           //misal inputnya 2 ==> fin =1 fiz dan fip = 0
  a1 = 120 - (FiN() * (fylambat[2] - fylambat[1])); //a1 = 0
  // if lumayan then sedang
  b1a = 25 + (FiZ() * (fysedang[1] - fysedang[0])); //b1a = 25
  b1b =  235 - (FiZ() * (fysedang[2] - fysedang[1])); //b1b = 235
  // if jauh then cepat
  c1 = 155 + (FiP() * (fycepat[1] - fycepat[0])); //c1 = 155
}

void luas_deffuzzy() {
  implikasi ();
  L1 = ((fylambat[2] - a1) * FiN()) / 2; //= 60
  L2 = (a1 - fylambat[0]) * FiN(); // = 0
  L3 = ((b1a - fysedang[0]) * FiZ()) / 2; // L3 = 0
  L4 = ((fysedang[2] - b1b) * FiZ()) / 2; // L4 = 0
  L5 = (b1b - b1a) * FiZ();               // L5 = 0
  L6 = ((c1 - fycepat[0]) * FiP()) / 2; // = 0
  L7 = (fycepat[2] - c1) * FiP(); // = 0
}

float f(float x) {

  if (B > 0 && sel_ == 0) {
    return ((x - A) / B) * x;
  }
  else if (B > 0 && sel_ == 1) {
    return ((A - x) / B) * x;
  }
  else {
    return A * x;
  }
}

/*Function definition to perform integration by Simpson's 1/3rd Rule */
float simpsons(float f(float x), float a, float b, float n) {
  float h, integral, x, sum = 0;
  int i;
  h = fabs(b - a) / n;
  for (i = 1; i < n; i++) {
    x = a + i * h;
    if (i % 2 == 0) {
      sum = sum + 2 * f(x);
    }
    else {
      sum = sum + 4 * f(x);
    }
  }
  integral = (h / 3) * (f(a) + f(b) + sum);
  return integral;
}

float fx(float limd, float limu, float a, float b, int sel) {
  int n, i = 2;
  float h, x, integral, eps = 0.1, integral_new;
  A = a;
  B = b;
  sel_ = sel;

  integral_new = simpsons(f, limd, limu, i);

  do {
    integral = integral_new;
    i = i + 2;
    integral_new = simpsons(f, limd, limu, i);
  } while (fabs(integral_new - integral) >= eps);

  return integral_new;
}

void moment() {
  luas_deffuzzy();
  M1 = fx(a1, fylambat[2], fylambat[2], (fylambat[2] - fylambat[0]), 1); //2400
  M2 = fx(fylambat[0], a1, FiN(), 0, 0); // 0
  M3 = fx(fysedang[0], b1a, fysedang[0], (fysedang[1] - fysedang[0]), 0); //0
  M4 = fx (b1b, fysedang[2], fysedang[2], (fysedang[2] - fysedang[1]), 1);//0
  M5 = fx (b1a, b1b, FiZ(), 0, 0); //0
  M6 = fx(fycepat[0], c1, fycepat[0], (fycepat[2] - fycepat[0]), 0);//0
  M7 = fx(c1, fycepat[2], FiP(), 0, 0);//0
}

float deffuzzyfikasi() {
  moment();
  return (M1 + M2 + M3 + M4 + M5 + M6 + M7) / (L1 + L2 + L3 + L4 + L5 + L6 + L7); //40
}
void setup() {
  Serial.begin(9600);
  pinMode(25, OUTPUT);
  pinMode(27, INPUT);
  pinMode(22, OUTPUT); // set pin PA0 sebagai output digital
  pinMode(23, OUTPUT); // set pin PA1 sebagai output digital
  pinMode(44, OUTPUT);
  Serial.println("---Program Fuzzy---");
}

// menghitung volume 
float volume(float t){
  float y = (((25 - t) /25) * 100);
  return y;
}

float dutycycle(){
  float duty ;
  duty = (int(deffuzzyfikasi())/255) * 100 ; 
  return duty ;
}

void loop() {
  long duration, distance_cm;
  int pwm ,t;
  //Code Motor L298N
  digitalWrite(22, HIGH); // set nilai "HIGH" pada pin PA0
  digitalWrite(23, LOW); // set nilai "LOW" pada pin PA1
  pwm = int((deffuzzyfikasi()+24));
  analogWrite(44, pwm ); // set nilai output PWM sebesar 255 (100%)

  //Code Ultrasonik
  digitalWrite(25, LOW); 
  delayMicroseconds(2); 
  digitalWrite(25, HIGH);
  delayMicroseconds(10); 
  digitalWrite(25, LOW);
  duration = pulseIn(27, HIGH);
  distance_cm = duration / 58; // konversi dari waktu ke jarak dalam satuan cm
  
  inputf = distance_cm;

      Serial.print("Input jarak :");
      Serial.print(inputf);
      Serial.println(" CM"); 
  
      Serial.print("volume :");
      Serial.print(volume(inputf));
      Serial.println("%");

      Serial.print("Output Deffuzyfikasi Mamdani :");
      Serial.println(pwm);
      
  delay(500);
  
}
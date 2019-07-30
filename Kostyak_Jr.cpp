#include <iostream>
#include <stdio.h>
#include <sstream>
#include <iterator>
#include "string.h"
#include <fstream>
#include <string>
#include <vector>
#include <chrono>
#include <time.h>
#include <ctime>
#include "Serial.h"
#include <cmath>

using namespace std;

/********************************************************Floats********************************************************/
 float CurrAz = 0, CurrAlt = 0, AzDiff, hours, AltDiff, rasumm , decsumm;
 float RA[100], AZ[100], DEC[100], ALT[100];
 float t0, a1, a2, a3, a4, s0, z;
 float pi = 3.14159265359;
 float dtor = pi/180;
 float LAT, LON;

/*********************************************************Ints*********************************************************/

 int v1, v2, md, i, result, choiceA, choiceB;
 int year, month, day, hrs, minutes, sec;
 int decdeg, decmm, decss;
 int rahh, ramm, rass;
 int strLineLen = 0;
 int listNums[100];
 int j = -1;

/*****************************************************CharsStrings*****************************************************/

 char* lastStrLine = new char[1024];
 char *port_num = "\\\\.\\COM8";
 char output[MAX_DATA_LENGTH];
 const char stop = '\n';
 vector<const char*> starNames;
 vector<const char*> starDec;
 vector<const char*> starRa;
 vector<string> names;
 vector<string> Dec;// Declination is degdeg:mm:ss
 vector<string> Ra; // Right ascension is hh:mm:ss

/**********************************************************************************************************************/
/*******************************************************FUNCTIONS******************************************************/
/**********************************************************************************************************************/

float day2000(int y, int m, int d, float h){
    float a, d1, b, c, greg;
    greg = y*10000 + m*100 + d;
    if (m == 1 || m == 2) {
        y = y - 1;
        m = m + 12;
    }

    if (greg > 15821004) {
        a = floor(y/100);
        b = 2 - a  + floor(a/4);
    }
    else {
        b = 0;
    }
    c = floor(365.25 * y);
    d1 = floor(30.6001 * (m + 1));
    return (b + c + d1 -730550.5 + d + h/24);
}
float JulianDay(int yy, int mm, int dd) {
    float tulos;
    tulos =367*yy-floor(7*(yy+floor((mm+9)/12))/4)
           -floor(3*(floor((yy+(mm-9)/7)/100)+1)/4)
           +floor(275*mm/9)+dd+1721028.5;
    return tulos;
}
float Greenwich (int yy, int mm, int dd, float time) {
    // calculates Greenwich sidereal time in hours
    float T,JD,GWT0,GW,GMST;
    JD = JulianDay(yy,mm,dd);
    // JD is Julian Day
    T = (JD - 2451545.0)/36525;
    // GWT0 is Greenwich Sidereal Time in hours
    GWT0 = (24110.54841+8640184.812866*T+0.093104*T*T-0.0000062*T*T*T)/3600;
    if (GWT0 >= 24)
        GW = GWT0 - 24 * floor(GWT0/24);
    if (GWT0 <= - 24)
        GW = 24 -(abs(GWT0) - 24 * floor(abs(GWT0)/24));
    GMST = GW + 1.00273790935 * time;
    if (GMST >= 24)
        GMST = GMST - 24;
    if (GMST < 0 )
        GMST = GMST + 24;
    return GMST;
}
float isleap(int y) {
    int a;
    //assume not a leap year...
    a = 0;
    //...flag leap year candidates...
    if (y % 4 == 0) a = 1;
    //...if year is a century year then not leap...
    if (y % 100 ==0 ) a = 0;
    //...except if century year divisible by 400...
    if (y % 400 == 0) a = 1;
    //...and so done according to Gregory's wishes
    return a;
}
float goodmonthday(int y, int m, int d) {
    int a, leap;
    leap = isleap(y);
    //	assume OK
    a = 1;
    //	first deal with zero day number!
    if (d == 0) a = 0;
    //	Sort Feburary next
    if ((m==2) && (leap ==1) && (d > 29)) a= 0;
    if ((m==2) && (d > 28) && (leap ==0))	a = 0;
    //	then the rest of the months - 30 days...
    if(((m==4) || (m == 6) || (m == 9) || (m==11)) && d > 30) a = 0;
    //	...31 days...
    if (d > 31) a = 0;
    //	...and so done.
    return a;
}
float round(float num, float dp) {
    //   dp = (!dp ? 2: dp);
    return round (num * pow(10, dp)) / pow(10, dp);
}
float datan(float x) {
    return 180/ pi * atan(x);
}
float datan2(float y, float x) {
    float a;
    if ((x == 0) && (y == 0)) {
        return 0;
    }
    else	{
        a = datan(y / x);
        if (x < 0) {
            a = a + 180;
        }
        if (y < 0 && x > 0) {
            a = a + 360;
        }
        return a;
    }
}
float dasin(float x) {
    return 180/ pi * asin(x);
}
float dacos(float x) {
    return 180/ pi * acos(x);
}
float ipart(float x) {
    float a;
    if (x> 0) {
        a = floor(x);
    }
    else {
        a = ceil(x);
    }
    return a;
}
float range(float x) {
    float a, b;
    b = x / 360;
    a = 360 * (b - ipart(b));
    if (a  < 0 ) {
        a = a + 360;
    }
    return a;
}
float dsin(float x) {
    return sin(pi / 180 * x);
}
float dcos(float x) {
    return cos(pi / 180 * x);
}
float dtan(float x) {
    return tan(pi / 180 * x);
}

/*********************************************************СOOL**********************************************************/
/*******************************************************FUNCTIONS*******************************************************/
/***********************************************************************************************************************/

float moon_position(int d, int m, int y, float h) {

    float g, days, t ,L1, M1, C1, V1, Ec1, R1, Th1, Om1, Lam1, Obl, Ra1, Dec1;
    float F, L2, Om2, M2, D, D2, R2, R3, Bm, Lm, HLm, HBm, Ra2, Dec2, EL, EB, W, X, Y, A;
    float Co, SLt, Psi, Il, K, P1, P2, bit, bk;

    days = day2000(y, m, d, h);
    t = days / 36525;


    L1 = range(280.466 + 36000.8 * t);
    M1 = range(357.529+35999*t - 0.0001536* t*t + t*t*t/24490000);

    C1 = (1.915 - 0.004817* t - 0.000014* t * t)* dsin(M1);
    C1 = C1 + (0.01999 - 0.000101 * t)* dsin(2*M1);
    C1 = C1 + 0.00029 * dsin(3*M1);
    V1 = M1 + C1;
    Ec1 = 0.01671 - 0.00004204 * t - 0.0000001236 * t*t;
    R1 = 0.99972 / (1 + Ec1 * dcos(V1));
    Th1 = L1 + C1;
    Om1 = range(125.04 - 1934.1 * t);
    Lam1 = Th1 - 0.00569 - 0.00478 * dsin(Om1);
    Obl = (84381.448 - 46.815 * t)/3600;
    Ra1 = datan2(dsin(Th1) * dcos(Obl) - dtan(0)* dsin(Obl), dcos(Th1));
    Dec1 = dasin(dsin(0)* dcos(Obl) + dcos(0)*dsin(Obl)*dsin(Th1));

    F = range(93.2721 + 483202 * t - 0.003403 * t* t - t * t * t/3526000);
    L2 = range(218.316 + 481268 * t);
    Om2 = range(125.045 - 1934.14 * t + 0.002071 * t * t + t * t * t/450000);
    M2 = range(134.963 + 477199 * t + 0.008997 * t * t + t * t * t/69700);
    D = range(297.85 + 445267 * t - 0.00163 * t * t + t * t * t/545900);
    D2 = 2*D;
    R2 = 1 + (-20954 * dcos(M2) - 3699 * dcos(D2 - M2) - 2956 * dcos(D2)) / 385000;
    R3 = (R2 / R1) / 379.168831168831;
    Bm = 5.128 * dsin(F) + 0.2806 * dsin(M2 + F);
    Bm = Bm + 0.2777 * dsin(M2 - F) + 0.1732 * dsin(D2 - F);
    Lm = 6.289 * dsin(M2) + 1.274 * dsin(D2 -M2) + 0.6583 * dsin(D2);
    Lm = Lm + 0.2136 * dsin(2*M2) - 0.1851 * dsin(M1) - 0.1143 * dsin(2 * F);
    Lm = Lm +0.0588 * dsin(D2 - 2*M2);
    Lm = Lm + 0.0572* dsin(D2 - M1 - M2) + 0.0533* dsin(D2 + M2);
    Lm = Lm + L2;
    Ra2 = datan2(dsin(Lm) * dcos(Obl) - dtan(Bm)* dsin(Obl), dcos(Lm));
    Dec2 = dasin(dsin(Bm)* dcos(Obl) + dcos(Bm)*dsin(Obl)*dsin(Lm));
    HLm = range(Lam1 + 180 + (180/pi) * R3 * dcos(Bm) * dsin(Lam1 - Lm));
    HBm = R3 * Bm;

    float I = 1.54242;
    W = Lm - Om2;
    Y = dcos(W) * dcos(Bm);
    X = dsin(W) * dcos(Bm) * dcos(I) - dsin(Bm) * dsin(I);
    A = datan2(X, Y);
    EL = A - F;
    EB = dasin(-dsin(W) * dcos(Bm) * dsin(I) - dsin(Bm) * dcos(I));

    W = range(HLm - Om2);
    Y = dcos(W) * dcos(HBm);
    X = dsin(W) * dcos(HBm) * dcos(I) - dsin(HBm) * dsin(I);
    A = datan2(X, Y);
    float SL = range(A - F);
    float SB = dasin(-dsin(W) * dcos(HBm) * dsin(I) - dsin(HBm) * dcos(I));

    if (SL < 90) {
        Co = 90 - SL;
    }
    else	{
        Co = 450 - SL;
    }

    if ((Co > 90) && (Co < 270)) {
        SLt = 180 - Co;
    }
    else	{
        if	(Co < 90) {
            SLt = 0 - Co;
        }
        else	{
            SLt = 360 - Co;
        }
    }

    A = dcos(Bm) * dcos(Lm - Lam1);
    Psi = 90 - datan(A / sqrt(1- A*A));
    X = R1 * dsin(Psi);
    Y = R3 - R1* A;
    Il = datan2(X, Y);
    K = (1 + dcos(Il))/2;

    X = dsin(Dec1) * dcos(Dec2) - dcos(Dec1) * dsin(Dec2) * dcos(Ra1 - Ra2);
    Y = dcos(Dec1) * dsin(Ra1 - Ra2);
    P1 = datan2(Y, X);

    X = dsin(I) * dsin(Om2);
    Y = dsin(I) * dcos(Om2) * dcos(Obl) - dcos(I) * dsin(Obl);
    W = datan2(X, Y);
    A = sqrt(X*X + Y*Y) * dcos(Ra2 - W);
    P2 = dasin(A / dcos(EB));

    float daynumber = round(days, 4);
    float julday = round(days + 2451545.0, 4);
    float SunDistance = round(R1, 4);
    float SunRa = round(Ra1 / 15, 3);
    float SunDec = round(Dec1, 2);

   //
   //	Write Moon numbers to form.
   //

    float MoonDist = round(R2 * 60.268511, 2);
    float MoonRa = round(Ra2 / 15, 3);
    float MoonDec = round(Dec2, 2);

   //
   //	Print the libration numbers.
   //

    float SelLatEarth = round(EB , 1);
    float SelLongEarth = round(EL, 1);

   //
   //	Print the Sub-solar numbers.
   //

    float SelLatSun = round(SB , 1);
    float SelLongSun = round(SL, 1);
    float SelColongSun = round(Co, 2);
    float SelLongTerm = round(SLt, 1);

   //
   //	Print the rest - position angles and illuminated fraction.
   //

    float SelIlum = round(K, 3);
    float SelPaBl = round(P1, 1);
    float SelPaPole = round(P2, 1);

    cout<< MoonRa <<stop;
    cout<< MoonDec <<stop;
}
void position_calculator(float lat, float lon, float ra, float dec) {

    if (month <= 2)
    {
        month += 12;
        year -=1;
    }
    v1 = (year/400) - (year/100) + (year/4);
    
    v2 = (365 * year) - 679004;

    md = v1 + v2 + 306001 * (month+1)/10000 + day;
    t0 = (md - 51544.5)/36525;

    a1 = 24110.54841;
    a2 = 8640184.812;
    a3 = 0.093104;
    a4 = 0.0000062;

    s0 = a1 + (a2 * t0) + (a3 * pow(t0,2)) - (a4 * pow(t0,3));

    float ut = hrs + minutes/60 + sec/3600;
    float Nsec = ut * 3600;

    float NsecS = Nsec * 366.2422 / 365.2422;

    float sg = (s0 + NsecS)/3600 * 15;
    float st = sg + lon; // Местное звёздное время
    float st1 = st - (int(st)/360)*360; // mestnoe zvezdnoe vremya

    AZ[i] = (st1 - ra);

    z = acos((sin(lat*dtor) * sin(dec*dtor)) + (cos(lat*dtor)
                                                * cos(dec*dtor)  *  cos(AZ[i]*dtor)));

    ALT[i] = 90 - (z/dtor);

    if (ALT[i]>0)
    {
        j+=1;
        listNums[j] = i;
        cout<< j+1 << ", " << starNames[i] << ", " << AZ[i] << ", " << ALT[i] <<stop;

    }
}
void diff_position_calculator() {

}
int after_parser_counter() {
    ifstream namesFile("C:\\Users\\User\\CLionProjects\\BigBoyTelescop\\StarNames.txt",
                       ios::in);                                 ///NAMES TEXT FILE START

    if (!namesFile) {
        std::cout << "File 1 is dead, long live the File 1 !" << '\n';
        return 0;
    }

    for (int i = 0; !namesFile.eof(); i++) {
        namesFile.getline(lastStrLine, 1024, stop); //считывание строк в распарсеном и записаном файле
        strLineLen++;
        names.emplace_back(lastStrLine);
    }
    delete[] lastStrLine;
    namesFile.close();                                                            ///NAMES TEXT FILE END

    ifstream raFile("C:\\Users\\User\\CLionProjects\\BigBoyTelescop\\RA.txt",
                    ios::in);                                           ///RA TEXT FILE START

    if (!raFile) {
        std::cout << "File 2 is dead, long live the File 2 !" << '\n';
        return 0;
    }

    for (int i = 0; !raFile.eof(); i++) {
        raFile.getline(lastStrLine, 1024, stop); //считывание строк в распарсеном и записаном файле
        strLineLen++;
        Ra.emplace_back(lastStrLine);
    }
    delete[] lastStrLine;
    raFile.close();                                                               ///RA TEXT FILE END

    ifstream decFile("C:\\Users\\User\\CLionProjects\\BigBoyTelescop\\Dec.txt",
                     ios::in);                                         ///DEC TEXT FILE START

    if (!decFile) {
        std::cout << "File 3 is dead, long live the File 3 !" << '\n';
        return 0;
    }

    for (int i = 0; !decFile.eof(); i++) {
        decFile.getline(lastStrLine, 1024, stop); //считывание строк в распарсеном и записаном файле
        strLineLen++;
        Dec.emplace_back(lastStrLine);
    }
    delete[] lastStrLine;
    decFile.close();                                                              ///DEC TEXT FILE END
 ///////////////////////////////////////////////////
    for (size_t i = 0; i < names.size(); ++i) {
        starNames.emplace_back(const_cast<char *>(names[i].c_str()));
    }
 ///////////////////////////////////////////////////
    for (size_t i = 0; i < Ra.size(); ++i) {
        starRa.emplace_back(const_cast<char *>(Ra[i].c_str()));
    }
 ///////////////////////////////////////////////////
    for (size_t i = 0; i < Dec.size(); ++i) {
        starDec.emplace_back(const_cast<const char *>(Dec[i].c_str()));
    }
 ///////////////////////////////////////////////////

    for (int i = 0; i < starRa.size(); i++) {
        result = sscanf(starRa[i], "%d:%d:%d", &rahh, &ramm, &rass);
        rasumm = rahh + float(ramm) / 60 + float(rass) / 3600;//std::cout << rasumm << '\n';
        RA[i] = rasumm;
    }

    for (int i = 0; i < starDec.size(); i++) {
        result = sscanf(starDec[i], "%d:%d:%d", &decdeg, &decmm, &decss);
        if (decdeg > 0) {
            decsumm = decdeg + float(decmm) / 60 + float(decss) / 3600;
            DEC[i] = decsumm;

        } else {
            decsumm = decdeg - float(decmm) / 60 - float(decss) / 3600;
            DEC[i] = decsumm;
        }
    }
}
int arduino_receiver() {
    string lonlat;
    char *incomingData = new char[256];
    string poop = "foo";
    char *ghg = new char[poop.size() + 1];
    std::copy(poop.begin(), poop.end(), ghg);
    //printf("%s\n",incomingData);
    int dataLength = 256;
    int readResult = 0;
    ghg[poop.size()] = '\n';

    SerialPort arduino(port_num);
    if (arduino.isConnected()) cout << "Connection Established" << endl;
    else cout << "ERROR, check port name";

    after_parser_counter();

    while (arduino.isConnected())
    {

        int read_result = arduino.readSerialPort(incomingData, MAX_DATA_LENGTH);
        //prints out dat
        //read_result = 1;
        if (read_result == 0){continue;}
        //cout<< read_result <<stop;
        int write_result = arduino.writeSerialPort(ghg, MAX_DATA_LENGTH);
        //cout<< incomingData <<stop;
        lonlat.assign(incomingData);

        vector<string> VecStr;
        istringstream ss(lonlat);
        string String;
        while(ss>>String)
            VecStr.push_back(String);

        LAT = stof(VecStr[0]);
        LON = stof(VecStr[1]);
        //LAT = 57.16;
        //LON = 60.20;

        cout<< "Choose a star, nibba: ";
        cin >> choiceA;
        cout<< "You've chosen: " << starNames[listNums[choiceA-1]] <<stop;
        cout<< AZ[listNums[choiceA-1]] <<stop;
        cout<< ALT[listNums[choiceA-1]] <<stop;

        cout<< LAT <<stop;
        cout<< LON <<stop;

        for (i = 0; i < starNames.size(); i++)
        {
            position_calculator(LAT, LON, RA[i], DEC[i]);
        }
        break;
    }
}
int time_date_parse() {
    time_t myTime = time(NULL);

    tm *tm = gmtime(&myTime);
    year = 1900 + tm->tm_year;
    month = 1 + tm->tm_mon;
    day = tm->tm_mday;
    hrs = 1 + tm->tm_hour;
    minutes = 1 + tm->tm_min;
    sec = 1 + tm->tm_sec;

    hours = hrs + float(minutes)/60 + float(sec) / 3600;
}
int arduino_sender() {
    {
        SerialPort arduino(port_num);
        if (arduino.isConnected()) cout << "Connection Established" << endl;
        else cout << "ERROR, check port name";

        while (arduino.isConnected()) {

            char *outputData = new char[256];
            //printf("%s\n",incomingData);
            int dataLength = 256;
            int readResult = 0;
            int writeResult = 0;

            std::ostringstream azal;
            azal << AzDiff;
            azal << " ";
            azal << AltDiff;
            string aZaL(azal.str());
            //cout << aZaL << stop;

            char *c_string = new char[aZaL.size() + 1];

            std::copy(aZaL.begin(), aZaL.end(), c_string);

            c_string[aZaL.size()] = '\n';
            //cout<< c_string <<stop;

            int pop = arduino.writeSerialPort(c_string, MAX_DATA_LENGTH);
            //Getting reply from arduino
            int poop = arduino.readSerialPort(output, MAX_DATA_LENGTH);
            cout<< pop <<stop;
            if(poop == 0){continue;}
            puts(output);
            break;

        }
    }
} ///cool stuff doe

/**********************************************************************************************************************/
/*********************************************************CODE*********************************************************/
/**********************************************************************************************************************/

int main()
{
    time_date_parse();

    //arduino_receiver();

    AzDiff = AZ[choiceB] - CurrAz; cout<< "AzDiff is: " << AzDiff << stop;
    AltDiff = ALT[choiceB] - CurrAlt; cout<< "AltDiff is: " << AltDiff << stop;

    CurrAlt = ALT[choiceB];
    CurrAz = AZ[choiceB];


    //arduino_sender();
}
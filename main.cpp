#include"lepkosc.h"

int main()
{

double laktyd = 6.762;
double kat1 = 0.002701; //0.04%
double kat2 = 0.003375; //0.05%
double kat3 = 0.002057; //0.03%
double kat4 = 0.004049; //0.06%
string nazwa1 = "stezenia1.txt";
string nazwa2 = "stezenia2.txt";
string nazwa3 = "stezenia3.txt";
string nazwa4 = "stezenia4.txt";
int n = 4000;
double kR = 300;
double kI = 10;
int czas = 200;
double krok = 0.001;


Mieszanina M1(n, laktyd, kat1, kR, kI, 0, nazwa1);
Mieszanina M2(n, laktyd, kat2, kR, kI, 0, nazwa2);
Mieszanina M3(n, laktyd, kat3, kR, kI, 0, nazwa3);
Mieszanina M4(n, laktyd, kat4, kR, kI, 0, nazwa4);
utworz_plik(M1, n);
utworz_plik(M2, n);
utworz_plik(M3, n);
utworz_plik(M4, n);

M1.symulacja(czas, krok);
M1.~Mieszanina();

M2.symulacja(czas, krok);
M2.~Mieszanina();

M3.symulacja(czas, krok);
M3.~Mieszanina();

M4.symulacja(czas, krok);
M4.~Mieszanina();

return 0;
}

#ifndef _lepkosc_h_
#define _lepkosc_h_

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <string>

using namespace std;

class Mieszanina{

	public:
	
	double * stezenia;
	double * stezenia_poczatkowe;
	double * ulamki_molowe;
	double * Mi;
	double * ulamki_wagowe;
	int n;
	double kR;
	float * K;
	double kI;
	double suma_polimerow;
	double Mn;
	double Mw;
	double Y;
	double lepkosc180, lepkosc190, lepkosc200;
	string plik;

	//destruktory, konstruktory...
	Mieszanina(int nn = 0, double laktyd = 0, double kat = 0, 
		   double kRR = 0, double kII = 0, double suma = 0, string nazwa = " ");
	~Mieszanina();	

	
	void krok_symulacji(double k);
	void krok_wstepny(double k);
	void symulacja(int t, double k);
	void aktualizuj_ulamki_mol();
	void aktualizuj_ulamki_wag();
	void licz_Mn();
	void licz_Mw();
	void zapisz_stezenia(double i);
	void licz_Y();
	void licz_lepkosc180(double lep);
	void licz_lepkosc190(double lep);
	void licz_lepkosc200(double lep);
	double licz_lep_polim180();
	double licz_lep_polim190();
	double licz_lep_polim200();

};


//void krok_symulacji(Mieszanina &M);
//void symulacja(Mieszanina &M);
void utworz_plik(Mieszanina &M, int n);

#endif

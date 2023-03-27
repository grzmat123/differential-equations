#include "lepkosc.h"

Mieszanina::Mieszanina(int nn, double laktyd, double kat, double kRR, double kII, double suma, string nazwa){

	n = nn;
	kR = kRR;
	kI = kII;
	suma_polimerow = suma;
	plik = nazwa;
	Mn = 0;
	Mw = 0;
	Y = 0;
	lepkosc180 = 0;
	lepkosc190 = 0;
	lepkosc200 = 0;

	stezenia = new double [n+2];
	stezenia_poczatkowe = new double[n+2];
	ulamki_molowe = new double [n];
	Mi = new double [n];
	ulamki_wagowe = new double [n];
	K = new float [n-1];
	stezenia[0] = laktyd;
	stezenia[n+1] = kat;
	stezenia_poczatkowe[0] = laktyd;
	stezenia_poczatkowe[n+1] = kat;

	for (int i = 1; i < n+1; i++)
	{
	stezenia[i] = 0;
	stezenia_poczatkowe[i] = 0;
	ulamki_molowe[i-1] = 0;
	ulamki_wagowe[i-1] = 0;
	Mi[i-1] = i*144 + 405;
	}

	for (int i = 0; i < n-1; i++)
	{
	//K[i] = (1 - i/n)*(1 - i/n)*(1 - i/n)*(1 - i/n)*(1 - i/n)*(1 - i/n)*(1 - i/n)*(1 - i/n)*kR;
	K[i] = exp(-70*i/n)*(1 - i/n)*kR;
	}
}

Mieszanina::~Mieszanina(){

	delete stezenia, stezenia_poczatkowe, ulamki_molowe, Mi, ulamki_wagowe, K;
}

void Mieszanina::zapisz_stezenia(double i){

	fstream output;
	const char * pomoc = plik.c_str();
     	output.open (pomoc); 
	output.seekg(0, ios::end); 
     	if( output.good() == true )
    	{
    	//cout << "Zapisuje wartosci stężeń do pliku stezenia!" << endl;
    	output<<i<<" "<<stezenia[0]<<" "<<stezenia[n+1]<<" "<<Mn<<" "
	<<Mw<<" "<<lepkosc180<<" "<<lepkosc190<<" "<<lepkosc200<<endl;
    	}
	else {cout<<"Brak dostepu do wyników!"<<endl;}
	output.close();
}

void utworz_plik(Mieszanina &M, int n){

	fstream output;
	const char * pomoc = M.plik.c_str();
	output.open (pomoc, ios::out);  
     	if( output.good() == true )
    	{
    	cout << "Utworzyłem plik do zapisu wartości stężeń!" << endl;
	output<<"czas "<<"laktyd "<<"kat "<<"Mn "<<"Mw "<<"Lepkosc180 "
	<<"Lepkosc190 "<<"Lepkosc200 "<<endl;    	
	output<<0<<" "<<M.stezenia[0]<<" "<<M.stezenia[n+1]<<" "<<M.Mn<<" "
	<<M.Mw<<" "<<M.lepkosc180<<" "<<M.lepkosc190<<" "<<M.lepkosc200<<endl;
    	}
	else {cout<<"Błąd utworzenia pliku!"<<endl;}
	output.close();	
}


void Mieszanina::krok_symulacji(double k){

		
	//tablice współczynników
	double tablica[n+2][4];
	double temp[n+2];	//pomocnicze tabela do stężeń


	//liczy wsp. zerowe
	tablica[1][0] = (kI*stezenia[n+1]-K[0]*stezenia[1])*k*stezenia[0];
	tablica[n+1][0] = -k*kI*stezenia[0]*stezenia[n+1];
	tablica[0][0] = (-kI*stezenia[n+1]-suma_polimerow)*k*stezenia[0];
	tablica[n][0] =	k*K[n-2]*stezenia[n-1]*stezenia[0];
	for (int i = 2; i < n; i++)
	{
	tablica[i][0] = k*(K[i-2]*stezenia[i-1]-K[i-1]*stezenia[i])*stezenia[0];
	}

	//aktualizacja stezen po 0. rzedzie
	for (int i = 0; i < n+2; i ++)
	{
	temp[i] = stezenia[i] + 0.5*tablica[i][0];
	}
	suma_polimerow = 0;
	for (int i = 1; i < n; i ++)
	{
	suma_polimerow = suma_polimerow + K[i-1]*temp[i];
	}	

	//liczy wsp. rzedu 1
	tablica[1][1] = (kI*temp[n+1]-K[0]*temp[1])*k*temp[0];
	tablica[n+1][1] = -k*kI*temp[0]*temp[n+1];
	tablica[0][1] = (-kI*temp[n+1]-suma_polimerow)*k*temp[0];
	tablica[n][1] =	k*K[n-2]*temp[n-1]*temp[0];
	for (int i = 2; i < n; i++)
	{
	tablica[i][1] = k*(K[i-2]*temp[i-1]-K[i-1]*temp[i])*temp[0];
	}

	//aktualizacja stezen po 1. rzedzie
	for (int i = 0; i < n+2; i ++)
	{
	temp[i] = stezenia[i] + 0.5*tablica[i][1];
	}
	suma_polimerow = 0;
	for (int i = 1; i < n; i ++)
	{
	suma_polimerow = suma_polimerow + K[i-1]*temp[i];
	}

	//liczy wsp. rzedu 2
	tablica[1][2] = (kI*temp[n+1]-K[0]*temp[1])*k*temp[0];
	tablica[n+1][2] = -k*kI*temp[0]*temp[n+1];
	tablica[0][2] = (-kI*temp[n+1]-suma_polimerow)*k*temp[0];
	tablica[n][2] =	k*K[n-2]*temp[n-1]*temp[0];
	for (int i = 2; i < n; i++)
	{
	tablica[i][2] = k*(K[i-2]*temp[i-1]-K[i-1]*temp[i])*temp[0];
	}

	//aktualizacja stezen po 2. rzedzie
	for (int i = 0; i < n+2; i ++)
	{
	temp[i] = stezenia[i] + tablica[i][2];
	}
	suma_polimerow = 0;
	for (int i = 1; i < n; i ++)
	{
	suma_polimerow = suma_polimerow + K[i-1]*temp[i];
	}
	
	//liczy wsp. rzedu 3
	tablica[1][3] = (kI*temp[n+1]-K[0]*temp[1])*k*temp[0];
	tablica[n+1][3] = -k*kI*temp[0]*temp[n+1];
	tablica[0][3] = (-kI*temp[n+1]-suma_polimerow)*k*temp[0];
	tablica[n][3] =	k*K[n-2]*temp[n-1]*temp[0];
	for (int i = 2; i < n; i++)
	{
	tablica[i][3] = k*(K[i-2]*temp[i-1]-K[i-1]*temp[i])*temp[0];
	}
	
	//liczy stezenia
	for (int i = 0; i < n+2; i++)
	{
	stezenia[i] = stezenia[i] + (tablica[i][0]+2*tablica[i][1]+2*tablica[i][2]+tablica[i][3])*1./6.;
	}
	suma_polimerow = 0;
	for (int i = 1; i < n; i ++)
	{
	suma_polimerow = suma_polimerow + K[i-1]*stezenia[i];
	}
}

void Mieszanina::krok_wstepny(double k){

	cout<<"Rozpoczynam krok wstępny."<<endl;
	//tablice współczynników
	double tablica[n+2][4];
	double temp[n+2];	//pomocnicze tabela do stężeń


	//liczy wsp. zerowe
	tablica[1][0] = kI*stezenia_poczatkowe[n+1]*k*stezenia_poczatkowe[0];
	tablica[n+1][0] = -k*kI*stezenia_poczatkowe[0]*stezenia_poczatkowe[n+1];
	tablica[0][0] = -kI*stezenia_poczatkowe[n+1]*k*stezenia_poczatkowe[0];
	tablica[n][0] =	0;
	for (int i = 2; i < n; i++)
	{
	tablica[i][0] = 0;
	}

	//aktualizacja stezen po 0. rzedzie
	for (int i = 2; i < n+1; i ++)
	{
	temp[i] = 0;
	}
	temp[1] = 0.5*tablica[1][0];
	temp[0] = stezenia_poczatkowe[0] + 0.5*tablica[0][0];
	temp[n+1] = stezenia_poczatkowe[n+1] + 0.5*tablica[n+1][0];
	suma_polimerow = K[0]*temp[1];	

	//liczy wsp. rzedu 1
	tablica[1][1] = (kI*temp[n+1]-K[0]*temp[1])*k*temp[0];
	tablica[n+1][1] = -k*kI*temp[0]*temp[n+1];
	tablica[0][1] = (-kI*temp[n+1]-suma_polimerow)*k*temp[0];
	tablica[n][1] =	0;
	tablica[2][1] = k*K[0]*temp[1]*temp[0];
	for (int i = 3; i < n; i++)
	{
	tablica[i][1] = 0;
	}

	//aktualizacja stezen po 1. rzedzie
	for (int i = 3; i < n+1; i ++)
	{
	temp[i] = 0;
	}
	temp[1] = 0.5*tablica[1][1];
	temp[2] = 0.5*tablica[2][1];
	temp[0] = stezenia_poczatkowe[0] + 0.5*tablica[0][1];
	temp[n+1] = stezenia_poczatkowe[n+1] + 0.5*tablica[n+1][1];
	suma_polimerow = K[0]*temp[1] + K[1]*temp[2];
	


	//liczy wsp. rzedu 2
	tablica[1][2] = (kI*temp[n+1]-K[0]*temp[1])*k*temp[0];
	tablica[n+1][2] = -k*kI*temp[0]*temp[n+1];
	tablica[0][2] = (-kI*temp[n+1]-suma_polimerow)*k*temp[0];
	tablica[n][2] =	0;
	tablica[2][2] = k*(K[0]*temp[1]-K[1]*temp[2])*temp[0];
	tablica[3][2] = k*K[1]*temp[2]*temp[0];
	for (int i = 4; i < n; i++)
	{
	tablica[i][2] = 0;
	}

	//aktualizacja stezen po 2. rzedzie
	for (int i = 4; i < n+1; i ++)
	{
	temp[i] = 0;
	}
	temp[1] = tablica[1][2];
	temp[2] = tablica[2][2];
	temp[3] = tablica[3][2];
	temp[0] = stezenia_poczatkowe[0] + tablica[0][2];
	temp[n+1] = stezenia_poczatkowe[n+1] + tablica[n+1][2];
	suma_polimerow = K[0]*temp[1] + K[1]*temp[2] + K[2]*temp[3];
	
	
	//liczy wsp. rzedu 3
	tablica[1][3] = (kI*temp[n+1]-K[0]*temp[1])*k*temp[0];
	tablica[n+1][3] = -k*kI*temp[0]*temp[n+1];
	tablica[0][3] = (-kI*temp[n+1]-suma_polimerow)*k*temp[0];
	tablica[n][3] =	0;
	tablica[2][3] = k*(K[0]*temp[1]-K[1]*temp[2])*temp[0];
	tablica[3][3] = k*(K[1]*temp[2]-K[2]*temp[3])*temp[0];
	tablica[4][3] = k*K[2]*temp[3]*temp[0];
	for (int i = 5; i < n; i++)
	{
	tablica[i][3] = 0;
	}
	
	//liczy stezenia
	for (int i = 0; i < 5; i++)
	{
	stezenia[i] = stezenia_poczatkowe[i] + (tablica[i][0]+2*tablica[i][1]+2*tablica[i][2]+tablica[i][3])*1./6.;
	}
	for (int i = 5; i < n+1; i++)
	{
	stezenia[i] = 0;
	}
	stezenia[n+1] = stezenia_poczatkowe[n+1] + (tablica[n+1][0]+2*tablica[n+1][1]+2*tablica[n+1][2]+tablica[n+1][3])*1./6.;  
	suma_polimerow = K[0]*stezenia[1] + K[1]*stezenia[2] + K[2]*stezenia[3] + K[3]*stezenia[4];
	cout<<"Zakończyłem krok wstępny"<<endl;
}

void Mieszanina::aktualizuj_ulamki_mol(){
	
	double temp = 0;
	for (int i = 1; i < n+1; i++)
	{
	temp = temp + stezenia[i];
	}
 
	for (int i = 1; i < n+1; i++)
	{
	ulamki_molowe[i-1] = stezenia[i]/temp;
	}
}

void Mieszanina::aktualizuj_ulamki_wag(){

	double temp = 0;
	for (int i = 1; i < n+1; i++)
	{
	temp = temp + stezenia[i]*Mi[i-1];
	}
 
	for (int i = 1; i < n+1; i++)
	{
	ulamki_wagowe[i-1] = stezenia[i]*Mi[i-1]/temp;
	}
}

void Mieszanina::licz_Mn(){
	
	double temp = 0;
	for (int i = 0; i < n; i++){
	temp = temp + ulamki_molowe[i]*Mi[i];
	}
	Mn = temp;
}

void Mieszanina::licz_Mw(){
	
	double temp = 0;
	for (int i = 0; i < n; i++){
	temp = temp + ulamki_wagowe[i]*Mi[i];
	}
	Mw = temp;
}

void Mieszanina::licz_Y(){

	Y = 1 - stezenia[0]/stezenia_poczatkowe[0];
}

double Mieszanina::licz_lep_polim180(){
	
	double lep = Mn*231.6-1.099*10000000;
	return lep;
}

double Mieszanina::licz_lep_polim190(){
	
	double lep = Mn*117.7-5.539*1000000;
	return lep;
}


double Mieszanina::licz_lep_polim200(){
	
	double lep = Mn*62.89-2.942*1000000;
	return lep;
}

void Mieszanina::licz_lepkosc180(double lep){
	
	double pom = 7.2511*Y ;
	if (lep < 0) lep = -lep;
	lepkosc180 = 716.84*exp(pom)*lep/1400448;
}

void Mieszanina::licz_lepkosc190(double lep){

	double pom = 6.74*Y ;
	if (lep < 0) lep = -lep;
	lepkosc190 = 911.67*exp(pom)*lep/952000;
}

void Mieszanina::licz_lepkosc200(double lep){

	double pom = 7.8794*Y ;
	if (lep < 0) lep = -lep;
	lepkosc200 = 228.69*exp(pom)*lep/422971;
}

void Mieszanina::symulacja(int t, double k){
	
	cout<<"Rozpoczynam symulacje"<<endl;
	int liczba_krokow = t/k;
	cout<<"Liczba krokow: "<<liczba_krokow<<endl;
	krok_wstepny(k);
	aktualizuj_ulamki_mol();
	aktualizuj_ulamki_wag();
	licz_Mn();
	double lep180 = licz_lep_polim180();
	double lep190 = licz_lep_polim190();
	double lep200 = licz_lep_polim200();
	licz_Mw();
	licz_Y();
	licz_lepkosc180(lep180);
	licz_lepkosc190(lep190);
	licz_lepkosc200(lep200);
	zapisz_stezenia(k);

	for (int i = 1; i < liczba_krokow ; i++)
	{
	double czas = (i+1)*k;
	krok_symulacji(k);
	aktualizuj_ulamki_mol();
	aktualizuj_ulamki_wag();
	licz_Mn();
	lep180 = licz_lep_polim180();
	lep190 = licz_lep_polim190();
	lep200 = licz_lep_polim200();
	licz_Mw();
	licz_Y();
	licz_lepkosc180(lep180);
	licz_lepkosc190(lep190);
	licz_lepkosc200(lep200);
	zapisz_stezenia(czas);
	//cout<<"\a";
	}
	cout<<"Zakonczylem symulacje"<<endl;
}


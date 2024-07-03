#include "TFile.h"
#include "TH2F.h"
#include "TRandom.h"
#include "TTree.h"
#include "TSpectrum.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TCutG.h"
#include "TProfile.h"
#include <iostream>
#include <fstream>
#include <istream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "TH3.h"
#include "TH2.h"
#include "TH1.h"
#include <TROOT.h>
#include <stdint.h>
#include "TMath.h"
#include <algorithm>
#include <string>
#include <vector>
#include <cstdlib>
#include "Kratta.h"
#include "Plastic.h"
#include "Paris.h"
#include "LaBr.h"
#include "Scaler.h"
#include "Silicon.h"


void rewrite_tree(TTree * tree){
//-------------------------------------------------------------------
// otwieranie pliku istniejacego, podwiazywanie drzewa i galezi
//-------------------------------------------------------------------


	struct {
     		int nrun;
     		int id;
     		int trig;
     		int mult;
     		} nfo;
         
    Plastic * plastic = new Plastic;
    Kratta *  kratta = new Kratta;
    Paris *  paris = new Paris;
    LaBr *  labr = new LaBr;
    Silicon * silicon = new Silicon;
	
    tree->SetBranchAddress("paris", &paris);
    tree->SetBranchAddress("labr", &labr);
    tree->SetBranchAddress("plastic", &plastic);
    tree->SetBranchAddress("kratta", &kratta);
    tree->SetBranchAddress("silicon", &silicon);
    tree->SetBranchAddress("info", &nfo);

//-------------------------------------------------------------------
// tworzenie nowego pliku, jego drzewa i galezi
//-------------------------------------------------------------------

TFile *f1 = new TFile("13C_data.root","recreate");
TTree *T1 = new TTree("tree2","tree2");

	struct 
	{
		int id;
    		float time;
    		float energyqshort;
    		float energyqlong;
    		float energyqshort_cal;
    	} 
	LaBr;


	struct 
	{
		int id;
    		float time;
    	} 
	Plastic;


	struct
	{
		int id;
		float time;
    		float energyqlong;
		float energyqshort;
    	} 
	Paris;


	struct 
	{
		int id;
    		float energyPD0;
    		float energyPD1;
    		float energyPD2;
    		float timePD0;
    		float timePD1;
    		float timePD2;
    		float energyPD2_cal;
    	} 
	Kratta;


	struct 
	{
		int id;
    		float energy;
		float energy_cal;
    		float time30;
		float time80;
		float difftime;
		float difftime_cal;
    	} 
	Silicon;
	
	
	T1->Branch("LaBr",&LaBr,"id/I:time/F:energyqshort/F:energyqlong/F:energyqshort_cal/F");
	T1->Branch("Plastic",&Plastic,"id/I:time/F");
	T1->Branch("Paris",&Paris,"id/I:time/F:energyqlong/F:energyqshort/F");
	T1->Branch("Kratta",&Kratta,"id/I:energyPD0/F:energyPD1/F:energyPD2/F:timePD0/F:timePD1/F:timePD2/F:energyPD2_cal/F");
	T1->Branch("Silicon",&Silicon,"id/I:energy/F:energy_cal/F:time30/F:time80/F:difftime/F:difftime_cal/F");
	T1->Branch("Info",&nfo,"run/I:id/I:trig/I:mult/I");
//T1->Branch("dowolna_nazwa",&wskaznik_do_tego_co_ma_byc_galezia,"lisc1/I:lisc2/D:lisc3/D");

//-------------------------------------------------------------------
// wczytanie współczynników kalibarcyjnych do energii Silicon
//-------------------------------------------------------------------

std::vector<double>sil_cal0, sil_cal1;

cout << "cal_silicon" << endl;

FILE *plik;
plik = fopen("calibrations/cal_silicon.txt", "r");
if (plik == NULL) {
        printf("Nie można otworzyć pliku.\n");
    }

double sil1[32], sil2[32];
for (int i = 0; i < 32; i++) {
	if (fscanf(plik, "%lf %lf", &sil1[i], &sil2[i]) != 2) {
		printf("Błąd wczytywania danych z wiersza %d.\n", i + 1);
		sil1[i] = 0;
		sil2[i] = 1;
        }
        if (std::isnan(sil1[i]) || std::isnan(sil2[i])) {
		printf("Błąd: Wartość nan w wierszu %d.\n", i + 1);
		sil1[i] = 0;
		sil2[i] = 1;
        }
        
        sil_cal0.push_back(sil1[i]);
	sil_cal1.push_back(sil2[i]);
    }

for (int i = 0; i < sil_cal0.size(); i++)
	cout << sil_cal0[i] << "\t" << sil_cal1[i] << endl;
cout << " " << endl;
fclose(plik);

//-------------------------------------------------------------------
// wczytanie współczynników kalibarcyjnych do energii Kratta
//-------------------------------------------------------------------

std::vector<double>kratta_cal0_I, kratta_cal1_I, kratta_cal2_I;
std::vector<double>kratta_cal0_II, kratta_cal1_II, kratta_cal2_II;
double kra11[24], kra12[24], kra13[24];
double kra21[24], kra22[24], kra23[24];

cout << "cal_kratta_I" << endl;

FILE *plik2;
plik2 = fopen("calibrations/cal_kratta_I.txt", "r");
if (plik2 == NULL) {
        printf("Nie można otworzyć pliku.\n");
    }
    
for (int i = 0; i < 24; i++) {
	if (fscanf(plik2, "%lf %lf %lf", &kra11[i], &kra12[i], &kra13[i]) != 3) {
		printf("Błąd wczytywania danych z wiersza %d.\n", i + 1);
		kra11[i] = 0;
		kra12[i] = 1;
		kra13[i] = 0;
        }
        if (std::isnan(kra11[i]) || std::isnan(kra12[i]) || std::isnan(kra13[i])) {
		printf("Błąd: Wartość nan w wierszu %d.\n", i + 1);
		kra11[i] = 0;
		kra12[i] = 1;
		kra13[i] = 0;
        }
        
        kratta_cal0_I.push_back(kra11[i]);
	kratta_cal1_I.push_back(kra12[i]);
	kratta_cal2_I.push_back(kra13[i]);
    }
    	
for (int i = 0; i < kratta_cal0_I.size(); i++)
	cout << kratta_cal0_I[i] << "\t" << kratta_cal1_I[i] << "\t"<< kratta_cal2_I[i] << endl;
cout << " " << endl;	
cout << "cal_kratta_II" << endl;

FILE *plik3;
plik3 = fopen("calibrations/cal_kratta_II.txt", "r");
if (plik3 == NULL) {
        printf("Nie można otworzyć pliku.\n");
    }
	
for (int i = 0; i < 24; i++) {
	if (fscanf(plik3, "%lf %lf %lf", &kra21[i], &kra22[i], &kra23[i]) != 3) {
		printf("Błąd wczytywania danych z wiersza %d.\n", i + 1);
		kra21[i] = 0;
		kra22[i] = 1;
		kra23[i] = 0;
        }
        if (std::isnan(kra21[i]) || std::isnan(kra22[i]) || std::isnan(kra23[i])) {
		printf("Błąd: Wartość nan w wierszu %d.\n", i + 1);
		kra21[i] = 0;
		kra22[i] = 1;
		kra23[i] = 0;
        }
        
        kratta_cal0_II.push_back(kra21[i]);
	kratta_cal1_II.push_back(kra22[i]);
	kratta_cal2_II.push_back(kra23[i]);
    }

for (int i = 0; i < kratta_cal0_II.size(); i++)
	cout << kratta_cal0_II[i] << "\t" << kratta_cal1_II[i] << "\t"<< kratta_cal2_II[i] << endl;
cout << " " << endl;
fclose(plik2);
fclose(plik3);

//-------------------------------------------------------------------
// wczytanie współczynników kalibarcyjnych do energii LaBr
//-------------------------------------------------------------------

std::vector<double>labr_cal0_I, labr_cal1_I, labr_cal2_I;
std::vector<double>labr_cal0_II_1, labr_cal1_II_1, labr_cal2_II_1;
std::vector<double>labr_cal0_II_2, labr_cal1_II_2, labr_cal2_II_2;
std::vector<double>labr_cal0_II_3, labr_cal1_II_3, labr_cal2_II_3;
std::vector<double>labr_cal0_III, labr_cal1_III, labr_cal2_III;
std::vector<double>labr_cal0_IV, labr_cal1_IV, labr_cal2_IV;
double labr11[4], labr12[4], labr13[4];
double labr21[4], labr22[4], labr23[4];
double labr31[4], labr32[4], labr33[4];
double labr41[4], labr42[4], labr43[4];
double labr51[4], labr52[4], labr53[4];
double labr61[4], labr62[4], labr63[4];

cout << "cal_labr_I" << endl;

FILE *plik4;
plik4 = fopen("calibrations/cal_labr_I.txt", "r");
if (plik4 == NULL) {
        printf("Nie można otworzyć pliku.\n");
    }

for (int i = 0; i < 4; i++) {
	if (fscanf(plik3, "%lf %lf %lf", &labr11[i], &labr12[i], &labr13[i]) != 3) {
		printf("Błąd wczytywania danych z wiersza %d.\n", i + 1);
		labr11[i] = 0;
		labr12[i] = 1;
		labr13[i] = 0;
        }
        if (std::isnan(labr11[i]) || std::isnan(labr12[i]) || std::isnan(labr13[i])) {
		printf("Błąd: Wartość nan w wierszu %d.\n", i + 1);
		labr11[i] = 0;
		labr12[i] = 1;
		labr13[i] = 0;
        }
        
        labr_cal0_I.push_back(labr11[i]);
	labr_cal1_I.push_back(labr12[i]);
	labr_cal2_I.push_back(labr13[i]);
    }
	
for (int i = 0; i < labr_cal0_I.size(); i++)
	cout << labr_cal0_I[i] << "\t" << labr_cal1_I[i] << "\t" << labr_cal2_I[i] << endl;
cout << " " << endl;
cout << "cal_labr_II_1" << endl;

FILE *plik5;
plik5 = fopen("calibrations/cal_labr_II_1.txt", "r");
if (plik5 == NULL) {
        printf("Nie można otworzyć pliku.\n");
    }
	
for (int i = 0; i < 4; i++) {
	if (fscanf(plik5, "%lf %lf %lf", &labr21[i], &labr22[i], &labr23[i]) != 3) {
		printf("Błąd wczytywania danych z wiersza %d.\n", i + 1);
		labr21[i] = 0;
		labr22[i] = 1;
		labr23[i] = 0;
        }
        if (std::isnan(labr21[i]) || std::isnan(labr22[i]) || std::isnan(labr23[i])) {
		printf("Błąd: Wartość nan w wierszu %d.\n", i + 1);
		labr21[i] = 0;
		labr22[i] = 1;
		labr23[i] = 0;
        }
        
        labr_cal0_II_1.push_back(labr21[i]);
	labr_cal1_II_1.push_back(labr22[i]);
	labr_cal2_II_1.push_back(labr23[i]);
    }
    
for (int i = 0; i < labr_cal0_II_1.size(); i++)
	cout << labr_cal0_II_1[i] << "\t" << labr_cal1_II_1[i] << "\t" << labr_cal2_II_1[i] << endl;
	
cout <<"\n"<< endl;
cout << "cal_labr_II_2" << endl;

FILE *plik6;
plik6 = fopen("calibrations/cal_labr_II_2.txt", "r");
if (plik6 == NULL) {
        printf("Nie można otworzyć pliku.\n");
    }
	
for (int i = 0; i < 4; i++) {
	if (fscanf(plik6, "%lf %lf %lf", &labr31[i], &labr32[i], &labr33[i]) != 3) {
		printf("Błąd wczytywania danych z wiersza %d.\n", i + 1);
		labr31[i] = 0;
		labr32[i] = 1;
		labr33[i] = 0;
        }
        if (std::isnan(labr31[i]) || std::isnan(labr32[i]) || std::isnan(labr33[i])) {
		printf("Błąd: Wartość nan w wierszu %d.\n", i + 1);
		labr31[i] = 0;
		labr32[i] = 1;
		labr33[i] = 0;
        }
        
        labr_cal0_II_2.push_back(labr31[i]);
	labr_cal1_II_2.push_back(labr32[i]);
	labr_cal2_II_2.push_back(labr33[i]);
    }

for (int i = 0; i < labr_cal0_II_2.size(); i++)
	cout << labr_cal0_II_2[i] << "\t" << labr_cal1_II_2[i] << "\t" << labr_cal2_II_2[i]  << endl;
cout << " " << endl;
cout << "cal_labr_II_3" << endl;
FILE *plik7;
plik7 = fopen("calibrations/cal_labr_II_3.txt", "r");
if (plik7 == NULL) {
        printf("Nie można otworzyć pliku.\n");
    }

for (int i = 0; i < 4; i++) {
	if (fscanf(plik7, "%lf %lf %lf", &labr41[i], &labr42[i], &labr43[i]) != 3) {
		printf("Błąd wczytywania danych z wiersza %d.\n", i + 1);
		labr41[i] = 0;
		labr42[i] = 1;
		labr43[i] = 0;
        }
        if (std::isnan(labr41[i]) || std::isnan(labr42[i]) || std::isnan(labr43[i])) {
		printf("Błąd: Wartość nan w wierszu %d.\n", i + 1);
		labr41[i] = 0;
		labr42[i] = 1;
		labr43[i] = 0;
        }
        
        labr_cal0_II_3.push_back(labr41[i]);
	labr_cal1_II_3.push_back(labr42[i]);
	labr_cal2_II_3.push_back(labr43[i]);
    }

for (int i = 0; i < labr_cal0_II_3.size(); i++)
	cout << labr_cal0_II_3[i] << "\t" << labr_cal1_II_3[i] << "\t" << labr_cal2_II_3[i] << endl;

cout << " " << endl;
cout << "cal_labr_III" << endl;
FILE *plik8;
plik8 = fopen("calibrations/cal_labr_III.txt", "r");
if (plik8 == NULL) {
        printf("Nie można otworzyć pliku.\n");
    }

for (int i = 0; i < 4; i++) {
	if (fscanf(plik8, "%lf %lf %lf", &labr51[i], &labr52[i], &labr53[i]) != 3) {
		printf("Błąd wczytywania danych z wiersza %d.\n", i + 1);
		labr51[i] = 0;
		labr52[i] = 1;
		labr53[i] = 0;
        }
        if (std::isnan(labr51[i]) || std::isnan(labr52[i]) || std::isnan(labr53[i])) {
		printf("Błąd: Wartość nan w wierszu %d.\n", i + 1);
		labr51[i] = 0;
		labr52[i] = 1;
		labr53[i] = 0;
        }
        
        labr_cal0_III.push_back(labr51[i]);
	labr_cal1_III.push_back(labr52[i]);
	labr_cal2_III.push_back(labr53[i]);
    }

for (int i = 0; i < labr_cal0_III.size(); i++)
	cout << labr_cal0_III[i] << "\t" << labr_cal1_III[i] << "\t" << labr_cal2_III[i] << endl;

cout << " " << endl;
cout << "cal_labr_IV" << endl;
FILE *plik9;
plik9 = fopen("calibrations/cal_labr_IV.txt", "r");
if (plik9 == NULL) {
        printf("Nie można otworzyć pliku.\n");
    }

for (int i = 0; i < 4; i++) {
	if (fscanf(plik9, "%lf %lf %lf", &labr61[i], &labr62[i], &labr63[i]) != 3) {
		printf("Błąd wczytywania danych z wiersza %d.\n", i + 1);
		labr61[i] = 0;
		labr62[i] = 1;
		labr63[i] = 0;
        }
        if (std::isnan(labr61[i]) || std::isnan(labr62[i]) || std::isnan(labr63[i])) {
		printf("Błąd: Wartość nan w wierszu %d.\n", i + 1);
		labr61[i] = 0;
		labr62[i] = 1;
		labr63[i] = 0;
        }
        
        labr_cal0_IV.push_back(labr61[i]);
	labr_cal1_IV.push_back(labr62[i]);
	labr_cal2_IV.push_back(labr63[i]);
    }

for (int i = 0; i < labr_cal0_IV.size(); i++)
	cout << labr_cal0_IV[i] << "\t" << labr_cal1_IV[i] << "\t" << labr_cal2_IV[i] << endl;

cout << " " << endl;

fclose(plik4);
fclose(plik5);
fclose(plik6);
fclose(plik7);
fclose(plik8);
fclose(plik9);
//-------------------------------------------------------------------
// wczytanie współczynników kalibarcyjnych do różnicy czasu Silicon
//-------------------------------------------------------------------
std::vector<double>diff_cal0, diff_cal1;

cout << "cal_difftime" << endl;

FILE *plik10;
plik10 = fopen("calibrations/cal_difftime.txt", "r");
if (plik10 == NULL) {
        printf("Nie można otworzyć pliku.\n");
    }

double diff1[32], diff2[32];
for (int i = 0; i < 32; i++) {
	if (fscanf(plik10, "%lf %lf", &diff1[i], &diff2[i]) != 2) {
		printf("Błąd wczytywania danych z wiersza %d.\n", i + 1);
		diff1[i] = 0;
		diff2[i] = 1;
        }
        if (std::isnan(diff1[i]) || std::isnan(diff2[i])) {
		printf("Błąd: Wartość nan w wierszu %d.\n", i + 1);
		diff1[i] = 0;
		diff2[i] = 1;
        }
        
        diff_cal0.push_back(diff1[i]);
	diff_cal1.push_back(diff2[i]);
    }

for (int i = 0; i < diff_cal0.size(); i++)
	cout << diff_cal0[i] << "\t" << diff_cal1[i] << endl;
	
fclose(plik10);
cout << " " << endl;
std::vector<double>diff_cal0_hist1, diff_cal1_hist1;

cout << "cal_difftime_hist1" << endl;

FILE *plik11;
plik11 = fopen("calibrations/cal_difftime_hist1.txt", "r");
if (plik11 == NULL) {
        printf("Nie można otworzyć pliku.\n");
    }

double diff1_hist1[4], diff2_hist1[4];
for (int i = 0; i < 4; i++) {
	if (fscanf(plik11, "%lf %lf", &diff1_hist1[i], &diff2_hist1[i]) != 2) {
		printf("Błąd wczytywania danych z wiersza %d.\n", i + 1);
		diff1_hist1[i] = 0;
		diff2_hist1[i] = 1;
        }
        if (std::isnan(diff1_hist1[i]) || std::isnan(diff2_hist1[i])) {
		printf("Błąd: Wartość nan w wierszu %d.\n", i + 1);
		diff1_hist1[i] = 0;
		diff2_hist1[i] = 1;
        }
        
        diff_cal0_hist1.push_back(diff1_hist1[i]);
	diff_cal1_hist1.push_back(diff2_hist1[i]);
    }

for (int i = 0; i < diff_cal0_hist1.size(); i++)
	cout << diff_cal0_hist1[i] << "\t" << diff_cal1_hist1[i] << endl;
	
fclose(plik11);
cout << " " << endl;
std::vector<double>diff_cal0_hist2, diff_cal1_hist2;

cout << "cal_difftime_hist2" << endl;

FILE *plik12;
plik12 = fopen("calibrations/cal_difftime_hist2.txt", "r");
if (plik12 == NULL) {
        printf("Nie można otworzyć pliku.\n");
    }

double diff1_hist2[4], diff2_hist2[4];
for (int i = 0; i < 4; i++) {
	if (fscanf(plik12, "%lf %lf", &diff1_hist2[i], &diff2_hist2[i]) != 2) {
		printf("Błąd wczytywania danych z wiersza %d.\n", i + 1);
		diff1_hist2[i] = 0;
		diff2_hist2[i] = 1;
        }
        if (std::isnan(diff1_hist2[i]) || std::isnan(diff2_hist2[i])) {
		printf("Błąd: Wartość nan w wierszu %d.\n", i + 1);
		diff1_hist2[i] = 0;
		diff2_hist2[i] = 1;
        }
        
        diff_cal0_hist2.push_back(diff1_hist2[i]);
	diff_cal1_hist2.push_back(diff2_hist2[i]);
    }

for (int i = 0; i < diff_cal0_hist2.size(); i++)
	cout << diff_cal0_hist2[i] << "\t" << diff_cal1_hist2[i] << endl;
	
fclose(plik12);

//-------------------------------------------------------------------
// wypelnianie galezi nowego drzewa 
//-------------------------------------------------------------------



Long64_t Nev = tree->GetEntries(); //liczba eventow
	std::cout << "----------------------------------" << std::endl;
	std::cout << " NUMBER OF EVENTS: " << Nev << std::endl;
	std::cout << "----------------------------------" << std::endl;

	for (Long64_t i=0;i<Nev;i++) 
		{ //Loop over events
		tree->GetEntry(i);
		std::cout<<"Processing entry no. "<<i<<" out of "<<Nev<<'\r'<<std::flush;
		if (i==Nev-1) std::cout<<'\n'<<std::endl;

		for (int l = 0;  l < labr->mult(); l++){
			LaBr.id = labr->det[l].id;
			LaBr.time = labr->det[l].time;
			LaBr.energyqshort = labr->det[l].qshort;
			LaBr.energyqlong = labr->det[l].qlong;
			if(nfo.nrun >= 4383 && nfo.nrun <= 4724){
				LaBr.energyqshort_cal = labr_cal2_I[LaBr.id]*LaBr.energyqshort*LaBr.energyqshort + labr_cal1_I[LaBr.id]*LaBr.energyqshort + labr_cal0_I[LaBr.id];
			}else if(nfo.nrun>=6984 && nfo.nrun<=7070){
				LaBr.energyqshort_cal = labr_cal2_II_1[LaBr.id]*LaBr.energyqshort*LaBr.energyqshort + labr_cal1_II_1[LaBr.id]*LaBr.energyqshort + labr_cal0_II_1[LaBr.id];
			}else if(nfo.nrun>=7071 && nfo.nrun<=7157){
				LaBr.energyqshort_cal = labr_cal2_II_2[LaBr.id]*LaBr.energyqshort*LaBr.energyqshort + labr_cal1_II_2[LaBr.id]*LaBr.energyqshort + labr_cal0_II_2[LaBr.id];
			}else if(nfo.nrun>=7158 && nfo.nrun<=7245){
				LaBr.energyqshort_cal = labr_cal2_II_3[LaBr.id]*LaBr.energyqshort*LaBr.energyqshort + labr_cal1_II_3[LaBr.id]*LaBr.energyqshort + labr_cal0_II_3[LaBr.id];
			}else if(nfo.nrun>=7368 && nfo.nrun<=7640){
				LaBr.energyqshort_cal = labr_cal2_III[LaBr.id]*LaBr.energyqshort*LaBr.energyqshort + labr_cal1_III[LaBr.id]*LaBr.energyqshort + labr_cal0_III[LaBr.id];
			}else if(nfo.nrun>=7641 && nfo.nrun<=7838){
				LaBr.energyqshort_cal = labr_cal2_IV[LaBr.id]*LaBr.energyqshort*LaBr.energyqshort + labr_cal1_IV[LaBr.id]*LaBr.energyqshort + labr_cal0_IV[LaBr.id];
			}else{
				LaBr.energyqshort_cal = LaBr.energyqshort;
			}
			
		T1->Fill();
		}
		
		for (int l = 0;  l < plastic->mult(); l++){
			Plastic.id = plastic->det[l].id;
			Plastic.time = plastic->det[l].time;
		T1->Fill();
		}

		for (int l = 0;  l < paris->mult(); l++){
			Paris.id = paris->det[l].id;
			Paris.time = paris->det[l].time;
			Paris.energyqshort = paris->det[l].qshort;
			Paris.energyqlong = paris->det[l].qlong;
		T1->Fill();
		}
		
		for (int l = 0; l < kratta->mult(); l++){
			Kratta.id = kratta->det[l].id;
			Kratta.energyPD0 = kratta->det[l].pd0.ampl;
			Kratta.energyPD1 = kratta->det[l].pd1.ampl;
			Kratta.energyPD2 = kratta->det[l].pd2.ampl;
			Kratta.timePD0 = kratta->det[l].pd0.time;
			Kratta.timePD1 = kratta->det[l].pd1.time;
			Kratta.timePD2 = kratta->det[l].pd2.time;
			if(nfo.nrun >= 4383 && nfo.nrun <= 4724){
				Kratta.energyPD2_cal = kratta_cal2_I[Kratta.id]*Kratta.energyPD2*Kratta.energyPD2 + kratta_cal1_I[Kratta.id]*Kratta.energyPD2 + kratta_cal0_I[Kratta.id];
			}else if(nfo.nrun >= 6984 && nfo.nrun <= 7838){
				Kratta.energyPD2_cal = kratta_cal2_II[Kratta.id]*Kratta.energyPD2*Kratta.energyPD2 + kratta_cal1_II[Kratta.id]*Kratta.energyPD2 + kratta_cal0_II[Kratta.id];
			}else {
				Kratta.energyPD2_cal = Kratta.energyPD2;
			}
		T1->Fill();
		
		}
		
		for (int l = 0;  l < silicon->mult(); l++){
			if(silicon->det[l].id != 26 && silicon->det[l].id != 27 && silicon->det[l].id != 28 && silicon->det[l].id != 29){
				
					Silicon.id = silicon->det[l].id;
					Silicon.energy = silicon->det[l].ampl;
					Silicon.energy_cal = sil_cal0[Silicon.id] + sil_cal1[Silicon.id] * Silicon.energy;
					Silicon.time30 = silicon->det[l].time30;
					Silicon.time80 = silicon->det[l].time80;
					Silicon.difftime = Silicon.time30 - Silicon.time80;
					if(silicon->det[l].id >=20 && silicon->det[l].id <=23){
						if (nfo.nrun >= 7635 && nfo.nrun <= 7838) {
            						Silicon.difftime_cal = diff_cal0_hist1[Silicon.id - 20] + diff_cal1_hist1[Silicon.id - 20] * Silicon.difftime;
        					} else if (nfo.nrun >= 6984 && nfo.nrun <= 7027) {
            						Silicon.difftime_cal = diff_cal0_hist1[Silicon.id - 20] + diff_cal1_hist1[Silicon.id - 20] * Silicon.difftime;
        					} else if (nfo.nrun >= 7028 && nfo.nrun <= 7634) {
            						Silicon.difftime_cal = diff_cal0_hist2[Silicon.id - 20] + diff_cal1_hist2[Silicon.id - 20] * Silicon.difftime;
        					}
        				} 
        				else {
            					Silicon.difftime_cal = diff_cal0[Silicon.id] + diff_cal1[Silicon.id] * Silicon.difftime;
        				}
    					
				
				}
			}
		T1->Fill();
		}
	
T1->Write();

}



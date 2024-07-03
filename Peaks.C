#include "TTree.h"
#include "TCutG.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TFile.h"
#include "TSpectrum.h"
#include "TCanvas.h"
#include "TLine.h"
#include <TH2.h>
#include <TF1.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <algorithm>
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

void Peaks()
{
Double_t times[2] = {1128.782994885779090,1236.521202591038673}; //ID = 25
Double_t timesprim[2];
Double_t cal[2];
Double_t par_g[6];
Double_t par_l[2];
Double_t par_e[2];
Double_t peaks[2];

Double_t div[2]={0,0};

int p=0;

// --------- OPENING THE FILE AND TAKING THE HISTOGRAMS: --------- //
    TFile *myFile = TFile::Open("./13C_data.root");
	
    TTreeReader myReader("tree2", myFile);
    TTreeReaderValue<Int_t> id(myReader, "Silicon.id");
    TTreeReaderValue<Float_t> difftime(myReader, "Silicon.difftime");
    TTreeReaderValue<Int_t> run(myReader, "Info.run");

    std::map<Int_t, TH1F*> histogram_1;
	std::map<Int_t, TH1F*> histogram_2;
	
	for (int i = 20; i <= 23; ++i) {
    	histogram_1[i] = new TH1F(Form("hist1_%02d", i), Form("Silicon.Id = %02d (Run 6984-7027, 7635-7838)", i), 500, 0, 2500);
    	histogram_2[i] = new TH1F(Form("hist2_%02d", i), Form("Silicon.Id = %02d (Run 7028-7634)", i), 500, 0, 2500);
		}
	
	std::map<Int_t, TH1F*> other_histograms;
	
    for (int i = 0; i <= 31; ++i) {
    	if (other_histograms.find(i) == other_histograms.end()) {
        	other_histograms[i] = new TH1F(Form("hist_other_%02d", i), Form("Silicon.Id = %02d", i), 500, 10, 2500);
    		}
		}
	
    while (myReader.Next())
    {
		Int_t currentId = *id;
    	Int_t currentRun = *run;
    	Float_t currentDifftime = *difftime;
    	    
        if (currentId >= 20 && currentId <= 23) {
        	if ((currentRun >= 6984 && currentRun <= 7027) || (currentRun >= 7635 && currentRun <= 7838)) {
            	histogram_1[currentId]->Fill(currentDifftime);
        		} 
        	else if (currentRun >= 7028 && currentRun <= 7634) {
            	histogram_2[currentId]->Fill(currentDifftime);
        		}
    		} 
   
        other_histograms[currentId]->Fill(currentDifftime);
    		   
    	}
printf("\n");
printf("Single histograms\n");
printf("\n");
// --------- NIEDZIELONE HISTOGRAMY: --------- //
	fstream plik;
    plik.open("cal_difftime.txt", std::fstream::out | std::fstream::trunc);
    TFile *fout = new TFile("./graph.root","recreate");
    Double_t FWHM[2][32];
	Double_t ChiSqr[2][32];
	Double_t cal_par[2][32];
// --------- SEARCHING FOR THE PEAKS: --------- //
	for(auto &pair : other_histograms)
	{
	
        TH1F *hist = pair.second;
        hist->Draw();

        printf("START for Silicon.id = %d\n", pair.first);
        printf("\n");

        Int_t nfound = 0;
        TSpectrum *s = new TSpectrum();
        nfound = s->Search(hist,4,"",0.06);
        hist->ShowPeaks(4,"",0.05);
        printf("nfound = %d \n",nfound);

        if (nfound == 0){
            cal[0] = 0;
            cal[1] = 1;
            plik << cal[0] << "               " << cal[1] << endl;
            continue;
        }

// --------- TAKING THE PEAK POSITIONS: --------- //

        Double_t *xpeaks = s->GetPositionX();

        int i=0;
        for (p = 0; p<2; p++){
            printf("Position of the %d peak: %3.2f\n",p+1,xpeaks[p]);

            if(xpeaks[p]>750){
                peaks[i]=xpeaks[p];
                i++;
            }
        }

        printf("\n");

        std::sort(peaks, peaks+2);

        for (p = 0; p<2; p++)
            printf("Position of the %d peak after swap: %3.2f\n",p+1,peaks[p]);

        printf("\n");

// --------- FITTING FUNCTIONS: --------- //

        for (p=0;p<2;p++) {
            Double_t r1 = peaks[p];
            Double_t r2 = peaks[p] + 30;
            Double_t r3 = peaks[p] - 50;
            Double_t r4 = peaks[p] - 15;
            
            TF1 *g = new TF1("g","gaus",r4,r2);
			hist->Fit(g,"NR+");
			g->GetParameters(&par_g[0]);
			peaks[p] = par_g[1];
		
			Double_t i1 = 0.76, i2 = 0.24;
			Double_t p0 = par_g[0];
			Double_t p2 = par_g[2];
			Double_t Mean2 = p0*i2/i1;
			Double_t pozycja = peaks[p] - 5;
		
			TF1 *g2 =new TF1("g2","gaus",r3,r4);
			g2->SetParameters(Mean2,pozycja,p2);
			hist->Fit(g2,"NR+");
			g2->GetParameters(&par_g[3]);
			Double_t pov = g->GetParameter(0);
			
			if (abs(pov)>100000){
					
				g->Delete();
				g2->Delete();
				TF1 *f1 = new TF1("ga","gaus",r4,r2);
        		TF1 *f2 = new TF1("ex","expo",r3,r4);

        		hist->Fit(f1,"NR+");
        		f1->GetParameters(&par_g[0]);
	
        		hist->Fit(f2,"NR+");
        		f2->GetParameters(&par_l[0]);
        		Double_t pov2 = f1->GetParameter(2);
        		
				if(abs(pov2)>98.9){
				
					f1->Delete();
					f2->Delete();
					TF1 *h = new TF1("h","gaus",r3,r2);
					h->SetLineColor(6);
					hist->Fit(h,"R+");
					h->GetParameters(&par_g[0]);
					ChiSqr[p][pair.first] = h->GetChisquare();
        			peaks[p] = h->GetParameter(1);
        			FWHM[p][pair.first] = 2.35*par_g[2];
        			
					}
					
				else{
        		TF1 *sum2 = new TF1("sum2","gaus(0)+[3]*exp((x-[1])/[2])*TMath::Erfc((x-[1])/[2])",r3,r2);
        		sum2->SetParameters(par_g[0],par_g[1],par_g[2],par_l[0]);
        		sum2->SetLineColor(9);
        		hist->Fit(sum2,"R+");
        		ChiSqr[p][pair.first] = sum2->GetChisquare();
        		peaks[p] = par_g[1];
        		FWHM[p][pair.first] = 2.35*par_g[2];
        		
        			}		
				} 
			else {	
						
				TF1 *total = new TF1("total","gaus(0)+gaus(3)",r3,r2);
				total->SetLineColor(3);
				total->SetParameters(par_g);
				hist->Fit(total,"R+");
				ChiSqr[p][pair.first] = total->GetChisquare();
				peaks[p] = par_g[1];
				FWHM[p][pair.first] = 2.35*par_g[2];
				
				}

            printf("Position of the %d peak from fit: %3.2f\n",p+1,peaks[p]);
            printf("FWHM: %3.2f\n",FWHM[p][pair.first]);
        }

// --------- ENERGY CALIBRATION: --------- //
	
        TGraph *gr = new TGraph((sizeof(peaks) / sizeof(Double_t)), peaks, times); 
		TCanvas *grp = new TCanvas(Form("Fit%d",pair.first));
		grp->Divide(1,2);
		gr->Fit("pol1");
		gr->SetMarkerStyle(20);
		gr->SetMarkerColor(4);
		gr->SetLineColor(0);
		gr->SetTitle(Form("Calibration Id = %d", pair.first));
		gr->GetYaxis()->SetTitle("Times");
		gr->GetXaxis()->SetTitle("Channel");
		grp->cd(1);
		gr->Draw();
	 
		TF1 *myfunc = gr->GetFunction("pol1");
	
		for(int i = 0; i<2 ; i++){
			cal[i] = myfunc->GetParameter(i);
			printf("%d cal. coef. for %d spectrum: %3.2f\n",i,pair.first,cal[i]);
			cal_par[i][pair.first] = cal[i];
		}
		plik << cal[0] << "               " << cal[1] << endl;
		
		for(p=0;p<2;p++){
			timesprim[p] = abs(cal[0]+cal[1]*peaks[p]);
			div[p]=times[p]-timesprim[p];
		}
		
		TGraph *gr2 = new TGraph((sizeof(peaks) / sizeof(Double_t)), peaks, div);
		gr2->SetMarkerStyle(20);
		gr2->SetMarkerColor(4);
		gr2->SetLineColor(0);
		gr2->SetTitle("Residuum");
		gr2->GetYaxis()->SetTitle("Delta Times");
		gr2->GetXaxis()->SetTitle("Channel");
	
		grp->cd(2);
		gr2->Draw();
		TF1 *l = new TF1("line","pol0");
		l->SetParameters(0);
		gr2->Fit(l);
	
		fout->cd();
		grp->Write();
    
	}   // end of loop on histograms
	
	plik.close();
    fout->Close();
printf("\n");
printf("Double histograms part 1\n");
printf("\n");
// --------- 1 CZESC PODZIELONEGO HISTOGRAMU: --------- //
	fstream plik2;
    plik2.open("cal_difftime_hist1.txt", std::fstream::out | std::fstream::trunc);
    TFile *fout2 = new TFile("./graph_hist1.root","recreate");
    Double_t FWHM2[2][32];
	Double_t ChiSqr2[2][32];
	Double_t cal_par2[2][32];
// --------- SEARCHING FOR THE PEAKS: --------- //
	for(auto &pair : histogram_1)
	{
	
        TH1F *hist2 = pair.second;
        hist2->Draw();

        printf("START for Silicon.id = %d\n", pair.first);
        printf("\n");

        Int_t nfound = 0;
        TSpectrum *s2 = new TSpectrum();
        nfound = s2->Search(hist2,4,"",0.06);
        hist2->ShowPeaks(4,"",0.05);
        printf("nfound = %d \n",nfound);

        if (nfound == 0){
            cal[0] = 0;
            cal[1] = 1;
            plik2 << cal[0] << "               " << cal[1] << endl;
            continue;
        }

// --------- TAKING THE PEAK POSITIONS: --------- //

        Double_t *xpeaks = s2->GetPositionX();

        int i=0;
        for (p = 0; p<2; p++){
            printf("Position of the %d peak: %3.2f\n",p+1,xpeaks[p]);

            if(xpeaks[p]>750){
                peaks[i]=xpeaks[p];
                i++;
            }
        }

        printf("\n");

        std::sort(peaks, peaks+2);

        for (p = 0; p<2; p++)
            printf("Position of the %d peak after swap: %3.2f\n",p+1,peaks[p]);

        printf("\n");

// --------- FITTING FUNCTIONS: --------- //

        for (p=0;p<2;p++) {
            Double_t r1 = peaks[p];
            Double_t r2 = peaks[p] + 30;
            Double_t r3 = peaks[p] - 50;
            Double_t r4 = peaks[p] - 15;
            
            TF1 *g = new TF1("g","gaus",r4,r2);
			hist2->Fit(g,"NR+");
			g->GetParameters(&par_g[0]);
			peaks[p] = par_g[1];
		
			Double_t i1 = 0.76, i2 = 0.24;
			Double_t p0 = par_g[0];
			Double_t p2 = par_g[2];
			Double_t Mean2 = p0*i2/i1;
			Double_t pozycja = peaks[p] - 5;
		
			TF1 *g2 =new TF1("g2","gaus",r3,r4);
			g2->SetParameters(Mean2,pozycja,p2);
			hist2->Fit(g2,"NR+");
			g2->GetParameters(&par_g[3]);
			Double_t pov = g->GetParameter(0);
			
			if (abs(pov)>100000){
					
				g->Delete();
				g2->Delete();
				TF1 *f1 = new TF1("ga","gaus",r4,r2);
        		TF1 *f2 = new TF1("ex","expo",r3,r4);

        		hist2->Fit(f1,"NR+");
        		f1->GetParameters(&par_g[0]);
	
        		hist2->Fit(f2,"NR+");
        		f2->GetParameters(&par_l[0]);
        		Double_t pov2 = f1->GetParameter(2);
        		
				if(abs(pov2)>98.9){
				
					f1->Delete();
					f2->Delete();
					TF1 *h = new TF1("h","gaus",r3,r2);
					h->SetLineColor(6);
					hist2->Fit(h,"R+");
					h->GetParameters(&par_g[0]);
					ChiSqr2[p][pair.first] = h->GetChisquare();
        			peaks[p] = h->GetParameter(1);
        			FWHM2[p][pair.first] = 2.35*par_g[2];
        			
					}
					
				else{
        		TF1 *sum2 = new TF1("sum2","gaus(0)+[3]*exp((x-[1])/[2])*TMath::Erfc((x-[1])/[2])",r3,r2);
        		sum2->SetParameters(par_g[0],par_g[1],par_g[2],par_l[0]);
        		sum2->SetLineColor(9);
        		hist2->Fit(sum2,"R+");
        		ChiSqr2[p][pair.first] = sum2->GetChisquare();
        		peaks[p] = par_g[1];
        		FWHM2[p][pair.first] = 2.35*par_g[2];
        		
        			}		
				} 
			else {	
						
				TF1 *total = new TF1("total","gaus(0)+gaus(3)",r3,r2);
				total->SetLineColor(3);
				total->SetParameters(par_g);
				hist2->Fit(total,"R+");
				ChiSqr2[p][pair.first] = total->GetChisquare();
				peaks[p] = par_g[1];
				FWHM2[p][pair.first] = 2.35*par_g[2];
				
				}

            printf("Position of the %d peak from fit: %3.2f\n",p+1,peaks[p]);
            printf("FWHM: %3.2f\n",FWHM2[p][pair.first]);
        }

// --------- ENERGY CALIBRATION: --------- //
	
        TGraph *gr = new TGraph((sizeof(peaks) / sizeof(Double_t)), peaks, times); 
		TCanvas *grp = new TCanvas(Form("Fit%d",pair.first));
		grp->Divide(1,2);
		gr->Fit("pol1");
		gr->SetMarkerStyle(20);
		gr->SetMarkerColor(4);
		gr->SetLineColor(0);
		gr->SetTitle(Form("Calibration Id = %d", pair.first));
		gr->GetYaxis()->SetTitle("Times");
		gr->GetXaxis()->SetTitle("Channel");
		grp->cd(1);
		gr->Draw();
	 
		TF1 *myfunc = gr->GetFunction("pol1");
	
		for(int i = 0; i<2 ; i++){
			cal[i] = myfunc->GetParameter(i);
			printf("%d cal. coef. for %d spectrum: %3.2f\n",i,pair.first,cal[i]);
			cal_par2[i][pair.first] = cal[i];
		}
		plik2 << cal[0] << "               " << cal[1] << endl;
		
		for(p=0;p<2;p++){
			timesprim[p] = abs(cal[0]+cal[1]*peaks[p]);
			div[p]=times[p]-timesprim[p];
		}
		
		TGraph *gr2 = new TGraph((sizeof(peaks) / sizeof(Double_t)), peaks, div);
		gr2->SetMarkerStyle(20);
		gr2->SetMarkerColor(4);
		gr2->SetLineColor(0);
		gr2->SetTitle("Residuum");
		gr2->GetYaxis()->SetTitle("Delta Times");
		gr2->GetXaxis()->SetTitle("Channel");
	
		grp->cd(2);
		gr2->Draw();
		TF1 *l = new TF1("line","pol0");
		l->SetParameters(0);
		gr2->Fit(l);
	
		fout2->cd();
		grp->Write();
    
	}   // end of loop on histograms
	
	plik2.close();
    fout2->Close();
printf("\n");
printf("Double histograms part 2\n");
printf("\n");
// ---------2 CZESC PODZIELONEGO HISTOGRAMU: --------- //
	fstream plik3;
    plik3.open("cal_difftime_hist2.txt", std::fstream::out | std::fstream::trunc);
    TFile *fout3 = new TFile("./graph_hist2.root","recreate");
    Double_t FWHM3[2][32];
	Double_t ChiSqr3[2][32];
	Double_t cal_par3[2][32];
// --------- SEARCHING FOR THE PEAKS: --------- //
	for(auto &pair : histogram_2)
	{
	
        TH1F *hist3 = pair.second;
        hist3->Draw();

        printf("START for Silicon.id = %d\n", pair.first);
        printf("\n");

        Int_t nfound = 0;
        TSpectrum *s3 = new TSpectrum();
        nfound = s3->Search(hist3,1.4,"",0.06);
        hist3->ShowPeaks(1.4,"",0.06);
        printf("nfound = %d \n",nfound);

        if (nfound == 0){
            cal[0] = 0;
            cal[1] = 1;
            plik3 << cal[0] << "               " << cal[1] << endl;
            continue;
        }

// --------- TAKING THE PEAK POSITIONS: --------- //

        Double_t *xpeaks = s3->GetPositionX();

        int i=0;
        for (p = 0; p<2; p++){
            printf("Position of the %d peak: %3.2f\n",p+1,xpeaks[p]);

            if(xpeaks[p]>750){
                peaks[i]=xpeaks[p];
                i++;
            }
        }

        printf("\n");

        std::sort(peaks, peaks+2);

        for (p = 0; p<2; p++)
            printf("Position of the %d peak after swap: %3.2f\n",p+1,peaks[p]);

        printf("\n");

// --------- FITTING FUNCTIONS: --------- //

        for (p=0;p<2;p++) {
            Double_t r1 = peaks[p];
            Double_t r2 = peaks[p] + 30;
            Double_t r3 = peaks[p] - 50;
            Double_t r4 = peaks[p] - 15;
            
            TF1 *g = new TF1("g","gaus",r4,r2);
			hist3->Fit(g,"NR+");
			g->GetParameters(&par_g[0]);
			peaks[p] = par_g[1];
		
			Double_t i1 = 0.76, i2 = 0.24;
			Double_t p0 = par_g[0];
			Double_t p2 = par_g[2];
			Double_t Mean2 = p0*i2/i1;
			Double_t pozycja = peaks[p] - 5;
		
			TF1 *g2 =new TF1("g2","gaus",r3,r4);
			g2->SetParameters(Mean2,pozycja,p2);
			hist3->Fit(g2,"NR+");
			g2->GetParameters(&par_g[3]);
			Double_t pov = g->GetParameter(0);
			
			if (abs(pov)>100000){
					
				g->Delete();
				g2->Delete();
				TF1 *f1 = new TF1("ga","gaus",r4,r2);
        		TF1 *f2 = new TF1("ex","expo",r3,r4);

        		hist3->Fit(f1,"NR+");
        		f1->GetParameters(&par_g[0]);
	
        		hist3->Fit(f2,"NR+");
        		f2->GetParameters(&par_l[0]);
        		Double_t pov2 = f1->GetParameter(2);
        		
				if(abs(pov2)>98.9){
				
					f1->Delete();
					f2->Delete();
					TF1 *h = new TF1("h","gaus",r3,r2);
					h->SetLineColor(6);
					hist3->Fit(h,"R+");
					h->GetParameters(&par_g[0]);
					ChiSqr3[p][pair.first] = h->GetChisquare();
        			peaks[p] = h->GetParameter(1);
        			FWHM3[p][pair.first] = 2.35*par_g[2];
        			
					}
					
				else{
        		TF1 *sum2 = new TF1("sum2","gaus(0)+[3]*exp((x-[1])/[2])*TMath::Erfc((x-[1])/[2])",r3,r2);
        		sum2->SetParameters(par_g[0],par_g[1],par_g[2],par_l[0]);
        		sum2->SetLineColor(9);
        		hist3->Fit(sum2,"R+");
        		ChiSqr3[p][pair.first] = sum2->GetChisquare();
        		peaks[p] = par_g[1];
        		FWHM3[p][pair.first] = 2.35*par_g[2];
        		
        			}		
				} 
			else {	
						
				TF1 *total = new TF1("total","gaus(0)+gaus(3)",r3,r2);
				total->SetLineColor(3);
				total->SetParameters(par_g);
				hist3->Fit(total,"R+");
				ChiSqr3[p][pair.first] = total->GetChisquare();
				peaks[p] = par_g[1];
				FWHM3[p][pair.first] = 2.35*par_g[2];
				
				}
			printf("\n");
            printf("Position of the %d peak from fit: %3.2f\n",p+1,peaks[p]);
            printf("FWHM: %3.2f\n",FWHM3[p][pair.first]);
            printf("\n");
        }

// --------- ENERGY CALIBRATION: --------- //
	
        TGraph *gr = new TGraph((sizeof(peaks) / sizeof(Double_t)), peaks, times); 
		TCanvas *grp = new TCanvas(Form("Fit%d",pair.first));
		grp->Divide(1,2);
		gr->Fit("pol1");
		gr->SetMarkerStyle(20);
		gr->SetMarkerColor(4);
		gr->SetLineColor(0);
		gr->SetTitle(Form("Calibration Id = %d", pair.first));
		gr->GetYaxis()->SetTitle("Times");
		gr->GetXaxis()->SetTitle("Channel");
		grp->cd(1);
		gr->Draw();
	 
		TF1 *myfunc = gr->GetFunction("pol1");
	    printf("\n");
		for(int i = 0; i<2 ; i++){
			cal[i] = myfunc->GetParameter(i);
			printf("%d cal. coef. for %d spectrum: %3.2f\n",i,pair.first,cal[i]);
			cal_par3[i][pair.first] = cal[i];
		}
		plik3 << cal[0] << "               " << cal[1] << endl;
		
		for(p=0;p<2;p++){
			timesprim[p] = abs(cal[0]+cal[1]*peaks[p]);
			div[p]=times[p]-timesprim[p];
		}
		printf("\n");
		TGraph *gr2 = new TGraph((sizeof(peaks) / sizeof(Double_t)), peaks, div);
		gr2->SetMarkerStyle(20);
		gr2->SetMarkerColor(4);
		gr2->SetLineColor(0);
		gr2->SetTitle("Residuum");
		gr2->GetYaxis()->SetTitle("Delta Times");
		gr2->GetXaxis()->SetTitle("Channel");
	
		grp->cd(2);
		gr2->Draw();
		TF1 *l = new TF1("line","pol0");
		l->SetParameters(0);
		gr2->Fit(l);
	
		fout3->cd();
		grp->Write();
    	printf("\n");
	}   // end of loop on histograms
	
	plik3.close();
    fout3->Close();
	

// --------- CREATING THE OUTPUT FILE WITH THE HISTOGRAMS: --------- //

    TFile *outputFile = new TFile("output.root", "RECREATE");
    for (auto &pair : histogram_1)
	{
    	pair.second->Write();
	}

	for (auto &pair : histogram_2)
	{
    	pair.second->Write();
	}

	for (auto &pair : other_histograms)
	{
    	pair.second->Write();
	}
	
    outputFile->Close();
    myFile->Close();
}


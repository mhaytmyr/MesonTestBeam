#include <iostream>
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TGraphErrors.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"

const int sampleSize = 1024;

int findMax(int n, double *a) 
{
  if (n <= 0 || !a) 
	return -1;

  double xmax = a[0];
  int loc = 0;
  for(int i = 0; i < n; i++) 
  {
  	if (xmax < a[i]) 
	{
      		xmax = a[i];
      		loc = i;
    	}
  }
  return loc;
}


int findMin(int n, double *a) 
{
  if (n <= 0 || !a) 
	return -1;

  double xmin = a[0];
  int loc = 0;
  for(int i = 0; i < n; i++) 
  {
  	if (a[i] < xmin ) 
	{
      		xmin = a[i];
      		loc = i;
    	}
  }
  return loc;
}

double getRiseTime(double p0, double p1, double halfMin)
{
	return (halfMin-p0)/p1;
}

void analyze_MCP(char *filename)
{
	//read input file
	TFile *fIn = new TFile(filename);
	TTree *t1 = (TTree*)fIn->Get("rec");


	//decalare variables to hold 
	Double_t b1_ch1[1024], b1_ch2[1024], b1_ch3[1024], b1_ch4[1024];
	Double_t b1_t1[1024], b1_t2[1024], b1_t3[1024], b1_t4[1024];
	//
	int ch1Min, ch2Min, ch3Min, ch4Min;
	double risingEdgeCh1X[4],risingEdgeCh1Y[4];
	double risingEdgeCh3X[4],risingEdgeCh3Y[4];
	double risingEdgeCh4X[4],risingEdgeCh4Y[4];
	int ch1Idx = 0, ch2Idx = 0, ch3Idx = 0, ch4Idx = 0;


	//get from branch
	t1->SetBranchAddress("b1_w1",&b1_ch1);
   	t1->SetBranchAddress("b1_w2",&b1_ch2);
   	t1->SetBranchAddress("b1_w3",&b1_ch3);
   	t1->SetBranchAddress("b1_w4",&b1_ch4);

	t1->SetBranchAddress("b1_t1",&b1_t1);
   	t1->SetBranchAddress("b1_t2",&b1_t2);
   	t1->SetBranchAddress("b1_t3",&b1_t3);
   	t1->SetBranchAddress("b1_t4",&b1_t4);


	TCanvas *c1 = new TCanvas("c1","Dynamic Filling Example",200,10,700,500);
	TH2F *Ch3Dynamic = new TH2F("Ch3Dynamic","",40000,-200,200,100,-0.6,0.2); //ns
	TH1F *deltaCh1Ch4 = new TH1F("deltaCh1Ch4","",50000,-1,1); //ns
	TH1F *deltaCh3Ch4 = new TH1F("deltaCh3Ch4","",50000,-1,1); //ns

	Ch3Dynamic->GetXaxis()->SetTitle("Pulse time (ns)");
	Ch3Dynamic->GetYaxis()->SetTitle("Pulse amplitude (V)");
	Ch3Dynamic->SetStats(0000);
	Ch3Dynamic->SetMarkerStyle(20);
	Ch3Dynamic->SetMarkerColor(2);
	Ch3Dynamic->SetMarkerSize(1);


	//read all entries and fill the histograms
   	Long64_t nentries = t1->GetEntries();
	for (Long64_t iEntry=0; iEntry<nentries; iEntry++) 
	//for (Long64_t iEntry=0; iEntry<5; iEntry++) 
	{
     		t1->GetEntry(iEntry);

		std::cout<<"Entry Number "<<iEntry<<std::endl;
		ch1Min = findMin(sampleSize,b1_ch1);
		ch2Min = findMin(sampleSize,b1_ch2);
		ch3Min = findMin(sampleSize,b1_ch3);
		ch4Min = findMin(sampleSize,b1_ch4);

		double ch1halfMin = b1_ch1[ch1Min]/2;
		double ch3halfMin = b1_ch3[ch3Min]/2;
		double ch4halfMin = b1_ch4[ch4Min]/2;
		ch1Idx = 0, ch3Idx = 0; ch4Idx = 0;

		for(int m=0; m<1024; m++)
                {
                        //Ch3Dynamic->Fill(b1_t3[m],b1_ch3[m]);
                        if(b1_ch1[ch1Min]<-0.1)
			Ch3Dynamic->Fill(b1_t1[m],b1_ch1[m]);

			//Do it for Channel 4 first
			if(b1_ch1[m]>ch1halfMin && ch1Min-m<3 && m<ch1Min)
                        {
                                risingEdgeCh1X[ch1Idx] = b1_t1[m];
                                risingEdgeCh1Y[ch1Idx] = b1_ch1[m];
                                ch1Idx+=1;
                        }
                        if(b1_ch1[m]<ch1halfMin && m<ch1Min)
                        {
                                risingEdgeCh1X[ch1Idx] = b1_t1[m];
                                risingEdgeCh1Y[ch1Idx] = b1_ch1[m];
                                ch1Idx+=1;
			}
			
			//Do it for Channel 3 first
			if(b1_ch3[m]>ch3halfMin && ch3Min-m<5 && m<ch3Min)
			{
				risingEdgeCh3X[ch3Idx] = b1_t3[m];
				risingEdgeCh3Y[ch3Idx] = b1_ch3[m];
				ch3Idx+=1;
			}
			if(b1_ch3[m]<ch3halfMin && m<ch3Min)
			{
				risingEdgeCh3X[ch3Idx] = b1_t3[m];
				risingEdgeCh3Y[ch3Idx] = b1_ch3[m];
				ch3Idx+=1;
			}

			//Do it for Channel 4 first
			if(b1_ch4[m]>ch4halfMin && ch4Min-m<5 && m<ch4Min)
                        {
                                risingEdgeCh4X[ch4Idx] = b1_t4[m];
                                risingEdgeCh4Y[ch4Idx] = b1_ch4[m];
                                ch4Idx+=1;
                        }
                        if(b1_ch4[m]<ch4halfMin && m<ch4Min)
                        {
                                risingEdgeCh4X[ch4Idx] = b1_t4[m];
                                risingEdgeCh4Y[ch4Idx] = b1_ch4[m];
                                ch4Idx+=1;
                        }

                }

		//Perform linear fit to Channel3
		TGraphErrors risingEdgeCh1(4,risingEdgeCh1X,risingEdgeCh1Y);
		TFitResultPtr ch1Fit = risingEdgeCh1.Fit("pol1","SFC");
		double ch1halfRise = getRiseTime(ch1Fit->Value(0),ch1Fit->Value(1),ch1halfMin);

		//Perform linear fit to Channel3
		TGraphErrors risingEdgeCh3(4,risingEdgeCh3X,risingEdgeCh3Y);
		TFitResultPtr ch3Fit = risingEdgeCh3.Fit("pol1","SFC");
		double ch3halfRise = getRiseTime(ch3Fit->Value(0),ch3Fit->Value(1),ch3halfMin);

		//Perform linear fit to Channel4
                TGraphErrors risingEdgeCh4(4,risingEdgeCh4X,risingEdgeCh4Y);
                TFitResultPtr ch4Fit = risingEdgeCh4.Fit("pol1","SFC");
                double ch4halfRise = getRiseTime(ch4Fit->Value(0),ch4Fit->Value(1),ch4halfMin);

		std::cout<<"Ch1 Rise Time "<<ch1halfRise<<std::endl;
		std::cout<<"Ch3 Rise Time "<<ch3halfRise<<std::endl;
		std::cout<<"Ch4 Rise Time "<<ch4halfRise<<std::endl;

		deltaCh1Ch4->Fill(ch1halfRise-ch4halfRise);
		deltaCh3Ch4->Fill(ch3halfRise-ch4halfRise);

		//myFit->Print("V");

                //Ch3Dynamic->GetXaxis()->SetRangeUser(b1_t3[ch3Min-50],b1_t3[ch3Min+30]);
                Ch3Dynamic->GetXaxis()->SetRangeUser(b1_t1[ch1Min-50],b1_t1[ch1Min+30]);
                Ch3Dynamic->Draw();
                c1->Modified();
                c1->Update();
                Ch3Dynamic->Reset();
                gSystem->Sleep(350);  //in mictroseconds

	}

	deltaCh3Ch4->Draw();

}

void plot_MCP(char *filename)
{
	//read input file
	TFile *fIn = new TFile(filename);
	TTree *t1 = (TTree*)fIn->Get("rec");

	//decalare variables to hold 
	Double_t b1_ch1[1024], b1_ch2[1024], b1_ch3[1024], b1_ch4[1024];
	Double_t b1_t1[1024], b1_t2[1024], b1_t3[1024], b1_t4[1024];

	//get from branch
	t1->SetBranchAddress("b1_w1",&b1_ch1);
   	t1->SetBranchAddress("b1_w2",&b1_ch2);
   	t1->SetBranchAddress("b1_w3",&b1_ch3);
   	t1->SetBranchAddress("b1_w4",&b1_ch4);

	t1->SetBranchAddress("b1_t1",&b1_t1);
   	t1->SetBranchAddress("b1_t2",&b1_t2);
   	t1->SetBranchAddress("b1_t3",&b1_t3);
   	t1->SetBranchAddress("b1_t4",&b1_t4);


	TCanvas *c1 = new TCanvas("c1","Dynamic Filling Example",200,10,700,500);
	//decalre histograms
	TH1F *delTMax = new TH1F("delTMax","Delta T btw maxima",10000,-10,10); //ns
	TH2F *Ch3Dynamic = new TH2F("Ch3Dynamic","",40000,-200,200,100,-0.6,0.2); //ns

	Ch3Dynamic->GetXaxis()->SetTitle("Pulse time (ns)");
	Ch3Dynamic->GetYaxis()->SetTitle("Pulse amplitude (V)");
	Ch3Dynamic->SetStats(0000);
	Ch3Dynamic->SetMarkerStyle(20);
	//Ch3Dynamic->SetMarkerColor(4);
	Ch3Dynamic->SetMarkerColor(2);
	Ch3Dynamic->SetMarkerSize(1);


	//read all entries and fill the histograms
   	Long64_t nentries = t1->GetEntries();
	for (Long64_t iEntry=0; iEntry<nentries; iEntry++) 
	//for (Long64_t iEntry=0; iEntry<5; iEntry++) 
	{
     		t1->GetEntry(iEntry);

		std::cout<<"Entry Number "<<iEntry<<std::endl;
		int ch1Max = findMax(sampleSize,b1_ch1);
		int ch2Max = findMax(sampleSize,b1_ch2);
		int ch3Max = findMax(sampleSize,b1_ch3);
		int ch4Max = findMax(sampleSize,b1_ch4);

		int ch1Min = findMin(sampleSize,b1_ch1);
		int ch2Min = findMin(sampleSize,b1_ch2);
		int ch3Min = findMin(sampleSize,b1_ch3);
		int ch4Min = findMin(sampleSize,b1_ch4);


		for(int m=0; m<1024; m++)
		{
			Ch3Dynamic->Fill(b1_t3[m],b1_ch3[m]);
			//Ch3Dynamic->Fill(b1_t2[m],b1_ch2[m]);
		}
		//Zoom to the +/-10 marks of the minimum value
		Ch3Dynamic->GetXaxis()->SetRangeUser(b1_t3[ch3Min-50],b1_t3[ch3Min+30]);
		std::cout<<"Min channel "<<b1_t3[ch3Min]<<std::endl;
		std::cout<<"Limits "<<b1_t3[ch3Min-10]<<" : "<<b1_t3[ch3Min+10]<<std::endl;

		Ch3Dynamic->Draw();
		c1->Modified();
		c1->Update();
		/*
		//Save good histograms
		if(b1_ch2[ch2Min]<-0.1)
		{
			char histo[20];
			//sprintf(histo,"b1_ch3_%d.gif",iEntry);
			sprintf(histo,"b1_ch2_%d.gif",iEntry);
			c1->SaveAs(histo);
		}
		*/

		Ch3Dynamic->Reset();
		gSystem->Sleep(500);  //in mictroseconds

  	}

	//delTMax->Draw();
	//Ch3Dynamic->SetMarkerStyle(6);
	//Ch3Dynamic->SetMarkerSize(1);
	//Ch3Dynamic->Draw();
}

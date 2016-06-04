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
#include "TProfile.h"

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


void analyze_MCP(char *filename)
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
	TH2F *Ch3Dynamic = new TH2F("Ch3Dynamic","",40000,-200,200,100,-0.6,0.2); //ns

	Ch3Dynamic->GetXaxis()->SetTitle("Pulse time (ns)");
	Ch3Dynamic->GetYaxis()->SetTitle("Pulse amplitude (V)");
	Ch3Dynamic->SetStats(0000);
	Ch3Dynamic->SetMarkerStyle(20);
	Ch3Dynamic->SetMarkerColor(2);
	Ch3Dynamic->SetMarkerSize(1);


	//read all entries and fill the histograms
   	Long64_t nentries = t1->GetEntries();
	//for (Long64_t iEntry=0; iEntry<nentries; iEntry++) 
	for (Long64_t iEntry=0; iEntry<5; iEntry++) 
	{
     		t1->GetEntry(iEntry);

		std::cout<<"Entry Number "<<iEntry<<std::endl;
		int ch1Min = findMin(sampleSize,b1_ch1);
		int ch2Min = findMin(sampleSize,b1_ch2);
		int ch3Min = findMin(sampleSize,b1_ch3);
		int ch4Min = findMin(sampleSize,b1_ch4);

		for(int m=0; m<1024; m++)
                {
                        Ch3Dynamic->Fill(b1_t3[m],b1_ch3[m]);
                }

		//TF1 *f1 = new TF1("f1","[0]+[1]*x",b1_t3[ch3Min-5],b1_t3[ch3Min]);
		TProfile *prof = Ch3Dynamic->ProfileX();
		prof->Fit("pol1","","",b1_t3[ch3Min-100],b1_t3[ch3Min+100]);		
		prof->SetMarkerColor(kBlue);

                //Zoom to the +/-10 marks of the minimum value
                Ch3Dynamic->GetXaxis()->SetRangeUser(b1_t3[ch3Min-50],b1_t3[ch3Min+30]);
                std::cout<<"Min channel "<<b1_t3[ch3Min]<<std::endl;
                std::cout<<"Limits "<<b1_t3[ch3Min-10]<<" : "<<b1_t3[ch3Min+10]<<std::endl;

                Ch3Dynamic->Draw();
		prof->Draw("same");
                c1->Modified();
                c1->Update();
                Ch3Dynamic->Reset();
		prof->Reset();
                gSystem->Sleep(1500);  //in mictroseconds


	}

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

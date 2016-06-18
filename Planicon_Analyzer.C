#include <iostream>
#include <iomanip>
#include <fstream>
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TLine.h"
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

double getAverage(std::vector<double> points)
{
	double avg = 0;
	for (unsigned int i=0; i<points.size(); i++)
	{
		avg += points[i];	
	}

	return avg/points.size();
}

std::pair<double,double> getFit(std::vector<double> x, std::vector<double> y)
{
	double avgX = getAverage(x);
	double avgY = getAverage(y);
	
	double num, denom;
	for (unsigned int i=0; i<x.size(); i++)
	{
		num+=(x[i]-avgX)*(y[i]-avgY);
		denom+=(x[i]-avgX)*(x[i]-avgX);
	}
	
	std::pair<double,double> myPair;
	myPair.first = num/denom;
	myPair.second = avgY-myPair.first*avgX;

	return myPair;
}

double getChi2(std::vector<double> x, std::vector<double> y, std::pair<double,double> fit)
{
	double chiSquare = 0;
	for (unsigned int i=0; i<x.size(); i++)
	{
		chiSquare += (fit.second+fit.first*x[i]-y[i])*(fit.second+fit.first*x[i]-y[i]);
	}

	return chiSquare/(2*x.size());
}

TLine getFitLine(double x1, double x2, std::pair<double,double> fit)
{

	double y1 = fit.second+fit.first*x1;
	double y2 = fit.second+fit.first*x2;

	TLine line = TLine(x1,y1,x2,y2);
	return line;
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
	std::ofstream outFile("PlaniconTimeDifference.txt",ios::out);
	//std::ofstream outFile("R1aR1bDifference.txt",ios::out);

	//decalare variables to hold 
	Double_t b1_ch1[1024], b1_ch2[1024], b1_ch3[1024], b1_ch4[1024];
	Double_t b1_t1[1024], b1_t2[1024], b1_t3[1024], b1_t4[1024];
	//
	int ch1Min, ch2Min, ch3Min, ch4Min, ch1Max, ch2Max;
	std::vector<double> risingEdgeCh1X, risingEdgeCh1Y;
	std::vector<double> risingEdgeCh2X, risingEdgeCh2Y;
	double risingEdgeCh3X[4],risingEdgeCh3Y[4];
	double risingEdgeCh4X[4],risingEdgeCh4Y[4];
	int ch1Idx = 0, ch2Idx = 0, ch3Idx = 0, ch4Idx = 0;
	double ch1halfRise, ch2halfRise, ch3halfRise, ch4halfRise;

	TGraphErrors risingEdgeCh1, risingEdgeCh2, risingEdgeCh3, risingEdgeCh4;
	TFitResultPtr ch3Fit, ch4Fit;
	std::pair<double,double> ch1Fit, ch2Fit;

	risingEdgeCh1.SetLineColor(kBlue);
	risingEdgeCh2.SetLineColor(kBlue);
	risingEdgeCh3.SetLineColor(kRed);
	risingEdgeCh4.SetLineColor(kRed);

	//get from branch
	t1->SetBranchAddress("b1_w1",&b1_ch1);
   	t1->SetBranchAddress("b1_w2",&b1_ch2);
   	t1->SetBranchAddress("b1_w3",&b1_ch3);
   	t1->SetBranchAddress("b1_w4",&b1_ch4);

	t1->SetBranchAddress("b1_t1",&b1_t1);
   	t1->SetBranchAddress("b1_t2",&b1_t2);
   	t1->SetBranchAddress("b1_t3",&b1_t3);
   	t1->SetBranchAddress("b1_t4",&b1_t4);

	TCanvas *c1 = new TCanvas("c1","Dynamic Filling Example",1400,600);
        c1->Divide(2);
        //decalre histograms
        TH2F *ch1Histo = new TH2F("ch2Histo","",40000,-200,200,100,-0.15,0.15); //ns
        TH2F *ch2Histo = new TH2F("ch2Histo","",40000,-200,200,100,-0.15,0.15); //ns
        TH2F *ch3Histo = new TH2F("ch3Histo","",40000,-200,200,100,-0.6,0.2); //ns
        TH2F *ch4Histo = new TH2F("ch4Histo","",40000,-200,200,100,-0.6,0.2); //ns

	ch1Histo->GetXaxis()->SetTitle("Pulse time (ns)");
        ch1Histo->GetYaxis()->SetTitle("Pulse amplitude (V)");
        ch1Histo->SetStats(0000);
        ch1Histo->GetYaxis()->SetTitleOffset(1.3);
        ch1Histo->SetMarkerStyle(20);
        ch1Histo->SetMarkerColor(kBlue);

	ch2Histo->GetXaxis()->SetTitle("Pulse time (ns)");
        ch2Histo->GetYaxis()->SetTitle("Pulse amplitude (V)");
        ch2Histo->SetStats(0000);
        ch2Histo->GetYaxis()->SetTitleOffset(1.3);
        ch2Histo->SetMarkerStyle(20);
        ch2Histo->SetMarkerColor(kBlue);

        ch3Histo->GetXaxis()->SetTitle("Pulse time (ns)");
        ch3Histo->GetYaxis()->SetTitle("Pulse amplitude (V)");
        ch3Histo->SetStats(0000);
        ch3Histo->SetMarkerStyle(20);
        ch3Histo->SetMarkerColor(kRed);
        ch3Histo->SetMarkerSize(1);

        ch4Histo->GetXaxis()->SetTitle("Pulse time (ns)");
        ch4Histo->GetYaxis()->SetTitle("Pulse amplitude (V)");
        ch4Histo->SetStats(0000);
        ch4Histo->SetMarkerStyle(20);
        ch4Histo->SetMarkerColor(kRed);

	TH1F *deltaCh1Ch4 = new TH1F("deltaCh1Ch4","",50000,-1,1); //ns
	TH1F *deltaCh3Ch4 = new TH1F("deltaCh3Ch4","",50000,-1,1); //ns
	TH1F *nPoint1 = new TH1F("nPoint1","",5,0,5); //ns
	TH1F *nPoint2 = new TH1F("nPoint2","",5,0,5); //ns

	//read all entries and fill the histograms
   	Long64_t nentries = t1->GetEntries();
	//for (Long64_t iEntry=0; iEntry<nentries; iEntry++) 
	for (int iEntry=0; iEntry<nentries; iEntry++) 
	//for (Long64_t iEntry=0; iEntry<5; iEntry++) 
	{
     		t1->GetEntry(iEntry);

		std::cout<<"Entry Number "<<iEntry<<std::endl;
		//for ch1 and ch2 fit between peak and min
		ch1Min = findMin(sampleSize,b1_ch1);
		ch2Min = findMin(sampleSize,b1_ch2);
		ch1Max = findMax(sampleSize,b1_ch1);
		ch2Max = findMax(sampleSize,b1_ch2);
		//for ch3 and ch4 find rising edge
		ch3Min = findMin(sampleSize,b1_ch3);
		ch4Min = findMin(sampleSize,b1_ch4);
	
		double ch1halfMin = b1_ch1[ch1Min]/2;
		double ch2halfMin = b1_ch2[ch2Min]/2;
		double ch3halfMin = b1_ch3[ch3Min]/2;
		double ch4halfMin = b1_ch4[ch4Min]/2;
		ch1Idx = 0, ch2Idx =0, ch3Idx = 0; ch4Idx = 0;

		for(int m=0; m<1024; m++)
                {
                        ch1Histo->Fill(b1_t1[m],b1_ch1[m]);
                        ch2Histo->Fill(b1_t2[m],b1_ch2[m]);
                        ch3Histo->Fill(b1_t3[m],b1_ch3[m]);
                        ch4Histo->Fill(b1_t4[m],b1_ch4[m]);


			//Do it for Channel 1 
			//if(m<=ch1Max && m>=ch1Min && ch1Idx<4 )
			if(m<=ch1Max && m>=ch1Min)
                        {
                                risingEdgeCh1X.push_back(b1_t1[m]);
                                risingEdgeCh1Y.push_back(b1_ch1[m]);
                        }

			//Do it for Channel 2
			//if(m<=ch2Max && m>=ch2Min && ch2Idx<4)
			if(m<=ch2Max && m>=ch2Min)
                        {
                                risingEdgeCh2X.push_back(b1_t2[m]);
                                risingEdgeCh2Y.push_back(b1_ch2[m]);
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
		nPoint1->Fill(ch1Idx);
		nPoint2->Fill(ch2Idx);

		//Perform linear fit to Channel1
		//risingEdgeCh1 = TGraphErrors(4,risingEdgeCh1X,risingEdgeCh1Y);
		//ch1Fit = risingEdgeCh1.Fit("pol1","SFC");
		//ch1halfRise = getRiseTime(ch1Fit->Value(0),ch1Fit->Value(1),0.00);


		//perform fit using my code
		std::pair<double,double> ch1Fit = getFit(risingEdgeCh1X,risingEdgeCh1Y);
		std::pair<double,double> ch2Fit = getFit(risingEdgeCh2X,risingEdgeCh2Y);

		ch1halfRise = getRiseTime(ch1Fit.second,ch1Fit.first,0.00);
		ch2halfRise = getRiseTime(ch2Fit.second,ch2Fit.first,0.00);

		double ch1 = getChi2(risingEdgeCh1X,risingEdgeCh1Y,ch1Fit);
		double ch2 = getChi2(risingEdgeCh2X,risingEdgeCh2Y,ch2Fit);
		std::cout<<"My Result slope "<<ch1Fit.first<<" ; "<<ch1Fit.second<<std::endl;
		std::cout<<"Chisquare "<<ch1<<" ; "<<ch2<<std::endl;

		TLine line1 = getFitLine(risingEdgeCh1X[0],risingEdgeCh1X.back(),ch1Fit);
		TLine line2 = getFitLine(risingEdgeCh2X[0],risingEdgeCh2X.back(),ch2Fit);
		line1.SetLineColor(kBlue);
		line1.SetLineWidth(2);
		line2.SetLineColor(kBlue);
		line2.SetLineWidth(2);


		/*
		//Perform linear fit to Channel3
		risingEdgeCh3 = TGraphErrors(4,risingEdgeCh3X,risingEdgeCh3Y);
		ch3Fit = risingEdgeCh3.Fit("pol1","SFC");
		ch3halfRise = getRiseTime(ch3Fit->Value(0),ch3Fit->Value(1),ch3halfMin);

		//Perform linear fit to Channel4
                risingEdgeCh4 = TGraphErrors(4,risingEdgeCh4X,risingEdgeCh4Y);
                ch4Fit = risingEdgeCh4.Fit("pol1","SFC");
                ch4halfRise = getRiseTime(ch4Fit->Value(0),ch4Fit->Value(1),ch4halfMin);


		//std::cout<<"Ch2 Array "<<sizeof(risingEdgeCh2X)/sizeof(risingEdgeCh2X[0])<<std::endl;
		//for (unsigned int i=0; i<sizeof(risingEdgeCh2X)/sizeof(risingEdgeCh2X[0]); i++)
		//for (unsigned int i=0; i<ch1TestX.size(); i++)
		//{
		//	std::cout<<"Points "<<ch1TestX[i]<<" : "<<ch1TestY[i]<<std::endl;
		//}

		outFile<<std::setprecision(std::numeric_limits<long double>::digits10 + 1)<<
		" "<<ch3halfRise-ch4halfRise<<
		" "<<ch1halfRise-ch2halfRise<<
		" "<<ch3halfRise-ch1halfRise<<
		" "<<ch3halfRise-ch2halfRise<<
		" "<<b1_ch1[ch1Min]<<" "<<b1_ch2[ch2Min]<<" "<<b1_ch3[ch3Min]<<" "<<b1_ch4[ch4Min]<<
		" "<<std::endl;

		//deltaCh1Ch4->Fill(ch1halfRise-ch4halfRise);
		//deltaCh3Ch4->Fill(ch3halfRise-ch4halfRise);
		*/

		/*
		//zoom in
                ch2Histo->GetXaxis()->SetRangeUser(b1_t2[ch2Min-30],b1_t2[ch2Min+30]);
                //ch3Histo->GetXaxis()->SetRangeUser(b1_t3[ch3Min-30],b1_t3[ch3Min+30]);
                ch1Histo->GetXaxis()->SetRangeUser(b1_t1[ch1Min-30],b1_t1[ch1Min+30]);
		//draw ch3 first
		c1->cd(1); ch1Histo->Draw(); line1.Draw("same");
		//draw ch2 next
		c1->cd(2); ch2Histo->Draw(); line2.Draw("same");

                c1->Modified();
                c1->Update();
		//char histo[50];
                //sprintf(histo,"R1aR2a_WithFit_%d.gif",iEntry);
                //c1->SaveAs(histo);

                ch2Histo->Reset(); //ch3Histo->Reset(); 
		ch1Histo->Reset();
                gSystem->Sleep(3500);  //in mictroseconds
		*/

		risingEdgeCh1X.clear(); risingEdgeCh1Y.clear();
		risingEdgeCh2X.clear(); risingEdgeCh2Y.clear();

	}

	//deltaCh3Ch4->Draw();
	c1->cd(1); nPoint1->Draw();
	c1->cd(2); nPoint2->Draw();

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


	TCanvas *c1 = new TCanvas("c1","Dynamic Filling Example",1400,600);
	c1->Divide(2,2);
	//c1->Divide(2);
	//decalre histograms
	TH2F *ch1Histo = new TH2F("ch2Histo","",40000,-200,200,100,-0.5,0.5); //ns
	TH2F *ch2Histo = new TH2F("ch2Histo","",40000,-200,200,100,-0.5,0.5); //ns
	TH2F *ch3Histo = new TH2F("ch3Histo","",40000,-200,200,100,-0.6,0.2); //ns
	TH2F *ch4Histo = new TH2F("ch4Histo","",40000,-200,200,100,-0.6,0.2); //ns

	ch3Histo->GetXaxis()->SetTitle("Pulse time (ns)");
	ch3Histo->GetYaxis()->SetTitle("Pulse amplitude (V)");
	ch3Histo->SetStats(0000);
	ch3Histo->SetMarkerStyle(20);
	ch3Histo->SetMarkerColor(kRed);

	ch4Histo->GetXaxis()->SetTitle("Pulse time (ns)");
	ch4Histo->GetYaxis()->SetTitle("Pulse amplitude (V)");
	ch4Histo->SetStats(0000);
	ch4Histo->SetMarkerStyle(20);
	ch4Histo->SetMarkerColor(kRed);

	ch2Histo->GetXaxis()->SetTitle("Pulse time (ns)");
	ch2Histo->GetYaxis()->SetTitle("Pulse amplitude (V)");
	ch2Histo->SetStats(0000);
        ch2Histo->GetYaxis()->SetTitleOffset(1.3);
	ch2Histo->SetMarkerStyle(20);
	ch2Histo->SetMarkerColor(kBlue);

	ch1Histo->GetXaxis()->SetTitle("Pulse time (ns)");
	ch1Histo->GetYaxis()->SetTitle("Pulse amplitude (V)");
	ch1Histo->SetStats(0000);
	ch1Histo->SetMarkerStyle(20);
	ch1Histo->SetMarkerColor(kBlue);


	//read all entries and fill the histograms
   	Long64_t nentries = t1->GetEntries();
	for (Long64_t iEntry=0; iEntry<nentries; iEntry++) 
	//for (Long64_t iEntry=0; iEntry<10; iEntry++) 
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

		//Draw only good events
		if(b1_ch4[ch4Min]>-0.2) continue;


		for(int m=0; m<1024; m++)
		{
			ch1Histo->Fill(b1_t1[m],b1_ch1[m]);
			ch2Histo->Fill(b1_t2[m],b1_ch2[m]);
			ch3Histo->Fill(b1_t3[m],b1_ch3[m]);
			ch4Histo->Fill(b1_t4[m],b1_ch4[m]);
		}
		//Zoom to the +/-10 marks of the minimum value
		ch1Histo->GetXaxis()->SetRangeUser(b1_t1[ch4Min-50],b1_t1[ch4Min+50]);
		ch2Histo->GetXaxis()->SetRangeUser(b1_t2[ch4Min-50],b1_t2[ch4Min+50]);
		ch3Histo->GetXaxis()->SetRangeUser(b1_t3[ch4Min-50],b1_t3[ch4Min+50]);
		ch4Histo->GetXaxis()->SetRangeUser(b1_t4[ch4Min-30],b1_t4[ch4Min+30]);

		std::cout<<"Min channel "<<b1_t3[ch3Min]<<std::endl;
		std::cout<<"Limits "<<b1_t3[ch3Min-10]<<" : "<<b1_t3[ch3Min+10]<<std::endl;

		c1->cd(1); ch4Histo->Draw();
		//c1->cd(2); ch4Histo->Draw();
		//c1->cd(3); ch2Histo->Draw();
		//c1->cd(4); ch1Histo->Draw();
		c1->Modified();
		c1->Update();
		ch1Histo->Reset(); ch2Histo->Reset(); ch3Histo->Reset(); ch4Histo->Reset();
		gSystem->Sleep(1500);  //in mictroseconds

  	}

	//ch3Histo->SetMarkerStyle(6);
	//ch3Histo->SetMarkerSize(1);
	//ch3Histo->Draw();
}

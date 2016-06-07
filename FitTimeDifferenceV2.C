//RooDataSet* dataSet(0);
//gSystem->Load("libRooFit") ;
using namespace RooFit;

void FitTimeDifferenceV2() 
{  
      
	double delTMin=-0.1, delTMax=0.1;

   	//Selection Variables
   	RooRealVar ch3ch4("ch3ch4","#Delta T (ch3,ch4) (ns)", delTMin, delTMax);
   	RooRealVar ch1ch2("ch1ch2","#Delta T (ch1,ch2) (ns)", -200, 200);
   	RooRealVar ch1ch4("ch1ch4","#Delta T (ch1,ch4) (ns)", -200, 200);
   	RooRealVar ch2ch4("ch2ch4","#Delta T (ch2,ch4) (ns)", -200, 200);
	RooRealVar ch1Min("ch1Min","",-1,1); 
	RooRealVar ch2Min("ch2Min","",-1,1);  
	RooRealVar ch3Min("ch3Min","",-1,1);  
	RooRealVar ch4Min("ch4Min","",-1,1);  

   	//Selection for New Variables
   	RooArgSet variables;
   	variables.add(ch3ch4);
   	variables.add(ch1ch2);
   	variables.add(ch1ch4);
   	variables.add(ch2ch4);
   	variables.add(ch1Min);
   	variables.add(ch2Min);
   	variables.add(ch3Min);
   	variables.add(ch4Min);

   	//Truth resolution model
	//Wide gaussian resolution model
        RooRealVar mean("mean","mean", 2.6,2.0,3.2);
        RooRealVar sigma("sigma","sigma", 0.05, 0, 0.12 );
        RooGaussModel gauss("gauss","gauss", ch2ch4, mean, sigma);
 
  	//Wide gaussian resolution model
  	RooRealVar mean1("mean1","mean", 0.01,0.0,0.02);
  	RooRealVar sigma1("sigma1","sigma", 0.005, 0.001,0.04 );
  	RooGaussModel gauss1("gauss1","gauss", ch3ch4, mean1, sigma1);

  	//Wide gaussian resolution model
  	RooRealVar mean2("mean2","mean", 0.01,-0.02,0.02);
  	RooRealVar sigma2("sigma2","sigma", 0.01, 0.001,0.08 );
  	RooGaussModel gauss2("gauss2","gauss", ch3ch4, mean2, sigma2);
  	//fraction
  	RooRealVar frac("frac","",0.3,0.,1.);

  	//Add two gaussians
  	RooAddModel totalPdf("totalPdf","",RooArgList(gauss1,gauss2),frac);

  	//construct a decay (x) gauss PDF
  	//RooDecay totalPdf("totalPdf","",ch3ch4,tau,gauss,RooDecay::DoubleSided);

  	data = RooDataSet::read("PlaniconTimeDifference.txt",variables,"Q");

  	//RooDataSet dataSet("dataSet","",variables,Import(*data)); 
  	//totalPdf.fitTo(*data,Minos(kTRUE),Save());
  	//gauss.fitTo(*data,Minos(kTRUE),Range(2,3.2));
 
  	//////////////////

  	TCanvas * cx=new TCanvas("cx","cx",800,600);
  
  	//RooPlot *framex = ch3ch4.frame(100);  
  	//RooPlot *framex = ch1ch2.frame(-5,5,1000);  
  	//RooPlot *framex = ch1ch4.frame(-2,7,1000);  
  	RooPlot *framex = ch2ch4.frame(-2,7,1000);  
  	framex->SetTitle("");
  	data->plotOn(framex,DataError(RooAbsData::SumW2),Name("dataX"));
  	//gauss.plotOn(framex,Name("Total_pdf"),LineColor(4));
	/* 
 	totalPdf.plotOn(framex,Name("Total_pdf"),LineColor(4));
  	totalPdf.plotOn(framex,Name("Gauss1"),Components(RooArgSet(gauss1)),LineColor(2),LineStyle(3));
  	totalPdf.plotOn(framex,Name("Gauss2"),Components(RooArgSet(gauss2)),LineColor(3),LineStyle(3));

  	// Legend
  	TLegend* leg = new TLegend(0.75,0.65,0.90,0.90,NULL,"brNDC"); //Standard format
  	//leg->SetFillColor(0);
  	leg->AddEntry(framex->findObject("dataX"),"Data","P");
  	leg->AddEntry(framex->findObject("Total_pdf"),"Total Fit","L");
  	leg->AddEntry(framex->findObject("Gauss1"),"Guassian 1","L");
  	leg->AddEntry(framex->findObject("Gauss2"),"Gaussain 2","L");
	*/
  	framex->Draw();
	//leg->Draw();


}

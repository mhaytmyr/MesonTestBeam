//RooDataSet* dataSet(0);
//gSystem->Load("libRooFit") ;
using namespace RooFit;

void FitTimeRef2aRef2b() 
{  
      
	double delTMin=-0.1, delTMax=0.1;

   	//Selection Variables
   	RooRealVar ch3ch4("ch3ch4","#Delta T (R1a, R1b) (ns)", delTMin, delTMax);
   	RooRealVar ch1ch2("ch1ch2","#Delta T (R2a, R2b) (ns)", -100,100);
   	RooRealVar ch1ch4("ch1ch4","#Delta T (R2a, R1b) (ns)", -200, 200);
   	RooRealVar ch2ch4("ch2ch4","#Delta T (R2b, R1b) (ns)", -200, 200);
	RooRealVar ch1Min("ch1Min","",-1,-0.05); 
	RooRealVar ch2Min("ch2Min","",-1,-0.05);  
	RooRealVar ch3Min("ch3Min","",-1,-0.1);  
	RooRealVar ch4Min("ch4Min","",-1,-0.1);  

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
        //RooRealVar mean("mean","mean", 1.5,2.2);
        RooRealVar mean("mean","mean", 0.1, 0.05, 0.15 );
        RooRealVar sigma("sigma","sigma", 0.05, 0, 0.15 );
        RooGaussModel gauss("gauss","gauss", ch1ch2, mean, sigma);
 
  	//Wide gaussian resolution model
  	RooRealVar mean1("mean1","mean", 0.1,0.0,0.2);
  	RooRealVar sigma1("sigma1","sigma", 0.005, 0.001,0.04 );
  	RooGaussModel gauss1("gauss1","gauss", ch1ch2, mean1, sigma1);

  	//Wide gaussian resolution model
  	RooRealVar mean2("mean2","mean", 0.01,0.0,0.2);
  	RooRealVar sigma2("sigma2","sigma", 0.01, 0.001,0.08 );
  	RooGaussModel gauss2("gauss2","gauss", ch1ch2, mean2, sigma2);
  	//fraction
  	RooRealVar frac("frac","",0.3,0.,1.);

  	//Add two gaussians
  	RooAddModel totalPdf("totalPdf","",RooArgList(gauss1,gauss2),frac);
  	data = RooDataSet::read("PlaniconTimeDifference.txt",variables,"Q");

  	//RooDataSet dataSet("dataSet","",variables,Import(*data)); 
  	totalPdf.fitTo(*data,Minos(kTRUE),Range(0.0,0.2));
 
  	//////////////////

  	TCanvas * cx=new TCanvas("cx","cx",800,600);
  
  	RooPlot *framex = ch1ch2.frame(-0.1,0.20,90);  
  	framex->SetTitle("");
  	data->plotOn(framex,DataError(RooAbsData::SumW2),Name("dataX"));
 	totalPdf.plotOn(framex,Name("Total_pdf"),LineColor(4));
  	//gauss.plotOn(framex,Name("Total_pdf"),LineColor(4));
	totalPdf.paramOn(framex);
 	totalPdf.plotOn(framex,Name("Total_pdf"),LineColor(4));
  	totalPdf.plotOn(framex,Name("Gauss1"),Components(RooArgSet(gauss1)),LineColor(2));
  	totalPdf.plotOn(framex,Name("Gauss2"),Components(RooArgSet(gauss2)),LineColor(3));

  	// Legend
  	TLegend* leg = new TLegend(0.75,0.65,0.90,0.90,NULL,"brNDC"); //Standard format
  	//leg->SetFillColor(0);
  	leg->AddEntry(framex->findObject("dataX"),"Data","P");
  	leg->AddEntry(framex->findObject("Total_pdf"),"Total Fit","L");
  	leg->AddEntry(framex->findObject("Gauss1"),"Guassian 1","L");
  	leg->AddEntry(framex->findObject("Gauss2"),"Gaussain 2","L");
  	framex->Draw();
	leg->Draw();


}

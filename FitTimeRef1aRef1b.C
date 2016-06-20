//RooDataSet* dataSet(0);
//gSystem->Load("libRooFit") ;
using namespace RooFit;

void FitTimeRef1aRef1b() 
{  
      
	double delTMin=-0.1, delTMax=0.1;

   	//Selection Variables
   	RooRealVar ch3ch4("ch3ch4","#Delta T (ch3,ch4) (ns)", delTMin, delTMax);
  
   	//Selection for New Variables
   	RooArgSet variables;
   	variables.add(ch3ch4);

   	//Truth resolution model
  	RooTruthModel tm("tm","Truth",ch3ch4);

  	//Unsmeared decat PDF
  	RooRealVar tau("tau","",-1,1);
  	RooDecay decay("decay","Decay Model",ch3ch4,tau,tm,RooDecay::DoubleSided); 
 
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
  	//data = RooDataSet::read("R1aR1bDifference.txt",variables,"Q");

  	//RooDataSet dataSet("dataSet","",variables,Import(*data)); 
  	totalPdf.fitTo(*data,Minos(kTRUE),Save());
  	//gauss1.fitTo(*data,Minos(kTRUE),Save());
 
  	//////////////////

  	TCanvas * cx=new TCanvas("cx","cx",800,600);
  
  	RooPlot *framex = ch3ch4.frame(100);  
  	framex->SetTitle("");
  	data->plotOn(framex,DataError(RooAbsData::SumW2),Name("dataX"));
  	//gauss1.plotOn(framex,Name("Total_pdf"),LineColor(4));
  	totalPdf.plotOn(framex,Name("Total_pdf"),LineColor(4));
  	totalPdf.plotOn(framex,Name("Gauss1"),Components(RooArgSet(gauss1)),LineColor(2));
  	totalPdf.plotOn(framex,Name("Gauss2"),Components(RooArgSet(gauss2)),LineColor(3));
	totalPdf.paramOn(framex);

  	// Legend
  	TLegend* leg = new TLegend(0.65,0.65,0.90,0.90,NULL,"brNDC"); //Standard format
  	//leg->SetFillColor(0);
  	leg->AddEntry(framex->findObject("dataX"),"Data","P");
  	leg->AddEntry(framex->findObject("Total_pdf"),"Total Fit","L");
  	leg->AddEntry(framex->findObject("Gauss1"),"Guassian 1","L");
  	leg->AddEntry(framex->findObject("Gauss2"),"Gaussain 2","L");

  	framex->Draw();
	leg->Draw();


}

/*
 
   Name:           ConvertBoard2Root.C
   Created by:     Maksat Haytmyradov <mhaytmyr@fnal.gov>
   Date:           May 10th, 2016
 
   Purpose:        Example program under ROOT to read a binary data file written 
                   by the DRSOsc program. Decode time and voltages from waveforms 
                   and display them as a graph. Put values into a ROOT Tree for 
                   further analysis.
 
                   To run it, do:
 
                   - Crate a file test.dat via the "Save" button in DRSOsc
                   - start ROOT
                   root [0] .L ConvertBoard2Root.C+
                   root [1] readDRS("test.dat");

		   - This will produce a file name. test.root
                   - The file contains:
                      TTree rec,
	              Under tree there are braches. Naming convention.
		      -------------------------------
                      Amplitude of board1 ch1: b1_w1
                      Times of board1 ch1: b1_t1 
		      -------------------------------
		      Amplitude of board1 ch2: b1_w2
                      Times of board1 ch12: b1_t2
		     --------------------------------
                     Amplitude of board1 ch3: b1_w3
                      Times of board1 ch3: b1_t3
                     ------------------------------
                     
                      

                     Amplitude of board2 ch1: b2_w1
                     Times of board2 ch1: b2_t1
		     ------------------------------
	             Amplitude of board3 ch1: b3_w1
                     Times of board3 ch1: b3_t1

			

*/



// include std libraries
#include <iostream>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// include ROOT libraries
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TTree.h"
#include "TChain.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TFolder.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TMath.h"
#include "TFile.h"
#include "TSystem.h"
#include "TProfile.h"

using namespace std;
class THEADER {
	public:
	char           bn[2]; //board_number
        unsigned short board_serial_number;
};
class TRIGGER {
        public:
        char           tc[2]; //board_number
        unsigned short trigger_cell;
};
class BINWIDTH{
        public:
                char    chn_header[4];
                float bin_width[1024];
};
class EHEADER {
	public:
		char           event_header[4];
        	unsigned int   event_serial_number;
        	unsigned short year;
        	unsigned short month;
        	unsigned short day;
        	unsigned short hour;
        	unsigned short minute;
        	unsigned short second;
        	unsigned short millisecond;
        	unsigned short reserved1;
};

//int main (int argc, char **argv)
void readDRS(char *filename)
{
	struct THEADER headers[3];
	struct BINWIDTH binwidths[12];

	THEADER th;
	TRIGGER trigger;
	EHEADER eh;
	BINWIDTH bw;
	unsigned short voltage[1024];
	//here first index is board number, second channel number, and third is point in axis;
	double waveform[3][4][1024], time[3][4][1024], t0[3][4][1024];
	//board serial numbers
	unsigned short boards[3];
	//timing calibration of boards
	float bin_width[3][4][1024];
	//channeld header
	char channelHdr[4];
	char rootfile[264];

  	ifstream file;		// read file directly
  	bool endoffile = false;


  	cout <<"********************************************************************" << endl;
  	cout <<"*****              Welcome to DRS4 data analysis               *****" <<endl;
  	cout <<"********************************************************************" << endl;
  	cout << endl;


      	file.open(filename, ios::in | ios::binary);
      	cout << ">> Opening file " << filename << " ......" << endl;
      	if (!file.is_open ())
	{// terminate if the file can't be opened
	 	cerr << "!! File open error:" << filename << endl;
	}


	//open the root file
        strcpy(rootfile, filename);
        if (strchr(rootfile, '.')) *strchr(rootfile, '.') = 0;

        strcat(rootfile, ".root");
        TFile *outfile = new TFile(rootfile, "RECREATE");

	// define the rec tree
        TTree *rec = new TTree("rec","rec");
        rec->Branch("b1_t1", time[0][0], "b1_t1[1024]/D");
        rec->Branch("b1_t2", time[0][1], "b1_t2[1024]/D");
        rec->Branch("b1_t3", time[0][2], "b1_t3[1024]/D");
        rec->Branch("b1_t4", time[0][3], "b1_t4[1024]/D");
        rec->Branch("b2_t1", time[1][0], "b2_t1[1024]/D");
        rec->Branch("b3_t1", time[2][0], "b3_t1[1024]/D");

	rec->Branch("b1_w1", waveform[0][0] ,"b1_w1[1024]/D");
        rec->Branch("b1_w2", waveform[0][1] ,"b1_w2[1024]/D");
        rec->Branch("b1_w3", waveform[0][2] ,"b1_w3[1024]/D");
        rec->Branch("b1_w4", waveform[0][3] ,"b1_w4[1024]/D");
        rec->Branch("b2_w1", waveform[1][0] ,"b2_w1[1024]/D");
        rec->Branch("b3_w1", waveform[2][0] ,"b3_w1[1024]/D");

	//reading header
	char TimeHeader[5];
	file.read((char*)&TimeHeader, 4);
	cout<<TimeHeader<<"\n";

	//loop over boards
	for(int i=0; i<3; i++)
	{
		file.read((char*) &th, sizeof(th));
		//if not board header break loop
                if(th.bn[0] != 'B')
		{
		file.seekg(-sizeof(th), ios::cur);
                break;
		}
		cout<<"Read Board "<<th.bn[0]<<"# "<<th.board_serial_number << "\n";
		//assign board serial number
		boards[i] = th.board_serial_number;

		//now loop over channels in each board
		for(int j=0; j<4; j++)
		{
			file.read((char*) &bw, sizeof(bw));
			bw.chn_header[4]='\0';

			//if event header break loop
			if(bw.chn_header[0]!='C')
			{
			file.seekg(-sizeof(bw), ios::cur);
			break;
			}
        		cout <<"Read Channel "<<bw.chn_header<<" "<< sizeof(bw.bin_width) <<"\n";

			//assign calibration to array
			for(int m=0; m<1024; m++)
			{
				bin_width[i][j][m] = bw.bin_width[m];	
			}

		}//end loop over channels

	}//end loop over boards


  	cout <<"********************************************************************" << endl;
  	cout <<"***********             Beginning event		      **************" <<endl;
  	cout <<"********************************************************************" << endl;
  	cout << endl;

	int totEvent = 0;
	while(!endoffile) //loop over all events
	//for(int n=0; n<1; n++) //loop over only 10 events
	{
		totEvent++;
		cout<<endl;
		cout<<"Processing Event # "<<totEvent<<endl;
		cout<<endl;

		if (file.eof())
		{
                        endoffile = true;
                        break;
		}

		//First Read Event Header
		file.read((char*) &eh, sizeof(eh));
        	eh.event_header[4] = '\0';
		/*
        	cout<<"Event: "<<eh.event_header<<" "<<eh.event_serial_number<<"\n";
        	cout<<"Date: "<<eh.year<<"/"<<eh.month<<"/"<<eh.day<<"\n";
        	cout<<"Time: "<<eh.hour<<":"<<eh.minute<<":"<<eh.second<<":"<<eh.millisecond<<"\n";
		*/

		for(int i=0; i<3; i++) //loop over boards
		{
			//print board serial number for each event
        		file.read((char*) &th, sizeof(th));
			if(th.bn[0] != 'B')
                	{
                		file.seekg(-sizeof(th), ios::cur);
                		break;
                	}
        		cout<<"Event Board #: "<<th.board_serial_number<<"\n";

        		//print trigger cell # for each event
        		file.read((char*) &trigger, sizeof(trigger));
        		cout<<"Event Trigger Cell: "<<trigger.trigger_cell<<"\n";

			//loop over channels in that board
			for(int j=0; j<4; j++)
			{
				//read channels 
                                file.read((char*)&channelHdr, sizeof(channelHdr));
                        	if(channelHdr[0] != 'C')
                        	{
                                	file.seekg(-sizeof(channelHdr), ios::cur);
                                	break;
                        	}
                                cout<<"Channel Header: "<<channelHdr<<"\n";

				//read saved voltages
                                file.read((char*)&voltage, sizeof(voltage));

				for(int k=0; k<1024; k++)
                                {                                               

                                	// convert data to volts
                                        waveform[i][j][k] = (voltage[k] / 65536. - 0.5);
                                        time[i][j][k] = k; 

					t0[i][j][k] = 0;                                        
                                        // calculate time for this cell
                                        for (int m=0; m<k ; m++)
					{
					//time[i][j][k] += bin_width[i][j][m];
                                        t0[i][j][k] += bin_width[i][j][(m+trigger.trigger_cell) % 1024];
					}

					//corrected data	
                                 	time[i][j][k] = t0[i][j][k]-t0[i][j][(1024-trigger.trigger_cell)%1024];
                                }

			}
			

		}//loop over all boards

	// fill root tree
      	rec->Fill();

	}//loop over all event

	// save and close root file
   	rec->Write();
   	outfile->Close();
	file.close();	
}

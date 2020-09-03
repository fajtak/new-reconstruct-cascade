#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"

#include<iostream>
#include<fstream>

TH1F* h_timeRes = new TH1F("h_timeRes","Overall time residuals; #Delta T [ns]; NoE [#]",200,-100,100);
TH2F* h_timeResOMID = new TH2F("h_timeResOMID","Time Residuals per OM;OMID [#]; #Delta T [ns]",300,0,300,200,-100,100);

int ReadFile(TString fileName)
{
	ifstream inputFile;
	inputFile.open(fileName);

	int OMID;
	double timeRes;

	while(!inputFile.eof())
	{
		inputFile >> OMID >> timeRes;
		h_timeRes->Fill(timeRes);
		h_timeResOMID->Fill(OMID,timeRes);
	}

	return 0;
}

void DrawResults()
{
	TCanvas* c_timeRes = new TCanvas("c_timeRes","TimeRes",800,600);
	h_timeRes->Draw();

	TCanvas* c_timeResOMID = new TCanvas("c_timeResOMID","TimeResOMID",800,600);
	h_timeResOMID->Draw("colz");
}

int timeRes()
{
	TString fileName = "/Data/BaikalData/exp16_barsv051/cluster0/0070/timeRes.txt";
	ReadFile(fileName);

	DrawResults();

	return 0;
}
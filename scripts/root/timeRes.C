#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TROOT.h"

#include<iostream>
#include<fstream>

TH1F* h_timeRes = new TH1F("h_timeRes","Overall time residuals; #Delta T [ns]; NoE [#]",200,-100,100);
TH1F* h_timeResOM271 = new TH1F("h_timeResOM271","Time residuals for OM 271; #Delta T [ns]; NoE [#]",200,-100,100);
TH1F* h_nHits = new TH1F("h_nHits","Number of hits per OM; OMID [#]; N_{hits} [#]",300,0,300);
TH2F* h_timeResOMID = new TH2F("h_timeResOMID","Time Residuals per OM;OMID [#]; #Delta T [ns]",300,0,300,100,-100,100);
TGraph* g_qSatCurve = new TGraph();

int ReadFile(TString fileName)
{
	ifstream inputFile;
	inputFile.open(fileName);

	int OMID;
	double timeRes, measQ, expQ;
	int nPoints = 0;

	while(!inputFile.eof())
	{
		inputFile >> OMID >> timeRes;
		// inputFile >> OMID >> timeRes >> measQ >> expQ;
		h_timeRes->Fill(timeRes);
		if (OMID == 271)
			h_timeResOM271->Fill(timeRes);
		h_timeResOMID->Fill(OMID,timeRes);
		h_nHits->Fill(OMID);
		// g_qSatCurve->SetPoint(nPoints,expQ,measQ/expQ);
		nPoints++;
	}

	return 0;
}

void DrawResults()
{
	TCanvas* c_timeRes = new TCanvas("c_timeRes","TimeRes",800,600);
	h_timeRes->Draw();

	TCanvas* c_timeResOM271 = new TCanvas("c_timeResOM271","TimeResOM271",800,600);
	h_timeResOM271->Draw();


	TCanvas* c_timeResOMID = new TCanvas("c_timeResOMID","TimeResOMID",800,600);
	h_timeResOMID->Draw("colz");
	h_timeResOMID->FitSlicesY();
	TH1D* h_timeResOMIDMmean = (TH1D*)gDirectory->Get("h_timeResOMID_1");

	TCanvas* c_timeResOMIDMean = new TCanvas("c_timeResOMIDMean","TimeResOMIDMeans",800,600);
	h_timeResOMIDMmean->Draw();

	TCanvas* c_nHits = new TCanvas("c_nHits","NHits",800,600);
	h_nHits->Draw();

	TCanvas* c_qSatCurve = new TCanvas("c_qSatCurve","QSatCurve",800,600);
	g_qSatCurve->Draw("AP");

}

int timeRes(int season, int cluster, TString fileFolder)
{
	TString fileName = Form("/%s/exp%d/cluster%d/timeRes_%d_%d.txt",fileFolder.Data(),season,cluster,season,cluster);
	// TString fileName = Form("/%s/exp%d_barsv051/cluster%d/0070/timeRes.txt",fileFolder.Data(),season,cluster,season,cluster);
	cout << fileName << endl;
	ReadFile(fileName);

	DrawResults();

	return 0;
}
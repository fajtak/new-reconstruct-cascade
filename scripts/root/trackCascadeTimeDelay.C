#include <vector>
#include "TMath.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TH2F.h"

const double c = 0.299792; // m/ns
const double n = 1.37;
const double v = c/n;
const double theta = TMath::ACos(1/n);

double ExpectedTimeCascade(double x, double y)
{
	return TMath::Sqrt(TMath::Power(x,2)+TMath::Power(y,2))/v;
}

double ExpectedTimeTrack(double x, double y)
{
	double b0 = TMath::Abs(y);
	double c0 = b0/TMath::Sin(theta);
	double a0 = c0*TMath::Cos(theta);
	return (x-a0)/c + c0/v;
}

vector<TGraph*> timeDelay;

int trackCascadeTimeDelay()
{
	double xMin = -200;
	double xMax = 200;
	int xSteps = 400;
	double yMin = -200;
	double yMax = 200;
	int ySteps = 400;
	vector<int> yValues {0,10,20,30,40,50,100};
	double yStepsMultiGraph = yValues.size();

	TMultiGraph* mg_timeDelay = new TMultiGraph("mg_timeDelay",";x [m]; #DeltaT [ns]");

	TH2F* h_timeDelay = new TH2F("h_timeDelay",";x [m]; y[m]; #DeltaT [ns]",400,-200,200,400,-200,200);

	for (int j = 0; j < yStepsMultiGraph; ++j)
	{
		timeDelay.push_back(new TGraph(xSteps));
		timeDelay.back()->SetMarkerColor(j+1);
		timeDelay.back()->SetTitle(Form("y = %d",yValues[j]));

		for (int i = 0; i <= xSteps; ++i)
		{
			timeDelay[j]->SetPoint(i,xMin+i,ExpectedTimeCascade(xMin+i,yValues[j])-ExpectedTimeTrack(xMin+i,yValues[j]));
		}
		mg_timeDelay->Add(timeDelay[j]);
	}

	for (int i = 0; i < xSteps; ++i)
	{
		for (int j = 0; j < ySteps; ++j)
		{
			double timeDelay = ExpectedTimeCascade(xMin+i,yMin+j)-ExpectedTimeTrack(xMin+i,yMin+j);
			if (TMath::Abs(timeDelay) < 50)
				h_timeDelay->Fill(xMin+i,yMin+j,timeDelay);
			// h_timeDelay->Fill(xMin+i,yMin+j,ExpectedTimeCascade(xMin+i,yMin+j));
			// h_timeDelay->Fill(xMin+i,yMin+j,ExpectedTimeTrack(xMin+i,yMin+j));
		}
	}

	TCanvas* c_timeDelay = new TCanvas("c_timeDelay","TimeDelay",800,600);
	gPad->SetGrid();
	mg_timeDelay->Draw("AP");
	gPad->BuildLegend(0.75,0.75,0.95,0.95,"");

	TCanvas* c_timeDelay2D = new TCanvas("c_timeDelay2D","TimeDelay2D",800,600);
	// gPad->SetGrid();
	h_timeDelay->Draw("colz");

	return 0;
}
#include "TCanvas.h"
#include "TString.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TVector3.h"
#include "TMath.h"
#include "TSystem.h"


#include <iostream>

TH1F* h_nHits = new TH1F("h_nHits","Number of hits per Event; N_{hits} [#]; NoE [#]",250,0,250);
TH1F* h_energy = new TH1F("h_energy","Reconstructed energy; E [TeV]; NoE [#]",1000,0,1000);
TH1F* h_theta = new TH1F("h_theta","Zenith angle (0 = up-going, 180 down-going); #theta [deg] ",180,0,180);

void DrawResults()
{
	TCanvas* c_nHits = new TCanvas("c_nHits","NHits",800,600);
	h_nHits->Draw();

	TCanvas* c_energy = new TCanvas("c_energy","Energy",800,600);
	h_energy->Draw();

	TCanvas* c_theta = new TCanvas("c_theta","Theta",800,600);
	h_theta->Draw();
}

int DatastudyRecCas(int year, int cluster = -1, bool upGoing = false, bool highEnergy = true)
{
	TChain reconstructedCascades("Tree/t_RecCasc");

	TString filesDir;

	int startID = cluster!=-1?cluster:0;
	int endID = cluster!=-1?cluster+1:10;

	for (int i = startID; i < endID; ++i)
	{
		filesDir = Form("/Data/BaikalData/dataVal/exp%d/cluster%d/",year,i);
		cout << filesDir << endl;

		auto dir = gSystem->OpenDirectory(filesDir.Data());
		while (auto f = gSystem->GetDirEntry(dir)) 
		{
		  	if (!strcmp(f,".") || !strcmp(f,"..")) continue;
		  	TString fullFilePath = filesDir + f + "/recCascResults.root";
		  	if (!gSystem->AccessPathName(fullFilePath))
		  		reconstructedCascades.Add(TString(filesDir) + f + "/recCascResults.root");
		}
		gSystem->FreeDirectory(dir);
	}

	int runID, eventID, nHits, nHitsAfterCaus, nHitsAfterTFilter, nStringsAfterCaus, nStringsAfterTFilter;
	double energy,theta,phi,mcEnergy,mcTheta,mcPhi;
	double energySigma,thetaSigma,phiSigma,directionSigma;
	double chi2AfterCaus, chi2AfterTFilter, time, likelihood, qTotal;
	TVector3* position = new TVector3();
	TVector3* mcPosition = new TVector3();

	reconstructedCascades.SetBranchAddress("runID", &runID);
	reconstructedCascades.SetBranchAddress("eventID", &eventID);
	reconstructedCascades.SetBranchAddress("nHits", &nHits);
	reconstructedCascades.SetBranchAddress("nHitsAfterCaus", &nHitsAfterCaus);
	reconstructedCascades.SetBranchAddress("nStringsAfterCaus", &nStringsAfterCaus);
	reconstructedCascades.SetBranchAddress("chi2AfterCaus", &chi2AfterCaus);
	reconstructedCascades.SetBranchAddress("nHitsAfterTFilter", &nHitsAfterTFilter);
	reconstructedCascades.SetBranchAddress("nStringsAfterTFilter", &nStringsAfterTFilter);
	reconstructedCascades.SetBranchAddress("chi2AfterTFilter", &chi2AfterTFilter);
	reconstructedCascades.SetBranchAddress("energy", &energy);
	reconstructedCascades.SetBranchAddress("energySigma", &energySigma);
	reconstructedCascades.SetBranchAddress("theta", &theta);
	reconstructedCascades.SetBranchAddress("thetaSigma", &thetaSigma);
	reconstructedCascades.SetBranchAddress("phi", &phi);
	reconstructedCascades.SetBranchAddress("phiSigma", &phiSigma);
	reconstructedCascades.SetBranchAddress("directionSigma", &directionSigma);
	reconstructedCascades.SetBranchAddress("position", &position);
	reconstructedCascades.SetBranchAddress("time", &time);
	reconstructedCascades.SetBranchAddress("mcEnergy", &mcEnergy);
	reconstructedCascades.SetBranchAddress("mcTheta", &mcTheta);
	reconstructedCascades.SetBranchAddress("mcPhi", &mcPhi);
	reconstructedCascades.SetBranchAddress("mcPosition", &mcPosition);
	reconstructedCascades.SetBranchAddress("likelihood", &likelihood);
	reconstructedCascades.SetBranchAddress("qTotal", &qTotal);

	// reconstructedCascades.Print();

	cout << reconstructedCascades.GetEntries() << endl;

	for (int i = 0; i < reconstructedCascades.GetEntries(); ++i)
	{
		reconstructedCascades.GetEntry(i);

		if (directionSigma > 10)
			continue;

		// if (directionSigma > 5 || energy < 10 || energySigma > 5)
			// continue;
		// if (mcEnergy < 20)
			// continue;

		TVector3 cascDirRec(0,0,1);
		cascDirRec.SetTheta(theta);
		cascDirRec.SetPhi(phi);

		if (energy > 100 && highEnergy)
		{
			cout << "Energy above 100 TeV - RunID: " << runID << " EventID: " << eventID << " E = " << energy << " L = " << likelihood << " S = " << directionSigma << " N = " << nHitsAfterTFilter << endl;
			cout << (*position).X() << " " << (*position).Y() << " " << (*position).Z() << endl;
		}

		if (theta/TMath::Pi()*180 < 80 && upGoing)
		{
			cout << "Up-going Event - RunID: " << runID << " EventID: " << eventID << " E = " << energy << " T = " << theta/TMath::Pi()*180 << " S = " << directionSigma << " N = " << nHitsAfterTFilter << endl;
			cout << (*position).X() << " " << (*position).Y() << " " << (*position).Z() << endl;
		}

		// cout << i << " " << runID << " " << eventID << " " << theta  << " " << phi << " " << nHitsAfterTFilter << " " << nStringsAfterTFilter << " " << likelihood << " " << qTotal << endl;
		// cout << energy << " " << energySigma << " " << directionSigma << endl;

		h_nHits->Fill(nHits);
		h_energy->Fill(energy);
		h_theta->Fill(theta/TMath::Pi()*180);
	}

	DrawResults();

	return 0;
}
#include "TCanvas.h"
#include "TString.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TVector3.h"
#include "TMath.h"

#include <iostream>

TH1F* h_mismatchAngle = new TH1F("h_mismatchAngle","Mismatch angle;Mismatch angle [deg];NoE [#]",360,0,360);
TH1F* h_mismatchTheta = new TH1F("h_mismatchTheta","Mismatch theta;Mismatch theta [deg];NoE [#]",360,-180,180);
TH1F* h_mismatchPhi = new TH1F("h_mismatchPhi","Mismatch phi;Mismatch phi [deg];NoE [#]",720,-360,360);
TH1F* h_mismatchPosition = new TH1F("h_mismatchPosition","Mismatch position; Mismatch position [m]; NoE [#]",100,0,100);

TH2F* h_mismatchAngleVsDirSig = new TH2F("h_mismatchAngleVsDirSig","Mismatch angle vs. direction sigma; Mismatch angle [deg]; Sigma [deg]",360,0,360,50,0,50);

void DrawResults()
{
	TCanvas* c_mismatchAngle = new TCanvas("c_mismatchAngle","MismatchAngle",800,600);
	h_mismatchAngle->Draw();

	TCanvas* c_mismatchTheta = new TCanvas("c_mismatchTheta","MismatchTheta",800,600);
	h_mismatchTheta->Draw();

	TCanvas* c_mismatchPhi = new TCanvas("c_mismatchPhi","Results",800,600);
	h_mismatchPhi->Draw();

	Double_t x, q;
	q = 0.5; // 0.5 for "median"
	h_mismatchAngle->ComputeIntegral(); // just a precaution
	h_mismatchAngle->GetQuantiles(1, &x, &q);
	std::cout << "Median mismatch angle = " << x << std::endl;

	// TCanvas* c_mismatchAngleEnergy = new TCanvas("c_mismatchAngleEnergy","Results",800,600);
	// h_mismatchAngleEnergy->Draw("colz");

	// TCanvas* c_mismatchEnergyLog = new TCanvas("c_mismatchEnergyLog","Results",800,600);
	// h_mismatchEnergyLog->Draw();

	// TCanvas* c_mismatchEnergy = new TCanvas("c_mismatchEnergy","Results",800,600);
	// h_mismatchEnergy->Draw();

	TCanvas* c_mismatchPosition = new TCanvas("c_mismatchPosition","MismatchPosition",800,600);	
	h_mismatchPosition->Draw();

	TCanvas* c_mismatchAngleVsDirSig = new TCanvas("c_mis","MismatchAngleVsDirSigma",800,600);
	h_mismatchAngleVsDirSig->Draw("colz");
}

int MCstudyRecCas(int inputFile = 0)
{
	TChain reconstructedCascades("Tree/t_RecCasc");

	TString filesDir;
	switch (inputFile) {
		case 0: 
			filesDir = "/Data/BaikalData/mc/nuatm_feb19/recCascResults.root";
			break;
		case 1:
			filesDir = "/Data/BaikalData/mc/2018may/recCascResults.root";
			break;
		case 2:
			filesDir = "/Data/BaikalData/mc/DZH_cascades/recCascResults.root";
		default:
			break;
	}

	reconstructedCascades.Add(filesDir);

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

		if (directionSigma > 5 || energy < 10 || energySigma > 5)
			continue;
		TVector3 cascDirTrue(0,0,1);
		cascDirTrue.SetTheta(mcTheta);
		cascDirTrue.SetPhi(mcPhi);
		TVector3 cascDirRec(0,0,1);
		cascDirRec.SetTheta(theta);
		cascDirRec.SetPhi(phi);
		cout << i << " " << eventID << " " << mcEnergy << " " << mcTheta << " " << theta << " " << mcPhi << " " << phi << " " << cascDirRec.Angle(cascDirTrue)/TMath::Pi()*180 << " " << nHitsAfterTFilter << " " << nStringsAfterTFilter << " " << likelihood << " " << qTotal << endl;
		cout << mcEnergy << " " << energy << " " << energySigma << " " << directionSigma << endl;
		h_mismatchPosition->Fill(((*position)-(*mcPosition)).Mag());
		h_mismatchAngle->Fill(cascDirRec.Angle(cascDirTrue)/TMath::Pi()*180);
		h_mismatchTheta->Fill((theta-mcTheta)/TMath::Pi()*180);
		h_mismatchPhi->Fill((phi-mcPhi)/TMath::Pi()*180);
		h_mismatchAngleVsDirSig->Fill(cascDirRec.Angle(cascDirTrue)/TMath::Pi()*180,directionSigma);

	}

	DrawResults();

	return 0;
}
#include "TCanvas.h"
#include "TString.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TVector3.h"
#include "TMath.h"

#include <iostream>

std::vector<double> stringXPositions = {5.22,52.13,57.54,25.17,-29.84,-53.6,-42.32,0};
std::vector<double> stringYPositions = {62.32,37.15,-13.92,-52.01,-52.36,-7.49,42.74,0};

TH1F* h_nHits = new TH1F("h_nHits","Number of hits per Event; N_{hits} [#]; NoE [#]",250,0,250);
TH1F* h_nHitsFull = new TH1F("h_nHitsFull","Number of hits per Event; N_{hits} [#]; NoE [#]",250,0,250);
TH1F* h_nHitsAfterTFilter = new TH1F("h_nHitsAfterTFilter","Number of hits per Event used for reconstruction; N_{hits} [#]; NoE [#]",250,0,250);
TH1F* h_nHitsAfterTFilterFull = new TH1F("h_nHitsAfterTFilterFull","Number of hits per Event used for reconstruction; N_{hits} [#]; NoE [#]",250,0,250);
TH1F* h_nStringsAfterTFilter = new TH1F("h_nStringsAfterTFilter","Number of strings per Event used for reconstruction; N_{strings} [#]; NoE [#]",10,0,10);
TH1F* h_nStringsAfterTFilterFull = new TH1F("h_nStringsAfterTFilterFull","Number of strings per Event used for reconstruction; N_{strings} [#]; NoE [#]",10,0,10);
TH1F* h_energy = new TH1F("h_energy","Reconstructed energy; E [TeV]; NoE [#]",1000,0,1000);
TH1F* h_energyFull = new TH1F("h_energyFull","Reconstructed energy; E [TeV]; NoE [#]",1000,0,1000);
TH1F* h_theta = new TH1F("h_theta","Zenith angle (0 = up-going, 180 = down-going); #theta [deg] ",180,0,180);
TH1F* h_thetaFull = new TH1F("h_thetaFull","Zenith angle (0 = up-going, 180 down-going); #theta [deg] ",180,0,180);
TH1F* h_phi = new TH1F("h_phi","Azimuth angle (0 = East, 90 = North); #phi [deg] ",360,0,360);
TH1F* h_phiFull = new TH1F("h_phiFull","Azimuth angle (0 = East, 90 = North); #phi [deg] ",360,0,360);
TH1F* h_thetaContained = new TH1F("h_thetaContained","Zenith angle, contained cascades (0 = up-going, 180 down-going); #theta [deg] ",180,0,180);
TH2F* h_energyNHits = new TH2F("h_energyNHits","Energy vs NHits; log_{10}(E [TeV]); N_{hits} [#]",50,0,5,100,0,100);
TH1F* h_qTotal = new TH1F("h_qTotal","Overall deposited charge;Q [p.e.];NoE [#]",1000,0,10000);
TH1F* h_qTotalFull = new TH1F("h_qTotalFull","Overall deposited charge;Q [p.e.];NoE [#]",1000,0,10000);
TH1F* h_likelihood = new TH1F("h_likelihood","Likelihood;L [#];NoE [#]",100,0,10);
TH1F* h_likelihoodFull = new TH1F("h_likelihoodFull","Likelihood;L [#];NoE [#]",100,0,10);
TGraph* g_cascadeXY = new TGraph();
TGraph* g_stringPositions = new TGraph(8,&stringXPositions[0],&stringYPositions[0]);


void SaveResults(int inputFile)
{
	TString suffix = "";
	if (inputFile == 0)
		suffix = "nuatm_feb19";
	if (inputFile == 1)
		suffix = "muatm_may19";
	TString outputFileName = Form("../../results/mcResults_%s.root",suffix.Data());
	TFile* outputFile = new TFile(outputFileName,"RECREATE");

	h_nHits->Write();
	h_nHitsAfterTFilter->Write();
	h_nStringsAfterTFilter->Write();
	h_energy->Write();
	h_theta->Write();
	h_phi->Write();
	h_qTotal->Write();
	h_likelihood->Write();

	h_nHitsFull->Write();
	h_nHitsAfterTFilterFull->Write();
	h_nStringsAfterTFilterFull->Write();
	h_energyFull->Write();
	h_thetaFull->Write();
	h_phiFull->Write();
	h_qTotalFull->Write();
	h_likelihoodFull->Write();

}

void DrawResults()
{
	TCanvas* c_nHits = new TCanvas("c_nHits","NHits",800,600);
	h_nHits->Draw();

	TCanvas* c_energy = new TCanvas("c_energy","Energy",800,600);
	h_energy->Draw();

	TCanvas* c_theta = new TCanvas("c_theta","Theta",800,600);
	h_theta->Draw();
	h_thetaContained->SetLineColor(kRed);
	h_thetaContained->Draw("same");

	TCanvas* c_energyNHits2 = new TCanvas("c_energyNHits2","EnergyVsNHits2",800,600);
	h_energyNHits->Draw("colz");

	TCanvas* c_cascadePositions = new TCanvas("c_cascadePositions","Results",800,600);

	g_stringPositions->SetMarkerStyle(20);
	g_stringPositions->SetMarkerColor(kRed);
	g_stringPositions->Draw("AP");
	g_stringPositions->SetTitle("Cascade XY positions; X [m]; Y [m]");
	g_cascadeXY->SetMarkerStyle(21);
	g_cascadeXY->SetTitle("Positions of reconstructed cascades;X [m];Y [m]");
	g_cascadeXY->Draw("PSame");

}

bool IsContained(TVector3* position)
{
	if (TMath::Sqrt(TMath::Power(position->X(),2)+TMath::Power(position->Y(),2)) < 60 && TMath::Abs(position->Z() < 265))
		return true;
	else
		return false;
}

int MCstudyRecCas(int inputFile = 0, bool upGoing = false, bool highEnergy = true)
{
	TChain reconstructedCascades("Tree/t_RecCasc");

	TString filesDir;
	const char* env_p = std::getenv("CEPH_MNT");

	switch (inputFile) {
		case 0:
			filesDir = Form("%s/mc/nuatm_feb19/recCascResults.root",env_p);
			break;
		case 1:
			filesDir = Form("%s/mc/muatm_may19/recCascResults.root",env_p);
			break;
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

	int nProcessedEvents = 0;
	int nHighEnergyEvents = 0;

	for (int i = 0; i < reconstructedCascades.GetEntries(); ++i)
	{
		reconstructedCascades.GetEntry(i);

		h_nHitsFull->Fill(nHits);
		h_nHitsAfterTFilterFull->Fill(nHitsAfterTFilter);
		h_nStringsAfterTFilterFull->Fill(nStringsAfterTFilter);
		h_energyFull->Fill(energy);
		h_thetaFull->Fill(theta/TMath::Pi()*180);
		h_phiFull->Fill(phi/TMath::Pi()*180);
		h_qTotalFull->Fill(qTotal);
		h_likelihoodFull->Fill(likelihood);

		if (directionSigma > 10 ||!IsContained(position) || nHitsAfterTFilter < 30)
			continue;

		// if (directionSigma > 5 || energy < 10 || energySigma > 5)
			// continue;
		// if (mcEnergy < 20)
			// continue;

		TVector3 cascDirRec(0,0,1);
		cascDirRec.SetTheta(theta);
		cascDirRec.SetPhi(phi);

		if (energy > 100 && highEnergy && IsContained(position))
		{
			cout << "Energy above 100 TeV - RunID: " << runID << " EventID: " << eventID << " E = " << energy << " L = " << likelihood << " S = " << directionSigma << " N = " << nHitsAfterTFilter << " T = " << theta/TMath::Pi()*180 << " P = " << phi/TMath::Pi()*180 << endl;
			cout << (*position).X() << " " << (*position).Y() << " " << (*position).Z() << endl;
			g_cascadeXY->SetPoint(nHighEnergyEvents,position->X(),position->Y());
			nHighEnergyEvents++;
		}

		if (theta/TMath::Pi()*180 < 80 && upGoing && IsContained(position))
		{
			cout << "Up-going Event - RunID: " << runID << " EventID: " << eventID << " E = " << energy << " T = " << theta/TMath::Pi()*180 << " S = " << directionSigma << " N = " << nHitsAfterTFilter << endl;
			cout << (*position).X() << " " << (*position).Y() << " " << (*position).Z() << endl;
		}

		// cout << i << " " << runID << " " << eventID << " " << theta  << " " << phi << " " << nHitsAfterTFilter << " " << nStringsAfterTFilter << " " << likelihood << " " << qTotal << endl;
		// cout << energy << " " << energySigma << " " << directionSigma << endl;

		h_nHits->Fill(nHits);
		h_nHitsAfterTFilter->Fill(nHitsAfterTFilter);
		h_nStringsAfterTFilter->Fill(nStringsAfterTFilter);
		h_energy->Fill(energy);
		h_theta->Fill(theta/TMath::Pi()*180);
		h_phi->Fill(phi/TMath::Pi()*180);
		h_qTotal->Fill(qTotal);
		h_likelihood->Fill(likelihood);

		if (IsContained(position))
		{
			h_thetaContained->Fill(theta/TMath::Pi()*180);
			h_energyNHits->Fill(TMath::Log10(energy),nHitsAfterTFilter);
		nProcessedEvents++;
			// position->Print();
		}
	}

	DrawResults();
	SaveResults(inputFile);

	cout << nProcessedEvents << endl;

	return 0;
}
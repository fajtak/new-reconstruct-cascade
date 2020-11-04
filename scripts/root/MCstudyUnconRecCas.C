#include "TCanvas.h"
#include "TString.h"
#include "TGraph.h"
#include "TChain.h"
#include "TProfile.h"
#include "TH1F.h"
#include "TFile.h"
#include "TH2F.h"
#include "TVector3.h"
#include "TMath.h"

#include <iostream>

std::vector<double> stringXPositions = {5.22,52.13,57.54,25.17,-29.84,-53.6,-42.32,0};
std::vector<double> stringYPositions = {62.32,37.15,-13.92,-52.01,-52.36,-7.49,42.74,0};

// TH1F* h_nHits = new TH1F("h_nHits","Number of hits per Event; N_{hits} [#]; NoE [#]",250,0,250);
// TH1F* h_nHitsFull = new TH1F("h_nHitsFull","Number of hits per Event; N_{hits} [#]; NoE [#]",250,0,250);
// TH1F* h_nHitsAfterTFilter = new TH1F("h_nHitsAfterTFilter","Number of hits per Event used for reconstruction; N_{hits} [#]; NoE [#]",250,0,250);
// TH1F* h_nHitsAfterTFilterFull = new TH1F("h_nHitsAfterTFilterFull","Number of hits per Event used for reconstruction; N_{hits} [#]; NoE [#]",250,0,250);
// TH1F* h_nStringsAfterTFilter = new TH1F("h_nStringsAfterTFilter","Number of strings per Event used for reconstruction; N_{strings} [#]; NoE [#]",10,0,10);
// TH1F* h_nStringsAfterTFilterFull = new TH1F("h_nStringsAfterTFilterFull","Number of strings per Event used for reconstruction; N_{strings} [#]; NoE [#]",10,0,10);
// TH1F* h_energy = new TH1F("h_energy","Reconstructed energy; E [TeV]; NoE [#]",1000,0,1000);
// TH1F* h_energyFull = new TH1F("h_energyFull","Reconstructed energy; E [TeV]; NoE [#]",1000,0,1000);
// TH1F* h_theta = new TH1F("h_theta","Zenith angle (0 = up-going, 180 = down-going); #theta [deg] ",180,0,180);
// TH1F* h_thetaFull = new TH1F("h_thetaFull","Zenith angle (0 = up-going, 180 down-going); #theta [deg] ",180,0,180);
// TH1F* h_phi = new TH1F("h_phi","Azimuth angle (0 = East, 90 = North); #phi [deg] ",360,0,360);
// TH1F* h_phiFull = new TH1F("h_phiFull","Azimuth angle (0 = East, 90 = North); #phi [deg] ",360,0,360);
// TH1F* h_thetaContained = new TH1F("h_thetaContained","Zenith angle, contained cascades (0 = up-going, 180 down-going); #theta [deg] ",180,0,180);
// TH2F* h_energyNHits = new TH2F("h_energyNHits","Energy vs NHits; log_{10}(E [TeV]); N_{hits} [#]",50,0,5,100,0,100);
// TH1F* h_qTotal = new TH1F("h_qTotal","Overall deposited charge;Q [p.e.];NoE [#]",1000,0,10000);
// TH1F* h_qTotalFull = new TH1F("h_qTotalFull","Overall deposited charge;Q [p.e.];NoE [#]",1000,0,10000);
// TH1F* h_likelihood = new TH1F("h_likelihood","Likelihood;L [#];NoE [#]",100,0,10);
// TH1F* h_likelihoodFull = new TH1F("h_likelihoodFull","Likelihood;L [#];NoE [#]",100,0,10);
// TH1F* h_z = new TH1F("h_z","Z position; z [m]; NoE [#]",60,-300,300);
// TGraph* g_cascadeXY = new TGraph();
// TGraph* g_stringPositions = new TGraph(8,&stringXPositions[0],&stringYPositions[0]);
// TH2F* h_dirError = new TH2F("h_dirError","Mismatch Angle vs. Estimated Error; Mismatch angle [deg.];Estimated Error [deg.]",100,0,100,50,0,25);
// TH1F* h_mismatchAngle = new TH1F("h_mismatchAngle",";Mismatch angle [deg.]; NoE [#]",180,0,180);

TH2F* h_horDistVsMisAngle = new TH2F("h_horDistVsMisAngle","Horizontal distance vs. mismatch angle",200,0,200,180,0,180);
TH2F* h_horDistVsMisEne = new TH2F("h_horDistVsMisEne","Horizontal distance vs. mismatch energy",50,0,200,200,0,100);
TH1F* h_misEne = new TH1F("h_misEne","Mismatch energy",100,0,10);
TH1F* h_horDist100TeV = new TH1F("h_horDist100TeV","Horizontal distance for cascades above 100 TeV",200,0,200);
TH1F* h_horDist100TeVTrue = new TH1F("h_horDist100TeVTrue","Horizontal distance for cascades above 100 TeV",200,0,200);
TH1F* h_horDist60TeV = new TH1F("h_horDist60TeV","Horizontal distance for cascades above 60 TeV",200,0,200);
TH1F* h_horDist60TeVTrue = new TH1F("h_horDist60TeVTrue","Horizontal distance for cascades above 60 TeV",200,0,200);
TH2F* h_nHitsVsMisAngle = new TH2F("h_nHitsVsMisAngle","Number of hits vs. mismatch angle",50,0,100,180,0,180);


void SaveResults(int inputFile)
{
	TString suffix = "";
	// if (inputFile == 0)
	// 	suffix = "nuatm_feb19";
	// if (inputFile == 1)
	// 	suffix = "muatm_may19";
	// if (inputFile == 2)
	// 	suffix = "muatm_jun20";
	// if (inputFile == 3)
	// 	suffix = "muatm_may19_nonHit";
	// if (inputFile == 4)
	// 	suffix = "muatm_jun20_nonHit";
	// TString outputFileName = Form("../../results/mcResults_%s.root",suffix.Data());
	// TFile* outputFile = new TFile(outputFileName,"RECREATE");

	// h_nHits->Write();
	// h_nHitsAfterTFilter->Write();
	// h_nStringsAfterTFilter->Write();
	// h_energy->Write();
	// h_theta->Write();
	// h_phi->Write();
	// h_qTotal->Write();
	// h_likelihood->Write();
	// h_z->Write();
	// h_dirError->Write();

	// h_nHitsFull->Write();
	// h_nHitsAfterTFilterFull->Write();
	// h_nStringsAfterTFilterFull->Write();
	// h_energyFull->Write();
	// h_thetaFull->Write();
	// h_phiFull->Write();
	// h_qTotalFull->Write();
	// h_likelihoodFull->Write();

}

void DrawResults()
{
	// TCanvas* c_nHits = new TCanvas("c_nHits","NHits",800,600);
	// h_nHits->Draw();

	// TCanvas* c_energy = new TCanvas("c_energy","Energy",800,600);
	// h_energy->Draw();

	// TCanvas* c_theta = new TCanvas("c_theta","Theta",800,600);
	// h_theta->Draw();
	// h_thetaContained->SetLineColor(kRed);
	// h_thetaContained->Draw("same");

	// TCanvas* c_energyNHits2 = new TCanvas("c_energyNHits2","EnergyVsNHits2",800,600);
	// h_energyNHits->Draw("colz");

	// TCanvas* c_z = new TCanvas("c_z","Z position",800,600);
	// h_z->Draw();

	// TCanvas* c_cascadePositions = new TCanvas("c_cascadePositions","Results",800,600);

	// g_stringPositions->SetMarkerStyle(20);
	// g_stringPositions->SetMarkerColor(kRed);
	// g_stringPositions->Draw("AP");
	// g_stringPositions->SetTitle("Cascade XY positions; X [m]; Y [m]");
	// g_cascadeXY->SetMarkerStyle(21);
	// g_cascadeXY->SetTitle("Positions of reconstructed cascades;X [m];Y [m]");
	// g_cascadeXY->Draw("PSame");

	// TCanvas* c_dirError = new TCanvas("c_dirError","DirError",800,600);
	// h_dirError->Draw("colz");

	// TCanvas* c_mismatchAngle = new TCanvas("c_mismatchAngle","MismatchAngle",800,600);
	// h_mismatchAngle->Draw();

	TCanvas* c_horDistVsMisAngle = new TCanvas("c_horDistVsMisAngle","HorDistMisAngle",800,600);
	h_horDistVsMisAngle->Draw("colz");

	TCanvas* c_horDistVsMisEne = new TCanvas("c_horDistVsMisEne","HorDistMisEne",800,600);
	h_horDistVsMisEne->Draw("colz");

	TCanvas* c_horDistVsMisEneProf = new TCanvas("c_horDistVsMisEneProf","HorDistMisEneProf",800,600);
	h_horDistVsMisEne->ProfileX()->Draw();

	TCanvas* c_mismatchEnergy = new TCanvas("c_mismatchEnergy","MisEne",800,600);
	h_misEne->Draw("colz");

	TCanvas* c_horDist100TeV = new TCanvas("c_horDist100TeV","HorDist100TeV",800,600);
	h_horDist100TeV->Draw();
	h_horDist100TeVTrue->Draw("same");
	h_horDist100TeVTrue->SetLineColor(kRed);

	TCanvas* c_horDist60TeV = new TCanvas("c_horDist60TeV","HorDist60TeV",800,600);
	h_horDist60TeV->Draw();
	h_horDist60TeVTrue->Draw("same");
	h_horDist60TeVTrue->SetLineColor(kRed);

	TCanvas* c_nHitsVsMisAngle = new TCanvas("c_nHitsVsMisAngle","NHitsMisAngle",800,600);
	h_nHitsVsMisAngle->Draw("colz");

	TCanvas* c_nHitsVsMisAngleProf = new TCanvas("c_nHitsVsMisAngleProf","NHitsMisAngleProf",800,600);
	h_nHitsVsMisAngle->ProfileX()->Draw();
}

bool IsContained(TVector3* position, double distFromCluster = 0)
{
	if (TMath::Sqrt(TMath::Power(position->X(),2)+TMath::Power(position->Y(),2)) < 60+distFromCluster && TMath::Abs(position->Z()) < 265+distFromCluster)
		return true;
	else
		return false;
}

double GetReconstructionError(double theta, double phi, double mcTheta, double mcPhi)
{
	TVector3 recDir(0,0,1);
	recDir.SetTheta(theta);
	recDir.SetPhi(phi);

	TVector3 mcDir(0,0,1);
	mcDir.SetTheta(mcTheta);
	mcDir.SetPhi(mcPhi);

	return recDir.Angle(mcDir);
}

int MCstudyUnconRecCas(int inputFile = 0)
{
	TChain reconstructedCascades("Tree/t_RecCasc");

	TString filesDir;
	const char* env_p = std::getenv("CEPH_MNT");

	switch (inputFile) {
		case 0:
			filesDir = Form("%s/mc/muatm_jun20_val/recCascResults.root",env_p);
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

		// h_nHitsFull->Fill(nHits);
		// h_nHitsAfterTFilterFull->Fill(nHitsAfterTFilter);
		// h_nStringsAfterTFilterFull->Fill(nStringsAfterTFilter);
		// h_energyFull->Fill(energy);
		// h_thetaFull->Fill(theta/TMath::Pi()*180);
		// h_phiFull->Fill(phi/TMath::Pi()*180);
		// h_qTotalFull->Fill(qTotal);
		// h_likelihoodFull->Fill(likelihood);


		double mismatchAngle = GetReconstructionError(theta,phi,mcTheta,mcPhi)/TMath::Pi()*180;
		double horDist = TMath::Sqrt(TMath::Power(position->X(),2)+TMath::Power(position->Y(),2));

		if (horDist > 200)
			continue;

		if (energy > 100)
			h_horDist100TeV->Fill(horDist);
		if (mcEnergy > 100)
			h_horDist100TeVTrue->Fill(horDist);
		if (energy > 60)
			h_horDist60TeV->Fill(horDist);
		if (mcEnergy > 60)
			h_horDist60TeVTrue->Fill(horDist);

		// if (mismatchAngle > 20 && horDist < 60)
		if (energy > 100)
			cout << i << " " << eventID << " D: " << horDist << " A: " << mismatchAngle << " E: " << energy << " MCE: " << mcEnergy << " Z: " << position->Z() << " Q: " << qTotal << endl;
		h_horDistVsMisAngle->Fill(horDist,mismatchAngle);
		h_horDistVsMisEne->Fill(horDist,energy/mcEnergy);
		h_nHitsVsMisAngle->Fill(nHitsAfterTFilter,mismatchAngle);
		h_misEne->Fill(energy/mcEnergy);
		nProcessedEvents++;
		// h_dirError->Fill(mismatchAngle,directionSigma);
		// h_mismatchAngle->Fill(mismatchAngle);

		// if (directionSigma > 20 ||!IsContained(position) || nHitsAfterTFilter < 20 || position->Z() > 240)
		// 	continue;



		// if (directionSigma > 10 || nHitsAfterTFilter < 30)
			// continue;

		// if (directionSigma > 5 || energy < 10 || energySigma > 5)
			// continue;
		// if (mcEnergy < 20)
			// continue;

		// TVector3 cascDirRec(0,0,1);
		// cascDirRec.SetTheta(theta);
		// cascDirRec.SetPhi(phi);

		// if (energy > 100 && highEnergy)
		// {
		// 	cout << "Energy above 100 TeV - RunID: " << runID << " EventID: " << eventID << " E = " << energy << " Et = " << mcEnergy << " L = " << likelihood << " S = " << directionSigma << " N = " << nHitsAfterTFilter << " T = " << theta/TMath::Pi()*180 << " P = " << phi/TMath::Pi()*180 << endl;
		// 	cout << (*position).X() << " " << (*position).Y() << " " << (*position).Z() << endl;
		// 	g_cascadeXY->SetPoint(nHighEnergyEvents,position->X(),position->Y());
		// 	nHighEnergyEvents++;
		// }

		// if (theta/TMath::Pi()*180 < 80 && upGoing)
		// {
		// 	cout << "Up-going Event - RunID: " << runID << " EventID: " << eventID << " E = " << energy << " T = " << theta/TMath::Pi()*180 << " S = " << directionSigma << " N = " << nHitsAfterTFilter << endl;
		// 	cout << (*position).X() << " " << (*position).Y() << " " << (*position).Z() << endl;
		// }

		// cout << i << " " << runID << " " << eventID << " " << theta  << " " << phi << " " << nHitsAfterTFilter << " " << nStringsAfterTFilter << " " << likelihood << " " << qTotal << endl;
		// cout << energy << " " << energySigma << " " << directionSigma << endl;

		// h_nHits->Fill(nHits);
		// h_nHitsAfterTFilter->Fill(nHitsAfterTFilter);
		// h_nStringsAfterTFilter->Fill(nStringsAfterTFilter);
		// h_energy->Fill(energy);
		// h_theta->Fill(theta/TMath::Pi()*180);
		// h_phi->Fill(phi/TMath::Pi()*180);
		// h_qTotal->Fill(qTotal);
		// h_likelihood->Fill(likelihood);
		// h_z->Fill(position->Z());

		// if (IsContained(position))
		// {
		// 	h_thetaContained->Fill(theta/TMath::Pi()*180);
		// 	h_energyNHits->Fill(TMath::Log10(energy),nHitsAfterTFilter);
		// 	nProcessedEvents++;
		// 	// position->Print();
		// }
	}

	DrawResults();
	// SaveResults(inputFile);

	cout << nProcessedEvents << endl;

	return 0;
}
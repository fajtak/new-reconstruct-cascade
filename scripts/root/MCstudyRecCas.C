#include "TCanvas.h"
#include "TString.h"
#include "TGraph.h"
#include "TChain.h"
#include "TH1F.h"
#include "TFile.h"
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
TH1F* h_cosTheta = new TH1F("h_cosTheta","Cosine of the zenith angle; cos(#theta) [1]; NoE [#]",180,-1,1);
TH1F* h_phi = new TH1F("h_phi","Azimuth angle (0 = East, 90 = North); #phi [deg] ",360,0,360);
TH1F* h_phiFull = new TH1F("h_phiFull","Azimuth angle (0 = East, 90 = North); #phi [deg] ",360,0,360);
TH1F* h_thetaContained = new TH1F("h_thetaContained","Zenith angle, contained cascades (0 = up-going, 180 down-going); #theta [deg] ",180,0,180);
TH2F* h_energyNHits = new TH2F("h_energyNHits","Energy vs NHits; log_{10}(E [TeV]); N_{hits} [#]",50,0,5,100,0,100);
TH1F* h_qTotal = new TH1F("h_qTotal","Overall deposited charge;Q [p.e.];NoE [#]",1000,0,10000);
TH1F* h_qTotalFull = new TH1F("h_qTotalFull","Overall deposited charge;Q [p.e.];NoE [#]",1000,0,10000);
TH1F* h_likelihood = new TH1F("h_likelihood","Likelihood;L [#];NoE [#]",1000,0,10);
TH1F* h_likelihoodHitOnly = new TH1F("h_likelihoodHitOnly","Likelihood;L [#];NoE [#]",100,0,10);
TH1F* h_chi2 = new TH1F("h_chi2","Chi2;#chi_{2} [1];NoE [#]",100,0,100);
TH1F* h_likelihoodFull = new TH1F("h_likelihoodFull","Likelihood;L [#];NoE [#]",100,0,10);
TH1F* h_z = new TH1F("h_z","Z position; z [m]; NoE [#]",60,-300,300);
TH1F* h_horDist = new TH1F("h_horDist","Horizontal distance; #rho [m]; NoE [#]",100,0,100);
TH1F* h_nTrackHits = new TH1F("h_nTrackHits","Number of Track Hits; #N_{trackHits} [#]; NoE [#]",30,0,30);
TH2F* h_distQTotal = new TH2F("h_distQTotal","Horizontal distance vs. Qtotal; #rho [m]; QTotal [p.e.]",100,0,100,100,0,3000);
TH2F* h_distEnergy = new TH2F("h_distEnergy","Horizontal distance vs. Energy; #rho [m]; E [TeV]",100,0,100,100,0,1000);
TH2F* h_zenithEnergy = new TH2F("h_zenithEnergy","Zenith vs. Energy; #theta [deq]; E [TeV]",180,0,180,100,0,1000);
TGraph* g_cascadeXY = new TGraph();
TGraph* g_stringPositions = new TGraph(8,&stringXPositions[0],&stringYPositions[0]);
TH2F* h_dirError = new TH2F("h_dirError","Mismatch Angle vs. Estimated Error; Mismatch angle [deg.];Estimated Error [deg.]",100,0,100,50,0,25);
TH1F* h_mismatchAngle = new TH1F("h_mismatchAngle",";Mismatch angle [deg.]; NoE [#]",180,0,180);
TH1F* h_energyRatio = new TH1F("h_energyRatio",";E_{reco}/E_{mc} [1]; NoE [#]",1000,0,100);
TH1F* h_nHitsUpGoing = new TH1F("h_nHitsUpGoing","N_{hits} for up-going events (#theta < 80 deg.); N_{hits} [#]; NoE [#]",100,0,100);
TH1F* h_energyUpGoing = new TH1F("h_energyUpGoing","Energy for up-going events (#theta < 80 deg.); E [TeV] ; NoE [#]",200,0,200);


void SaveResults(int inputFile, int clusterID)
{
	TString suffix = "";
	if (inputFile == 0)
		suffix = "nuatm_feb19";
	if (inputFile == 1)
		suffix = "muatm_may19";
	if (inputFile == 2)
		suffix = "muatm_jun20";
	if (inputFile == 3)
		suffix = "muatm_may19_nonHit";
	if (inputFile == 4)
		suffix = "muatm_jun20_nonHit";
	if (inputFile == 8)
		suffix = "muatm_jun20_val";
	if (inputFile == 9)
		suffix = "muatm_sep20";
	if (inputFile == 10)
		suffix = "nuatm_sep20";
	if (inputFile == 11)
		suffix = "nuatm_sep20";
	if (inputFile == 12)
		suffix = "muatm_sep20";
	if (clusterID != -2)
		suffix += Form("_cluster%d",clusterID);
	TString outputFileName = Form("../../results/mcResults_%s.root",suffix.Data());
	TFile* outputFile = new TFile(outputFileName,"RECREATE");

	h_nHits->Write();
	h_nHitsAfterTFilter->Write();
	h_nStringsAfterTFilter->Write();
	h_energy->Write();
	h_theta->Write();
	h_cosTheta->Write();
	h_phi->Write();
	h_qTotal->Write();
	h_likelihood->Write();
	h_likelihoodHitOnly->Write();
	h_chi2->Write();
	h_z->Write();
	h_horDist->Write();
	h_nTrackHits->Write();
	h_distQTotal->Write();
	h_distEnergy->Write();
	h_zenithEnergy->Write();
	h_dirError->Write();

	h_nHitsFull->Write();
	h_nHitsAfterTFilterFull->Write();
	h_nStringsAfterTFilterFull->Write();
	h_energyFull->Write();
	h_thetaFull->Write();
	h_phiFull->Write();
	h_qTotalFull->Write();
	h_likelihoodFull->Write();

	h_nHitsUpGoing->Write();

}

void DrawResults()
{
	TCanvas* c_nHits = new TCanvas("c_nHits","NHits",800,600);
	h_nHits->Draw();

	TCanvas* c_nHitsAfterTFilter = new TCanvas("c_nHitsAfterTFilter","NHitsAfterTFilter",800,600);
	h_nHitsAfterTFilter->Draw();

	TCanvas* c_energy = new TCanvas("c_energy","Energy",800,600);
	h_energy->Draw();

	TCanvas* c_theta = new TCanvas("c_theta","Theta",800,600);
	h_theta->Draw();
	h_thetaContained->SetLineColor(kRed);
	h_thetaContained->Draw("same");

	TCanvas* c_cosTheta = new TCanvas("c_cosTheta","CosTheta",800,600);
	h_cosTheta->Draw();

	TCanvas* c_energyNHits2 = new TCanvas("c_energyNHits2","EnergyVsNHits2",800,600);
	h_energyNHits->Draw("colz");

	TCanvas* c_z = new TCanvas("c_z","Z position",800,600);
	h_z->Draw();

	TCanvas* c_horDist = new TCanvas("c_horDist","HorizontalDistance",800,600);
	h_horDist->Draw();

	TCanvas* c_distQTotal = new TCanvas("c_distQTotal","DistQTotal",800,600);
	h_distQTotal->Draw("colz");

	TCanvas* c_distEnergy = new TCanvas("c_distEnergy","DistEnergy",800,600);
	h_distEnergy->Draw("colz");

	TCanvas* c_zenithEnergy = new TCanvas("c_zenithEnergy","ZenithEnergy",800,600);
	h_zenithEnergy->Draw("colz");

	TCanvas* c_cascadePositions = new TCanvas("c_cascadePositions","Results",800,600);

	g_stringPositions->SetMarkerStyle(20);
	g_stringPositions->SetMarkerColor(kRed);
	g_stringPositions->Draw("AP");
	g_stringPositions->SetTitle("Cascade XY positions; X [m]; Y [m]");
	g_cascadeXY->SetMarkerStyle(21);
	g_cascadeXY->SetTitle("Positions of reconstructed cascades;X [m];Y [m]");
	g_cascadeXY->Draw("PSame");

	TCanvas* c_dirError = new TCanvas("c_dirError","DirError",800,600);
	h_dirError->Draw("colz");

	TCanvas* c_mismatchAngle = new TCanvas("c_mismatchAngle","MismatchAngle",800,600);
	h_mismatchAngle->Draw();

	TCanvas* c_energyRatio = new TCanvas("c_energyRatio","EnergyRatio",800,600);
	h_energyRatio->Draw();

	TCanvas* c_nHitsUpGoing = new TCanvas("c_nHitsUpGoing","NHitsUpGoing",800,600);
	h_nHitsUpGoing->Draw();

	TCanvas* c_energyUpGoing = new TCanvas("c_energyUpGoing","EnergyUpGoing",800,600);
	h_energyUpGoing->Draw();

}

bool IsContained(TVector3* position, double distFromCluster = 0)
{
	if (TMath::Sqrt(TMath::Power(position->X(),2)+TMath::Power(position->Y(),2)) < 60+distFromCluster && TMath::Abs(position->Z()) < 265+distFromCluster)
		return true;
	else
		return false;
}

bool IsUncontained(TVector3* position, double near, double far)
{
	double horizontalDist = TMath::Sqrt(TMath::Power(position->X(),2)+TMath::Power(position->Y(),2));
	double verticalDist = TMath::Abs(position->Z());
	if ((horizontalDist < far && horizontalDist > near && verticalDist < 263) || (horizontalDist < far && verticalDist < 263+(far-60) && verticalDist > 263+(near-60)))
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

int MCstudyRecCas(int inputFile = 0, int clusterID = -2, bool upGoing = false, bool highEnergy = true)
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
		case 2:
			filesDir = Form("%s/mc/muatm_jun20/recCascResults.root",env_p);
			break;
		case 3:
			filesDir = Form("%s/mc/muatm_may19_nonHit/recCascResults.root",env_p);
			break;
		case 4:
			filesDir = Form("%s/mc/muatm_jun20_nonHit/recCascResults.root",env_p);
			break;
		case 5:
			filesDir = Form("%s/mc/cascades/recCascResults.root",env_p);
			break;
		case 6:
			filesDir = Form("%s/mc/cascades/recCascResults_grid.root",env_p);
			break;
		case 7:
			filesDir = Form("%s/mc/cascades/recCascResults_multiFit.root",env_p);
			break;
		case 8:
			filesDir = Form("%s/mc/muatm_jun20_val/recCascResults.root",env_p);
			break;
		case 9:
			filesDir = Form("%s/mc/muatm_sep20/recCascResults.root",env_p);
			break;
		case 10:
			filesDir = Form("%s/mc/nuatm_sep20_root/recCascResults.root",env_p);
			break;
		case 11:
			filesDir = Form("/Data/BaikalData/mc/nuatm_sep20_root/recCascResults.root",env_p);
			break;
		case 12:
			filesDir = Form("/Data/BaikalData/mc/muatm_sep20/recCascResults.root",env_p);
			break;
		default:
			break;
	}


	if (clusterID != -2)
	{
		int startID = clusterID!=-1?clusterID:0;
		int endID = clusterID!=-1?clusterID+1:10;

		int index = filesDir.Index("recCasc");
		for (int i = startID; i < endID; ++i)
		{
			TString tempDir = filesDir;
			tempDir.Insert(index-1,Form("/cluster%d",i));
			reconstructedCascades.Add(tempDir);
		}
	}else
	{
		reconstructedCascades.Add(filesDir);
	}

	int runID, eventID, nHits, nHitsAfterCaus, nHitsAfterTFilter, nStringsAfterCaus, nStringsAfterTFilter, nTrackHits;
	double energy,theta,phi,mcEnergy,mcTheta,mcPhi;
	double energySigma,thetaSigma,phiSigma,directionSigma;
	double chi2AfterCaus, chi2AfterTFilter, time, likelihood, likelihoodHitOnly, qTotal;
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
	reconstructedCascades.SetBranchAddress("likelihoodHitOnly", &likelihoodHitOnly);
	reconstructedCascades.SetBranchAddress("qTotal", &qTotal);
	reconstructedCascades.SetBranchAddress("nTrackHits", &nTrackHits);


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


		double mismatchAngle = GetReconstructionError(theta,phi,mcTheta,mcPhi)/TMath::Pi()*180;
		h_dirError->Fill(mismatchAngle,directionSigma);
		// if (mismatchAngle > 120 && mismatchAngle < 160)
		// {
		// 	cout << "High mismatch angle " << runID << " EventID: " << eventID << " E = " << energy << " Et = " << mcEnergy << " L = " << likelihood << " S = " << directionSigma << " N = " << nHitsAfterTFilter << " T = " << theta/TMath::Pi()*180 << " Tt = " << mcTheta/TMath::Pi()*180 << " P = " << phi/TMath::Pi()*180 << " Pt = " << mcPhi/TMath::Pi()*180 << " Q = " << qTotal << endl;
		// 	cout << (*position).X() << " " << (*position).Y() << " " << (*position).Z() << endl;
		// 	cout << (*mcPosition).X() << " " << (*mcPosition).Y() << " " << (*mcPosition).Z() << endl;
		// }

		if (!IsContained(position) || likelihoodHitOnly > 1.5 || theta/TMath::Pi()*180 > 80)
		// if (!IsContained(position) || likelihoodHitOnly > 3)
		// if (!IsContained(position,40) || likelihoodHitOnly > 3 || position->Z() > 200)
		// if (!IsUncontained(position,60,100) || likelihoodHitOnly > 3)
			continue;

		h_mismatchAngle->Fill(mismatchAngle);
		h_energyRatio->Fill(energy/mcEnergy);


		// if (directionSigma > 10 || nHitsAfterTFilter < 30)
			// continue;

		// if (directionSigma > 5 || energy < 10 || energySigma > 5)
			// continue;
		// if (mcEnergy < 20)
			// continue;

		TVector3 cascDirRec(0,0,1);
		cascDirRec.SetTheta(theta);
		cascDirRec.SetPhi(phi);

		if (energy > 100 && highEnergy)
		{
			cout << "Energy above 100 TeV - RunID: " << runID << " EventID: " << eventID << " E = " << energy << " Et = " << mcEnergy << " L = " << likelihood << " S = " << directionSigma << " N = " << nHitsAfterTFilter << " T = " << theta/TMath::Pi()*180 << " P = " << phi/TMath::Pi()*180 << " Q = " << qTotal << endl;
			cout << (*position).X() << " " << (*position).Y() << " " << (*position).Z() << endl;
			cout << (*mcPosition).X() << " " << (*mcPosition).Y() << " " << (*mcPosition).Z() << endl;
			g_cascadeXY->SetPoint(nHighEnergyEvents,position->X(),position->Y());
			nHighEnergyEvents++;
		}

		if (theta/TMath::Pi()*180 < 80 && upGoing && likelihoodHitOnly < 1.5)
		{
			cout << "Up-going Event - RunID: " << runID << " EventID: " << eventID << " E = " << energy << " T = " << theta/TMath::Pi()*180 << " S = " << directionSigma << " N = " << nHitsAfterTFilter << endl;
			cout << (*position).X() << " " << (*position).Y() << " " << (*position).Z() << endl;
			h_nHitsUpGoing->Fill(nHitsAfterTFilter);
			h_energyUpGoing->Fill(energy);
		}

		// cout << i << " " << runID << " " << eventID << " " << theta  << " " << phi << " " << nHitsAfterTFilter << " " << nStringsAfterTFilter << " " << likelihood << " " << qTotal << endl;
		// cout << energy << " " << energySigma << " " << directionSigma << endl;

		double horizontalDist = TMath::Sqrt(TMath::Power(position->X(),2)+TMath::Power(position->Y(),2));

		h_nHits->Fill(nHits);
		h_nHitsAfterTFilter->Fill(nHitsAfterTFilter);
		h_nStringsAfterTFilter->Fill(nStringsAfterTFilter);
		h_energy->Fill(energy);
		h_theta->Fill(theta/TMath::Pi()*180);
		h_cosTheta->Fill(TMath::Cos(theta)*(-1));
		h_phi->Fill(phi/TMath::Pi()*180);
		h_qTotal->Fill(qTotal);
		h_likelihood->Fill(likelihood);
		h_likelihoodHitOnly->Fill(likelihoodHitOnly);
		h_chi2->Fill(chi2AfterTFilter);
		h_z->Fill(position->Z());
		h_horDist->Fill(horizontalDist);
		h_nTrackHits->Fill(nTrackHits);
		h_distQTotal->Fill(horizontalDist,qTotal);
		h_distEnergy->Fill(horizontalDist,energy);
		h_zenithEnergy->Fill(theta/TMath::Pi()*180,energy);

		// if (IsContained(position))
		{
			h_thetaContained->Fill(theta/TMath::Pi()*180);
			h_energyNHits->Fill(TMath::Log10(energy),nHitsAfterTFilter);
			nProcessedEvents++;
			// position->Print();
		}
	}

	DrawResults();
	SaveResults(inputFile,clusterID);

	cout << nProcessedEvents << endl;

	return 0;
}
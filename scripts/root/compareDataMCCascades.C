#include "TFile.h"
#include "TH1F.h"
#include "THStack.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TRatioPlot.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TPaveStats.h"
#include "TLegend.h"

#include <iostream>

int compareDataMCCascades(int clusterID)
{
	bool showFull = false;
	bool showTwoMC = false;
	bool showSignalMC = true;
	double mcDataDuration = 0;
	TString mcFile2 = "../../results/mcResults_muatm_jun20_val.root";
	double realDataDuration = 0;
	double newMcDataDuration = 0;
	double signalMcDataDuration = 0;
	TString dataFile = "";
	TString newMcFile = "";
	TString signalMcFile = "";

	switch(clusterID)
	{
		case -1:
			// All Clusters
			realDataDuration = 306.48;
			newMcDataDuration = 167.10/1.47;
			signalMcDataDuration = 10.5*60*TMath::Pi()*1e7*5/3600/24;
			dataFile = "../../results/mcResults_data_y19c-1.root";
			newMcFile = "../../results/mcResults_muatm_sep20_cluster-1.root";
			signalMcFile = "../../results/mcResults_nuatm_sep20_cluster-1.root";
			break;
		case 0:
			// Cluster1
			realDataDuration = 64.064;
			newMcDataDuration = 33.36/1.086;
			signalMcDataDuration = 10.5*60*TMath::Pi()*1e7/3600/24;
			dataFile = "../../results/mcResults_data_y19c0.root";
			newMcFile = "../../results/mcResults_muatm_sep20_cluster0.root";
			signalMcFile = "../../results/mcResults_nuatm_sep20_cluster0.root";
			break;
		case 1:
			// Cluster2
			realDataDuration = 68.50;
			newMcDataDuration = 33.46/1.359;
			signalMcDataDuration = 10.5*60*TMath::Pi()*1e7/3600/24;
			dataFile = "../../results/mcResults_data_y19c1.root";
			newMcFile = "../../results/mcResults_muatm_sep20_cluster1.root";
			signalMcFile = "../../results/mcResults_nuatm_sep20_cluster1.root";
			break;
		case 2:
			// Cluster3
			realDataDuration = 65.12;
			newMcDataDuration = 33.36/1.69;
			signalMcDataDuration = 10.5*60*TMath::Pi()*1e7/3600/24;
			dataFile = "../../results/mcResults_data_y19c2.root";
			newMcFile = "../../results/mcResults_muatm_sep20_cluster2.root";
			signalMcFile = "../../results/mcResults_nuatm_sep20_cluster2.root";
			break;
		case 3:
			// Cluster4
			realDataDuration = 61.26;
			newMcDataDuration = 33.39/1.53;
			signalMcDataDuration = 10.5*60*TMath::Pi()*1e7/3600/24;
			dataFile = "../../results/mcResults_data_y19c3.root";
			newMcFile = "../../results/mcResults_muatm_sep20_cluster3.root";
			signalMcFile = "../../results/mcResults_nuatm_sep20_cluster3.root";
			break;
		case 4:
			// Cluster5
			realDataDuration = 47.5416;
			newMcDataDuration = 33.53/2.00;
			signalMcDataDuration = 10.5*60*TMath::Pi()*1e7/3600/24;
			dataFile = "../../results/mcResults_data_y19c4.root";
			newMcFile = "../../results/mcResults_muatm_sep20_cluster4.root";
			signalMcFile = "../../results/mcResults_nuatm_sep20_cluster4.root";
			break;
	}

	TFile* mcData = new TFile(mcFile2,"READ");
	TH1F* mcNHits  = (TH1F*)mcData->Get("h_nHits");
	mcNHits->SetTitle("Muon group MC");
	mcNHits->SetLineColor(kRed);
	mcNHits->Rebin(5);
	mcNHits->Scale(1/(mcDataDuration*24*3600));
	TH1F* mcNHitsTFilter  = (TH1F*)mcData->Get("h_nHitsAfterTFilter");
	mcNHitsTFilter->SetTitle("Muon group MC");
	mcNHitsTFilter->SetLineColor(kRed);
	mcNHitsTFilter->Rebin(2);
	mcNHitsTFilter->Scale(1/(mcDataDuration*24*3600));
	TH1F* mcNStringsTFilter  = (TH1F*)mcData->Get("h_nStringsAfterTFilter");
	mcNStringsTFilter->SetTitle("Muon group MC");
	mcNStringsTFilter->SetLineColor(kRed);
	// mcNStringsTFilter->Rebin(2);
	mcNStringsTFilter->Scale(1/(mcDataDuration*24*3600));
	TH1F* mcEnergy  = (TH1F*)mcData->Get("h_energy");
	mcEnergy->SetTitle("Muon group MC");
	mcEnergy->SetLineColor(kRed);
	mcEnergy->Rebin(5);
	mcEnergy->Scale(1/(mcDataDuration*24*3600));
	TH1F* mcTheta  = (TH1F*)mcData->Get("h_theta");
	mcTheta->SetTitle("Muon group MC");
	mcTheta->SetLineColor(kRed);
	mcTheta->Rebin(5);
	mcTheta->Scale(1/(mcDataDuration*24*3600));
	TH1F* mcPhi  = (TH1F*)mcData->Get("h_phi");
	mcPhi->SetTitle("Muon group MC");
	mcPhi->SetLineColor(kRed);
	mcPhi->Rebin(10);
	mcPhi->Scale(1/(mcDataDuration*24*3600));
	TH1F* mcQTotal  = (TH1F*)mcData->Get("h_qTotal");
	mcQTotal->SetTitle("Muon group MC");
	mcQTotal->SetLineColor(kRed);
	mcQTotal->Rebin(20);
	mcQTotal->Scale(1/(mcDataDuration*24*3600));
	TH1F* mclikelihood  = (TH1F*)mcData->Get("h_likelihood");
	mclikelihood->SetTitle("Muon group MC");
	mclikelihood->SetLineColor(kRed);
	// mclikelihood->Rebin(2);
	mclikelihood->Scale(1/(mcDataDuration*24*3600));
	TH1F* mcNHitsFull  = (TH1F*)mcData->Get("h_nHitsFull");
	mcNHitsFull->SetTitle("Muon group MC");
	mcNHitsFull->SetLineColor(kRed);
	mcNHitsFull->Rebin(5);
	mcNHitsFull->Scale(1/(mcDataDuration*24*3600));
	TH1F* mcNHitsTFilterFull  = (TH1F*)mcData->Get("h_nHitsAfterTFilterFull");
	mcNHitsTFilterFull->SetTitle("Muon group MC");
	mcNHitsTFilterFull->SetLineColor(kRed);
	mcNHitsTFilterFull->Rebin(2);
	mcNHitsTFilterFull->Scale(1/(mcDataDuration*24*3600));
	TH1F* mcNStringsTFilterFull  = (TH1F*)mcData->Get("h_nStringsAfterTFilterFull");
	mcNStringsTFilterFull->SetTitle("Muon group MC");
	mcNStringsTFilterFull->SetLineColor(kRed);
	// mcNStringsTFilterFull->Rebin(2);
	mcNStringsTFilterFull->Scale(1/(mcDataDuration*24*3600));
	TH1F* mcEnergyFull  = (TH1F*)mcData->Get("h_energyFull");
	mcEnergyFull->SetTitle("Muon group MC");
	mcEnergyFull->SetLineColor(kRed);
	mcEnergyFull->Rebin(5);
	mcEnergyFull->Scale(1/(mcDataDuration*24*3600));
	TH1F* mcThetaFull  = (TH1F*)mcData->Get("h_thetaFull");
	mcThetaFull->SetTitle("Muon group MC");
	mcThetaFull->SetLineColor(kRed);
	mcThetaFull->Rebin(5);
	mcThetaFull->Scale(1/(mcDataDuration*24*3600));
	TH1F* mcPhiFull  = (TH1F*)mcData->Get("h_phiFull");
	mcPhiFull->SetTitle("Muon group MC");
	mcPhiFull->SetLineColor(kRed);
	mcPhiFull->Rebin(10);
	mcPhiFull->Scale(1/(mcDataDuration*24*3600));
	TH1F* mcQTotalFull  = (TH1F*)mcData->Get("h_qTotalFull");
	mcQTotalFull->SetTitle("Muon group MC");
	mcQTotalFull->SetLineColor(kRed);
	mcQTotalFull->Rebin(20);
	mcQTotalFull->Scale(1/(mcDataDuration*24*3600));
	TH1F* mclikelihoodFull  = (TH1F*)mcData->Get("h_likelihoodFull");
	mclikelihoodFull->SetTitle("Muon group MC");
	mclikelihoodFull->SetLineColor(kRed);
	// mclikelihoodFull->Rebin(2);
	mclikelihoodFull->Scale(1/(mcDataDuration*24*3600));


	// TFile* newMcData = new TFile("../../results/mcResults_muatm_jun20.root","READ");
	TFile* newMcData = new TFile(newMcFile,"READ");
	TH1F* newMcNHits  = (TH1F*)newMcData->Get("h_nHits");
	newMcNHits->SetTitle("muatm_sep20");
	newMcNHits->SetLineColor(kBlue);
	newMcNHits->Rebin(5);
	newMcNHits->Scale(1/(newMcDataDuration*24*3600));
	TH1F* newMcNHitsTFilter  = (TH1F*)newMcData->Get("h_nHitsAfterTFilter");
	newMcNHitsTFilter->SetTitle("muatm_sep20");
	newMcNHitsTFilter->SetLineColor(kBlue);
	newMcNHitsTFilter->Rebin(2);
	newMcNHitsTFilter->Scale(1/(newMcDataDuration*24*3600));
	TH1F* newMcNStringsTFilter  = (TH1F*)newMcData->Get("h_nStringsAfterTFilter");
	newMcNStringsTFilter->SetTitle("muatm_sep20");
	newMcNStringsTFilter->SetLineColor(kBlue);
	// newMcNStringsTFilter->Rebin(2);
	newMcNStringsTFilter->Scale(1/(newMcDataDuration*24*3600));
	TH1F* newMcEnergy  = (TH1F*)newMcData->Get("h_energy");
	newMcEnergy->SetTitle("muatm_sep20");
	newMcEnergy->SetLineColor(kBlue);
	newMcEnergy->Rebin(5);
	newMcEnergy->Scale(1/(newMcDataDuration*24*3600));
	TH1F* newMcTheta  = (TH1F*)newMcData->Get("h_theta");
	newMcTheta->SetTitle("muatm_sep20");
	newMcTheta->SetLineColor(kBlue);
	newMcTheta->Rebin(5);
	newMcTheta->Scale(1/(newMcDataDuration*24*3600));
	TH1F* newMcPhi  = (TH1F*)newMcData->Get("h_phi");
	newMcPhi->SetTitle("muatm_sep20");
	newMcPhi->SetLineColor(kBlue);
	newMcPhi->Rebin(10);
	newMcPhi->Scale(1/(newMcDataDuration*24*3600));
	TH1F* newMcQTotal  = (TH1F*)newMcData->Get("h_qTotal");
	newMcQTotal->SetTitle("muatm_sep20");
	newMcQTotal->SetLineColor(kBlue);
	newMcQTotal->Rebin(20);
	newMcQTotal->Scale(1/(newMcDataDuration*24*3600));
	TH1F* newMclikelihood  = (TH1F*)newMcData->Get("h_likelihood");
	newMclikelihood->SetTitle("muatm_sep20");
	newMclikelihood->SetLineColor(kBlue);
	// newMclikelihood->Rebin(2);
	newMclikelihood->Scale(1/(newMcDataDuration*24*3600));
	TH1F* newMclikelihoodHitOnly  = (TH1F*)newMcData->Get("h_likelihoodHitOnly");
	newMclikelihoodHitOnly->SetTitle("muatm_sep20");
	newMclikelihoodHitOnly->SetLineColor(kBlue);
	// newMclikelihoodHitOnly->Rebin(2);
	newMclikelihoodHitOnly->Scale(1/(newMcDataDuration*24*3600));
	TH1F* newMcchi2  = (TH1F*)newMcData->Get("h_chi2");
	newMcchi2->SetTitle("muatm_sep20");
	newMcchi2->SetLineColor(kBlue);
	// newMcchi2->Rebin(2);
	newMcchi2->Scale(1/(newMcDataDuration*24*3600));
	TH1F* newMcZ  = (TH1F*)newMcData->Get("h_z");
	newMcZ->SetTitle("muatm_sep20");
	newMcZ->SetLineColor(kBlue);
	newMcZ->Rebin(2);
	newMcZ->Scale(1/(newMcDataDuration*24*3600));
	TH1F* newMcHorDist  = (TH1F*)newMcData->Get("h_horDist");
	newMcHorDist->SetTitle("muatm_sep20");
	newMcHorDist->SetLineColor(kBlue);
	newMcHorDist->Rebin(2);
	newMcHorDist->Scale(1/(newMcDataDuration*24*3600));
	// TH1F* mcNHitsFull  = (TH1F*)newMcData->Get("h_nHitsFull");
	// mcNHitsFull->SetTitle("NEW Muon group MC");
	// mcNHitsFull->SetLineColor(kBlue);
	// mcNHitsFull->Rebin(5);
	// mcNHitsFull->Scale(1/(newMcDataDuration*24*3600));
	// TH1F* mcNHitsTFilterFull  = (TH1F*)newMcData->Get("h_nHitsAfterTFilterFull");
	// mcNHitsTFilterFull->SetTitle("NEW Muon group MC");
	// mcNHitsTFilterFull->SetLineColor(kBlue);
	// mcNHitsTFilterFull->Rebin(2);
	// mcNHitsTFilterFull->Scale(1/(newMcDataDuration*24*3600));
	// TH1F* mcNStringsTFilterFull  = (TH1F*)newMcData->Get("h_nStringsAfterTFilterFull");
	// mcNStringsTFilterFull->SetTitle("NEW Muon group MC");
	// mcNStringsTFilterFull->SetLineColor(kBlue);
	// // mcNStringsTFilterFull->Rebin(2);
	// mcNStringsTFilterFull->Scale(1/(newMcDataDuration*24*3600));
	// TH1F* mcEnergyFull  = (TH1F*)newMcData->Get("h_energyFull");
	// mcEnergyFull->SetTitle("NEW Muon group MC");
	// mcEnergyFull->SetLineColor(kBlue);
	// mcEnergyFull->Rebin(5);
	// mcEnergyFull->Scale(1/(newMcDataDuration*24*3600));
	// TH1F* mcThetaFull  = (TH1F*)newMcData->Get("h_thetaFull");
	// mcThetaFull->SetTitle("NEW Muon group MC");
	// mcThetaFull->SetLineColor(kBlue);
	// mcThetaFull->Rebin(5);
	// mcThetaFull->Scale(1/(newMcDataDuration*24*3600));
	// TH1F* mcPhiFull  = (TH1F*)newMcData->Get("h_phiFull");
	// mcPhiFull->SetTitle("NEW Muon group MC");
	// mcPhiFull->SetLineColor(kBlue);
	// mcPhiFull->Rebin(10);
	// mcPhiFull->Scale(1/(newMcDataDuration*24*3600));
	// TH1F* mcQTotalFull  = (TH1F*)newMcData->Get("h_qTotalFull");
	// mcQTotalFull->SetTitle("NEW Muon group MC");
	// mcQTotalFull->SetLineColor(kBlue);
	// mcQTotalFull->Rebin(20);
	// mcQTotalFull->Scale(1/(newMcDataDuration*24*3600));
	// TH1F* mclikelihoodFull  = (TH1F*)newMcData->Get("h_likelihoodFull");
	// mclikelihoodFull->SetTitle("NEW Muon group MC");
	// mclikelihoodFull->SetLineColor(kBlue);
	// mclikelihoodFull->Rebin(2);
	// mclikelihoodFull->Scale(1/(newMcDataDuration*24*3600));

	TFile* signalMcData = new TFile(signalMcFile,"READ");
	TH1F* signalMcNHits  = (TH1F*)signalMcData->Get("h_nHits");
	signalMcNHits->SetTitle("nuatm_sep20");
	signalMcNHits->SetLineColor(kRed);
	signalMcNHits->Rebin(5);
	signalMcNHits->Scale(1/(signalMcDataDuration*24*3600));
	TH1F* signalMcNHitsTFilter  = (TH1F*)signalMcData->Get("h_nHitsAfterTFilter");
	signalMcNHitsTFilter->SetTitle("nuatm_sep20");
	signalMcNHitsTFilter->SetLineColor(kRed);
	signalMcNHitsTFilter->Rebin(2);
	signalMcNHitsTFilter->Scale(1/(signalMcDataDuration*24*3600));
	TH1F* signalMcNStringsTFilter  = (TH1F*)signalMcData->Get("h_nStringsAfterTFilter");
	signalMcNStringsTFilter->SetTitle("nuatm_sep20");
	signalMcNStringsTFilter->SetLineColor(kRed);
	// signalMcNStringsTFilter->Rebin(2);
	signalMcNStringsTFilter->Scale(1/(signalMcDataDuration*24*3600));
	TH1F* signalMcEnergy  = (TH1F*)signalMcData->Get("h_energy");
	signalMcEnergy->SetTitle("nuatm_sep20");
	signalMcEnergy->SetLineColor(kRed);
	signalMcEnergy->Rebin(5);
	signalMcEnergy->Scale(1/(signalMcDataDuration*24*3600));
	TH1F* signalMcTheta  = (TH1F*)signalMcData->Get("h_theta");
	signalMcTheta->SetTitle("nuatm_sep20");
	signalMcTheta->SetLineColor(kRed);
	signalMcTheta->Rebin(5);
	signalMcTheta->Scale(1/(signalMcDataDuration*24*3600));
	TH1F* signalMcPhi  = (TH1F*)signalMcData->Get("h_phi");
	signalMcPhi->SetTitle("nuatm_sep20");
	signalMcPhi->SetLineColor(kRed);
	signalMcPhi->Rebin(10);
	signalMcPhi->Scale(1/(signalMcDataDuration*24*3600));
	TH1F* signalMcQTotal  = (TH1F*)signalMcData->Get("h_qTotal");
	signalMcQTotal->SetTitle("nuatm_sep20");
	signalMcQTotal->SetLineColor(kRed);
	signalMcQTotal->Rebin(20);
	signalMcQTotal->Scale(1/(signalMcDataDuration*24*3600));
	TH1F* signalMclikelihood  = (TH1F*)signalMcData->Get("h_likelihood");
	signalMclikelihood->SetTitle("nuatm_sep20");
	signalMclikelihood->SetLineColor(kRed);
	// signalMclikelihood->Rebin(2);
	signalMclikelihood->Scale(1/(signalMcDataDuration*24*3600));
	TH1F* signalMclikelihoodHitOnly  = (TH1F*)signalMcData->Get("h_likelihoodHitOnly");
	signalMclikelihoodHitOnly->SetTitle("nuatm_sep20");
	signalMclikelihoodHitOnly->SetLineColor(kRed);
	// signalMclikelihoodHitOnly->Rebin(2);
	signalMclikelihoodHitOnly->Scale(1/(signalMcDataDuration*24*3600));
	TH1F* signalMcchi2  = (TH1F*)signalMcData->Get("h_chi2");
	signalMcchi2->SetTitle("nuatm_sep20");
	signalMcchi2->SetLineColor(kRed);
	// signalMcchi2->Rebin(2);
	signalMcchi2->Scale(1/(signalMcDataDuration*24*3600));
	TH1F* signalMcZ  = (TH1F*)signalMcData->Get("h_z");
	signalMcZ->SetTitle("nuatm_sep20");
	signalMcZ->SetLineColor(kRed);
	signalMcZ->Rebin(2);
	signalMcZ->Scale(1/(signalMcDataDuration*24*3600));
	TH1F* signalMcHorDist  = (TH1F*)signalMcData->Get("h_horDist");
	signalMcHorDist->SetTitle("nuatm_sep20");
	signalMcHorDist->SetLineColor(kRed);
	signalMcHorDist->Rebin(2);
	signalMcHorDist->Scale(1/(signalMcDataDuration*24*3600));


	TFile* realData = new TFile(dataFile,"READ");
	// TFile* realData = new TFile("../../results/mcResults_data_y18c-1.root","READ");
	// TFile* realData = new TFile("../../results/mcResults_data_y19c-1.root","READ");
	TH1F* realNHits  = (TH1F*)realData->Get("h_nHits");
	realNHits->SetTitle("Data");
	realNHits->SetMarkerStyle(kFullCircle);
	realNHits->Rebin(5);
	realNHits->Scale(1/(realDataDuration*24*3600));
	TH1F* realNHitsTFilter  = (TH1F*)realData->Get("h_nHitsAfterTFilter");
	realNHitsTFilter->SetTitle("Data");
	realNHitsTFilter->SetMarkerStyle(kFullCircle);
	realNHitsTFilter->Rebin(2);
	realNHitsTFilter->Scale(1/(realDataDuration*24*3600));
	TH1F* realNStringsTFilter  = (TH1F*)realData->Get("h_nStringsAfterTFilter");
	realNStringsTFilter->SetTitle("Data");
	realNStringsTFilter->SetMarkerStyle(kFullCircle);
	// realNStringsTFilter->Rebin(2);
	realNStringsTFilter->Scale(1/(realDataDuration*24*3600));
	TH1F* realEnergy  = (TH1F*)realData->Get("h_energy");
	realEnergy->SetTitle("Data");
	realEnergy->SetMarkerStyle(kFullCircle);
	realEnergy->Rebin(5);
	realEnergy->Scale(1/(realDataDuration*24*3600));
	TH1F* realTheta  = (TH1F*)realData->Get("h_theta");
	realTheta->SetTitle("Data");
	realTheta->SetMarkerStyle(kFullCircle);
	realTheta->Rebin(5);
	realTheta->Scale(1/(realDataDuration*24*3600));
	TH1F* realPhi  = (TH1F*)realData->Get("h_phi");
	realPhi->SetTitle("Data");
	realPhi->SetMarkerStyle(kFullCircle);
	realPhi->Rebin(10);
	realPhi->Scale(1/(realDataDuration*24*3600));
	TH1F* realQTotal  = (TH1F*)realData->Get("h_qTotal");
	realQTotal->SetTitle("Data");
	realQTotal->SetMarkerStyle(kFullCircle);
	realQTotal->Rebin(20);
	realQTotal->Scale(1/(realDataDuration*24*3600));
	TH1F* reallikelihood  = (TH1F*)realData->Get("h_likelihood");
	reallikelihood->SetTitle("Data");
	reallikelihood->SetMarkerStyle(kFullCircle);
	// reallikelihood->Rebin(2);
	reallikelihood->Scale(1/(realDataDuration*24*3600));
	TH1F* reallikelihoodHitOnly  = (TH1F*)realData->Get("h_likelihoodHitOnly");
	reallikelihoodHitOnly->SetTitle("Data");
	reallikelihoodHitOnly->SetMarkerStyle(kFullCircle);
	// reallikelihoodHitOnly->Rebin(2);
	reallikelihoodHitOnly->Scale(1/(realDataDuration*24*3600));
	TH1F* realchi2  = (TH1F*)realData->Get("h_chi2");
	realchi2->SetTitle("Data");
	realchi2->SetMarkerStyle(kFullCircle);
	// realchi2->Rebin(2);
	realchi2->Scale(1/(realDataDuration*24*3600));
	TH1F* realZ  = (TH1F*)realData->Get("h_z");
	realZ->SetTitle("Data");
	realZ->SetMarkerStyle(kFullCircle);
	realZ->Rebin(2);
	realZ->Scale(1/(realDataDuration*24*3600));
	TH1F* realHorDist  = (TH1F*)realData->Get("h_horDist");
	realHorDist->SetTitle("Data");
	realHorDist->SetMarkerStyle(kFullCircle);
	realHorDist->Rebin(2);
	realHorDist->Scale(1/(realDataDuration*24*3600));
	TH1F* realNHitsFull  = (TH1F*)realData->Get("h_nHitsFull");
	realNHitsFull->SetTitle("Data");
	realNHitsFull->SetMarkerStyle(kFullCircle);
	realNHitsFull->Rebin(5);
	realNHitsFull->Scale(1/(realDataDuration*24*3600));
	TH1F* realNHitsTFilterFull  = (TH1F*)realData->Get("h_nHitsAfterTFilterFull");
	realNHitsTFilterFull->SetTitle("Data");
	realNHitsTFilterFull->SetMarkerStyle(kFullCircle);
	realNHitsTFilterFull->Rebin(2);
	realNHitsTFilterFull->Scale(1/(realDataDuration*24*3600));
	TH1F* realNStringsTFilterFull  = (TH1F*)realData->Get("h_nStringsAfterTFilterFull");
	realNStringsTFilterFull->SetTitle("Data");
	realNStringsTFilterFull->SetMarkerStyle(kFullCircle);
	// realNStringsTFilterFull->Rebin(2);
	realNStringsTFilterFull->Scale(1/(realDataDuration*24*3600));
	TH1F* realEnergyFull  = (TH1F*)realData->Get("h_energyFull");
	realEnergyFull->SetTitle("Data");
	realEnergyFull->SetMarkerStyle(kFullCircle);
	realEnergyFull->Rebin(5);
	realEnergyFull->Scale(1/(realDataDuration*24*3600));
	TH1F* realThetaFull  = (TH1F*)realData->Get("h_thetaFull");
	realThetaFull->SetTitle("Data");
	realThetaFull->SetMarkerStyle(kFullCircle);
	realThetaFull->Rebin(5);
	realThetaFull->Scale(1/(realDataDuration*24*3600));
	TH1F* realPhiFull  = (TH1F*)realData->Get("h_phiFull");
	realPhiFull->SetTitle("Data");
	realPhiFull->SetMarkerStyle(kFullCircle);
	realPhiFull->Rebin(10);
	realPhiFull->Scale(1/(realDataDuration*24*3600));
	TH1F* realQTotalFull  = (TH1F*)realData->Get("h_qTotalFull");
	realQTotalFull->SetTitle("Data");
	realQTotalFull->SetMarkerStyle(kFullCircle);
	realQTotalFull->Rebin(20);
	realQTotalFull->Scale(1/(realDataDuration*24*3600));
	TH1F* reallikelihoodFull  = (TH1F*)realData->Get("h_likelihoodFull");
	reallikelihoodFull->SetTitle("Data");
	reallikelihoodFull->SetMarkerStyle(kFullCircle);
	// reallikelihoodFull->Rebin(2);
	reallikelihoodFull->Scale(1/(realDataDuration*24*3600));

	cout << "Background expected: " << mcTheta->Integral()*realDataDuration*24*3600 << endl;
	cout << "New Background expected: " << newMcTheta->Integral()*realDataDuration*24*3600 << endl;
	cout << "Signal expected: " << signalMcTheta->Integral()*realDataDuration*24*3600 << endl;
	cout << "Data: " << realTheta->Integral()*realDataDuration*24*3600 << endl;

	THStack* s_nHits = new THStack("s_nHits","; N_{hits} [#];dN/dN_{hits} [Hz / 5 hits]");
	if (showTwoMC)
		s_nHits->Add(mcNHits,"HIST");
	if (showSignalMC)
		s_nHits->Add(signalMcNHits,"HIST");
	s_nHits->Add(newMcNHits,"HIST");
	s_nHits->Add(realNHits,"");

	TCanvas* c_nHits = new TCanvas("c_nHits","NHits",800,600);
	gPad->SetGrid();
	s_nHits->Draw("nostack");
	gPad->BuildLegend(0.75,0.75,0.95,0.95,"");

	THStack* s_nHitsTFilter = new THStack("s_nHitsTFilter","; N_{hits}^{reco} [#];dN/dN_{hits}^{reco} [Hz / 2 hits]");
	if (showTwoMC)
		s_nHitsTFilter->Add(mcNHitsTFilter,"HIST");
	if (showSignalMC)
		s_nHitsTFilter->Add(signalMcNHitsTFilter,"HIST");
	s_nHitsTFilter->Add(newMcNHitsTFilter,"HIST");
	s_nHitsTFilter->Add(realNHitsTFilter,"");

	TCanvas* c_nHitsTFilter = new TCanvas("c_nHitsTFilter","NHitsTFilter",800,600);
	gPad->SetGrid();
	s_nHitsTFilter->Draw("nostack");
	gPad->BuildLegend(0.75,0.75,0.95,0.95,"");

	THStack* s_nStringsTFilter = new THStack("s_nStringsTFilter","; N_{strings}^{reco} [#];dN/dN_{strings}^{reco} [Hz]");
	if (showTwoMC)
		s_nStringsTFilter->Add(mcNStringsTFilter,"HIST");
	if (showSignalMC)
		s_nStringsTFilter->Add(signalMcNStringsTFilter,"HIST");
	s_nStringsTFilter->Add(newMcNStringsTFilter,"HIST");
	s_nStringsTFilter->Add(realNStringsTFilter,"");

	TCanvas* c_nStringsTFilter = new TCanvas("c_nStringsTFilter","NStringsTFilter",800,600);
	gPad->SetGrid();
	s_nStringsTFilter->Draw("nostack");
	gPad->BuildLegend(0.75,0.75,0.95,0.95,"");

	THStack* s_energy = new THStack("s_energy",";E_{rec} [TeV];dN/dE_{rec} [Hz / 5 TeV]");
	if (showTwoMC)
		s_energy->Add(mcEnergy,"HIST");
	if (showSignalMC)
		s_energy->Add(signalMcEnergy,"HIST");
	s_energy->Add(newMcEnergy,"HIST");
	s_energy->Add(realEnergy,"");

	TCanvas* c_energy = new TCanvas("c_energy","Energy",800,600);
	gPad->SetGrid();
	s_energy->Draw("nostack");
	gPad->BuildLegend(0.75,0.75,0.95,0.95,"");

	THStack* s_theta = new THStack("s_theta",";#theta_{rec} [deg.];dN/d#theta_{rec} [Hz / 5 deg]");
	if (showTwoMC)
		s_theta->Add(mcTheta,"HIST");
	if (showSignalMC)
		s_theta->Add(signalMcTheta,"HIST");
	s_theta->Add(newMcTheta,"HIST");
	s_theta->Add(realTheta,"");

	TCanvas* c_theta = new TCanvas("c_theta","Theta",800,600);
	gPad->SetGrid();
	s_theta->Draw("nostack");
	gPad->BuildLegend(0.75,0.75,0.95,0.95,"");

	TH1F* newMcThetaScaled  = (TH1F*)newMcTheta->Clone("h_thetaScaled");
	newMcThetaScaled->Scale(realTheta->Integral()/newMcThetaScaled->Integral());
	newMcThetaScaled->SetLineStyle(5);
	newMcThetaScaled->SetLineWidth(2);
	TH1F* mcThetaScaled  = (TH1F*)mcTheta->Clone("h_thetaScaled");
	mcThetaScaled->Scale(realTheta->Integral()/mcThetaScaled->Integral());
	mcThetaScaled->SetLineStyle(5);
	mcThetaScaled->SetLineWidth(2);

	THStack* s_thetaScaled = new THStack("s_thetaScaled",";#theta_{rec} [deg.];dN/d#theta_{rec} [Hz / 5 deg]");
	if (showTwoMC)
		s_thetaScaled->Add(mcThetaScaled,"HIST");
	s_thetaScaled->Add(newMcThetaScaled,"HIST");
	s_thetaScaled->Add(realTheta,"");

	TCanvas* c_thetaScaled = new TCanvas("c_thetaScaled","ThetaScaled",800,600);
	gPad->SetGrid();
	s_thetaScaled->Draw("nostack");
	gPad->BuildLegend(0.75,0.75,0.95,0.95,"");

	THStack* s_phi = new THStack("s_phi",";#phi_{rec} [deg.];dN/d#phi_{rec} [Hz / 10 deg]");
	if (showTwoMC)
		s_phi->Add(mcPhi,"HIST");
	if (showSignalMC)
		s_phi->Add(signalMcPhi,"HIST");
	s_phi->Add(newMcPhi,"HIST");
	s_phi->Add(realPhi,"");

	TCanvas* c_phi = new TCanvas("c_phi","Phi",800,600);
	gPad->SetGrid();
	s_phi->Draw("nostack");
	gPad->BuildLegend(0.75,0.75,0.95,0.95,"");

	THStack* s_qTotal = new THStack("s_qTotal",";Q [p.e.];dN/dQ [Hz / 200 p.e.]");
	if (showTwoMC)
		s_qTotal->Add(mcQTotal,"HIST");
	if (showSignalMC)
		s_qTotal->Add(signalMcQTotal,"HIST");
	s_qTotal->Add(newMcQTotal,"HIST");
	s_qTotal->Add(realQTotal,"");

	TCanvas* c_qTotal = new TCanvas("c_qTotal","QTotal",800,600);
	gPad->SetGrid();
	s_qTotal->Draw("nostack");
	gPad->BuildLegend(0.75,0.75,0.95,0.95,"");

	THStack* s_likelihood = new THStack("s_likelihood",";L [#];dN/dL [Hz / 0.1]");
	if (showTwoMC)
		s_likelihood->Add(mclikelihood,"HIST");
	if (showSignalMC)
		s_likelihood->Add(signalMclikelihood,"HIST");
	s_likelihood->Add(newMclikelihood,"HIST");
	s_likelihood->Add(reallikelihood,"");

	TCanvas* c_likelihood = new TCanvas("c_likelihood","Likelihood",800,600);
	gPad->SetGrid();
	s_likelihood->Draw("nostack");
	gPad->BuildLegend(0.75,0.75,0.95,0.95,"");

	THStack* s_likelihoodHitOnly = new THStack("s_likelihoodHitOnly",";L [#];dN/dL [Hz / 0.01]");
	// if (showTwoMC)
		// s_likelihoodHitOnly->Add(mclikelihoodHitOnly,"HIST");
	if (showSignalMC)
		s_likelihoodHitOnly->Add(signalMclikelihoodHitOnly,"HIST");
	s_likelihoodHitOnly->Add(newMclikelihoodHitOnly,"HIST");
	s_likelihoodHitOnly->Add(reallikelihoodHitOnly,"");

	TCanvas* c_likelihoodHitOnly = new TCanvas("c_likelihoodHitOnly","LikelihoodHitOnly",800,600);
	gPad->SetGrid();
	s_likelihoodHitOnly->Draw("nostack");
	gPad->BuildLegend(0.75,0.75,0.95,0.95,"");

	THStack* s_chi2 = new THStack("s_chi2",";#chi^{2} [1];dN/d#chi^{2} [Hz / 1]");
	// if (showTwoMC)
		// s_chi2->Add(mclikelihoodHitOnly,"HIST");
	if (showSignalMC)
		s_chi2->Add(signalMcchi2,"HIST");
	s_chi2->Add(newMcchi2,"HIST");
	s_chi2->Add(realchi2,"");

	TCanvas* c_chi2 = new TCanvas("c_chi2","Chi2",800,600);
	gPad->SetGrid();
	s_chi2->Draw("nostack");
	gPad->BuildLegend(0.75,0.75,0.95,0.95,"");

	THStack* s_z = new THStack("s_z",";Z [m];dN/dZ [Hz / 20 m]");
	if (showSignalMC)
		s_z->Add(signalMcZ,"HIST");
	s_z->Add(newMcZ,"HIST");
	s_z->Add(realZ,"");

	TCanvas* c_z = new TCanvas("c_z","Z",800,600);
	gPad->SetGrid();
	s_z->Draw("nostack");
	gPad->BuildLegend(0.75,0.75,0.95,0.95,"");

	THStack* s_horDist = new THStack("s_horDist",";#rho [m];dN/d#rho [Hz / 2 m]");
	if (showSignalMC)
		s_horDist->Add(signalMcHorDist,"HIST");
	s_horDist->Add(newMcHorDist,"HIST");
	s_horDist->Add(realHorDist,"");

	TCanvas* c_horDist = new TCanvas("c_horDist","HorizontalDistance",800,600);
	gPad->SetGrid();
	s_horDist->Draw("nostack");
	gPad->BuildLegend(0.75,0.75,0.95,0.95,"");

	if (showFull)
	{
		THStack* s_nHitsFull = new THStack("s_nHitsFull","; N_{hits} [#];dN/dN_{hits} [#]");
		s_nHitsFull->Add(mcNHitsFull,"HIST");
		s_nHitsFull->Add(realNHitsFull,"");

		TCanvas* c_nHitsFull = new TCanvas("c_nHitsFull","NHitsFull",800,600);
		gPad->SetGrid();
		s_nHitsFull->Draw("nostack");
		gPad->BuildLegend(0.75,0.75,0.95,0.95,"");

		THStack* s_nHitsTFilterFull = new THStack("s_nHitsTFilterFull","; N_{hits}^{reco} [#];dN/dN_{hits}^{reco} [#]");
		s_nHitsTFilterFull->Add(mcNHitsTFilterFull,"HIST");
		s_nHitsTFilterFull->Add(realNHitsTFilterFull,"");

		TCanvas* c_nHitsTFilterFull = new TCanvas("c_nHitsTFilterFull","NHitsTFilterFull",800,600);
		gPad->SetGrid();
		s_nHitsTFilterFull->Draw("nostack");
		gPad->BuildLegend(0.75,0.75,0.95,0.95,"");

		THStack* s_nStringsTFilterFull = new THStack("s_nStringsTFilterFull","; N_{strings}^{reco} [#];dN/dN_{strings}^{reco} [#]");
		s_nStringsTFilterFull->Add(mcNStringsTFilterFull,"HIST");
		s_nStringsTFilterFull->Add(realNStringsTFilterFull,"");

		TCanvas* c_nStringsTFilterFull = new TCanvas("c_nStringsTFilterFull","NStringsTFilterFull",800,600);
		gPad->SetGrid();
		s_nStringsTFilterFull->Draw("nostack");
		gPad->BuildLegend(0.75,0.75,0.95,0.95,"");

		THStack* s_energyFull = new THStack("s_energyFull",";E_{rec} [TeV];dN/dE_{rec} [TeV]");
		s_energyFull->Add(mcEnergyFull,"HIST");
		s_energyFull->Add(realEnergyFull,"");

		TCanvas* c_energyFull = new TCanvas("c_energyFull","EnergyFull",800,600);
		gPad->SetGrid();
		s_energyFull->Draw("nostack");
		gPad->BuildLegend(0.75,0.75,0.95,0.95,"");

		THStack* s_thetaFull = new THStack("s_thetaFull",";#theta_{rec} [deg.];dN/d#theta_{rec} [Hz / 5 deg]");
		s_thetaFull->Add(mcThetaFull,"HIST");
		s_thetaFull->Add(realThetaFull,"");

		TCanvas* c_thetaFull = new TCanvas("c_thetaFull","ThetaFull",800,600);
		gPad->SetGrid();
		s_thetaFull->Draw("nostack");
		gPad->BuildLegend(0.75,0.75,0.95,0.95,"");

		THStack* s_phiFul = new THStack("s_phiFul",";#phi_{rec} [deg.];dN/d#phi_{rec} [Hz / 5 deg]");
		s_phiFul->Add(mcPhiFull,"HIST");
		s_phiFul->Add(realPhiFull,"");

		TCanvas* c_phiFull = new TCanvas("c_phiFull","PhiFull",800,600);
		gPad->SetGrid();
		s_phiFul->Draw("nostack");
		gPad->BuildLegend(0.75,0.75,0.95,0.95,"");

		THStack* s_qTotalFull = new THStack("s_qTotalFull",";Q [p.e.];dN/dQ [Hz]");
		s_qTotalFull->Add(mcQTotalFull,"HIST");
		s_qTotalFull->Add(realQTotalFull,"");

		TCanvas* c_qTotalFull = new TCanvas("c_qTotalFull","QTotalFull",800,600);
		gPad->SetGrid();
		s_qTotalFull->Draw("nostack");
		gPad->BuildLegend(0.75,0.75,0.95,0.95,"");

		THStack* s_likelihoodFull = new THStack("s_likelihoodFull",";L [#];dN/dL [Hz]");
		s_likelihoodFull->Add(mclikelihoodFull,"HIST");
		s_likelihoodFull->Add(reallikelihoodFull,"");

		TCanvas* c_likelihoodFull = new TCanvas("c_likelihoodFull","LikelihoodFull",800,600);
		gPad->SetGrid();
		s_likelihoodFull->Draw("nostack");
		gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
	}
	return 0;
}
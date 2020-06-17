#include "TFile.h"
#include "TH1F.h"
#include "THStack.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TPaveStats.h"
#include "TLegend.h"

int compareDataMCCascades(double realDataDuration = 237.632, double mcDataDuration = 76.78)
{
	bool showFull = false;

	TFile* mcData = new TFile("../../results/mcResults_muatm_may19.root","READ");
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
	mclikelihood->Rebin(2);
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
	mclikelihoodFull->Rebin(2);
	mclikelihoodFull->Scale(1/(mcDataDuration*24*3600));



	TFile* realData = new TFile("../../results/mcResults_data_y16c0.root","READ");
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
	reallikelihood->Rebin(2);
	reallikelihood->Scale(1/(realDataDuration*24*3600));

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
	reallikelihoodFull->Rebin(2);
	reallikelihoodFull->Scale(1/(realDataDuration*24*3600));

	cout << "Background expected: " << mcNHits->Integral()*realDataDuration*24*3600 << endl;
	cout << "Data: " << realNHits->Integral()*realDataDuration*24*3600 << endl;

	THStack* s_nHits = new THStack("s_nHits","; N_{hits} [#];dN/dN_{hits} [#]");
	s_nHits->Add(mcNHits,"HIST");
	s_nHits->Add(realNHits,"");

	TCanvas* c_nHits = new TCanvas("c_nHits","NHits",800,600);
	gPad->SetGrid();
	s_nHits->Draw("nostack");
	gPad->BuildLegend(0.75,0.75,0.95,0.95,"");

	THStack* s_nHitsTFilter = new THStack("s_nHitsTFilter","; N_{hits}^{reco} [#];dN/dN_{hits}^{reco} [#]");
	s_nHitsTFilter->Add(mcNHitsTFilter,"HIST");
	s_nHitsTFilter->Add(realNHitsTFilter,"");

	TCanvas* c_nHitsTFilter = new TCanvas("c_nHitsTFilter","NHitsTFilter",800,600);
	gPad->SetGrid();
	s_nHitsTFilter->Draw("nostack");
	gPad->BuildLegend(0.75,0.75,0.95,0.95,"");

	THStack* s_nStringsTFilter = new THStack("s_nStringsTFilter","; N_{strings}^{reco} [#];dN/dN_{strings}^{reco} [#]");
	s_nStringsTFilter->Add(mcNStringsTFilter,"HIST");
	s_nStringsTFilter->Add(realNStringsTFilter,"");

	TCanvas* c_nStringsTFilter = new TCanvas("c_nStringsTFilter","NStringsTFilter",800,600);
	gPad->SetGrid();
	s_nStringsTFilter->Draw("nostack");
	gPad->BuildLegend(0.75,0.75,0.95,0.95,"");

	THStack* s_energy = new THStack("s_energy",";E_{rec} [TeV];dN/dE_{rec} [TeV]");
	s_energy->Add(mcEnergy,"HIST");
	s_energy->Add(realEnergy,"");

	TCanvas* c_energy = new TCanvas("c_energy","Energy",800,600);
	gPad->SetGrid();
	s_energy->Draw("nostack");
	gPad->BuildLegend(0.75,0.75,0.95,0.95,"");

	THStack* s_theta = new THStack("s_theta",";#theta_{rec} [deg.];dN/d#theta_{rec} [Hz / 5 deg]");
	s_theta->Add(mcTheta,"HIST");
	s_theta->Add(realTheta,"");

	TCanvas* c_theta = new TCanvas("c_theta","Theta",800,600);
	gPad->SetGrid();
	s_theta->Draw("nostack");
	gPad->BuildLegend(0.75,0.75,0.95,0.95,"");

	THStack* s_phi = new THStack("s_phi",";#phi_{rec} [deg.];dN/d#phi_{rec} [Hz / 5 deg]");
	s_phi->Add(mcPhi,"HIST");
	s_phi->Add(realPhi,"");

	TCanvas* c_phi = new TCanvas("c_phi","Phi",800,600);
	gPad->SetGrid();
	s_phi->Draw("nostack");
	gPad->BuildLegend(0.75,0.75,0.95,0.95,"");

	THStack* s_qTotal = new THStack("s_qTotal",";Q [p.e.];dN/dQ [Hz]");
	s_qTotal->Add(mcQTotal,"HIST");
	s_qTotal->Add(realQTotal,"");

	TCanvas* c_qTotal = new TCanvas("c_qTotal","QTotal",800,600);
	gPad->SetGrid();
	s_qTotal->Draw("nostack");
	gPad->BuildLegend(0.75,0.75,0.95,0.95,"");

	THStack* s_likelihood = new THStack("s_likelihood",";L [#];dN/dL [Hz]");
	s_likelihood->Add(mclikelihood,"HIST");
	s_likelihood->Add(reallikelihood,"");

	TCanvas* c_likelihood = new TCanvas("c_likelihood","Likelihood",800,600);
	gPad->SetGrid();
	s_likelihood->Draw("nostack");
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

		TCanvas* c_likelihoodFull = new TCanvas("c_likelihoodFull","Likelihood",800,600);
		gPad->SetGrid();
		s_likelihoodFull->Draw("nostack");
		gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
	}
	return 0;
}
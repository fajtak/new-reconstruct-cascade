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

int compareLikelihoods()
{
	TString combiLike20 = "../../results/mcResults_combiLikelihood_20steps.root";
	TString combiLike36 = "../../results/mcResults_combiLikelihood_36steps.root";
	TString multiLike = "../../results/mcResults_multiLikelihood.root";

	TFile* combi20 = new TFile(combiLike20,"READ");
	TH1F* combi20MismatchAngle  = (TH1F*)combi20->Get("h_mismatchAngle");
	combi20MismatchAngle->SetTitle("Combined likelihood - 20 steps");
	combi20MismatchAngle->SetLineColor(kRed);
	// combiMismatchAngle->Rebin(5);
	TH1F* combi20MismatchTheta  = (TH1F*)combi20->Get("h_mismatchTheta");
	combi20MismatchTheta->SetTitle("Combined likelihood - 20 steps");
	combi20MismatchTheta->SetLineColor(kRed);
	TH1F* combi20Likelihood  = (TH1F*)combi20->Get("h_likelihood");
	combi20Likelihood->SetTitle("Combined likelihood - 20 steps");
	combi20Likelihood->SetLineColor(kRed);

	TFile* combi36 = new TFile(combiLike36,"READ");
	TH1F* combi36MismatchAngle  = (TH1F*)combi36->Get("h_mismatchAngle");
	combi36MismatchAngle->SetTitle("Combined likelihood - 36 steps");
	combi36MismatchAngle->SetLineColor(kGreen);
	// combiMismatchAngle->Rebin(5);
	TH1F* combi36MismatchTheta  = (TH1F*)combi36->Get("h_mismatchTheta");
	combi36MismatchTheta->SetTitle("Combined likelihood - 36 steps");
	combi36MismatchTheta->SetLineColor(kGreen);
	TH1F* combi36Likelihood  = (TH1F*)combi36->Get("h_likelihood");
	combi36Likelihood->SetTitle("Combined likelihood - 36 steps");
	combi36Likelihood->SetLineColor(kGreen);

	TFile* multi = new TFile(multiLike,"READ");
	TH1F* multiMismatchAngle  = (TH1F*)multi->Get("h_mismatchAngle");
	multiMismatchAngle->SetTitle("Multi likelihood");
	multiMismatchAngle->SetLineColor(kBlue);

	TH1F* multiMismatchTheta  = (TH1F*)multi->Get("h_mismatchTheta");
	multiMismatchTheta->SetTitle("Multi likelihood");
	multiMismatchTheta->SetLineColor(kBlue);

	TH1F* multiLikelihood  = (TH1F*)multi->Get("h_likelihood");
	multiLikelihood->SetTitle("Multi likelihood");
	multiLikelihood->SetLineColor(kBlue);

	THStack* s_mismatchAngle = new THStack("s_mismatchAngle","; Mismatch angle [deg];NoE [#]");
	s_mismatchAngle->Add(combi20MismatchAngle,"HIST");
	s_mismatchAngle->Add(combi36MismatchAngle,"HIST");
	s_mismatchAngle->Add(multiMismatchAngle,"HIST");

	TCanvas* c_mismatchAngle = new TCanvas("c_mismatchAngle","MismatchAngle",800,600);
	gPad->SetGrid();
	s_mismatchAngle->Draw("nostack");
	gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
	// s_mismatchAngle->Write();

	THStack* s_mismatchTheta = new THStack("s_mismatchTheta","; Mismatch theta [deg];NoE [#]");
	s_mismatchTheta->Add(combi20MismatchTheta,"HIST");
	s_mismatchTheta->Add(combi36MismatchTheta,"HIST");
	s_mismatchTheta->Add(multiMismatchTheta,"HIST");

	TCanvas* c_mismatchTheta = new TCanvas("c_mismatchTheta","MismatchTheta",800,600);
	gPad->SetGrid();
	s_mismatchTheta->Draw("nostack");
	gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
	// s_mismatchTheta->Write();

	THStack* s_likelihood = new THStack("s_likelihood","; Likelihood [1];NoE [#]");
	s_likelihood->Add(combi20Likelihood,"HIST");
	s_likelihood->Add(combi36Likelihood,"HIST");
	s_likelihood->Add(multiLikelihood,"HIST");

	TCanvas* c_likelihood = new TCanvas("c_likelihood","Likelihood",800,600);
	gPad->SetGrid();
	s_likelihood->Draw("nostack");
	gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
	// s_likelihood->Write();

	return 0;
}
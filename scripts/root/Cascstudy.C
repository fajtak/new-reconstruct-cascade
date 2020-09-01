#include "TString.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TMath.h"
#include "TChain.h"
#include "THStack.h"

#include<iostream>

TH1F* h_cost = new TH1F("h_cost","cos(Theta); cos(Theta); NoE [#]",200,-1,1);
TH1F* h_astroCascEnergy = new TH1F("h_astroCascEnergy","Astrophysical neutrinos; E [TeV]; NoE [#]",200,0,1000);
TH1F* h_atmCascEnergy = new TH1F("h_atmCascEnergy","Atmospheric neutrinos; E [TeV]; NoE [#]",200,0,1000);
TH1F* h_zenithAngle = new TH1F("h_zenithAngle","Zenith Angle; #theta [deg]; NoE [#]",36,0,180);
TH1F* h_energy = new TH1F("h_energy","",500,0,5);


// structure necessary for the reading of the Zhan-Arys cascade simulations
struct mcCascade
{
	UShort_t nHitct;
	float cosTheta, position[3], showerEnergy, neutrinoEnergy, pweight, probint, s0, sp246, spat;
};

void SetBranches(mcCascade* cascade, TChain* mcFiles)
{
	mcFiles->SetBranchAddress("cost",&cascade->cosTheta);
	mcFiles->SetBranchAddress("xtr",cascade->position);
	mcFiles->SetBranchAddress("Esh",&cascade->showerEnergy);
	mcFiles->SetBranchAddress("Ein",&cascade->neutrinoEnergy);
	mcFiles->SetBranchAddress("pwait",&cascade->pweight);
	mcFiles->SetBranchAddress("probint",&cascade->probint);
	// mcFiles->SetBranchAddress("spwe2",&cascade->showerEnergy);
	mcFiles->SetBranchAddress("s0",&cascade->s0);
	mcFiles->SetBranchAddress("spat",&cascade->spat);
	// mcFiles->SetBranchAddress("sp23",&cascade->showerEnergy);
	mcFiles->SetBranchAddress("sp246",&cascade->sp246);
	// mcFiles->SetBranchAddress("twait",&cascade->showerEnergy);
	mcFiles->SetBranchAddress("Nhitct",&cascade->nHitct);
	// mcFiles->SetBranchAddress("ocht",&cascade->showerEnergy);
	// mcFiles->SetBranchAddress("oxt",&cascade->showerEnergy);
	// mcFiles->SetBranchAddress("ocha",&cascade->showerEnergy);
	// mcFiles->SetBranchAddress("oxa",&cascade->showerEnergy);
	// mcFiles->SetBranchAddress("ohnh",&cascade->showerEnergy);
	// mcFiles->SetBranchAddress("ohma",&cascade->showerEnergy);
	// mcFiles->SetBranchAddress("odang",&cascade->showerEnergy);
}

void DrawResults()
{
	TCanvas* c_cost = new TCanvas("c_cost","CosTheta",800,600);
	h_cost->Draw();

	TCanvas* c_CascEnergy = new TCanvas("c_CascEnergy","Cascade Energy",800,600);
	THStack* hs_cascadeEnergy = new THStack("hs_cascadeEnergy","Cascade Energy; E [TeV];dN/dE [Hz / 5 TeV]");
	hs_cascadeEnergy->Add(h_astroCascEnergy);
	hs_cascadeEnergy->Add(h_atmCascEnergy);
	h_atmCascEnergy->SetLineColor(kRed);
	hs_cascadeEnergy->Draw("nostack");
	gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
	// h_astroCascEnergy->Draw();
	// h_atmCascEnergy->Draw();

	TCanvas* c_zenithAngle = new TCanvas("c_zenithAngle","Zenith Angle",800,600);
	h_zenithAngle->Draw();

	TCanvas* c_energy = new TCanvas("c_energy","Energy",800,600);
	h_energy->Draw();

}

int Cascstudy()
{
	TString inputFolder = "/Data/BaikalData/mc/newCascades/nu_e_MC/";
	// TString inputFolder = "/Data/BaikalData/mc/newCascades/nu_tau_MC/";
	TString filePath = Form("%s/a*_hitne16_c*.root",inputFolder.Data());
	// TString filePath = Form("%s/a*_hitntau16_c*.root",inputFolder.Data());
	TChain* mcFiles = new TChain("h11");
	mcFiles->Add(filePath.Data());
	if (mcFiles->GetEntries() == 0)
	{
		std::cout << "Files: " << filePath << " were not found!" << endl;
    	return -2;
	}

	// Sets necessary pointers to access data through TChain
	mcCascade* cascade = new mcCascade;
	SetBranches(cascade,mcFiles);

	// double time = 1;
	double time = 3.15e7/365*230;
	int hitCut = 20;
	// double ICflux = 4.1e-9;
	double ICFlux = 1.7e-10;
	double electronAtmFlux = 1.1e-8;
	double muonAtmFlux = 1.2e-7;

	for (int i = 0; i < mcFiles->GetEntries(); ++i)
	{
		mcFiles->GetEntry(i);
		h_cost->Fill(cascade->cosTheta);
		h_zenithAngle->Fill(TMath::ACos(cascade->cosTheta)/TMath::Pi()*180,ICFlux*time*cascade->s0*(4*3.14/20/20)/10000*cascade->pweight*cascade->probint*cascade->sp246*(cascade->nHitct>hitCut));

		h_astroCascEnergy->Fill(cascade->showerEnergy,ICFlux*time*cascade->s0*(4*3.14/20/20)/10000*cascade->pweight*cascade->probint*cascade->sp246*(cascade->nHitct>hitCut));
		h_atmCascEnergy->Fill(cascade->showerEnergy,muonAtmFlux*time*cascade->s0*(4*3.14/20/20)/10000*cascade->pweight*cascade->probint*cascade->spat*(cascade->nHitct>hitCut));
		h_energy->Fill(TMath::Log10(cascade->showerEnergy));
	}

	DrawResults();
	return 0;
}

// root> h23.Draw("log10(Ein)","1.7e-10*3.15e7*s0*(4*3.14/20/20)/10000*pwait*probint*sp246*(Nhitct>19&...)")
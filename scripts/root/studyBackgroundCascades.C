#include "TChain.h"
#include "TString.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TMath.h"
#include "BEvent.h"
#include "BMCEvent.h"

#include <iostream>

TH1F* h_nMuons = new TH1F("h_nMuons","Number of muons; N_{muons} [#]; NoE [#]", 20,0,20);
TH1F* h_nCascades = new TH1F("h_nCascades","Number of cascades; N_{cascades} [#]; NoE [#]", 20,0,20);
TH1F* h_cascadeEnergies = new TH1F("h_cascadeEnergies","Cascade energies; E [TeV]; NoE [#]", 100,0,100);
TH1F* h_highestCascadeEnergies = new TH1F("h_highestCascadeEnergies","Highest Cascade energies; E [TeV]; NoE [#]", 100,0,100);
TH2F* h_cascadeEnergiesDistance = new TH2F("h_cascadeEnergiesDistance","Cascade energies vs. Distance; E [TeV]; Distance [m]", 100,0,100,100,0,100);
TH2F* h_horizontalPos = new TH2F("h_horizontalPos","Horizontal position; X [m]; Y [m]", 300 ,-150,150,300,-150,150);
TH1F* h_verticalPos = new TH1F("h_verticalPos","Vertical position; Z [m]; NoE [#]", 600 ,-300,300);

int DrawResults()
{
	TCanvas* c_nMuons = new TCanvas("c_nMuons","NMuons",800,600);
	h_nMuons->Draw();

	TCanvas* c_nCascades = new TCanvas("c_nCascades","NCascades",800,600);
	h_nCascades->Draw();

	TCanvas* c_cascadeEnergies = new TCanvas("c_cascadeEnergies","CascadeEnergies",800,600);
	h_cascadeEnergies->Draw();
	h_highestCascadeEnergies->Draw("same");
	h_highestCascadeEnergies->SetLineColor(kRed);

	TCanvas* c_cascadeEnergiesDistance = new TCanvas("c_cascadeEnergiesDistance","CascadeEnergiesDistance",800,600);
	h_cascadeEnergiesDistance->Draw("colz");

	TCanvas* c_horizontalPos = new TCanvas("c_horizontalPos","HorizontalPosition",800,600);
	h_horizontalPos->Draw("colz");

	TCanvas* c_verticalPos = new TCanvas("c_verticalPos","VerticalPos",800,600);
	h_verticalPos->Draw();

	return 0;
}

double CalculateDistance(BMCInteraction* interaction)
{
	// return TMath::Sqrt(TMath::Power(interaction->GetX(),2)+TMath::Power(interaction->GetY(),2)+TMath::Power(interaction->GetZ(),2));
	return TMath::Sqrt(TMath::Power(interaction->GetX(),2)+TMath::Power(interaction->GetY(),2));
}

int studyBackgroundCascades()
{
	// TString gFileInputFolder = "/Data/BaikalData/mc/nuatm_jun20/";
	TString gFileInputFolder = "/Data/BaikalData/mc/nuatm_sep20_root/";

	// const char* filePath = Form("%s/n_nuatm_gs_n2m_cl2016_x*.root",gFileInputFolder.Data());
	const char* filePath = Form("%s/cluster0/n_nuatm_gs_n2m_cl2016_x*.root",gFileInputFolder.Data());
	// const char* filePath = Form("%s/n_nuatm_gs_n2m_cl2016_x1001.root",gFileInputFolder.Data());

	TChain* mcFiles = new TChain("Events");
	mcFiles->Add(filePath);
	if (mcFiles->GetEntries() == 0)
	{
		std::cout << "Files: " << filePath << " were not found!" << endl;
    	return -2;
	}

	// Sets necessary pointers to access data through TChain
	BEvent* event = NULL;
    mcFiles->SetBranchAddress("BEvent.",&event);
    BEventMaskMC* eventMask = NULL;
    mcFiles->SetBranchAddress("MCEventMask.",&eventMask);
    // BSourceEAS* sourceEAS = NULL;
    // mcFiles->SetBranchAddress("MCEventSource.",&sourceEAS);
    BMCEvent* mcEvent = NULL;
    mcFiles->SetBranchAddress("BMCEvent.",&mcEvent);

    for (int i = 0; i < mcFiles->GetEntries(); ++i)
    {
    	mcFiles->GetEntry(i);

    	double highestEnergy = 0;
    	int muonID = -1;
    	int cascadeID = -1;

    	h_nMuons->Fill(mcEvent->GetResponseMuonsN());

		for (int k = 0; k < mcEvent->GetResponseMuonsN(); ++k)
		{
			h_nCascades->Fill(mcEvent->GetTrack(k)->GetInteractionN());
			// cout << mcEvent->GetTrack(k)->GetMuonEnergy() << endl;
			for (int j = 0; j < mcEvent->GetTrack(k)->GetInteractionN(); ++j)
			{
				h_cascadeEnergies->Fill(mcEvent->GetTrack(k)->GetInteraction(j)->GetEnergy());
				h_cascadeEnergiesDistance->Fill(mcEvent->GetTrack(k)->GetInteraction(j)->GetEnergy(),CalculateDistance(mcEvent->GetTrack(k)->GetInteraction(j)));
				// cout << mcEvent->GetTrack(k)->GetInteraction(i)->GetEnergy() << endl;
				if (mcEvent->GetTrack(k)->GetInteraction(j)->GetEnergy() > highestEnergy)
				{
					highestEnergy = mcEvent->GetTrack(k)->GetInteraction(j)->GetEnergy();
					muonID = k;
					cascadeID = j;
				}
			};
    	}
    	if (muonID != -1 && cascadeID != -1)
		{
			h_highestCascadeEnergies->Fill(mcEvent->GetTrack(muonID)->GetInteraction(cascadeID)->GetEnergy());
			h_verticalPos->Fill(mcEvent->GetTrack(muonID)->GetInteraction(cascadeID)->GetZ());
			h_horizontalPos->Fill(mcEvent->GetTrack(muonID)->GetInteraction(cascadeID)->GetX(),mcEvent->GetTrack(muonID)->GetInteraction(cascadeID)->GetY());
			if (mcEvent->GetTrack(muonID)->GetInteraction(cascadeID)->GetEnergy() > 20)
				cout << i << " " << mcEvent->GetTrack(muonID)->GetInteraction(cascadeID)->GetEnergy() << " " << CalculateDistance(mcEvent->GetTrack(muonID)->GetInteraction(cascadeID)) << endl;
			// unifiedEvent.mcEnergy = mcEvent->GetTrack(muonID)->GetInteraction(cascadeID)->GetEnergy();
			// unifiedEvent.mcPosition.SetX(mcEvent->GetTrack(muonID)->GetInteraction(cascadeID)->GetX());
			// unifiedEvent.mcPosition.SetY(mcEvent->GetTrack(muonID)->GetInteraction(cascadeID)->GetY());
			// unifiedEvent.mcPosition.SetZ(mcEvent->GetTrack(muonID)->GetInteraction(cascadeID)->GetZ());
			// unifiedEvent.mcFlagID = (cascadeID+1)*1000+(muonID+1);
			// cout << muonID << " " << cascadeID << " " << mcEvent->GetTrack(muonID)->GetInteraction(cascadeID)->GetEnergy() << endl;
		}
	}

	DrawResults();
	return 0;
}



	// // cout << "New Event" << endl;
	// for (int k = 0; k < mcEvent->GetResponseMuonsN(); ++k)
	// {
	// 	// cout << mcEvent->GetTrack(k)->GetMuonEnergy() << endl;
	// 	for (int i = 0; i < mcEvent->GetTrack(k)->GetInteractionN(); ++i)
	// 	{
	// 		// cout << mcEvent->GetTrack(k)->GetInteraction(i)->GetEnergy() << endl;
	// 		if (mcEvent->GetTrack(k)->GetInteraction(i)->GetEnergy() > highestEnergy)
	// 		{
	// 			highestEnergy = mcEvent->GetTrack(k)->GetInteraction(i)->GetEnergy();
	// 			muonID = k;
	// 			cascadeID = i;
	// 		}
	// 	};
	// 	// for (int j = 0; j < mcEvent->GetMuonTrack(k)->GetNumShowers(); ++j)
	// 	// {
	// 	// 	if (mcEvent->GetMuonTrack(k)->GetShower(j)->GetE() > highestEnergy)
	// 	// 	{
	// 	// 		highestEnergy = mcEvent->GetMuonTrack(k)->GetShower(j)->GetE();
	// 	// 		muonID = k;
	// 	// 		cascadeID = j;
	// 	// 	}
	// 	// 	mcEvent->GetMuonTrack(k)->GetShower(j)->GetDirection().Print();
	// 	// 	mcEvent->GetMuonTrack(k)->GetShower(j)->GetPosition().Print();
	// 	// 	cout << k << " " << j << " " << mcEvent->GetMuonTrack(k)->GetShower(j)->GetX() << " " << mcEvent->GetMuonTrack(k)->GetShower(j)->GetY() << " " << mcEvent->GetMuonTrack(k)->GetShower(j)->GetZ() << " " << mcEvent->GetMuonTrack(k)->GetShower(j)->GetE() << endl;
	// 	// }
	// }
	// if (muonID != -1 && cascadeID != -1)
	// {
	// 	unifiedEvent.mcEnergy = mcEvent->GetTrack(muonID)->GetInteraction(cascadeID)->GetEnergy();
	// 	unifiedEvent.mcPosition.SetX(mcEvent->GetTrack(muonID)->GetInteraction(cascadeID)->GetX());
	// 	unifiedEvent.mcPosition.SetY(mcEvent->GetTrack(muonID)->GetInteraction(cascadeID)->GetY());
	// 	unifiedEvent.mcPosition.SetZ(mcEvent->GetTrack(muonID)->GetInteraction(cascadeID)->GetZ());
	// 	unifiedEvent.mcFlagID = (cascadeID+1)*1000+(muonID+1);
	// 	// cout << muonID << " " << cascadeID << " " << unifiedEvent.mcEnergy << " " << unifiedEvent.mcPosition.X() << " " << unifiedEvent.mcPosition.Y() << " " << unifiedEvent.mcPosition.Z() << endl;
	// }
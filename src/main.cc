//program headers
#include "opts.h"

//system headers
#include <iostream>
#include <fstream>

//ROOT dependencies
#include "TTree.h"
#include "TFile.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TH1D.h"
#include "TChain.h"
#include "TGraph2D.h"

//BARS dependencies
#include "BARS.h"
#include "BExtractedImpulseTel.h"
#include "BExtractedHeader.h"
#include "BGeomTel.h"
#include "BDynamicGeometry.h"
#include "BEventMask.h"
#include "BEvent.h"
#include "BMCEvent.h"
#include "BSource.h"
#include "BGeomTelMC.h"


using namespace BARS;
using namespace std;

TH1D* h_nHits = new TH1D("h_nHits","Number of hits per event;N_{hits} [#];NoE [#]",300,0,300);
TH1D* h_totalQ = new TH1D("h_totalQ","Total charge per event;Q [p.e.];NoE [#]",500,0,5000);
TH1D* h_nHitsCaus = new TH1D("h_nHitsCaus","Number of hits per event after causality filter;N_{hits} [#];NoE [#]",300,0,300);
TH1D* h_nStringsCaus = new TH1D("h_nStringsCaus","Number of hit strings per event after causality filter;N_{hits} [#];NoE [#]",10,0,10);
TH1D* h_chi2Caus = new TH1D("h_chi2Caus","#chi^{2} of the position fit after causality filter;#chi^{2} [#];NoE [#]",500,0,500);
TH1D* h_nHitsTFilter = new TH1D("h_nHitsTFilter","Number of hits per event after TFilter;N_{hits} [#];NoE [#]",300,0,300);
TH1D* h_nStringsTFilter = new TH1D("h_nStringsTFilter","Number of hit strings per event after TFilter;N_{hits} [#];NoE [#]",10,0,10);
TH1D* h_chi2TFilter = new TH1D("h_chi2TFilter","#chi^{2} of the position fit after TFilter;#chi^{2} [#];NoE [#]",500,0,500);
TH1D* h_nHitsChange = new TH1D("h_nHitsChange","Change in the number of hits after TFilter;#delta N_{hits} [#];NoE [#]",200,-50,150);
TH1D* h_exitStatus = new TH1D("h_exitStatus","Exit status of the event;#Status [#];NoE [#]",20,-10,10);
TH1D* h_nHitsTrack = new TH1D("h_nHitsTrack","Number of track hits; N_{hits} [#]; NoE [#]",100,0,100);

TH1F* h_likelihood = new TH1F("h_likelihood","Likelihood value; #mathcal{L}; NoE [#]",100,0,100);
TH1F* h_mcEnergy = new TH1F("h_mcEnergy","MC energy; E_{mc} [TeV]; NoE [#]",1000,0,1000);
TH1F* h_mcTheta = new TH1F("h_mcTheta","MC theta; #theta [deg]; NoE [#]",400,0,4);
TH1F* h_mcPhi = new TH1F("h_mcPhi","Likelihood value; #phi [deg]; NoE [#]",800,0,8);


void SaveHistograms()
{
	h_nHits->Write();
	h_totalQ->Write();
	h_nHitsCaus->Write();
	h_nStringsCaus->Write();
	h_chi2Caus->Write();
	h_nHitsTFilter->Write();
	h_nStringsTFilter->Write();
	h_chi2TFilter->Write();
	h_nHitsChange->Write();
	h_likelihood->Write();
	h_mcEnergy->Write();
	h_mcTheta->Write();
	h_mcPhi->Write();
	h_exitStatus->Write();
	h_nHitsTrack->Write();
}

int SaveCascadeJSON(int eventID, UnifiedEvent& event)
{
	if (gInputType != 0)
		return 1;

	TString outputFileName;
	if (App::Output == "" || App::Output == "a")
		outputFileName = BARS::Data::Directory(BARS::Data::JOINT, BARS::App::Season, BARS::App::Cluster, BARS::App::Run, gProductionID.c_str());
	else
		outputFileName = Form("%s/exp%d/cluster%d/%04d/",App::Output.Data(),BARS::App::Season,BARS::App::Cluster,BARS::App::Run);
	char fname[100];
	std::sprintf(fname,"cascade_season%d_cluster%d_run%d_evt%d.json",BARS::App::Season, BARS::App::Cluster, BARS::App::Run, eventID);
	outputFileName += fname;
	std::ofstream fOutputFile;
	fOutputFile.open(outputFileName);
	fOutputFile<<"{"<<std::endl;
	fOutputFile<<"\t\"eventID\": "<<eventID<<","<<std::endl;
	fOutputFile<<"\t\"season\": "<<BARS::App::Season<<","<<std::endl;
	fOutputFile<<"\t\"cluster\": "<<BARS::App::Cluster<<","<<std::endl;
	fOutputFile<<"\t\"run\": "<<BARS::App::Run<<","<<std::endl;
	fOutputFile<<"\t\"pulses\": {"<<std::endl;

	Int_t nImpulse = gPulses.size();
	for (int i = 0; i < nImpulse-1; i++)
	{
		fOutputFile<<"\t\t\"" << i << "\": \{ \"mask\": "<< 1 <<
		  ", \"amplitude\": "<< gPulses[i].charge <<
		  ", \"charge\": "<< gPulses[i].charge <<
		  ", \"time\": "<< gPulses[i].time <<
		  ", \"channelID\": "<< gPulses[i].OMID <<
		  " },"<<std::endl;
	}
	fOutputFile<<"\t\t\"" << nImpulse-1 << "\": \{ \"mask\": "<< 1 <<
	  ", \"amplitude\": "<< gPulses[nImpulse-1].charge <<
	  ", \"charge\": "<< gPulses[nImpulse-1].charge <<
	  ", \"time\": "<< gPulses[nImpulse-1].time <<
	  ", \"channelID\": "<< gPulses[nImpulse-1].OMID <<
	  " }"<<std::endl;
	fOutputFile<<"\t},"<<std::endl;
	fOutputFile<<"\t\"geometry\": ["<<std::endl;
	for (int i=0; i<gNOMs-1; i++){
	 fOutputFile<<"\t\t{\"channelID\": "<<i<<
	            ", \"x\": "<<gOMpositions[i].X()<<
	            ", \"y\": "<<gOMpositions[i].Y()<<
	            ", \"z\": "<<gOMpositions[i].Z()<<"},"<<std::endl;
	}
	fOutputFile<<"\t\t{\"channelID\": "<<gNOMs-1<<
	            ", \"x\": "<<gOMpositions[gNOMs-1].X()<<
	            ", \"y\": "<<gOMpositions[gNOMs-1].Y()<<
	            ", \"z\": "<<gOMpositions[gNOMs-1].Z()<<"}"<<std::endl;
	fOutputFile<<"\t],"<<std::endl;
	fOutputFile<<"\t\"origins\": {"<<std::endl;
	fOutputFile<<"\t\t\"cascades\": [{"<<std::endl;
	fOutputFile<<"\t\t\t\"mc\": false,"<<std::endl;
	fOutputFile<<"\t\t\t\"title\": \"cascadeFit\","<<std::endl;
	fOutputFile<<"\t\t\t\"direction\": {"<<std::endl;
	fOutputFile<<"\t\t\t\t\"theta\": "<<std::right<<event.theta<<","<<std::endl;
	fOutputFile<<"\t\t\t\t\"phi\": "<<std::right<<event.phi<<","<<std::endl;
	fOutputFile<<"\t\t\t\t\"x\": "<<std::right<<event.position.X()<<","<<std::endl;
	fOutputFile<<"\t\t\t\t\"y\": "<<std::right<<event.position.Y()<<","<<std::endl;
	fOutputFile<<"\t\t\t\t\"z\": "<<std::right<<event.position.Z()<<""<<std::endl;
	// fOutputFile<<"\t\t\t\t\"time\": "<<std::right<<time<<std::endl;
	fOutputFile<<"\t\t\t}"<<std::endl;
	fOutputFile<<"\t\t}]"<<std::endl;
	fOutputFile<<"\t}"<<std::endl;
	fOutputFile<<"}"<<std::endl;
	fOutputFile.close();

	return 0;
}

void TransformToUnifiedEvent(BExtractedImpulseTel* impulseTel, UnifiedEvent &unifiedEvent)
{
	unifiedEvent.hits.clear();
	unifiedEvent.nHits = 0;
	unifiedEvent.qTotal = 0;
	for (int i = 0; i < impulseTel->GetNimpulse(); ++i)
	{
		int OMID = impulseTel->GetNch(i);
		if (gOMqCal[OMID != -1] && gOMtimeCal[OMID] != -1 && gOMpositions[OMID].Mag() !=0 && impulseTel->GetQ(i) > 0)
		{
			unifiedEvent.nHits++;
			unifiedEvent.qTotal += impulseTel->GetQ(i)/gOMqCal[OMID];
			unifiedEvent.hits.push_back(UnifiedHit{impulseTel->GetNch(i),5*(impulseTel->GetT(i)-512)-gOMtimeCal[OMID],impulseTel->GetQ(i)/gOMqCal[OMID],-1,false,-1});
		}
	}
}

void TransformToUnifiedEvent(BEvent* event, BMCEvent* mcEvent, BEventMaskMC* eventMask, UnifiedEvent &unifiedEvent)
{

	// unifiedEvent.mcTheta = mcEvent->GetMuonTrack(0)->GetPolarAngle();
	unifiedEvent.mcTheta = mcEvent->GetPrimaryParticlePolar()/180*TMath::Pi();
	// unifiedEvent.mcPhi = mcEvent->GetMuonTrack(0)->GetAzimuthAngle();
	unifiedEvent.mcPhi = mcEvent->GetPrimaryParticleAzimuth()/180*TMath::Pi();

	double highestEnergy = 0;
	int muonID = -1;
	int cascadeID = -1;

	// cout << "New Event" << endl;
	for (int k = 0; k < mcEvent->GetResponseMuonsN(); ++k)
	{
		// cout << mcEvent->GetTrack(k)->GetMuonEnergy() << endl;
		for (int i = 0; i < mcEvent->GetTrack(k)->GetInteractionN(); ++i)
		{
			// cout << mcEvent->GetTrack(k)->GetInteraction(i)->GetEnergy() << endl;
			if (mcEvent->GetTrack(k)->GetInteraction(i)->GetEnergy() > highestEnergy)
			{
				highestEnergy = mcEvent->GetTrack(k)->GetInteraction(i)->GetEnergy();
				muonID = k;
				cascadeID = i;
			}
		};
		// for (int j = 0; j < mcEvent->GetMuonTrack(k)->GetNumShowers(); ++j)
		// {
		// 	if (mcEvent->GetMuonTrack(k)->GetShower(j)->GetE() > highestEnergy)
		// 	{
		// 		highestEnergy = mcEvent->GetMuonTrack(k)->GetShower(j)->GetE();
		// 		muonID = k;
		// 		cascadeID = j;
		// 	}
		// 	mcEvent->GetMuonTrack(k)->GetShower(j)->GetDirection().Print();
		// 	mcEvent->GetMuonTrack(k)->GetShower(j)->GetPosition().Print();
		// 	cout << k << " " << j << " " << mcEvent->GetMuonTrack(k)->GetShower(j)->GetX() << " " << mcEvent->GetMuonTrack(k)->GetShower(j)->GetY() << " " << mcEvent->GetMuonTrack(k)->GetShower(j)->GetZ() << " " << mcEvent->GetMuonTrack(k)->GetShower(j)->GetE() << endl;
		// }
	}
	if (muonID != -1 && cascadeID != -1)
	{
		unifiedEvent.mcEnergy = mcEvent->GetTrack(muonID)->GetInteraction(cascadeID)->GetEnergy();
		unifiedEvent.mcPosition.SetX(mcEvent->GetTrack(muonID)->GetInteraction(cascadeID)->GetX());
		unifiedEvent.mcPosition.SetY(mcEvent->GetTrack(muonID)->GetInteraction(cascadeID)->GetY());
		unifiedEvent.mcPosition.SetZ(mcEvent->GetTrack(muonID)->GetInteraction(cascadeID)->GetZ());
		unifiedEvent.mcFlagID = (cascadeID+1)*1000+(muonID+1);
		// cout << muonID << " " << cascadeID << " " << unifiedEvent.mcEnergy << " " << unifiedEvent.mcPosition.X() << " " << unifiedEvent.mcPosition.Y() << " " << unifiedEvent.mcPosition.Z() << endl;
	}
	unifiedEvent.hits.clear();
	int nHits = 0;
	int nNoiseHits = 0;
	unifiedEvent.qTotal = 0;
	for (int i = 0; i < event->NHits(); ++i)
	{
		// if (event->Q(i) > 0 && eventMask->GetFlag(i) != 0 && eventMask->GetFlag(i) == unifiedEvent.mcFlagID)
		if (event->Q(i) > 0)
		{
			nHits++;
			unifiedEvent.qTotal += event->Q(i);
			if (eventMask->GetFlag(i) == 0)
			{
				nNoiseHits++;
				unifiedEvent.hits.push_back(UnifiedHit{event->HitChannel(i),event->T(i),event->Q(i),-1,true,0});
			}
			else
				unifiedEvent.hits.push_back(UnifiedHit{event->HitChannel(i),event->T(i),event->Q(i),-1,false,eventMask->GetFlag(i)});
		}
	}
	unifiedEvent.nHits = nHits;
	unifiedEvent.nNoiseHits = nNoiseHits;
	unifiedEvent.nSignalHits = nHits-nNoiseHits;
}

void TransformToUnifiedEvent(mcCascade* cascade, UnifiedEvent &unifiedEvent)
{
	unifiedEvent.hits.clear();
	unifiedEvent.nHits = 0;
	unifiedEvent.qTotal = 0;
	int nHits = 0;
	int nNoiseHits = 0;
	for (int i = 0; i < cascade->nHits; ++i)
	{
		nHits++;
		unifiedEvent.hits.push_back(UnifiedHit{cascade->chID[i]-1,cascade->time[cascade->chID[i]-1],cascade->charge[cascade->chID[i]-1],-1,false,1});
		unifiedEvent.qTotal += cascade->charge[cascade->chID[i]-1];
	}

	unifiedEvent.nHits = nHits;
	unifiedEvent.nNoiseHits = nNoiseHits;
	unifiedEvent.nSignalHits = nHits-nNoiseHits;

	unifiedEvent.mcTheta = TMath::Pi()-TMath::ACos(cascade->cosTheta);
	unifiedEvent.mcPhi = cascade->phi+TMath::Pi();
	if (unifiedEvent.mcPhi > 2*TMath::Pi())
		unifiedEvent.mcPhi -= 2*TMath::Pi();

	unifiedEvent.mcEnergy = cascade->showerEnergy;
	unifiedEvent.mcPosition.SetX(cascade->position[0]);
	unifiedEvent.mcPosition.SetY(cascade->position[1]);
	unifiedEvent.mcPosition.SetZ(cascade->position[2]);
}

int GenerateNoise(UnifiedEvent &event)
{
	int nGeneratedNoiseHits = 0;
	for (int i = 0; i < gNOMs; ++i)
	{
		int noisePulses = gRanGen.Poisson(gMCNoiseRateInkHz/200.0);
		// cout << "OMID: " << i << " nNoisePulses: " << noisePulses << endl;
		for (int j = 0; j < noisePulses; ++j)
		{
			double noiseCharge = 0;//gRanGen.Landau(1,2);
			double ranNumber = gRanGen.Uniform(1);
			for (unsigned int j = 0; j < gNoiseTable.size(); ++j)
			{
				if (gNoiseTable[j] > ranNumber)
				{
					noiseCharge = j*0.1-0.05;
					break;
				}
			}
			double noiseTime = gRanGen.Uniform(5120)-2560;
			event.hits.push_back(UnifiedHit{i,noiseTime,noiseCharge,-1,true,0});
			nGeneratedNoiseHits++;
			// cout << "\t" << j << " " << noiseTime << " " << noiseCharge << endl;
		}
	}
	event.nHits += nGeneratedNoiseHits;
	event.nNoiseHits += nGeneratedNoiseHits;
	return 0;
}

// structure holding processing and filtration statistics
void PrintEventStats(EventStats* es)
{
	std::cout << std::string(81,'*') << std::endl;
	cout << "Nentries: \t" << es->nEntries << "\t(" << es->nEntries/(double)es->nEntries*100 << "%)" << endl;
	cout << "After NFilter: \t" << es->nNFilter << "\t(" << es->nNFilter/(double)es->nEntries*100 << "%)" << endl;
	cout << "After SixThreeFilter: \t" << es->nSixThrees << "\t(" << es->nSixThrees/(double)es->nEntries*100 << "%)" << endl;
	// cout << "After QFilter: \t" << es->nQFilter << endl;
	cout << "After QFilterChi2: \t" << es->nQFilterChi2 << "\t(" << es->nQFilterChi2/(double)es->nEntries*100 << "%)" << endl;
	cout << "After TFilter: \t" << es->nTFilter << "\t(" << es->nTFilter/(double)es->nEntries*100 << "%)" << endl;
	cout << "After TFilterChi2: \t" << es->nTFilterChi2 << "\t(" << es->nTFilterChi2/(double)es->nEntries*100 << "%)" << endl;
	// cout << "After ZFilter: \t" << es->nZFilter << endl;
	// cout << "After TDelayFilter: \t" << es->nTDelayFilter << endl;
	// cout << "After QRatioFilter: \t" << es->nQRatioFilter << endl;
	// cout << "After BranchFilter: \t" << es->nBranchFilter << endl;
	// cout << "After CloseHitsFilter: \t" << es->nCloseHitsFilter << endl;
	cout << "After LikelihoodFiter: \t" << es->nLikelihoodFilter  << "\t(" << es->nLikelihoodFilter/(double)es->nEntries*100 << "%)" << endl;
	std::cout << std::string(81,'*') << std::endl;
}

// Prints values of all important variables
void PrintConfig(void)
{
	// std::cout << std::string(81,'*') << std::endl;
	std::cout << "Reconstruction configuration: " << std::endl;
	std::cout << std::string(40,'-') << std::endl;
	std::cout << "NCut: " << gNCut << endl;
	std::cout << "QTotalCut: " << gQTotalCut << endl;
	std::cout << "CausQCut: " << gCausQCut << endl;
	std::cout << "CausFriendCut: " << gCausFriendCut << endl;
	std::cout << "CausChi2Cut: " << gQCutChi2 << endl;
	std::cout << "TCutTimeWindowNs: " << gTCutTimeWindowNs << endl;
	std::cout << "NCutT: " << gNCutT << endl;
	std::cout << "TFilterChi2Cut: " << gTCutChi2 << endl;
	std::cout << "LikelihoodCut: " << gLikelihoodCut << endl;
	std::cout << "UseMultiDirFit: " << gUseMultiDirFit << endl;
	std::cout << std::string(81,'*') << std::endl;
}

// Prints Welcome screen and calls PrintConfig()
void PrintHeader(void)
{
	std::cout << std::string(81,'*') << std::endl;
	std::cout << std::string(1,'*') << std::string(79,' ') << std::string(1,'*') << std::endl;
    std::cout << std::string(1,'*') <<  "\t\tWelcome to the cascade reconstruction App ver2.0!\t\t" << std::string(1,'*') << std::endl;
	std::cout << std::string(1,'*') << std::string(79,' ') << std::string(1,'*') << std::endl;
    std::cout << std::string(81,'*') << std::endl;

    PrintConfig();
}

// Prints important information about experimental data run
void PrintRunInfo(TTree* tree, BExtractedHeader* header)
{
	tree->GetEntry(0);
	long startTime = header->GetTime().GetSec();
	tree->GetEntry(tree->GetEntries()-1);
	long endTime = header->GetTime().GetSec();

	cout << "RunInfo (Number of entries, RunTime [hours], runTime [days])" << endl;
	cout << "Experimental Data" << endl;
	cout << "Season: " << BARS::App::Season << " Cluster: " << BARS::App::Cluster << " Run: " <<  BARS::App::Run << endl;
	cout << "! " << tree->GetEntries() << " " << (endTime-startTime)/3600.0 << " " << (endTime-startTime)/3600.0/24.0 << endl;
	std::cout << std::string(81,'*') << std::endl;
}

void PrintRunInfo(const char* filePath, TChain* events)
{
	cout << "RunInfo (Number of entries, RunTime [hours], runTime [days])" << endl;
	cout << "MC Data: up-going single muons (" << (gInputType == 2 ? 1 : 0) << ") down-going muon bundles (" << (gInputType == 3 ? 1 : 0) << ")" << endl;
	cout << "Files: " << filePath << " Number of files: " << events->GetListOfFiles()->GetEntries() << " Time constant: " <<  gMCTimeConstant << endl;
    cout << "! " << events->GetEntries() << " " << (events->GetListOfFiles()->GetEntries()*gMCTimeConstant)/3600.0 << "  " << (events->GetListOfFiles()->GetEntries()*gMCTimeConstant)/3600.0/24.0 << endl;

	std::cout << std::string(81,'*') << std::endl;
}

void PrintRunInfoMCCascades(const char* filePath, TChain* events)
{
	cout << "RunInfo (Number of entries, RunTime [hours], runTime [days])" << endl;
	cout << "MC Data: Zhan-Arys' cascades" << endl;
	cout << "Files: " << filePath << " Number of files: " << events->GetListOfFiles()->GetEntries() << " Time constant: " <<  0 << endl;
    cout << "! " << events->GetEntries() << " " << 0 << "  " << 0 << endl;

	std::cout << std::string(81,'*') << std::endl;
}

void PrintGPulses()
{
	for (int i = 0; i < gPulses.size(); ++i)
	{
		cout << i << " " << gPulses[i].OMID << " " << gPulses[i].charge << " " << gPulses[i].time << " " << gPulses[i].MCflag << " " << gOMpositions[gPulses[i].OMID].X() << " " << gOMpositions[gPulses[i].OMID].Y() << " " << gOMpositions[gPulses[i].OMID].Z() << endl;
	}
}

void PrintUnifiedEvent(UnifiedEvent &event)
{
	cout << "Nhits:" << event.nHits << endl;
	for (int i = 0; i < event.nHits; ++i)
	{
		cout << "\t" << event.hits[i].OMID << " " << event.hits[i].time << " " << event.hits[i].charge << endl;
	}
}

double ExpectedTime(TVector3 cascPos, double cascadeTime,int OMID)
{
  	double distanceToCascade = (cascPos-gOMpositions[OMID]).Mag();
  	double expectedTime = cascadeTime + distanceToCascade*gRecCinWater;// + scattering_correction; //v nanosekundach
  	return  expectedTime;
}

// minimization function
void chi2(Int_t &npar, Double_t* gin, Double_t &f, Double_t* par, Int_t iflag) //keep all these nonsense parameters
{
	double theChi2 = 0;
	double constant = 1.0/(gPulses.size() - 4);
  	double error = 1.0/4; //error is 3 ns - photomultiplier
  	double theChi;

	for (unsigned int i = 0; i < gPulses.size(); ++i)
	{
		// chi calculation (measured - theoretical)
		TVector3 cascPos(par[0], par[1], par[2]);
		theChi = error*(gPulses[i].time - ExpectedTime(cascPos,par[3],gPulses[i].OMID));
		theChi2 += theChi*theChi;
	}
	f = constant*theChi2; // function returns calculated chi2, it is global function which is used in SetFCN()
}

// Input: calculated parameters R,Z,phi,cosTheta Output: given lower indexes in 4D array
int GetIndexes(double* inputs, int* outputs)
{
	if (inputs[0] < 0)
		return -1;
	if (inputs[2] < 0 || inputs[2] > TMath::Pi())
		return -1;
	if (inputs[3] < -1 || inputs[3] > 1)
		return -1;

	if (inputs[0] < 199 )
		outputs[0] = (int)floor(inputs[0]);
	else
	{
		outputs[0] = 198;
		inputs[0] = 199;
	}
	if (inputs[1] < 150 && inputs[1] >= -200)
		outputs[1] = (int)floor(inputs[1]+200);
	else
	{
		if (inputs[1] < -200)
		{
			inputs[1] = -200;
			outputs[1] = 0;
		}else{
			outputs[1] = 349;
			inputs[1] = 150;
		}
	}
	outputs[3] = (int)floor((1-inputs[3])/0.1);
	outputs[2] = (int)floor(inputs[2]/(TMath::Pi()/20));

	return 0;
}

double interpolate(double x, double x1, double x2, double f1, double f2)
{
	return (x-x1)/(x2-x1)*f2 + (x2-x)/(x2-x1)*f1;
}

double interpolate4D(double* inputs, int* indexes)
{
	double R0Z0Phi0 = interpolate(inputs[3],1-(indexes[3])*0.1,1-(indexes[3]+1)*0.1,gLogTable4D[indexes[0]][indexes[1]][indexes[2]][indexes[3]],gLogTable4D[indexes[0]][indexes[1]][indexes[2]][indexes[3]+1]);
	double R0Z0Phi1 = interpolate(inputs[3],1-(indexes[3])*0.1,1-(indexes[3]+1)*0.1,gLogTable4D[indexes[0]][indexes[1]][indexes[2]+1][indexes[3]],gLogTable4D[indexes[0]][indexes[1]][indexes[2]+1][indexes[3]+1]);

	double R0Z1Phi0 = interpolate(inputs[3],1-(indexes[3])*0.1,1-(indexes[3]+1)*0.1,gLogTable4D[indexes[0]][indexes[1]+1][indexes[2]][indexes[3]],gLogTable4D[indexes[0]][indexes[1]][indexes[2]][indexes[3]+1]);
	double R0Z1Phi1 = interpolate(inputs[3],1-(indexes[3])*0.1,1-(indexes[3]+1)*0.1,gLogTable4D[indexes[0]][indexes[1]+1][indexes[2]+1][indexes[3]],gLogTable4D[indexes[0]][indexes[1]][indexes[2]][indexes[3]+1]);

	double R1Z0Phi0 = interpolate(inputs[3],1-(indexes[3])*0.1,1-(indexes[3]+1)*0.1,gLogTable4D[indexes[0]+1][indexes[1]][indexes[2]][indexes[3]],gLogTable4D[indexes[0]][indexes[1]][indexes[2]][indexes[3]+1]);
	double R1Z0Phi1 = interpolate(inputs[3],1-(indexes[3])*0.1,1-(indexes[3]+1)*0.1,gLogTable4D[indexes[0]+1][indexes[1]][indexes[2]+1][indexes[3]],gLogTable4D[indexes[0]][indexes[1]][indexes[2]][indexes[3]+1]);

	double R1Z1Phi0 = interpolate(inputs[3],1-(indexes[3])*0.1,1-(indexes[3]+1)*0.1,gLogTable4D[indexes[0]+1][indexes[1]+1][indexes[2]][indexes[3]],gLogTable4D[indexes[0]][indexes[1]][indexes[2]][indexes[3]+1]);
	double R1Z1Phi1 = interpolate(inputs[3],1-(indexes[3])*0.1,1-(indexes[3]+1)*0.1,gLogTable4D[indexes[0]+1][indexes[1]+1][indexes[2]+1][indexes[3]],gLogTable4D[indexes[0]][indexes[1]][indexes[2]][indexes[3]+1]);

	double R0Z0 = interpolate(inputs[2],(indexes[2])*(TMath::Pi()/20),(indexes[2]+1)*(TMath::Pi()/20),R0Z0Phi0,R0Z0Phi1);
	double R0Z1 = interpolate(inputs[2],(indexes[2])*(TMath::Pi()/20),(indexes[2]+1)*(TMath::Pi()/20),R0Z1Phi0,R0Z1Phi1);

	double R1Z0 = interpolate(inputs[2],(indexes[2])*(TMath::Pi()/20),(indexes[2]+1)*(TMath::Pi()/20),R1Z0Phi0,R1Z0Phi1);
	double R1Z1 = interpolate(inputs[2],(indexes[2])*(TMath::Pi()/20),(indexes[2]+1)*(TMath::Pi()/20),R1Z1Phi0,R1Z1Phi1);

	double R0 = interpolate(inputs[1],(indexes[1]-200),(indexes[1]+1-200),R0Z0,R0Z1);
	double R1 = interpolate(inputs[1],(indexes[1]-200),(indexes[1]+1-200),R1Z0,R1Z1);

	return interpolate(inputs[0],(indexes[0]),(indexes[0]+1),R0,R1);
}

// interpolation of the values from 4D log likelihood table. Input: 4D array, R,Z,phi,cosTheta
double GetInterpolatedValue(double* in)
{
	int out[4]{0};

	int returnValue = GetIndexes(in,out);
	if (returnValue == -1)
	{
		cerr << "Table parameters ERROR" << endl;
		exit(1);
	}
	// cout << "IN values " << in[0] << " " << in[1] << " " << in[2] << " " << in[3] << endl;
	// cout << "Out values " << out[0] << " " << out[1] << " " << out[2] << " " << out[3] << endl;
	// cout << interpolate4D(in,out) << endl;
	return interpolate4D(in,out);
}

// cascadeParameters[7]: X,Y,Z,Time,Energy,Theta,Phi
// tableParameters[4]: R,Z,Phi',CosTheta'
void GetParameters(const Double_t* cascadeParameters,const int OMID, double* tableParameters)
{
	TVector3 OMpos(gOMpositions[OMID].X()-cascadeParameters[0],gOMpositions[OMID].Y()-cascadeParameters[1],gOMpositions[OMID].Z()-cascadeParameters[2]);
	// cout << gOMpositions[OMID].X() << " " << gOMpositions[OMID].Y() << " " << gOMpositions[OMID].Z() << endl;
	// cout << cascadeParameters[0] << " " << cascadeParameters[1] << " " << cascadeParameters[2] << endl;
	// OMpos.Print();

	TVector3 showerRef(0,0,1);
	showerRef.SetTheta(cascadeParameters[5] - TMath::Pi());
	showerRef.SetPhi((cascadeParameters[6] + TMath::Pi()));
	// showerRef.Print();
	// OMpos.Print();

	tableParameters[0] = OMpos.Perp(showerRef);
	tableParameters[1] = OMpos*showerRef;

	// cout << "New" << endl;
	TVector3 y = OMpos.Cross(showerRef);
	// y.Print();
	y.SetMag(1.0);
	TVector3 x = OMpos.Cross(y);
	// x.Print();
	x.SetMag(1.0);

	// tableParameters[3] = rhoProjection.Angle(omegaProjection);
	TVector3 OMorien(0,0,-1);
	tableParameters[3] = OMpos*OMorien/OMpos.Mag();

	TVector3 par = OMpos;
	// par.Print();
	par.SetMag(OMpos*OMorien/OMpos.Mag());
	TVector3 phi = OMorien - par;

	tableParameters[2] = phi.Angle(x);

}

bool NotInGPulses(int OMID)
{
	bool inGPulses = false;
	for (unsigned int i = 0; i < gPulses.size(); ++i)
	{
		if (gPulses[i].OMID == OMID)
		{
			inGPulses = true;
			break;
		}
	}
	return !inGPulses;
}

double GetNoiseProbability(double measuredCharge)
{
	if(measuredCharge >= 50.0 || gUseNoiseHitLikelihoodTerm == false) // Noise Probability for charge > 50 p.e. is not stored in gNoiseProbability vector
		return 0;
	else{
		int roundedCharge = int(measuredCharge); // e.g measured charge = 14.8 p.e. is rounded to 14 and it is used for finding noise probability in gNoiseProbability vector
		
		return gNoiseProbability[roundedCharge];
	}	
}

void logLikelihood(Int_t &npar, Double_t* gin, Double_t &f, Double_t* par, Int_t iflag)
{
	// cout << "In log" << endl;
	double logLike = 0;
	double tableParameters[4]{0};
	// cout << "Calculating logLike" << endl;
	// cout << par[0] << " " << par[1] << " " << par[2] << " " << par[3] << " " << par[4] << " " << par[5] << " " << par[6] << endl;
	for (unsigned int i = 0; i < gPulses.size(); ++i)
	{
		// cout << "GPulses: " << gPulses[i].OMID << endl;
		GetParameters(par,gPulses[i].OMID,tableParameters);
		// gOMpositions[gPulses[i].OMID].Print();
		// cout << i << " " << tableParameters[0] << " " << tableParameters[1] << " " << tableParameters[2] << " " << tableParameters[3] << endl;
		double expectedNPE = GetInterpolatedValue(tableParameters);
		// cout << expectedNPE << " " << expectedNPE*100000000*par[4] << " " << gPulses[i].charge << " " << TMath::PoissonI(gPulses[i].charge,expectedNPE*100000000*par[4]) << endl;

		if ((TMath::Poisson(gPulses[i].charge,expectedNPE*110000000*par[4]) + GetNoiseProbability(gPulses[i].charge)) > 10e-320)
		{
			logLike -= TMath::Log10(TMath::Poisson(gPulses[i].charge,expectedNPE*110000000*par[4]) + GetNoiseProbability(gPulses[i].charge));
			// cout << TMath::Log10(TMath::PoissonI(gPulses[i].charge,expectedNPE*100000000*par[4])) << endl;
		}
		else
		{
			logLike -= -320;
			// cout << -350 << endl;
		}
	}

	if(gUseNonHitLikelihoodTerm)
	{
		for (int i = 0; i < gNOMs; ++i)
		{
			if (NotInGPulses(i))
			{
				GetParameters(par,i,tableParameters);
				double expectedNPE = GetInterpolatedValue(tableParameters);
				// cout << i << " " << TMath::Poisson(0,expectedNPE*100000000*par[4]) << endl;
	
				if (TMath::Poisson(0,expectedNPE*110000000*par[4]) > 10e-320)
				{
					logLike -= TMath::Log10(TMath::Poisson(0,expectedNPE*110000000*par[4]));
					// cout << TMath::Log10(TMath::PoissonI(gPulses[i].charge,expectedNPE*100000000*par[4])) << endl;
					// cout << i << " " << TMath::Log10(TMath::PoissonI(0,expectedNPE*100000000*par[4])) << endl;
				}
				else
				{
					logLike -= -320;
					// cout << i << " " << -350 << endl;
				}
			}
		}
		// cout << "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF: " << logLike << endl;
	}


	f = logLike;
}


// Function that checks if all arguments necessary for Experimental Data processing have been set
bool CheckInputParamsExpData()
{
	if (BARS::App::Season == -1)
    {
    	std::cout << "Set season with -s" << std::endl;
    	return false;
    }

    if (BARS::App::Cluster == -1)
    {
    	std::cout << "Set cluster with -c" << std::endl;
    	return false;
    }

    if (BARS::App::Run == -1)
    {
    	std::cout << "Set run with -r" << std::endl;
    	return false;
    }
    return true;
}

bool CheckInputParamsMCData()
{
	// if (gMCNu && gMCMu)
	// {
	// 	std::cout << "ERROR: You can't process both MC neutrinos and MC muons at once!" << std::endl;
	// 	return false;
	// }
    return true;
}

int SetOMsDynamic(BGeomTel* bgeom) //dynamic posiions
{
	int nOKOMs = 0;

	for (int i = 0; i < bgeom->GetNumOMs(); ++i)
	{
		gReferencePosition = (bgeom->At(270)->GetPosition() + bgeom->At(271)->GetPosition());
		gReferencePosition *= 0.5;
		if (bgeom->At(i)->HasData())
		{
			gOMpositions[i] = bgeom->At(i)->GetPosition() - gReferencePosition;
			// cout << i << " " << gOMpositions[i].X() << " " << gOMpositions[i].Y() << endl;
			if (gVerbose == 2)
				cout << i << " " << gOMpositions[i].X() << " " << gOMpositions[i].Y() << " " << gOMpositions[i].Z() << endl;
			nOKOMs++;
		}else
		{
			gOMpositions[i].SetXYZ(0,0,0);
		}
	}
	if (nOKOMs != gNOMs)
		cerr << "Number of OMs with known geometry: " << nOKOMs << endl;
	return nOKOMs;
}

int ReadGeometry(TTree* tree, BExtractedHeader* header) // read dynamic geometry
{
	const char* geometryFileName;
	if (BARS::App::Season != 2020)
		geometryFileName = BARS::Geom::File(BARS::App::Cluster, BARS::App::Season, BARS::Geom::OM_EXPORT_LINEAR);
	else
		geometryFileName = Form("/home/fajtak/geometry-tmp/2019/cluster%d/geometry.export-linear.root",BARS::App::Cluster);

	TTree* geometryTree = nullptr;
	TFile* geomFile = new TFile(geometryFileName,"READ");
	geometryTree = (TTree*)geomFile->Get("Events");

 	// check if the file is opened
	if (!geomFile || !geometryTree)
	{
		cout<<"No TTree called Event in the geometry was found!"<<endl;
		return -1;
	}

	int entryID = 0;
	Int_t startTime = 0;

	do
	{
		tree->GetEntry(entryID);
		startTime = header->GetTime().GetSec();
		entryID++;
	}while (startTime == 0 && entryID < tree->GetEntries());

	BDynamicGeometry* telGeometry = NULL;
	geometryTree->SetBranchAddress("BGeomTel.", &telGeometry);

    //extract time of the first and the last geometry record
	geometryTree->GetEntry(0);
	Int_t geometryStartTime = telGeometry->GetTime().GetSec();
	geometryTree->GetEntry(geometryTree->GetEntries()-1);
	Int_t geometryEndTime = telGeometry->GetTime().GetSec();


	if (startTime == 0)
	{
		cerr << "The run startTime could not be extracted from BExtractedHeader" << endl;
		cerr << "Set Geometry (the first geometryEntry) can be completely wrong!!!" << endl;
	}
  	//check if the geometry file covers the whole run
  	if (geometryStartTime > startTime)
  	{
  		geometryTree->GetEntry(0);
  		int nOKOMs = SetOMsDynamic(telGeometry);
    	cerr<< "The precise dynamic geometry for this run was not available (startGeometry > startRun)" << endl;
    	cerr << "StartGeometry: " << geometryStartTime << " startRun: " << startTime << endl;
    	cerr<< "The first accessible detector geometry is used. The time difference: " << (geometryStartTime-startTime)/3600.0/24.0 << " days." << endl;
    	return 0;
  	}
  	if (geometryEndTime < startTime)
  	{
  		geometryTree->GetEntry(geometryTree->GetEntries()-1);
  		int nOKOMs = SetOMsDynamic(telGeometry);
    	cerr<< "The precise dynamic geometry for this run was not available (endGeometry < startRun)" << endl;
    	cerr<< "The last accessible detector geometry is used. The time difference: " << (startTime-geometryEndTime)/3600.0/24.0 << " days." << endl;
    	return 0;
  	}

  	Int_t geometryTime = 0;
  	// iterate through all the geometry entries
	for (int i = 0; i < geometryTree->GetEntries(); ++i)
	{
		geometryTree->GetEntry(i);

	    if(telGeometry->GetTime().GetSec() >= startTime)
	    {
	    	geometryTime = telGeometry->GetTime().GetSec();
	    	// cout << "SubRun start: " << telGeometry->GetTime().AsDouble()<<endl;
	    	int nOKOMs = SetOMsDynamic(telGeometry);
	    	if (gVerbose == 1)
	    	{
	    		cout << endl;
	    		cout << "Run startTime: " << startTime << " GeometryTime: " << geometryTime << " First geometry: " << geometryStartTime << " Last geometry: " << geometryEndTime << endl;
	    	}
	    	return 0;
	    }
	}
	return -1;
}

int ReadGeometryMC(TChain* event)
{
	TFile* file = event->GetFile();
	// cout << file->GetName() << endl;
	TTree* tree = (TTree*)file->Get("ArrayConfig");

	BGeomTelMC* geomMC = NULL;
	tree->SetBranchAddress("BGeomTel.",&geomMC);

	int nOKOMs = 0;

	tree->GetEntry(0);
	for (int i = 0; i < geomMC->GetNumOMs(); ++i)
	{
		gOMpositions[i] = geomMC->At(i)->GetPosition();
		nOKOMs++;
		if ((i > 35 && i < 60) || (i > 71 && i < 84) || i == 130 || i == 131 || i == 245 || i == 246 || i == 247 || i == 256)
			gOMqCal[i] = -1;
	}
	return nOKOMs;
}

int ReadGeometryMCCascades()
{
	TString fileName = "/media/zuzana/Data/BaikalData/MC_cascadesDZH/MC_cascadeDZH/array2016_phys_kebkal.dat";

	ifstream inputFile;
    inputFile.open(fileName);

    if (!inputFile)
    {
    	cerr << "Geometry file for MC cascades: " << fileName << " was NOT found. Program termination!" << endl;
    	return -1;
  	}

    int OMID,dummyI;
    double x,y,z,dummyD;

    for (int i = 0; i < gNOMs; ++i)
    {
		inputFile >> OMID >> x >> y >> z >> dummyD >> dummyI;
		gOMpositions[i] = TVector3(x,y,z);
    }
    inputFile.close();

	return 0;
}

//  Reads time calibration files (offsets, timecalib_dzh) with the given structure
int ReadTimeCalFile(const char* fileName, double multConst)
{
	ifstream inputFile;
    inputFile.open(fileName);

    if (!inputFile)
    {
    	cerr << "Calibration file: " << fileName << " was NOT found. Program termination!" << endl;
    	return -1;
  	}

    string dummyLine;

    getline(inputFile,dummyLine);
    getline(inputFile,dummyLine);

    double readValue = 0;

    for (int i = 0; i < gNOMs; ++i)
    {
		if (i%36 == 0)
		{
		    getline(inputFile,dummyLine);
		}
		inputFile >> readValue;
		if (inputFile.fail())
		{
			cerr << "Wrong value in the time calibration file: " << fileName << " Try to check it!" << endl;
			return -2;
		}
		if (readValue != -10000 && readValue != -1000)
			gOMtimeCal[i] += readValue*multConst;
		else
			gOMtimeCal[i] = -1;
		// cout << i << " " << gOMtimeCal[i] << endl;
		if ((i+1)%36 == 0)
		{
		    getline(inputFile,dummyLine);
		    getline(inputFile,dummyLine);
		}
    }
    inputFile.close();
    return 0;
}

void NormalizeTimeCal()
{
	double meanValue = 0;
	int nValues = 0;
	for (int i = 0; i < gNOMs; ++i)
	{
		if (gOMtimeCal[i] != -1)
		{
			if (BARS::App::Cluster == 0 && BARS::App::Season == 2016 && BARS::App::Run > 389 && i > 251)
				gOMtimeCal[i] -= 2500;
			nValues++;
			meanValue += gOMtimeCal[i];
		}
	}
	meanValue /= nValues;
	for (int i = 0; i < gNOMs; ++i)
	{
		if (gOMtimeCal[i] != -1)
		{
			gOMtimeCal[i] -= meanValue;
		}
	}
}

// Both time calibration files (offsets, timecalib_dzh) are read and multiplied by corresponding constant
int ReadTimeCal()
{
	const char* offsetFileName = BARS::Calib::File(BARS::App::Cluster, BARS::App::Season, BARS::Calib::OFFSET, "offsets");
	if (ReadTimeCalFile(offsetFileName,-2.5) < 0)
		return -1;
	const char* timeCalName = BARS::Calib::File(BARS::App::Cluster, BARS::App::Season, BARS::Calib::TIME, "timecalib_dzh");
	if (ReadTimeCalFile(timeCalName,1) < 0)
		return -1;

	NormalizeTimeCal();

    if (gVerbose == 2)
    {
    	cout << endl;
    	for (int i = 0; i < gNOMs; ++i)
    	{
    		cout << gOMtimeCal[i] << "\t";
    		if ((i+1)%12 == 0)
    			cout << endl;
    	}
    }

    return 0;
}

// The charge calibration parameters are read directly from qcalib file created by Zhenya's DQM in the same folder
// It also let us know about operating OMs. If there is -1 in the qcalib for OM we assume it is dead.
int ReadQCal(void)
{
	const char* filePath = BARS::Data::File(BARS::Data::JOINT, BARS::App::Season, BARS::App::Cluster, BARS::App::Run, gProductionID.c_str());
	TString jointFileName(filePath);
	TString chargeFileName(jointFileName(0,jointFileName.Length()-17));
	chargeFileName += "qcalib";

	ifstream inputFile;
    inputFile.open(chargeFileName);

    if (!inputFile)
    {
    	cerr << "Calibration file: " << chargeFileName << " was NOT found. Program termination!" << endl;
    	return -1;
  	}

    string dummyLine;
    double readValue = 0;

    for (int i = 0; i < gNOMs; ++i)
    {
		inputFile >> readValue;
		if (inputFile.fail())
		{
			cerr << "Wrong value in the charge calibration file: " << filePath << " Try to check it!" << endl;
			return -2;
		}
		gOMqCal[i] = readValue;
		
		if ((i+1)%36 == 0)
		{
		    getline(inputFile,dummyLine);
		    getline(inputFile,dummyLine);
		}
    }

    if (gVerbose == 2)
    {
    	cout << endl;
    	for (int i = 0; i < gNOMs; ++i)
    	{
    		cout << gOMqCal[i] << "\t";
    		if ((i+1)%12 == 0)
    			cout << endl;
    	}
    }

    inputFile.close();
    return 0;
}

// Function reading the logLikelihood tables for direction and energy reconstruction
// gLogTable4D structure: Rho, Z, Phi, cosTheta
int ReadLogTable()
{
	// cout << "4D LogTable reading starts" << endl;
	ifstream fTab;
	if (App::Output == ""){
		fTab.open("/home/zuzana/software/barsNew/bars/programs/reconstruct-cascade/inputFiles/hq001200_double.dqho2011", ios::in|ios::binary|ios::ate);
		// fTab.open("/home/zubardac/showerTable/hq001200_double.dqho2011", ios::in|ios::binary|ios::ate);
	}
	else{
		fTab.open("/home/zubardac/showerTable/hq001200_double.dqho2011", ios::in|ios::binary|ios::ate);		
		// fTab.open("./hq001200_double.dqho2011", ios::in|ios::binary|ios::ate);
	}

	if (!fTab.is_open())
		return -1;

	streampos size = 8;
	char * memblock = new char[size];

	fTab.seekg (0, ios::beg);
	fTab.read(memblock,size);
	int dummy;

	for(int i = 0; i < 200; i++){           // step of R in meters
	  for(int j = 0; j < 351; j++){         // step of Z in meters
	      for(int k = 0; k < 21; k++){      // step of phi
	        for(int m = 0; m < 21; m++){     // step of cos(theta)
				fTab.read(memblock,size);
				double* value = (double*)memblock;
				gLogTable4D[i][j][k][m] = *value;
				// cout << i << " " << j << " " << k << " " << m << " " << gLogTable4D[i][j][k][m] << endl;
	        }
	        // cin >> dummy;
	      }
	  }
	}
	fTab.close();
	// cout << "LogTable ends" << endl;
	return 0;
}

int ReadNoiseChargeTable()
{
	ifstream inputFile("/home/zuzana/software/barsNew/bars/programs/reconstruct-cascade/inputFiles/mc-noise-charge-2016.dat", ios::in);

   	if (!inputFile)
    {
    	cerr << "Calibration file: " << "/home/zuzana/software/barsNew/bars/programs/reconstruct-cascade/inputFiles/mc-noise-charge-2016.dat" << " was NOT found. Program termination!" << endl;
    	return -1;
  	}

  	gNoiseTable.clear();

    string dummyLine;
    double readValue = 0;

    getline(inputFile,dummyLine);
    vector<double> noiseTable;

    for (int i = 0; i < 500; ++i)
    {
		inputFile >> readValue;
		gNoiseTable.push_back(readValue);
		// cout << i << " " << readValue << endl;
    }
    inputFile.close();

	return 0;
}

int ReadNoiseProbability()
{
	ifstream inputFile("/home/zuzana/software/barsNew/bars/programs/reconstruct-cascade/inputFiles/noiseProbability_Binwidth_1pe.dat", ios::in);

	if (!inputFile)
    {
    	cerr << "Noise Probability File: " << "/home/zuzana/software/barsNew/bars/programs/reconstruct-cascade/inputFiles/noiseProbability_Binwidth_1pe.dat" << " was NOT found. Program termination!" << endl;
    	return -1;
  	}
  	gNoiseProbability.clear();

	double readValue = 0;	

 	while (inputFile >> readValue){
   		gNoiseProbability.push_back(readValue);  
    }
    inputFile.close();

    return 0;
}

int ReadInputParamFiles(TTree* tree, BExtractedHeader* header)
{
	cout << "Reading Geometry ... ";
	cout << std::flush;
	if (ReadGeometry(tree,header) == -1)
	{
		std::cout << "Problem with geometry file!" << std::endl;
		return -1;
	}
	cout << "DONE!" << endl;
	cout << "Reading Time Calibration ... ";
	cout << std::flush;
	if (ReadTimeCal() == -1)
	{
		std::cout << "Problem with time calibration files!" << std::endl;
		return -1;
	}
	cout << "DONE!" << endl;
	cout << "Reading Charge Calibration ... ";
	cout << std::flush;
	if (ReadQCal() < 0)
	{
		std::cout << "Problem with charge calibration files!" << std::endl;
		return -1;
	}
	cout << "DONE!" << endl;
	cout << "Reading Log Table ... ";
	cout << std::flush;
	if (ReadLogTable() == -1)
	{
		std::cout << "Problem with 4D LogLikelihood file!" << std::endl;
		return -1;
	}
	cout << "DONE!" << endl;
	cout << "Reading Noise Probability File ... ";
	cout << std::flush;
    if (ReadNoiseProbability() == -1)
	{
		std::cout << "Problem with NoiseProbability File file!" << std::endl;
		return -1;
	}
	cout << "DONE!" << endl;
    return 0;
}

int ReadInputParamFiles(TChain* events)
{
	cout << "Reading Geometry ... ";
	cout << std::flush;
    if (ReadGeometryMC(events) == -1)
	{
		std::cout << "Problem with geometry file!" << std::endl;
		return -1;
	}
	cout << "DONE!" << endl;
	cout << "Reading Log Table ... ";
	cout << std::flush;
	if (ReadLogTable() == -1)
	{
		std::cout << "Problem with 4D LogLikelihood file!" << std::endl;
		return -1;
	}
	cout << "DONE!" << endl;
    return 0;
}

int ReadInputParamFiles()
{
	cout << "Reading Geometry ... ";
	cout << std::flush;
    if (ReadGeometryMCCascades() == -1)
	{
		std::cout << "Problem with geometry file!" << std::endl;
		return -1;
	}
	cout << "DONE!" << endl;
	cout << "Reading Noise Charge Table ... ";
	cout << std::flush;
    if (ReadNoiseChargeTable() == -1)
	{
		std::cout << "Problem with NoiseChargeTable file!" << std::endl;
		return -1;
	}
	cout << "DONE!" << endl;
	cout << "Reading Noise Probability File ... ";
	cout << std::flush;
    if (ReadNoiseProbability() == -1)
	{
		std::cout << "Problem with NoiseProbability File file!" << std::endl;
		return -1;
	}
	cout << "DONE!" << endl;
	cout << "Reading Log Table ... ";
	cout << std::flush;
	if (ReadLogTable() == -1)
    {
    	std::cout << "Problem with 4D LogLikelihood file!" << std::endl;
    	return -1;
    }
    cout << "DONE!" << endl;
    return 0;
}

bool OMIDAlreadyInGPulses(UnifiedHit &pulse)
{
	bool OMIDAlreadyIn = false;
	for (int k = 0; k < gPulses.size(); ++k)
	{
		if (gPulses[k].OMID == pulse.OMID)
		{
			if (gPulses[k].charge < pulse.charge)
			{
				gPulses[k].charge = pulse.charge;
				gPulses[k].time = pulse.time;
				gPulses[k].expectedCharge = pulse.expectedCharge;
				gPulses[k].noise = pulse.noise;
				gPulses[k].MCflag = pulse.MCflag;
			}
			OMIDAlreadyIn = true;
			break;
		}
	}
	return OMIDAlreadyIn;
}

bool bigger_charge(const UnifiedHit& x, const UnifiedHit& y) { return x.charge > y.charge; }

int CaussalityFilter(UnifiedEvent &event)
{
	gPulses.clear();

	sort(event.hits.begin(), event.hits.end(), bigger_charge);

	for (int i = 0; i < event.nHits; ++i)
	{
		if (event.hits[i].charge < gCausQCut)
			break;
		bool addPulse = true;
		for (unsigned int j = 0; j < gPulses.size(); ++j)
		{
			double distance = (gOMpositions[event.hits[i].OMID] - gOMpositions[gPulses[j].OMID]).Mag();

			if(abs(gPulses[j].time - event.hits[i].time) >= (distance*gRecCinWater + 50))
			{
				addPulse = false;
				break;
			}
		}
		bool friendFound = false;
		for (int j = 0; j < event.nHits; ++j)
		{
			if (i != j && abs(event.hits[i].OMID-event.hits[j].OMID) <= gCausFriendCut && abs(event.hits[i].time-event.hits[j].time) < abs(event.hits[i].OMID-event.hits[j].OMID)*70+50)
			{
				friendFound = true;
				break;
			}
		}
		if (addPulse && friendFound && !OMIDAlreadyInGPulses(event.hits[i]))
			gPulses.push_back(UnifiedHit{event.hits[i].OMID,event.hits[i].time,event.hits[i].charge,0,event.hits[i].noise,event.hits[i].MCflag});
	}
	return 0;
}

int TFilter(UnifiedEvent &event, TVector3& cascPos, double& cascTime)
{
	gPulses.clear();
	int nPulses = 0;

	for(int i = 0; i < event.nHits; i++)
	{
		double distanceToCascade = (cascPos - gOMpositions[event.hits[i].OMID]).Mag();
        double expectedTime = cascTime + distanceToCascade*gRecCinWater;// + scattering_correction;

        if(TMath::Abs(event.hits[i].time-expectedTime) <= 50)
        {
        	if (OMIDAlreadyInGPulses(event.hits[i]))
        		continue;
        	// if (gOMpositions[cascade->chID[i]-1].Z()-cascPos.Z() < -140 || gOMpositions[cascade->chID[i]-1].Z()-cascPos.Z() > 180)
        		// continue;
        	gPulses.push_back(UnifiedHit{event.hits[i].OMID,event.hits[i].time,event.hits[i].charge,0,event.hits[i].noise,event.hits[i].MCflag});
        	nPulses++;
        }
	}
	return nPulses;
}

int TrackFilter(UnifiedEvent &event, TVector3 cascPos, double cascTime)
{
	int nTrackPulses = 0;

	for(int i = 0; i < event.nHits; i++)
	{
		double distanceToCascade = (cascPos - gOMpositions[event.hits[i].OMID]).Mag();
        double expectedTime = cascTime + distanceToCascade*gRecCinWater;// + scattering_correction;

        if(TMath::Abs(event.hits[i].time-expectedTime) > 50 && TMath::Abs(event.hits[i].time-expectedTime) <= 100)
        {
        	nTrackPulses++;
        }
	}

	return nTrackPulses;
}

int GetNStrings()
{
	int nHitStrings = 0;
	int hitStrings[8]{0};

	for (unsigned int i = 0; i < gPulses.size(); ++i)
	{
		hitStrings[(gPulses[i].OMID)/36] = 1;
	}

	for (int i = 0; i < 8; ++i)
	{
		nHitStrings += hitStrings[i];
	}
	return nHitStrings;
}

int GetNTrackHits(UnifiedEvent &event)
{
	int nTrackHits = 0;
	for (unsigned int i = 0; i < gPulses.size(); ++i)
	{
		if (gPulses[i].MCflag != event.mcFlagID && !gPulses[i].noise)
			nTrackHits++;
	}
	return nTrackHits;
}

bool SixThreeFilterPassed()
{
	int nHitStrings = GetNStrings();
	// cout << nHitStrings << endl;
	// PrintCascade(cascade);

	h_nStringsCaus->Fill(nHitStrings);
	h_nHitsCaus->Fill(gPulses.size());

	if (nHitStrings > 2 && gPulses.size() > 5)
		return true;
	else
		return false;
}

bool NFilterPassed(UnifiedEvent &event)
{
	h_nHits->Fill(event.nHits);
	h_totalQ->Fill(event.qTotal);
	if (event.nHits >= gNCut && event.qTotal > gQTotalCut)
	{
		return true;
	}
	return false;
}

double EstimateInitialPosMatrix(TVector3 &cascPos, double &cascTime)
{
	TMatrixD A(gPulses.size()-1,4);
	TVectorD b(gPulses.size()-1);

	for (unsigned int i = 0; i < gPulses.size()-1; ++i)
	{
		double temp = TMath::Power(gOMpositions[gPulses[i+1].OMID].X(),2) - TMath::Power(gOMpositions[gPulses[i].OMID].X(),2);
		temp += TMath::Power(gOMpositions[gPulses[i+1].OMID].Y(),2) - TMath::Power(gOMpositions[gPulses[i].OMID].Y(),2);
		temp += TMath::Power(gOMpositions[gPulses[i+1].OMID].Z(),2) - TMath::Power(gOMpositions[gPulses[i].OMID].Z(),2);
		temp -= (TMath::Power(gPulses[i+1].time,2) - TMath::Power(gPulses[i].time,2))/TMath::Power(gRecCinWater,2);
		b[i] = temp;

		// cout << cascade->chID[i+1]-1 << " " << cascade->chID[i]-1 << endl;
		A[i][0] = 2*(gOMpositions[gPulses[i+1].OMID].X() - gOMpositions[gPulses[i].OMID].X());
		A[i][1] = 2*(gOMpositions[gPulses[i+1].OMID].Y() - gOMpositions[gPulses[i].OMID].Y());
		A[i][2] = 2*(gOMpositions[gPulses[i+1].OMID].Z() - gOMpositions[gPulses[i].OMID].Z());
		A[i][3] = -2*(gPulses[i+1].time - gPulses[i].time)/TMath::Power(gRecCinWater,2);
	}

	TMatrixD B(4,gPulses.size()-1);
	B.Transpose(A);
	TMatrixD C = (B*A).Invert();
	TVectorD D = A.T()*b;
	TVectorD X = C*D;

	// cout << X[0] << " " << X[1] << " " << X[2] << " " << X[3] << endl;

	cascPos[0] = X[0];
	cascPos[1] = X[1];
	cascPos[2] = X[2];
	cascTime = X[3];

	return 0;
}

double FitCascPos(TVector3 &cascPos, double &cascTime)
{
	gMinuit->ReleaseParameter(0);
	gMinuit->ReleaseParameter(1);
	gMinuit->ReleaseParameter(2);
	gMinuit->ReleaseParameter(3);

	//cout<<"fitEVent "<<endl;
	gMinuit->SetFCN(chi2);
	gMinuit->SetParameter(0,"casc X",cascPos.X(),1,-750,750); //estimation of start point
	gMinuit->SetParameter(1,"casc Y",cascPos.Y(),1,-750,750);
	gMinuit->SetParameter(2,"casc Z",cascPos.Z(),1,-750,750);
	gMinuit->SetParameter(3,"cast T",cascTime,1,-10000,50000);

	gMinuit->FixParameter(4);
	gMinuit->FixParameter(5);
	gMinuit->FixParameter(6);

	double arglist[10] = {500,0.01}; //500 means 500 iterations and then MIGRAD stops
	gMinuit->ExecuteCommand("MIGRAD",arglist,2); //this line performs minimalisation

	double chi2, edm, errdef;
	int nvpar, nparx;

	gMinuit->GetStats(chi2,edm,errdef,nvpar,nparx);

	cascPos.SetX(gMinuit->GetParameter(0));
	cascPos.SetY(gMinuit->GetParameter(1));
	cascPos.SetZ(gMinuit->GetParameter(2));
	cascTime = gMinuit->GetParameter(3);

	return chi2;
}

double FitCascDirection(UnifiedEvent &event, double &energy, double &theta, double &phi, double &eSig, double &tSig, double &pSig)
{
	// cout << "Fit Matrix direction" << endl;
	gMinuit->ReleaseParameter(4);
	gMinuit->ReleaseParameter(5);
	gMinuit->ReleaseParameter(6);

	// cout << "StartFitting" << endl;
	gMinuit->SetFCN(logLikelihood);
	gMinuit->SetParameter(0,"LED X",event.position.X(),0.01,-1000,1000); //estimation of start point
	gMinuit->SetParameter(1,"LED Y",event.position.Y(),0.01,-1000,1000);
	gMinuit->SetParameter(2,"LED Z",event.position.Z(),0.01,0,1000);
	gMinuit->SetParameter(3,"Time",event.time,1,0,20000);
	gMinuit->SetParameter(4,"Energy",energy,1,0,10000);
	gMinuit->SetParameter(5,"Theta",theta,0.1,0,TMath::Pi());
	gMinuit->SetParameter(6,"Phi",phi,0.1,0,2*TMath::Pi());

	gMinuit->FixParameter(0);
	gMinuit->FixParameter(1);
	gMinuit->FixParameter(2);
	gMinuit->FixParameter(3);

	// cout << "FIT" << endl;
	double arglist[10] = {500,0.01}; //500 means 500 iterations and then MIGRAD stops
	gMinuit->ExecuteCommand("MIGRAD",arglist,2); //this line performs minimalisation

	double chi2, edm, errdef;
	int nvpar, nparx;

	gMinuit->GetStats(chi2,edm,errdef,nvpar,nparx);
	energy = gMinuit->GetParameter(4);
	theta = gMinuit->GetParameter(5);
	phi = gMinuit->GetParameter(6);

	eSig = gMinuit->GetParError(4);
	tSig = gMinuit->GetParError(5);
	pSig = gMinuit->GetParError(6);

	// cout << chi2 << " " << gPulses.size() << " " <<  chi2/gPulses.size() << " " << fMinuit->GetParameter(5) << " " << fMinuit->GetParameter(6) << endl;
	// return chi2/gPulses.size();
	// cout<<"chi2/gNOMs "<<chi2/gNOMs<<endl;
	if(gUseNonHitLikelihoodTerm){
		return chi2/gNOMs;
	}else{
		return chi2/gPulses.size();
	}

}

double EstimateInitialDirection(TVector3& cascPos, double& cascTime, double& energy, double &theta, double &phi)
{
	int nPar = 0;
	double* gin = new double(0);
	int iflag = 0;
	double likelihoodValue = 0;
	double cascadeParameters[7];

	double lowestLog = 1000000;

	cascadeParameters[0] = cascPos[0];
	cascadeParameters[1] = cascPos[1];
	cascadeParameters[2] = cascPos[2];
	cascadeParameters[3] = cascTime;


	for (int i = 0; i < 17; ++i)
	{
		for (int j = 0; j < 36; ++j)
		{
			for (int k = 0; k < 1; ++k)
			{
				// Int_t &npar, Double_t* gin, Double_t &f, Double_t* par, Int_t iflag
				cascadeParameters[4] = energy;
				// cascadeParameters[4] = TMath::Power(10,k);
				// cascadeParameters[4] = cascade->showerEnergy;
				// cascadeParameters[4] = 1+k*50;
				cascadeParameters[5] = TMath::Pi()/18*(i+1);
				cascadeParameters[6] = TMath::Pi()*2/36*j;
				logLikelihood(nPar,gin,likelihoodValue,cascadeParameters,iflag);
				// cout << likelihoodValue << " " << cascadeParameters[0] << " " << cascadeParameters[1] << " " << cascadeParameters[2] << endl;
				if (likelihoodValue < lowestLog)
				{
					lowestLog = likelihoodValue;
					energy = cascadeParameters[4];
					theta = cascadeParameters[5];
					phi = cascadeParameters[6];
				}
			}
		}
	}
	return 0;
}

double LikelihoodFilterPassed(UnifiedEvent &event)
{
	event.likelihood = FitCascDirection(event,event.energy,event.theta,event.phi,event.energySigma,event.thetaSigma,event.phiSigma);
	return event.likelihood;
}

double LikelihoodFilterPassedGrid(UnifiedEvent &event)
{
	// cout << "In likelihood" << endl;
	int nThetaSteps = 4;
	int nPhiSteps = 6;
	int nEnergySteps = 4;
	double lowestLog = 10000;
	for (int k = 0; k < nThetaSteps-1; ++k)
	{
		for (int l = 0; l < nPhiSteps-1; ++l)
		{
			for (int i = 0; i < nEnergySteps; ++i)
			{
				double cascadeEnergy = TMath::Power(10,1);
				double cascadeTheta = TMath::Pi()/(nThetaSteps)*(k+1);
				double cascadePhi = TMath::Pi()*2/nPhiSteps*(l+1);
				double cascadeEnergySigma = 0;
				double cascadeThetaSigma = 0;
				double cascadePhiSigma = 0;
				double recentLog = FitCascDirection(event,cascadeEnergy,cascadeTheta,cascadePhi,cascadeEnergySigma,cascadeThetaSigma,cascadePhiSigma);
				if (recentLog < lowestLog)
				{
					lowestLog = recentLog;
					event.energy = cascadeEnergy;
					event.theta = cascadeTheta;
					event.phi = cascadePhi;
					event.energySigma = cascadeEnergySigma;
					event.thetaSigma = cascadeThetaSigma;
					event.phiSigma = cascadePhiSigma;
				}
				// cout << k << " " << l << " " << recentLog  << " " << cascadeEnergy << " " << cascadeTheta << " " << cascadePhi << endl;
			}
			// cout << k << " " << l << endl;
		}
	}
	event.likelihood = lowestLog;
	// cout << "End likelihood" << endl;
	cout <<"likelihood "<< lowestLog << " energy " << event.energy << " theta " << event.theta << " phi " << event.phi << endl;
	return lowestLog;
}

int inRange(int roundedEnergy)
{
    int bin = 0;
    int arraySize = sizeof(gEnergyCorrectionArray)/sizeof(gEnergyCorrectionArray[0]);
    int energyBin = 30; //xlow used in the rebinned histogram of h_mismatchEnergyvsEnergy_pfx

    if(roundedEnergy < energyBin) //if energy is lower than xlow in the histogram, the first element of array is used (i.e energy is lower than 1 TeV)
    	bin = 0;

    for(int i = 0; i < arraySize; i ++) 
    {
        if(energyBin <= roundedEnergy && roundedEnergy < (energyBin+2)) //check if the roundedEnegy is in one of the energy bins, which corresponds to gEnergyCorrectionArray
            bin = i;
        
        energyBin += 2; //increment by 2 due to binning in the rebinned histogram of h_mismatchEnergyvsEnergy_pfx
    }

    if(roundedEnergy >= energyBin) //if energy is higher than xhigh in the histogram, the last element of array is used (i.e energy is higher than 10 PeV)
    	bin = arraySize - 1;

    return bin;   
}

double GetCorrectedEnergy(double energy)
{
    double logEnergy = TMath::Log10(energy*1000)/0.1; // energy (in TeV) is converted to log_10(E [GeV]) to get a number on the xaxis in the histograms Energy vs MismatchEnergy
	int roundedEnergy = (int)logEnergy;
    int bin = inRange(roundedEnergy);
    
	return energy/gEnergyCorrectionArray[bin];     
}

int EventVisualization(int eventID, UnifiedEvent &event, TVector3& cascPos, double& cascTime)
{
	int nHitsPerString[8]{0};
	int nNoiseHitsPerString[8]{0};
	int nTrackHitsPerString[8]{0};
	TGraph* g_hits[gNStrings];
	TGraph* g_noiseHits[gNStrings];
	TGraph* g_trackHits[gNStrings];
	TGraph* g_lowerLimit[gNStrings];
	TGraph* g_upperLimit[gNStrings];
	TGraph* g_trackLimit[gNStrings];
	TMultiGraph* mg_hitsMatrix[gNStrings];
	TGraph* g_ledMatrix[gNStrings];
	TGraph* g_QvsL = new TGraph(gPulses.size());
	g_QvsL->SetMarkerStyle(20);

	int nOMsPerString = gNOMs/gNStrings;

	for(int k = 0; k < event.nHits; k++)
	{
		// PrintPulse(event.hits[k]);
		if (event.hits[k].noise)
			nNoiseHitsPerString[event.hits[k].OMID/36]++;
		else{
			if (event.hits[k].MCflag != event.mcFlagID)
				nTrackHitsPerString[event.hits[k].OMID/36]++;
			else
				nHitsPerString[event.hits[k].OMID/36]++;
		}
	}
	for (int i = 0; i < gNStrings; ++i)
	{
		g_hits[i] = new TGraph(nHitsPerString[i]);
		g_noiseHits[i] = new TGraph(nNoiseHitsPerString[i]);
		g_trackHits[i] = new TGraph(nTrackHitsPerString[i]);
		g_lowerLimit[i] = new TGraph(nOMsPerString);
		g_upperLimit[i] = new TGraph(nOMsPerString);
		g_trackLimit[i] = new TGraph();
		mg_hitsMatrix[i] = new TMultiGraph(Form("mg_%d",i),Form("String_%d;Calibrated time [ns]; OM Z position [m]",i+1));
		nHitsPerString[i] = 0;
		nNoiseHitsPerString[i] = 0;
		nTrackHitsPerString[i] = 0;
		g_ledMatrix[i] = new TGraph(1);
		g_ledMatrix[i]->SetPoint(0,cascTime,cascPos.Z());
	}
	for(int k = 0; k < event.nHits; k++)
	{
		int stringID = event.hits[k].OMID/36;

		if (event.hits[k].noise)
		{
			g_noiseHits[stringID]->SetPoint(nNoiseHitsPerString[stringID],event.hits[k].time,gOMpositions[event.hits[k].OMID].Z());
			nNoiseHitsPerString[stringID]++;
		}
		else{
			if (event.hits[k].MCflag != event.mcFlagID)
			{
				g_trackHits[stringID]->SetPoint(nTrackHitsPerString[stringID],event.hits[k].time,gOMpositions[event.hits[k].OMID].Z());
				nTrackHitsPerString[stringID]++;
			}
			else
			{
				g_hits[stringID]->SetPoint(nHitsPerString[stringID],event.hits[k].time,gOMpositions[event.hits[k].OMID].Z());
				nHitsPerString[stringID]++;
			}
		}

		// g_hits[stringID]->SetPoint(nHitsPerString[stringID],event.hits[k].time,gOMpositions[event.hits[k].OMID].Z());
		// nHitsPerString[stringID]++;
	}

	TVector3 cascDir(0,0,1);
	cascDir.SetTheta(event.theta);
	cascDir.SetPhi(event.phi);
	int nPoints = 0;

	for (int j = 0; j < gNOMs; ++j)
	{
		if (j%36 == 0)
			nPoints = 0;
		double distanceToCascade = (cascPos - gOMpositions[j]).Mag();
	    double expectedTime = cascTime + distanceToCascade*gRecCinWater;
		g_lowerLimit[j/nOMsPerString]->SetPoint(j%nOMsPerString,expectedTime-gTCutTimeWindowNs,gOMpositions[j].Z());
		g_upperLimit[j/nOMsPerString]->SetPoint(j%nOMsPerString,expectedTime+gTCutTimeWindowNs,gOMpositions[j].Z());

		if ((gOMpositions[j]-cascPos).Mag() > 100)
			continue;
		double sPerp = (gOMpositions[j]-cascPos).Perp(cascDir);
		double sLong = (gOMpositions[j]-cascPos)*(cascDir);
		double lLong = sPerp/TMath::Tan(0.719887);
		double expTime = cascTime + (sLong-lLong)*gRecC + TMath::Sqrt(TMath::Power(sPerp,2)+TMath::Power(lLong,2))*gRecCinWater;

		g_trackLimit[j/nOMsPerString]->SetPoint(nPoints,expTime,gOMpositions[j].Z());
		nPoints++;
	}
	for (int i = 0; i < gPulses.size(); ++i)
	{
		double distanceToCascade = (cascPos - gOMpositions[gPulses[i].OMID]).Mag();
		g_QvsL->SetPoint(i,distanceToCascade,gPulses[i].charge);
	}

	TCanvas *cEvent = new TCanvas(Form("cEvent_%d",eventID),Form("cEvent_%d",eventID),200,10,600,400);
	cEvent->Divide(3,3);
	for (int i = 0; i < gNStrings; ++i)
	{
		cEvent->cd(i+1);
		mg_hitsMatrix[i]->Add(g_hits[i],"P");
		mg_hitsMatrix[i]->Add(g_noiseHits[i],"P");
		mg_hitsMatrix[i]->Add(g_trackHits[i],"P");
		mg_hitsMatrix[i]->Add(g_ledMatrix[i],"P");
		mg_hitsMatrix[i]->Add(g_lowerLimit[i],"L");
		mg_hitsMatrix[i]->Add(g_upperLimit[i],"L");
		mg_hitsMatrix[i]->Add(g_trackLimit[i],"L");
		mg_hitsMatrix[i]->Draw("AP");
		// mg_hitsMatrix[i]->SetBit(kCanDelete);
		g_hits[i]->SetMarkerStyle(20);
		g_trackHits[i]->SetMarkerStyle(20);
		g_trackHits[i]->SetMarkerColor(kBlue);
		g_noiseHits[i]->SetMarkerStyle(20);
		g_noiseHits[i]->SetMarkerColor(kMagenta);
		g_ledMatrix[i]->SetMarkerStyle(20);
		g_ledMatrix[i]->SetMarkerColor(kRed);
		g_lowerLimit[i]->SetLineColor(kGreen);
		g_lowerLimit[i]->SetLineWidth(3);
		g_upperLimit[i]->SetLineColor(kGreen);
		g_upperLimit[i]->SetLineWidth(3);
		g_trackLimit[i]->SetLineColor(kOrange);
		g_trackLimit[i]->SetLineWidth(3);
	}
	cEvent->cd(9);
	g_QvsL->Draw("AP");
	g_QvsL->SetTitle("Charge vs. Distance;Distance from cascade [m]; Charge [p.e.]");
	cEvent->Write();

	delete cEvent;
	delete g_QvsL;

	for (int i = 0; i < gNStrings; ++i)
	{
		delete mg_hitsMatrix[i];
	}
	return 0;
}

int ChargeVisualization(int eventID, TVector3 cascPos, double energy, double theta, double phi)
{
	double tableParameters[4]{0};
	double par[7]{cascPos.X(),cascPos.Y(),cascPos.Z(),0,energy,theta,phi};

	int nHitsPerString[8]{0};
	TGraph* g_MeasQ[gNStrings];
	TGraph* g_DeadOM[gNStrings];
	TGraph* g_ExpQ[gNStrings];
	TMultiGraph* mg_QSum[gNStrings];

	int nOMsPerString = gNOMs/gNStrings;

	for(int k = 0; k < gPulses.size(); k++)
	{
		nHitsPerString[gPulses[k].OMID/36]++;
	}
	for (int i = 0; i < gNStrings; ++i)
	{
		g_MeasQ[i] = new TGraph(nHitsPerString[i]);
		g_ExpQ[i] = new TGraph(nHitsPerString[i]);
		g_DeadOM[i] = new TGraph();
		mg_QSum[i] = new TMultiGraph(Form("mg_%d",i),Form("String_%d; OM ID [#]; Charge [p.e.]",i+1));
		nHitsPerString[i] = 0;
	}
	for(int k = 0; k < gPulses.size(); k++)
	{
		int stringID = gPulses[k].OMID/36;
		g_MeasQ[stringID]->SetPoint(nHitsPerString[stringID],gPulses[k].OMID,gPulses[k].charge);
		GetParameters(par,gPulses[k].OMID,tableParameters);
		gPulses[k].expectedCharge = GetInterpolatedValue(tableParameters)*110000000*energy;
		nHitsPerString[stringID]++;
	}

	for (int i = 0; i < gNOMs; ++i)
	{
		int stringID = i/36;
		GetParameters(par,i,tableParameters);
		g_ExpQ[stringID]->SetPoint(i%36,i,GetInterpolatedValue(tableParameters)*110000000*energy);
		if (gOMqCal[i] == -1)
			g_DeadOM[stringID]->SetPoint(g_DeadOM[stringID]->GetN(),i,0);
	}

	TCanvas *cCharge = new TCanvas(Form("cCharge_%d",eventID),Form("cCharge_%d",eventID),200,10,600,400);
	cCharge->Divide(3,3);
	for (int i = 0; i < gNStrings; ++i)
	{
		cCharge->cd(i+1);
		mg_QSum[i]->Add(g_MeasQ[i],"P");
		mg_QSum[i]->Add(g_ExpQ[i],"P");
		mg_QSum[i]->Add(g_DeadOM[i],"P");
		mg_QSum[i]->Draw("AP");
		g_MeasQ[i]->SetMarkerStyle(20);
		g_MeasQ[i]->SetMarkerColor(kBlue);
		g_ExpQ[i]->SetMarkerStyle(20);
		g_ExpQ[i]->SetMarkerColor(kGreen);
		g_DeadOM[i]->SetMarkerStyle(20);
		g_DeadOM[i]->SetMarkerColor(kRed);

	}
	cCharge->Write();

	delete cCharge;

	for (int i = 0; i < gNStrings; ++i)
	{
		delete mg_QSum[i];
	}
	return 0;
}

void ScanLogLikelihoodEnergy(int eventID,UnifiedEvent &event)
{
	int nPar = 0;
	double* gin = new double(0);
	double likelihoodValue = 0;
	double cascadeParameters[7];
	int iflag = 0;
	int nPoints = 40;

	cascadeParameters[0] = event.position.X();
	cascadeParameters[1] = event.position.Y();
	cascadeParameters[2] = event.position.Z();
	cascadeParameters[3] = event.time;
	cascadeParameters[4] = event.energy;
	cascadeParameters[5] = event.theta;
	cascadeParameters[6] = event.phi;

	TCanvas* c_energyLikelihoodScan = new TCanvas(Form("cEnergyScan_%d",eventID),Form("cEnergyScan_%d",eventID));
	TGraph* g_energyLikelihoodScan = new TGraph();
	TGraph* g_trueEnergy = new TGraph(1);
	g_trueEnergy->SetPoint(0,event.mcEnergy,0);
	g_trueEnergy->SetMarkerStyle(5);
	g_trueEnergy->SetMarkerSize(2);
	g_trueEnergy->SetMarkerColor(kRed);
	g_energyLikelihoodScan->SetTitle("Energy LogLikelihood Scan;Energy [TeV];Likelihood(E_{var})-Likelihood(E_{H0})");
	g_energyLikelihoodScan->SetMarkerStyle(5);
	g_energyLikelihoodScan->SetMarkerSize(2);

	int nRealPoints = 0;
	for (int i = 0; i <= nPoints; ++i)
	{
		cascadeParameters[4] = event.energy + i - 20;
		if (cascadeParameters[4] <= 0)
			continue;
		logLikelihood(nPar,gin,likelihoodValue,cascadeParameters,iflag);
		g_energyLikelihoodScan->SetPoint(nRealPoints,cascadeParameters[4],likelihoodValue-event.likelihood*(gPulses.size()));
		nRealPoints++;
	}
	g_energyLikelihoodScan->Draw("AP");
	g_trueEnergy->Draw("PSAME");
	c_energyLikelihoodScan->Write();

	delete c_energyLikelihoodScan;
}

void ScanLogLikelihoodDirection(int eventID, UnifiedEvent &event)
{
	int nPar = 0;
	double* gin = new double(0);
	double likelihoodValue = 0;
	double cascadeParameters[7];
	int iflag = 0;
	int nPoints = 150;
	double degInRad = 0.001745;

	cascadeParameters[0] = event.position.X();
	cascadeParameters[1] = event.position.Y();
	cascadeParameters[2] = event.position.Z();
	cascadeParameters[3] = event.time;
	cascadeParameters[4] = event.energy;
	cascadeParameters[5] = event.theta;
	cascadeParameters[6] = event.phi;

	TCanvas* c_positionLikelihoodScan = new TCanvas(Form("cDirectionScanZoom_%d",eventID),Form("cDirectionScanZoom_%d",eventID));
	TGraph2D* g_positionLikelihoodScan = new TGraph2D();
	g_positionLikelihoodScan->SetTitle("Direction LogLikelihood Scan;Theta [degree];Phi [degree];Likelihood(theta_{var},phi_{var})-Likelihood(E_{H0})");
	TGraph* g_truePosition = new TGraph(1);
	g_truePosition->SetTitle("true");
	g_truePosition->SetName("true");
	g_truePosition->SetPoint(0,event.mcTheta/TMath::Pi()*180,event.mcPhi/TMath::Pi()*180);
	TGraph* g_fitPosition = new TGraph(1);
	g_fitPosition->SetTitle("fit");
	g_fitPosition->SetName("fit");
	g_fitPosition->SetPoint(0,event.theta/TMath::Pi()*180,event.phi/TMath::Pi()*180);
	// cout << theta/TMath::Pi()*180 << " " << phi/TMath::Pi()*180 << " " << event.mcTheta/TMath::Pi()*180 << " " << event.mcPhi/TMath::Pi()*180 << endl;


	int pointID = 0;
	for (int i = 0; i <= nPoints; ++i)
	{
		for (int j = 0; j <= nPoints; ++j)
		{
			cascadeParameters[5] = event.theta + (i - nPoints/2)*degInRad;
			cascadeParameters[6] = event.phi + (j - nPoints/2)*degInRad;
			logLikelihood(nPar,gin,likelihoodValue,cascadeParameters,iflag);
			// if (likelihoodValue-event.likelihood*gPulses.size() < 0.5)
			{
				g_positionLikelihoodScan->SetPoint(pointID,cascadeParameters[5]/TMath::Pi()*180,cascadeParameters[6]/TMath::Pi()*180,likelihoodValue-event.likelihood*gPulses.size());
				pointID++;
			}
		}
	}
	// g_positionLikelihoodScan->Draw("surf7");
	g_positionLikelihoodScan->Draw("cont1z");
	g_truePosition->SetMarkerStyle(5);
	g_truePosition->SetMarkerSize(5);
	g_truePosition->SetMarkerColor(kRed);
	g_truePosition->Draw("P");
	g_fitPosition->SetMarkerStyle(5);
	g_fitPosition->SetMarkerSize(5);
	g_fitPosition->SetMarkerColor(kGreen);
	g_fitPosition->Draw("same P");
	// g_positionLikelihoodScan->Draw("cont1z");
	c_positionLikelihoodScan->Write();

	delete c_positionLikelihoodScan;
	delete g_truePosition;
	delete g_fitPosition;
}

void ScanLogLikelihoodDirectionCircular(int eventID, UnifiedEvent &event)
{
	int nPar = 0;
	double* gin = new double(0);
	double likelihoodValue = 0;
	double cascadeParameters[7];
	int iflag = 0;
	double degInRad = 0.001745;


	cascadeParameters[0] = event.position.X();
	cascadeParameters[1] = event.position.Y();
	cascadeParameters[2] = event.position.Z();
	cascadeParameters[3] = event.time;
	cascadeParameters[4] = event.energy;
	cascadeParameters[5] = event.theta;
	cascadeParameters[6] = event.phi;

	TCanvas* c_positionLikelihoodScan = new TCanvas(Form("cDirectionScan_%d",eventID),Form("cDirectionScan_%d",eventID));
	TGraph2D* g_positionLikelihoodScan = new TGraph2D();
	g_positionLikelihoodScan->SetTitle("Direction LogLikelihood Scan;Theta [degree];Phi [degree];Likelihood(theta_{var},phi_{var})-Likelihood(E_{H0})");
	TGraph2D* g_truePosition = new TGraph2D(1);
	g_truePosition->SetTitle("true");
	g_truePosition->SetName("true");
	g_truePosition->SetPoint(0,event.mcTheta/TMath::Pi()*180,event.mcPhi/TMath::Pi()*180,0.3);
	TGraph2D* g_fitPosition = new TGraph2D(1);
	g_fitPosition->SetTitle("fit");
	g_fitPosition->SetName("fit");
	g_fitPosition->SetPoint(0,event.theta/TMath::Pi()*180,event.phi/TMath::Pi()*180,0.3);
	// cout << theta/TMath::Pi()*180 << " " << phi/TMath::Pi()*180 << " " << event.mcTheta/TMath::Pi()*180 << " " << event.mcPhi/TMath::Pi()*180 << endl;


	int pointID = 0;

	for (int i = 0; i <= 180; ++i)
	{
		for (int j = 0; j <= 360; ++j)
		{
			cascadeParameters[5] = (i)*TMath::Pi()/180.0;
			cascadeParameters[6] = (j)*2*TMath::Pi()/360.0;
			logLikelihood(nPar,gin,likelihoodValue,cascadeParameters,iflag);
			// if (likelihoodValue-event.likelihood*gPulses.size() < 0.5)
			{
				g_positionLikelihoodScan->SetPoint(pointID,cascadeParameters[5]/TMath::Pi()*180,cascadeParameters[6]/TMath::Pi()*180,likelihoodValue-event.likelihood*gPulses.size());
				pointID++;
			}
		}
	}
	// g_positionLikelihoodScan->Draw("surf7");
	g_positionLikelihoodScan->Draw("surf1");
	g_truePosition->SetMarkerStyle(5);
	g_truePosition->SetMarkerSize(5);
	g_truePosition->SetMarkerColor(kRed);
	g_truePosition->Draw("same p");
	g_fitPosition->SetMarkerStyle(5);
	g_fitPosition->SetMarkerSize(5);
	g_fitPosition->SetMarkerColor(kGreen);
	g_fitPosition->Draw("same p");
	// g_positionLikelihoodScan->Draw("cont1z");
	c_positionLikelihoodScan->Write();

	delete c_positionLikelihoodScan;
	delete g_truePosition;
	delete g_fitPosition;
}

double CalculateDirectionError(UnifiedEvent &event)
{
	int nPar = 0;
	double* gin = new double(0);
	double likelihoodValue = 0;
	double cascadeParameters[7];
	int iflag = 0;

	cascadeParameters[0] = event.position.X();
	cascadeParameters[1] = event.position.Y();
	cascadeParameters[2] = event.position.Z();
	cascadeParameters[3] = event.time;
	cascadeParameters[4] = event.energy;
	cascadeParameters[5] = event.theta;
	cascadeParameters[6] = event.phi;

	double likelihoodLimit = 0.5;
	double minIterLikelihood = 0;
	double stepInDegree = 0.5;
	double step = stepInDegree/180*TMath::Pi();

	int iteration = 1;
	while(minIterLikelihood < likelihoodLimit)
	{
    	minIterLikelihood = 10000;
	    for (int i = -iteration; i <= iteration; ++i)
	    {
	    	for (int j = -iteration; j <= iteration; ++j)
	    	{
	    		if (abs(i)-iteration == 0 || abs(j)-iteration == 0)
    			{
    				cascadeParameters[5] = event.theta + step*i;
					cascadeParameters[6] = event.phi + step*j;
					logLikelihood(nPar,gin,likelihoodValue,cascadeParameters,iflag);
					// cout << i << " " << j << " " << cascadeParameters[5]/TMath::Pi()*180 << " " << cascadeParameters[6]/TMath::Pi()*180 << " "  << likelihoodValue-event.likelihood*gPulses.size() << endl;
					if (likelihoodValue-event.likelihood*gPulses.size() < minIterLikelihood)
					{
						minIterLikelihood = likelihoodValue-event.likelihood*gPulses.size();
					}
    			}
	    	}
	    }
	    // cout << iteration << " " << minIterLikelihood << endl;
	    iteration++;
	    if (iteration == 50)
	    	break;
	}
	return (iteration-1)*stepInDegree;
}

int DoTheMagicUnified(int i, UnifiedEvent &event, EventStats* eventStats)
{
	if (!NFilterPassed(event))
		return -1;
	h_mcEnergy->Fill(event.mcEnergy);
	h_mcTheta->Fill(event.mcTheta);
	h_mcPhi->Fill(event.mcPhi);
	eventStats->nNFilter++;
	CaussalityFilter(event);

	event.nHitsAfterCaus = gPulses.size();
	event.nStringsAfterCaus = GetNStrings();

	if (!SixThreeFilterPassed())
		return -2;
	eventStats->nSixThrees++;

	// TVector3 cascPos(0,0,0);
	// double cascTime = 0;

	EstimateInitialPosMatrix(event.position,event.time);

	event.chi2AfterCaus = FitCascPos(event.position,event.time);
	h_chi2Caus->Fill(event.chi2AfterCaus);
	if (event.chi2AfterCaus > gQCutChi2)
		return -3;
	eventStats->nQFilterChi2++;

	event.nHitsAfterTFilter = TFilter(event,event.position,event.time);

	event.nStringsAfterTFilter = GetNStrings();
	event.mcNTrackHitsAfterTFilter = GetNTrackHits(event);
	h_nHitsTFilter->Fill(event.nHitsAfterTFilter);
	h_nStringsTFilter->Fill(GetNStrings());
	h_nHitsChange->Fill(event.nHitsAfterTFilter-event.nHitsAfterCaus);
	if (event.nHitsAfterTFilter-event.nHitsAfterCaus < 0 || event.nHitsAfterTFilter < gNCutT)
		return -4;
	eventStats->nTFilter++;

	event.chi2AfterTFilter = FitCascPos(event.position,event.time);
	// cout<<"TFilterChi2 "<<event.chi2AfterTFilter<<endl;
	h_chi2TFilter->Fill(event.chi2AfterTFilter);
	if (event.chi2AfterTFilter > gTCutChi2)
		return -5;
	eventStats->nTFilterChi2++;

	h_nHitsTrack->Fill(TrackFilter(event,event.position,event.time));

	event.energy = TMath::Power(10,3.30123+0.0447574*gPulses.size()-0.000135729*gPulses.size()*gPulses.size())/1000;

	if (gUseMultiDirFit)
		LikelihoodFilterPassedGrid(event);
	else
	{
		EstimateInitialDirection(event.position,event.time,event.energy,event.theta,event.phi);
		LikelihoodFilterPassed(event);
	}

	h_likelihood->Fill(event.likelihood);


	if (event.likelihood > gLikelihoodCut)
		return -6;
	eventStats->nLikelihoodFilter++;
	EventVisualization(i,event,event.position,event.time);
	ChargeVisualization(i,event.position,event.energy,event.theta,event.phi);
	event.correctedEnergy = GetCorrectedEnergy(event.energy);
	SaveCascadeJSON(i,event);
	ScanLogLikelihoodEnergy(i,event);
	ScanLogLikelihoodDirectionCircular(i,event);
	ScanLogLikelihoodDirection(i,event);
	event.directionSigma = CalculateDirectionError(event);
	
	//}
	return 0;
}

void InitializeOutputTTree(TTree* outputTree, UnifiedEvent &event)
{
	outputTree->Branch("runID",&event.runID);
	outputTree->Branch("eventID",&event.eventID);
	outputTree->Branch("nHits",&event.nHits);
	outputTree->Branch("nHitsAfterCaus",&event.nHitsAfterCaus);
	outputTree->Branch("nStringsAfterCaus",&event.nStringsAfterCaus);
	outputTree->Branch("chi2AfterCaus",&event.chi2AfterCaus);
	outputTree->Branch("nHitsAfterTFilter",&event.nHitsAfterTFilter);
	outputTree->Branch("nStringsAfterTFilter",&event.nStringsAfterTFilter);
	outputTree->Branch("chi2AfterTFilter",&event.chi2AfterTFilter);
	outputTree->Branch("energy",&event.energy);
	outputTree->Branch("energySigma",&event.energySigma);
	outputTree->Branch("correctedEnergy",&event.correctedEnergy);
	outputTree->Branch("theta",&event.theta);
	outputTree->Branch("thetaSigma",&event.thetaSigma);
	outputTree->Branch("directionSigma",&event.directionSigma);
	outputTree->Branch("phi",&event.phi);
	outputTree->Branch("phiSigma",&event.phiSigma);
	outputTree->Branch("position","TVector3",&event.position);
	outputTree->Branch("time",&event.time);
	outputTree->Branch("likelihood",&event.likelihood);
	outputTree->Branch("mcEnergy",&event.mcEnergy);
	outputTree->Branch("mcTheta",&event.mcTheta);
	outputTree->Branch("mcPhi",&event.mcPhi);
	outputTree->Branch("mcPosition","TVector3",&event.mcPosition);
	outputTree->Branch("qTotal",&event.qTotal);
	outputTree->Branch("mcNTrackHitsAfterTFilter",&event.mcNTrackHitsAfterTFilter);
}

// Main function responsible for the processing of the real data
// required input parameters are seasonID, clusterID, RunID
// joint.events.root files are found automatically based on the enviroment variables
int ProcessExperimentalData()
{
	if (!CheckInputParamsExpData()) // Check input parameters
	{
		return -1;
	}

	if (gProductionID == "") // Sets default productionID (can be set with "-t" switch )
		gProductionID = "barsv051";


	// Sets the file path and checks its existence
    const char* filePath = "";
    if (!gUseNewFolderStructure)
    	filePath = BARS::Data::File(BARS::Data::JOINT, BARS::App::Season, BARS::App::Cluster, BARS::App::Run, gProductionID.c_str());
    else
    	filePath = Form("/eos/baikalgvd/processed/%d/cluster%d/exp/joint/j01/%04d/%s",BARS::App::Season,BARS::App::Cluster,BARS::App::Run,BARS::Data::Filename(BARS::Data::JOINT_MARKED));
    if (!BARS::App::FileExists(filePath))
    {
    	std::cout << "File: " << filePath << " was not found!" << endl;
    	return -2;
    }

	const char* eosName = "root://eos.jinr.ru//";
    std::string buf(gUseEOSRead?eosName:"");
    buf.append(filePath);

    // Sets necessary pointers to access data through TTree
    //TFile* file = new TFile(filePath,"READ");
    TFile* file = new TFile(buf.c_str());
	TTree* tree = (TTree*)file->Get("Events");

	// Sets necessary pointers to access data through TTree
	// TFile* file = new TFile(filePath,"READ");
    // TTree* tree = (TTree*)file->Get("Events");
    if (!tree || !tree->GetBranch("BJointImpulseTel.") || !tree->GetBranch("BJointHeader."))
    {
    	std::cout << "No Events TTree or BJointImpulseTel/BJointHeader branches found!" << endl;
    	return -2;
    }
	BExtractedImpulseTel* impulseTel = NULL;
	BExtractedHeader* header = NULL;
	tree->SetBranchAddress("BJointImpulseTel.",&impulseTel);
	tree->SetBranchAddress("BJointHeader.",&header);

	PrintRunInfo(tree,header);

	if (ReadInputParamFiles(tree,header) == -1)
		return -3;

	TString outputFileName = "";
	if (App::Output == "" || App::Output == "a"){
    	outputFileName = BARS::Data::Directory(BARS::Data::JOINT, BARS::App::Season, BARS::App::Cluster, BARS::App::Run, gProductionID.c_str());
    	cout<<"outputFile "<<outputFileName<<endl;
	}
    else
    	outputFileName =  Form("%s/exp%d/cluster%d/%04d/",App::Output.Data(),BARS::App::Season,BARS::App::Cluster,BARS::App::Run);
	if (gEventID == -1)
		outputFileName += "recCascResults.root";
	else
		outputFileName += Form("singleRecCasc_%d.root",gEventID);
	TFile* outputFile = new TFile(outputFileName,"RECREATE");
	if (!outputFile->IsOpen())
	{
		std::cout << "Output File: " << outputFileName << " can not be recreated!" << endl;
    	return -4;
	}
	TDirectory *cdTree = outputFile->mkdir("Tree");
	UnifiedEvent unifiedEvent;
	TTree* t_RecCasc = new TTree("t_RecCasc","Reconstructed Cascades");
	InitializeOutputTTree(t_RecCasc,unifiedEvent);
	TDirectory *cdHist = outputFile->mkdir("Histograms");
	TDirectory *cdVis = outputFile->mkdir("Visualizations");
   	cdVis->cd();

	EventStats* eventStats = new EventStats();
	eventStats->nEntries = (gNEventsProcessed == -1)?tree->GetEntries():gNEventsProcessed;

	int startEventID = (gEventID == -1)?0:gEventID;
	int endEventID = (gEventID == -1)?eventStats->nEntries:gEventID+1;

	for (int i = startEventID; i < endEventID; ++i)
	{
		if (i%(eventStats->nEntries/10) == 0)
		{
			cout << round((double)(i)/eventStats->nEntries*100) << "% ";
			cout << std::flush;
		}
		tree->GetEntry(i);
		unifiedEvent.runID = BARS::App::Run;
		unifiedEvent.eventID = i;
		TransformToUnifiedEvent(impulseTel,unifiedEvent);
		GenerateNoise(unifiedEvent);
		int status = DoTheMagicUnified(i,unifiedEvent,eventStats);
		if (status == 0)
			t_RecCasc->Fill();
		h_exitStatus->Fill(status);
	}
	cout << endl;

	PrintEventStats(eventStats);
	cdHist->cd();
	SaveHistograms();
	cdTree->cd();
	t_RecCasc->Write();
	delete t_RecCasc;
	outputFile->Close();

    return 0;
}

int ProcessMCCascades()
{
	cout << "Processing MC Cascades Data" << endl;

	const char* filePath = "/media/zuzana/Data/BaikalData/MC_cascadesDZH/MC_cascadeDZH/ne16_tin_c*_00*.root";
	TChain* mcFiles = new TChain("h11");
	mcFiles->Add(filePath);
	if (mcFiles->GetEntries() == 0)
	{
		std::cout << "Files: " << filePath << " were not found!" << endl;
    	return -2;
	}

	// Sets necessary pointers to access data through TChain
	mcCascade* cascade = new mcCascade;
	mcFiles->SetBranchAddress("L0",&cascade->eventID);
	mcFiles->SetBranchAddress("Esh",&cascade->showerEnergy);
	mcFiles->SetBranchAddress("cost",&cascade->cosTheta);
	mcFiles->SetBranchAddress("fj",&cascade->phi);
	mcFiles->SetBranchAddress("jch",&cascade->nHits);
	mcFiles->SetBranchAddress("xtr",cascade->position);
	mcFiles->SetBranchAddress("Npmt",cascade->chID);
	mcFiles->SetBranchAddress("tre",cascade->time);
	mcFiles->SetBranchAddress("are",cascade->charge);

	PrintRunInfoMCCascades(filePath,mcFiles);

	if (ReadInputParamFiles() == -1)
		return -3;

	TString outputFileName = "/media/zuzana/Backup/mcCascadesResults/skuska/";
	outputFileName += "recCascResults.root";
	TFile* outputFile = new TFile(outputFileName,"RECREATE");
	TDirectory *cdTree = outputFile->mkdir("Tree");
	UnifiedEvent unifiedEvent;
	TTree* t_RecCasc = new TTree("t_RecCasc","Reconstructed Cascades");
	InitializeOutputTTree(t_RecCasc,unifiedEvent);
	TDirectory *cdHist = outputFile->mkdir("Histograms");
	TDirectory *cdVis = outputFile->mkdir("Visualizations");
   	cdVis->cd();

	EventStats* eventStats = new EventStats();
	eventStats->nEntries = (gNEventsProcessed == -1)?mcFiles->GetEntries():gNEventsProcessed;

	int nProcessed = 0;

	for (int i = 0; i < eventStats->nEntries; ++i)
	{
		if (i%(eventStats->nEntries/10) == 0)
		{
			cout << round((double)(i)/eventStats->nEntries*100) << "% ";
			cout << std::flush;
		}
		mcFiles->GetEntry(i);
		if (cascade->showerEnergy > 1000 || i % 1000 != 0)
			continue;

		nProcessed++;
		unifiedEvent.eventID = i;
		TransformToUnifiedEvent(cascade,unifiedEvent);
		GenerateNoise(unifiedEvent);
		int status = DoTheMagicUnified(i,unifiedEvent,eventStats);
		if (status == 0)
			t_RecCasc->Fill();
		h_exitStatus->Fill(status);
	}
	cout << endl;
	eventStats->nEntries = nProcessed;

	PrintEventStats(eventStats);
	cdHist->cd();
	SaveHistograms();
	cdTree->cd();
	t_RecCasc->Write();
	outputFile->Close();

    return 0;
}

int ProcessMCData()
{
	cout << "Processing MC Data" << endl ;

	if (!CheckInputParamsMCData()) // Check input parameters
	{
		return -1;
	}

	const char* filePath = " ";

	if (gFileInputFolder == "")
	{
		if (gInputType == 3)
			filePath = "/media/zuzana/Backup/atmBundle/MC_atmMuBundle/n_cors_n2m_cl2016_x*.root";
		if (gInputType == 2)
			filePath = "/media/zuzana/Backup/nuatm_feb19/n_nuatm_gs_n2m_cl2016_x*.root";
	}else
	{
		if (gInputType == 3)
			filePath = Form("%s/n_cors_n2m_cl2016_x*.root",gFileInputFolder.c_str());
		if (gInputType == 2)
			filePath = Form("%s/n_nuatm_gs_n2m_cl2016_x*.root",gFileInputFolder.c_str());
	}


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

    PrintRunInfo(filePath,mcFiles);

    if (ReadInputParamFiles(mcFiles) == -1)
		return -3;

	TString outputFileName = "";
	if (gFileInputFolder == "")
	{
		if (gInputType == 3)
			outputFileName = "/media/zuzana/Backup/atmBundle/normalL3/";
		if (gInputType == 2)
			outputFileName = "/media/zuzana/Backup/nuatm_feb19/newResults/";
	}else
	{
		outputFileName = gFileInputFolder;
	}
	outputFileName += "recCascResults.root";
	TFile* outputFile = new TFile(outputFileName,"RECREATE");
	TDirectory *cdTree = outputFile->mkdir("Tree");
	UnifiedEvent unifiedEvent;
	TTree* t_RecCasc = new TTree("t_RecCasc","Reconstructed Cascades");
	InitializeOutputTTree(t_RecCasc,unifiedEvent);
	TDirectory *cdHist = outputFile->mkdir("Histograms");
	TDirectory *cdVis = outputFile->mkdir("Visualizations");
   	cdVis->cd();

	EventStats* eventStats = new EventStats();
	eventStats->nEntries = (gNEventsProcessed == -1)?mcFiles->GetEntries():gNEventsProcessed;

	for (int i = 0; i < eventStats->nEntries; ++i)
	{
		if (i%(eventStats->nEntries/10) == 0)
		{
			cout << round((double)(i)/eventStats->nEntries*100) << "% ";
			cout << std::flush;
		}
		mcFiles->GetEntry(i);
		unifiedEvent.eventID = i;
		TransformToUnifiedEvent(event,mcEvent,eventMask,unifiedEvent);
		int status = DoTheMagicUnified(i,unifiedEvent,eventStats);
		if (status == 0)
			t_RecCasc->Fill();
		h_exitStatus->Fill(status);
	}
	cout << endl;

	PrintEventStats(eventStats);
	cdHist->cd();
	SaveHistograms();
	cdTree->cd();
	t_RecCasc->Write();
	outputFile->Close();
	return 0;
}

// Fitter setting
void SetFitter(void)
{
	gMinuit = TVirtualFitter::Fitter(0,7); // the second number is number of parameters
	double arg = -1;
	gMinuit->ExecuteCommand("SET PRINTOUT",&arg,1); // these two lines means that it wont be able to print results on the screen
	gMinuit->ExecuteCommand("SET NOW", &arg ,1);
	gMinuit->SetFCN(chi2);
	// fMinuit->SetFCN(MEstimator);
}

// Main processing function
// It sets things that are common for all data input types (like Fitter)
// and starts processing of the given input type
int main(int argc, char** argv)
{
	clock_t begin = clock();
    // Init should be called at the beggining of all BARS programms
    App::Init(argc, argv, 0, parseOpts, readRC, checkParams);

    PrintHeader();
    SetFitter();

    switch (gInputType) {
    	case 0:
    		ProcessExperimentalData();
    		break;
    	case 1:
    		ProcessMCCascades();
    		break;
    	case 2:
    	case 3:
    		ProcessMCData();
    		break;
    	default:
    		break;
    }

    clock_t end = clock();
    cout << endl << "Elapsed time: " << double(end-begin)/CLOCKS_PER_SEC << endl;

    return 0;
}

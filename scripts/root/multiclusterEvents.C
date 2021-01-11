#include <iostream>
#include <fstream>
#include <vector>

#include "TString.h"
#include "TChain.h"
#include "TH1F.h"
#include "TCanvas.h"

#include "BExtractedImpulseTel.h"
#include "BMultiJointHeader.h"

TH1F* h_nClusters = new TH1F("h_nClusters","Number of clusters in the coincidence; N_{clusters} [#]; NoE [#]",10,0,10);
TH1F* h_clusterIDs = new TH1F("h_clusterIDs","Cluster IDs in the coincidence; Cluster ID [1]; NoE [#]",10,0,10);
TH1F* h_eventIDs = new TH1F("h_eventIDs","Events IDs in the coincidence; Event ID [1]; NoE [#]",1000,0,10000000);

struct InputCascade
{
	int season;
	int cluster;
	int run;
	int event;
	double theta;
	double phi;
	double x;
	double y;
	double z;
	int nHits;
	double energy;
	double likelihood;

	// equality comparison. doesn't modify object. therefore const.
    bool operator==(const InputCascade& a) const
    {
        return (season == a.season && cluster == a.cluster && run == a.run && event == a.event);
        // return (season == a.season && cluster == a.cluster && run == a.run);
    }
};

vector<InputCascade> inputCascades;

const int nClusters = 5;
const int nStringsPerCluster = 8;
const int nOMsPerCluster = 288;
double xPos[nClusters*nStringsPerCluster] = {-13.76,32.14,45.06,5.13,-45.03,-76.21,-59.85,-14.47,-195.19,-164.79,-180.08,-227.51,-276.24,-279.59,-248.17,-222.70,-270.25,-228.58,-220.89,-261.89,-309.86,-337.48,-319.74,-282.27,65.85,108.73,113.87,74.19,25.1,-2.48,16.08,58.37,-163.91,-119.26,-113.90,-152.28,-202.59,-230.83,-213.25,-170.30};
double yPos[nClusters*nStringsPerCluster] = {-211.35,-235.88,-285.45,-325.83,-319.82,-281.63,-231.37,-270.17,-340.62,-384.09,-435.13,-450.13,-424.31,-372.59,-337.03,-391.09,-37.36,-65.26,-117.78,-153.57,-146.26,-101.43,-55.24,-96.82,-435.47,-462.39,-514.68,-549.90,-544.25,-500.53,-453,-491.97,-628.26,-656.49,-707.52,-744.24,-738.58,-694.13,-645.06,-685.35};

void PrintCascade(InputCascade &cascade)
{
	cout << "Season: " << cascade.season << " Cluster: " << cascade.cluster << " Run: " << cascade.run << " Event: " << cascade.event << endl;
}

int DrawResults()
{
	TCanvas* c_nClusters = new TCanvas("c_nClusters","NClusters",800,600);
	h_nClusters->Draw();

	TCanvas* c_clusterIDs = new TCanvas("c_clusterID","ClusterIDs",800,600);
	h_clusterIDs->Draw();

	TCanvas* c_eventIDs = new TCanvas("c_eventIDs","EventIDs",800,600);
	h_eventIDs->Draw();

	return 0;
}

int ReadInputCascades(TString fileName)
{
	ifstream inputFile;
    inputFile.open(fileName);

    if (!inputFile)
    {
    	cerr << " File: " << fileName << " with reconstructed cascades was NOT found. Program termination!" << endl;
    	return -1;
  	}

  	int season, cluster, run, event, nHits = 0;
  	double theta, phi, x, y, z, energy, likelihood = 0;

  	while(!inputFile.eof())
  	{
  		inputFile >> season >> cluster >> run >> event >> theta >> phi >> x >> y >> z >> nHits >> energy >> likelihood;
  		if (inputFile.eof())
  			break;
  		inputCascades.push_back(InputCascade{season,cluster,run,event,theta,phi,x,y,z,nHits,energy,likelihood});
  		PrintCascade(inputCascades.back());
  	}

	inputFile.close();

	cout << "N = " << inputCascades.size() << " has been read!" << endl;

	return 0;

}

int PrintHits(BExtractedImpulseTel* impulseTel)
{
	cout << "NImpulses: " << impulseTel->GetNimpulse() << endl;
	for (int i = 0; i < impulseTel->GetNimpulse(); ++i)
	{
		cout << i << "\t" << impulseTel->GetNch(i) << "\t" << impulseTel->GetNch(i)/288 << "\t" << impulseTel->GetQ(i) << "\t" << impulseTel->GetT(i) << endl;
	}
	return 0;
}

int SaveGeometry(std::ofstream &fOutputFile, int OMID)
{
	fOutputFile << "\t\t{" << std::endl;
	fOutputFile << "\t\t\t\"cluster\": " << OMID/288 << "," << std::endl;
	fOutputFile << "\t\t\t\"string\": " << (OMID%288)/36 << "," << std::endl;
	fOutputFile << "\t\t\t\"channelID\": " << OMID << "," << std::endl;
	fOutputFile << "\t\t\t\"x\": " << xPos[(OMID/288)*nStringsPerCluster+(OMID%288)/36] << "," << std::endl;
	fOutputFile << "\t\t\t\"y\": " << yPos[(OMID/288)*nStringsPerCluster+(OMID%288)/36] << "," << std::endl;
	fOutputFile << "\t\t\t\"z\": " << ((OMID%288)%36)*15 << std::endl;
	if (OMID != nOMsPerCluster*nClusters-1)
		fOutputFile << "\t\t}," << std::endl;
	else
		fOutputFile << "\t\t}" << std::endl;

	return 0;
}

int SavePulse(std::ofstream &fOutputFile, int pulseID, int maxPulseID, BExtractedImpulseTel* impulseTel)
{
	fOutputFile << "\t\t{" << std::endl;
	fOutputFile << "\t\t\t\"amplitude\": " << impulseTel->GetA(pulseID)/25 << "," << std::endl;
	fOutputFile << "\t\t\t\"channelID\": " << impulseTel->GetNch(pulseID) << "," << std::endl;
	fOutputFile << "\t\t\t\"charge\": " << impulseTel->GetQ(pulseID)/150 << "," << std::endl;
	if (impulseTel->GetQ(pulseID)/150 > 1.5)
		fOutputFile << "\t\t\t\"mask\": " << 1 << "," << std::endl;
	else
		fOutputFile << "\t\t\t\"mask\": " << 0 << "," << std::endl;
	fOutputFile << "\t\t\t\"time\": " << impulseTel->GetT(pulseID) << std::endl;
	if (pulseID != maxPulseID-1)
		fOutputFile << "\t\t}," << std::endl;
	else
		fOutputFile << "\t\t}" << std::endl;

	return 0;
}

int SaveJSON(BExtractedImpulseTel* impulseTel, int season, int cluster, int run, int event, int cascadeID)
{
	TString jsonFile = Form("../../results/s%d_c%d_r%d_evt%d.json",season,cluster,run,event);
	std::ofstream fOutputFile;
	fOutputFile.open(jsonFile);

	fOutputFile << "{" << std::endl;

	fOutputFile << "\t\"season\": " << season << "," <<std::endl;
	fOutputFile << "\t\"cluster\": " << cluster << "," <<std::endl;
	fOutputFile << "\t\"run\": " << run << "," <<std::endl;
	fOutputFile << "\t\"eventID\": " << event << "," <<std::endl;
	fOutputFile << "\t\"geometry\":[" <<std::endl;
	for (int i = 0; i < nOMsPerCluster*nClusters; ++i)
	{
		SaveGeometry(fOutputFile,i);
	}
	fOutputFile << "\t]," <<std::endl;
	fOutputFile << "\t\"pulses\":[" <<std::endl;
	for (int i = 0; i < impulseTel->GetNimpulse(); ++i)
	{
		SavePulse(fOutputFile,i,impulseTel->GetNimpulse(),impulseTel);
	}
	fOutputFile << "\t]," <<std::endl;
	fOutputFile<<"\t\"origins\": {"<<std::endl;
	fOutputFile<<"\t\t\"cascades\": [{"<<std::endl;
	fOutputFile<<"\t\t\t\"mc\": false,"<<std::endl;
	fOutputFile<<"\t\t\t\"title\": \"cascadeFit\","<<std::endl;
	fOutputFile<<"\t\t\t\"direction\": {"<<std::endl;
	// fOutputFile<<"\t\t\t\t\"theta\": "<<std::right<<TMath::Pi()-inputCascades[cascadeID].theta<<","<<std::endl;
	fOutputFile<<"\t\t\t\t\"theta\": "<<std::right<<inputCascades[cascadeID].theta<<","<<std::endl;
	// fOutputFile<<"\t\t\t\t\"phi\": "<<std::right<<TMath::Pi()+inputCascades[cascadeID].phi<<","<<std::endl;
	fOutputFile<<"\t\t\t\t\"phi\": "<<std::right<<inputCascades[cascadeID].phi<<","<<std::endl;
	fOutputFile<<"\t\t\t\t\"x\": "<<std::right<<inputCascades[cascadeID].x+xPos[8*cluster-1]<<","<<std::endl;
	fOutputFile<<"\t\t\t\t\"y\": "<<std::right<<inputCascades[cascadeID].y+yPos[8*cluster-1]<<","<<std::endl;
	fOutputFile<<"\t\t\t\t\"z\": "<<std::right<<inputCascades[cascadeID].z+262.5<<","<<std::endl;
	fOutputFile<<"\t\t\t\t\"time\": "<<std::right<<0<<""<<std::endl;
	// fOutputFile<<"\t\t\t\t\"time\": "<<std::right<<time<<std::endl;
	fOutputFile<<"\t\t\t}"<<std::endl;
	fOutputFile<<"\t\t}]"<<std::endl;
	fOutputFile<<"\t}"<<std::endl;
	fOutputFile << "}" << std::endl;

	fOutputFile.close();
	return 0;
}

int multiclusterEvents()
{
	TString filePath = "/Data/BaikalData/multicluster/imulticluster.*.root";
	// TString filePath = "/media/fajtak/Alpha/BaikalData/multicluster/imulticluster.*.root";
	// TString inputCascadesFile = "../../results/recCasc.txt";
	TString inputCascadesFile = "../../results/recCasc_y19c-1.txt";

	if (ReadInputCascades(inputCascadesFile) != 0)
		return -1;

	TChain* multiclusterFiles = new TChain("Events");
	multiclusterFiles->Add(filePath);
	if (multiclusterFiles->GetEntries() == 0)
	{
		std::cout << "Files: " << filePath << " were not found!" << endl;
    	return -2;
	}

	BExtractedImpulseTel* impulseTel = NULL;
    multiclusterFiles->SetBranchAddress("BExtractedImpulseTel",&impulseTel);
    BMultiJointHeader* jointHeader = NULL;
    multiclusterFiles->SetBranchAddress("BMultiJointHeader",&jointHeader);

    cout << "Number of entries: " << multiclusterFiles->GetEntries() << endl;

    int nMulticlusterEvents = 0;
    int n3clusterEvents = 0;
    int nBigOnes = 0;

    for (int i = 0; i < multiclusterFiles->GetEntries(); ++i)
    // for (int i = 0; i < 100000; ++i)
    {
    	multiclusterFiles->GetEntry(i);

    	h_nClusters->Fill(jointHeader->GetClusters());

    	for (int j = 0; j < jointHeader->GetClusters(); ++j)
    	{
    		h_clusterIDs->Fill(jointHeader->GetCluster(j));
    		h_eventIDs->Fill(jointHeader->GetEventIDCC(j));
    		InputCascade readEvent{jointHeader->GetSeason(j),jointHeader->GetCluster(j),jointHeader->GetRun(j),(int)jointHeader->GetEventIDCC(j)};
    		// PrintCascade(readEvent);
   //  		if (jointHeader->GetClusters() > 4)
			// {
			// 	cout << "FOUND BIG ONE!!! " << jointHeader->GetClusters() << endl;
			// 	PrintCascade(readEvent);
			// 	SaveJSON(impulseTel,jointHeader->GetSeason(j),jointHeader->GetCluster(j),jointHeader->GetRun(j),(int)jointHeader->GetEventIDCC(j),0);
			// 	nBigOnes++;
			// }
    		for (int k = 0; k < inputCascades.size(); ++k)
    		{
    			if (inputCascades[k] == readEvent)
    			{
    				cout << "FOUND!!! " << jointHeader->GetClusters() << endl;
    				if (jointHeader->GetClusters() == 3)
    					n3clusterEvents++;
    				PrintCascade(inputCascades[k]);
    				PrintCascade(readEvent);
    				// PrintHits(impulseTel);
    				SaveJSON(impulseTel,jointHeader->GetSeason(j),jointHeader->GetCluster(j),jointHeader->GetRun(j),(int)jointHeader->GetEventIDCC(j),k);
    				nMulticlusterEvents++;
    			}
    		}
    	}
    	// cout << jointHeader->GetClusters() << endl;
    }

    DrawResults();
    cout << "Number of multicluster events: " << nMulticlusterEvents << endl;
    cout << "Number of three cluster events: " << n3clusterEvents << endl;
    cout << "Number of big events: " << nBigOnes << endl;

    return 0;
}
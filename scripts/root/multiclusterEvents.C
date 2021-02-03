#include <iostream>
#include <fstream>
#include <vector>

#include "TString.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TFile.h"

#include "BExtractedImpulseTel.h"
#include "BMultiJointHeader.h"


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
	bool passesCluster;
	int expCoincSize;
	bool coincidenceFound;
	int coincidenceSize;

	// equality comparison. doesn't modify object. therefore const.
    bool operator==(const InputCascade& a) const
    {
        return (season == a.season && cluster == a.cluster && run == a.run && event == a.event);
        // return (season == a.season && cluster == a.cluster && run == a.run);
    }
};

vector<InputCascade> inputCascades;

const int nClusters = 7;
const int nStringsPerCluster = 8;
const int nOMsPerCluster = 288;
double xPos[nClusters*nStringsPerCluster] = {-13.76,32.14,45.06,5.13,-45.03,-76.21,-59.85,-14.47,-195.19,-164.79,-180.08,-227.51,-276.24,-279.59,-248.17,-222.70,-270.25,-228.58,-220.89,-261.89,-309.86,-337.48,-319.74,-282.27,65.85,108.73,113.87,74.19,25.1,-2.48,16.08,58.37,-163.91,-119.26,-113.90,-152.28,-202.59,-230.83,-213.25,-170.30,-500.14,-454.30,-430.70,-456.29,-507.00,-546.56,-544.70,-488.82,-494.67,-455.80,-451.87,-485.72,-538.53,-564.08,-546.77,-503.92};
double yPos[nClusters*nStringsPerCluster] = {-211.35,-235.88,-285.45,-325.83,-319.82,-281.63,-231.37,-270.17,-340.62,-384.09,-435.13,-450.13,-424.31,-372.59,-337.03,-391.09,-37.36,-65.26,-117.78,-153.57,-146.26,-101.43,-55.24,-96.82,-435.47,-462.39,-514.68,-549.90,-544.25,-500.53,-453,-491.97,-628.26,-656.49,-707.52,-744.24,-738.58,-694.13,-645.06,-685.35,-466.59,-479.59,-526.30,-571.16,-581.96,-547.72,-496.46,-525.08,-181.00,-208.44,-260.83,-296.66,-288.98,-243.64,-195.04,-240.24};

TH1F* h_nClusters = new TH1F("h_nClusters","Number of clusters in the coincidence; N_{clusters} [#]; NoE [#]",10,0,10);
TH1F* h_clusterIDs = new TH1F("h_clusterIDs","Cluster IDs in the coincidence; Cluster ID [1]; NoE [#]",10,0,10);
TH1F* h_eventIDs = new TH1F("h_eventIDs","Events IDs in the coincidence; Event ID [1]; NoE [#]",1000,0,10000000);
TH2F* h_coincidences = new TH2F("h_coincidences","Coincidences;Cluster ID [1]; Cluster ID [2]",10,0,10,10,0,10);
TH1F* h_coincMultiplicities = new TH1F("h_coincMultiplicities","Multiplicities;N_{clusters} [#]; NoE [#]",10,0,10);

void PrintCascade(InputCascade &cascade)
{
	cout << "Size: " << cascade.coincidenceSize << "/" << cascade.expCoincSize << " Season: " << cascade.season << " Cluster: " << cascade.cluster << " Run: " << cascade.run << " Event: " << cascade.event << "NHits: " << cascade.nHits << " Theta: " << cascade.theta << "/" << cascade.theta/TMath::Pi()*180 << " Phi: " << cascade.phi << "/" << cascade.phi/TMath::Pi()*180 << " X: " << cascade.x << " Y: " << cascade.y << " Z: " << cascade.z  << endl;
}

int DrawResults()
{
	TCanvas* c_nClusters = new TCanvas("c_nClusters","NClusters",800,600);
	h_nClusters->Draw();

	TCanvas* c_clusterIDs = new TCanvas("c_clusterID","ClusterIDs",800,600);
	h_clusterIDs->Draw();

	TCanvas* c_eventIDs = new TCanvas("c_eventIDs","EventIDs",800,600);
	h_eventIDs->Draw();

	TCanvas* c_coincidences = new TCanvas("c_coincidences","Coincidences",800,600);
	h_coincidences->Draw("colz");

	TCanvas* c_coincMultiplicities = new TCanvas("c_coincMultiplicities","Multiplicities",800,600);
	h_coincMultiplicities->Draw();

	return 0;
}

void SaveResults(int season)
{
	TString outputFileName = Form("../../results/multicluster_%d.root",season);
	TFile* outputFile = new TFile(outputFileName,"RECREATE");

	h_eventIDs->Write();
	h_clusterIDs->Write();
	h_nClusters->Write();
	h_coincidences->Write();
	h_coincMultiplicities->Write();
}

void PrintResults()
{
	cout << endl;
	cout << std::string(80,'*') << endl;
	cout << "Results: " << endl;
	for (int i = 3; i < h_coincMultiplicities->GetNbinsX(); ++i)
	{
		cout << i-1 << "\t" << h_coincMultiplicities->GetBinContent(i) << endl;
	}

	int nCoincidences = 0;
	cout << "Coincidences: " << endl;
	for (int i = 0; i < inputCascades.size(); ++i)
	{
		if (inputCascades[i].coincidenceFound)
		{
			nCoincidences++;
			PrintCascade(inputCascades[i]);
		}
	}
	cout << "Ncoincidences: " << nCoincidences << endl;

	int nUpgoing = 0;
	cout << "Up-going: " << endl;
	for (int i = 0; i < inputCascades.size(); ++i)
	{
		if (inputCascades[i].coincidenceFound && inputCascades[i].theta < TMath::Pi()/2)
		{
			nUpgoing++;
			PrintCascade(inputCascades[i]);
		}
	}
	cout << "Nupgoing: " << nUpgoing << endl;

	int nNeutrinos = 0;
	cout << "Neutrinos: " << endl;
	for (int i = 0; i < inputCascades.size(); ++i)
	{
		if (inputCascades[i].passesCluster && !inputCascades[i].coincidenceFound)
		{
			nNeutrinos++;
			PrintCascade(inputCascades[i]);
		}
	}
	cout << "Nneutrinos: " << nNeutrinos << endl;

	// cout << "Strange: " << endl;
	// for (int i = 0; i < inputCascades.size(); ++i)
	// {
	// 	if (!inputCascades[i].passesCluster && inputCascades[i].coincidenceFound)
	// 		PrintCascade(inputCascades[i]);
	// }
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
  		// PrintCascade(inputCascades.back());
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

int SaveCoincidences(int cascadeID, BMultiJointHeader* jointHeader,ofstream &multiCoincData)
{
	cout << "Saving" << endl;
	for (int i = 0; i < jointHeader->GetClusters(); ++i)
	{
		multiCoincData << cascadeID << "\t" << i << "\t" << jointHeader->GetSeason(i) << "\t" << jointHeader->GetCluster(i) << "\t" << jointHeader->GetRun(i) << "\t" << jointHeader->GetEventIDCC(i) << endl;
	}
	return 0;
}

int HorizontalIntersectionExists(double x0, double y0, double ux, double uy, double a0, double b0, double R)
{
	// cout<< x0 << "\t" << y0 << "\t" << ux << "\t" << uy << "\t" << a0 << "\t" << b0 << "\t" << R << endl;
	double a = TMath::Power(ux,2)+TMath::Power(uy,2);
	double b = (2)*((x0-a0)*ux + (y0-b0)*uy);
	double c = TMath::Power(x0-a0,2)+TMath::Power(y0-b0,2)-TMath::Power(R,2);
	double discriminant = TMath::Power(b,2)-4*a*c;
	double t1 = (-b+TMath::Sqrt(discriminant))/2/a;
	double t2 = (-b-TMath::Sqrt(discriminant))/2/a;
	double x1 = x0+t1*ux;
	double x2 = x0+t2*ux;
	double y1 = y0+t1*uy;
	double y2 = y0+t2*uy;
	if (discriminant < 0)
	{
		return 0;
	}
	else
	{
		// cout<< x0 << "\t" << y0 << "\t" << ux << "\t" << uy << "\t" << a0 << "\t" << b0 << "\t" << R << endl;
		// cout << a << "\t" << b << "\t" << c << "\t" << discriminant << endl;
		// cout << t1 << "\t" << t2 << "\t" << x1 << "\t" << y1 << "\t" << x2 << "\t" << y2 << endl;
		if (t1 > 0)
			return 1;
		else
			return -1;
	}
}

bool VerticalIntersectionExists(double z0, double theta, double distance, int returnValue)
{
	double z = 0;
	if (returnValue > 0)
	{
		z = z0 + distance/TMath::Tan(theta);
	}else{
		z = z0 + distance/TMath::Tan(TMath::Pi()-theta);
	}
	// cout << z0 << "\t" << z << "\t" << theta << "\t" << distance << "\t" << returnValue << endl;
	// if (z > -10 && z < 535)
	if (z > -10 && z < 535)
		return true;
	else
		return false;
}

bool PassesThroughOtherCluster(int cascadeID)
{
	bool passes = false;
	int returnValue = 0;
	int clusterUnderStudy = inputCascades[cascadeID].cluster;
	double x0 = inputCascades[cascadeID].x+xPos[8*clusterUnderStudy-1];
	double y0 = inputCascades[cascadeID].y+yPos[8*clusterUnderStudy-1];
	double z0 = inputCascades[cascadeID].z+262.5;
	double ux = TMath::Cos(inputCascades[cascadeID].phi);
	double uy = TMath::Sin(inputCascades[cascadeID].phi);
	int nPasses = 0;
	for (int i = 0; i < nClusters; ++i)
	{
		if (i+1 == clusterUnderStudy)
			continue;
		double a0 = xPos[8*(i+1)-1];
		double b0 = yPos[8*(i+1)-1];
		double R = 70;
		double dist = TMath::Sqrt(TMath::Power(x0-a0,2)+TMath::Power(y0-b0,2));
		// cout<< x0 << "\t" << y0 << "\t" << ux << "\t" << uy << "\t" << a0 << "\t" << b0 << "\t" << R << endl;

		if ((returnValue = HorizontalIntersectionExists(x0,y0,ux,uy,a0,b0,R)) != 0)
		{
			if (VerticalIntersectionExists(z0,inputCascades[cascadeID].theta,dist,returnValue))
			{
				// cout << i+1 << "\t" << returnValue << endl;
				nPasses++;
				passes = true;
			}
		}
	}
	inputCascades[cascadeID].expCoincSize = nPasses+1;
	return passes;
}


int multiclusterEvents(int season, TString multiclusterPath, TString inputCascadesFile)
{
	TString filePath = Form("%s/*multicluster.*.root",multiclusterPath.Data());
	// TString filePath = "/media/fajtak/Alpha/BaikalData/multicluster/imulticluster.*.root";
	// TString inputCascadesFile = "../../results/recCasc.txt";
	// TString inputCascadesFile = "../../results/recCasc_y19c-1.txt";
	// TString inputCascadesFile = "../../results/recCasc_y19c-1_Zuzka.txt";
	// TString outputCoincidenceFile = "../../results/multiCoinc_y19c-1.txt";

	// ofstream multiCoincData;
	// multiCoincData.open(outputCoincidenceFile);

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

    int nPassingCascasdes = 0;
    for (unsigned int i = 0; i < inputCascades.size(); ++i)
    {
    	// PrintCascade(inputCascades[i]);
    	bool studyMore = PassesThroughOtherCluster(i);
    	if (studyMore)
    	{
    		inputCascades[i].passesCluster = true;
    		nPassingCascasdes++;
    	}
    }

    for (int i = 0; i < multiclusterFiles->GetEntries(); ++i)
    {
    	if (multiclusterFiles->GetEntries() > 10 && i%(multiclusterFiles->GetEntries()/10) == 0)
		{
			cout << round((double)(i)/multiclusterFiles->GetEntries()*100) << "% ";
			cout << std::flush;
		}
    	multiclusterFiles->GetEntry(i);

    	h_nClusters->Fill(jointHeader->GetClusters());

    	for (int j = 0; j < jointHeader->GetClusters(); ++j)
    	{
    		h_clusterIDs->Fill(jointHeader->GetCluster(j));
    		h_eventIDs->Fill(jointHeader->GetEventIDCC(j));
    		InputCascade readEvent{jointHeader->GetSeason(j),jointHeader->GetCluster(j),jointHeader->GetRun(j),(int)jointHeader->GetEventIDCC(j),0,0,0,0,0,0,0,0,false,false};
    		// PrintCascade(readEvent);
   //  		if (jointHeader->GetClusters() > 4)
			// {
			// 	cout << "FOUND BIG ONE!!! " << jointHeader->GetClusters() << endl;
			// 	PrintCascade(readEvent);
			// 	SaveJSON(impulseTel,jointHeader->GetSeason(j),jointHeader->GetCluster(j),jointHeader->GetRun(j),(int)jointHeader->GetEventIDCC(j),0);
			// 	nBigOnes++;
			// }
    		for (unsigned int k = 0; k < inputCascades.size(); ++k)
    		{
    			if (inputCascades[k] == readEvent)
    			{
    				// cout << "FOUND!!! " << jointHeader->GetClusters() << endl;
    				for (int l = 0; l < jointHeader->GetClusters(); ++l)
    				{
    					if (jointHeader->GetCluster(l) != jointHeader->GetCluster(j) )
	    					h_coincidences->Fill(jointHeader->GetCluster(j),jointHeader->GetCluster(l));
    				}
    				h_coincMultiplicities->Fill(jointHeader->GetClusters());
    				// PrintCascade(inputCascades[k]);
    				// PrintCascade(readEvent);
    				inputCascades[k].coincidenceFound = true;
    				inputCascades[k].coincidenceSize = jointHeader->GetClusters();
    				SaveJSON(impulseTel,jointHeader->GetSeason(j),jointHeader->GetCluster(j),jointHeader->GetRun(j),(int)jointHeader->GetEventIDCC(j),k);
    				// SaveCoincidences(k,jointHeader,multiCoincData);
    			}
    		}
    	}
    	// cout << jointHeader->GetClusters() << endl;
    }

    DrawResults();
    SaveResults(season);
    PrintResults();
    // multiCoincData.close();

    return 0;
}
#include <iostream>
#include <vector>
#include <fstream>
#include <string>

#include "TString.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"

using namespace std;

TH1F* h_energyDif = new TH1F("h_energyDif","Energy difference;E_{zh} - E_{luk} [TeV];NoE [#]",1000,-500,500);
TH1F* h_thetaDif = new TH1F("h_thetaDif","Theta difference;#theta_{zh} - #theta_{luk} [deg.];NoE [#]",360,-180,180);
TH1F* h_phiDif = new TH1F("h_phiDif","Phi difference;#phi_{zh} - #phi_{luk} [deg.];NoE [#]",360,-180,180);
TH1F* h_decDif = new TH1F("h_decDif","Declination difference;#dec_{zh} - #dec_{luk} [deg.];NoE [#]",360,-180,180);
TH1F* h_raDif = new TH1F("h_raDif","Right Ascension difference;#Ra_{zh} - #Ra_{luk} [deg.];NoE [#]",360,-180,180);
TH2F* h_energyScat = new TH2F("h_energyScat","Energy scatter;E_{zh}[TeV];E_{luk} [Tev]; NoE [#]",100,0,1000,100,0,1000);
TH1F* h_nHitsAll = new TH1F("h_nHitsAll","NHits Zhan-Arys;NHits_{zh} [#];NoE [#]",50,0,50);
TH1F* h_nHitsCoinc = new TH1F("h_nHitsCoinc","NHits Zhan-Arys in Coincidence;NHits_{zh} [#];NoE [#]",50,0,50);
TH1F* h_xDif = new TH1F("h_xDif","X difference;x_{zh} - x_{luk} [m];NoE [#]",50,-25,25);
TH1F* h_yDif = new TH1F("h_yDif","Y difference;y_{zh} - y_{luk} [m];NoE [#]",50,-25,25);
TH1F* h_zDif = new TH1F("h_zDif","Z difference;z_{zh} - z_{luk} [m];NoE [#]",50,-25,25);


struct InputCascade
{
	int season;
	int cluster;
	int run;
	int event;
	double theta;
	double phi;
	double declination;
	double ra;
	double x;
	double y;
	double z;
	int nHits;
	double energy;
	bool inCoincidence;

	// equality comparison. doesn't modify object. therefore const.
    bool operator==(const InputCascade& a) const
    {
        return (season == a.season && cluster == a.cluster && run == a.run && event == a.event);
        // return (season == a.season && cluster == a.cluster && run == a.run);
    }
};

void PrintCascade(InputCascade &cascade)
{
	cout << "Season: " << cascade.season << " Cluster: " << cascade.cluster << " Run: " << cascade.run << " Event: " << cascade.event << " NHits: " << cascade.nHits << " Theta: " << cascade.theta << "/" << cascade.theta/TMath::Pi()*180 << " Phi: " << cascade.phi << "/" << cascade.phi/TMath::Pi()*180 << " X: " << cascade.x << " Y: " << cascade.y << " Z: " << cascade.z << " E: " << cascade.energy  << endl;
}

int DrawResults()
{
	TCanvas* c_energyDif = new TCanvas("c_energyDif","EnergyDifference",800,600);
	h_energyDif->Draw();

	TCanvas* c_energyScat = new TCanvas("c_energyScat","EnergyScatter",800,600);
	h_energyScat->Draw("colz");

	TCanvas* c_thetaDif = new TCanvas("c_thetaDif","ThetaDifference",800,600);
	h_thetaDif->Draw();

	TCanvas* c_phiDif = new TCanvas("c_phiDif","PhiDifference",800,600);
	h_phiDif->Draw();

	TCanvas* c_decDif = new TCanvas("c_decDif","DeclinationDifference",800,600);
	h_decDif->Draw();

	TCanvas* c_raDif = new TCanvas("c_raDif","RightAscensionDifference",800,600);
	h_raDif->Draw();

	TCanvas* c_nHitsZhan = new TCanvas("c_nHitsZhan","NHitsZhanArys",800,600);
	h_nHitsAll->Draw();
	h_nHitsCoinc->SetLineColor(kRed);
	h_nHitsCoinc->Draw("same");

	TCanvas* c_xDif = new TCanvas("c_xDif","XDifference",800,600);
	h_xDif->Draw();

	TCanvas* c_yDif = new TCanvas("c_yDif","YDifference",800,600);
	h_yDif->Draw();

	TCanvas* c_zDif = new TCanvas("c_zDif","ZDifference",800,600);
	h_zDif->Draw();

	return 0;
}

vector<InputCascade> ReadOlgasCascades(TString fileName, int season)
{
	vector<InputCascade> cascades;
	ifstream inputFile;
    inputFile.open(fileName);

    if (!inputFile)
    {
    	cerr << " File: " << fileName << " with reconstructed cascades was NOT found. Program termination!" << endl;
    	return cascades;
  	}

  	int cluster, run, event, nHits = 0;
  	double theta, phi, x, y, z, energy, likelihood = 0;
  	long tastro = 0;
  	double datjul,dist,dec,ra = 0;

  	string dummy;

  	getline(inputFile,dummy);

  	while(!inputFile.eof())
  	{
  		inputFile >> cluster >> run >> event >> tastro >> datjul >> x >> y >> z >> theta >> phi >> dec >> ra >> nHits >> energy;
  		if (inputFile.eof())
  			break;
  		// phi = (90-phi) > 0? (90-phi) : 360+(90-phi);
  		phi = (phi-90) > 0? (phi-90) : 360+(phi-90);
  		cascades.push_back(InputCascade{season,cluster,run,event-1,(180-theta)/180*TMath::Pi(),(phi)/180*TMath::Pi(),dec/180*TMath::Pi(),ra/180*TMath::Pi(),-y,x,z,nHits,energy,false});
  		// PrintCascade(cascades.back());
  	}

	inputFile.close();

	cout << "N = " << cascades.size() << " has been read!" << endl;

	return cascades;
}

vector<InputCascade> ReadLukasCascades(TString fileName)
{
	vector<InputCascade> cascades;
	ifstream inputFile;
    inputFile.open(fileName);

    if (!inputFile)
    {
    	cerr << " File: " << fileName << " with reconstructed cascades was NOT found. Program termination!" << endl;
    	return cascades;
  	}

  	int season, cluster, run, event, nHits = 0;
  	double theta, phi, x, y, z, energy, likelihood = 0;
  	long tastro = 0;
  	double datjul,dist,dec,ra = 0;

  	while(!inputFile.eof())
  	{
  		inputFile >> season >> cluster >> run >> event >> theta >> phi >> dec >> ra >> x >> y >> z >> nHits >> energy >> likelihood;
  		if (inputFile.eof())
  			break;
  		cascades.push_back(InputCascade{season,cluster,run,event,theta,phi,dec,ra,x,y,z+15,nHits,energy,false});
  		// PrintCascade(cascades.back());
  	}

	inputFile.close();

	cout << "N = " << cascades.size() << " has been read!" << endl;

	return cascades;
}

int CompareDatasets(vector<InputCascade> &largerSet, vector<InputCascade> &smallerSet)
{
	int nFoundCoincidences = 0;
	for (unsigned int i = 0; i < smallerSet.size(); ++i)
	{
		h_nHitsAll->Fill(smallerSet[i].nHits);
		// PrintCascade(smallerSet[i]);
		for (unsigned int j = 0; j < largerSet.size(); ++j)
		{
			if (smallerSet[i] == largerSet[j])
			{
				PrintCascade(smallerSet[i]);
				PrintCascade(largerSet[j]);
				nFoundCoincidences++;
				// cout << "Coincidence found" << endl;
				h_energyDif->Fill(smallerSet[i].energy - largerSet[j].energy);
				h_energyScat->Fill(smallerSet[i].energy,largerSet[j].energy);
				h_thetaDif->Fill((smallerSet[i].theta-largerSet[j].theta)/TMath::Pi()*180);
				h_phiDif->Fill((smallerSet[i].phi-largerSet[j].phi)/TMath::Pi()*180);
				h_decDif->Fill((smallerSet[i].declination-largerSet[j].declination)/TMath::Pi()*180);
				h_raDif->Fill((smallerSet[i].ra-largerSet[j].ra)/TMath::Pi()*180);
				h_nHitsCoinc->Fill(smallerSet[i].nHits);
				h_xDif->Fill(smallerSet[i].x - largerSet[j].x);
				h_yDif->Fill(smallerSet[i].y - largerSet[j].y);
				h_zDif->Fill(smallerSet[i].z - largerSet[j].z);

				smallerSet[i].inCoincidence = true;
				largerSet[j].inCoincidence = true;
				break;
			}
		}
	}
	return nFoundCoincidences;
}

int StudyNotInCoincidence(vector<InputCascade> &cascades)
{
	cout << "No coincidence found for: " << endl;
	for (unsigned int i = 0; i < cascades.size(); ++i)
	{
		if (!cascades[i].inCoincidence)
		{
			;// PrintCascade(cascades[i]);
		}
	}

	return 0;
}

int findCascOverlap()
{
	TString lukasFile = "/Data/BaikalData/cascadeOverlap/recCasc_y19c-1.txt";
	vector<InputCascade> lukasCascades = ReadLukasCascades(lukasFile);

	// TString olgaFile = "/Data/BaikalData/cascadeOverlap/forLukasZuzana_sample2019_Nh12_Esh10_Esh100.txt";
	// TString olgaFile = "/Data/BaikalData/cascadeOverlap/forLukas_casvades2019_Eless100_Nhit12.txt";
	TString olgaFile = "/Data/BaikalData/cascadeOverlap/forLukas_HEcscd2019.txt";
	vector<InputCascade> olgasCascades = ReadOlgasCascades(olgaFile,2019);

	int nFoundCoincidences = CompareDatasets(lukasCascades,olgasCascades);
	cout << "Found coincidences: " << nFoundCoincidences << endl;

	StudyNotInCoincidence(olgasCascades);

	DrawResults();

	return 0;
}
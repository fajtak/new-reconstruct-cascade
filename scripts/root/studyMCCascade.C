#include <iostream>
#include <vector>
#include <fstream>

#include "TChain.h"
#include "TH2F.h"
#include "TFile.h"
#include "TMath.h"
#include "TVector3.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TProfile.h"
#include "TRandom2.h"
#include "TCanvas.h"
#include "TVirtualFitter.h"
#include "TRatioPlot.h"

#include "BEvent.h"
#include "BEventMask.h"
#include "BSource.h"
#include "BMCEvent.h"
#include "BGeomTelMC.h"



struct mcCascade
{
	int eventID, nHits, nNoiseHits = 0;
	float showerEnergy, cosTheta, phi, position[3], charge[288], time[288], noiseCharge[1000], noiseTime[1000];
	unsigned short chID[288], noiseChID[1000];
};

struct PulsesVariables
{
	int OMID;
	double time;
	double charge;
	bool noise;
};

vector<TVector3> gOMpositions(288);
const int gNOMs = 288;
const double ReciprocalSpeedOfLightinWater = 4.57;
vector<double> gNoiseTable;
vector<double> gOMQCal(288,0);
TRandom2 ranGen;

double gLogTable4D[200][351][21][21]{0};

vector<PulsesVariables> g_pulses;

TVirtualFitter* fMinuit;


TH2F* h_xyPosition = new TH2F("h_xyPosition","XY position;X [m];Y [m]",100,-500,500,100,-500,500);
TH1F* h_overallCharge = new TH1F("h_overallCharge","Overall charge; Q [p.e.]; NoE [#]",10000,0,100000);
TH1F* h_overallChargeWithNoise = new TH1F("h_overallChargeWithNoise","Overall charge (noise included); Q [p.e.]; NoE [#]",10000,0,100000);
TH2F* h_overallChargeVsEnergy = new TH2F("h_overallChargeVsEnergy","Overall charge vs Cascade energy; log(Q[p.e]);log_{10}(E_{sh}[GeV])",700,0,7,100,0,10);
TH2F* h_nHitsVsEnergy = new TH2F("h_nHitsVsEnergy","Number of hits after TFilter vs Cascade energy; Number of hits [#];log_{10}(E_{sh}[GeV])",200,0,200,100,0,10);
TH1F* h_nHitStrings = new TH1F("h_nHitStrings","Number of hit strings",10,0,10);
TH1F* h_nHits = new TH1F("h_nHits","Number of hits",300,0,300);
TH1F* h_nHitsWithNoise = new TH1F("h_nHitsWithNoise","Number of hits (noise included)",300,0,300);
TH1F* h_nHitsQover3p5 = new TH1F("h_nHitsQover3p5","Number of Hits with Q > 3.5 p.e.",100,0,100);
TH1F* h_initDist = new TH1F("h_initDist","Distance of the initial position estimation",100,0,100);
TH1F* h_fitDist = new TH1F("h_fitDist","Distance of the first fit",100,0,100);
TH1F* h_fitDistWell = new TH1F("h_fitDistWell","Distance of the first fit of the well fitted events",100,0,100);
TH1F* h_fitTDist = new TH1F("h_fitTDist","Distance of the second fit",100,0,100);
TH2F* h_initDistvsEnergy = new TH2F("h_initDistvsEnergy","Distance of the initial position estimation",100,0,10,100,0,100);
TH1F* h_xDiffCoor = new TH1F("h_xDiffCoor","X coordinate difference",100,-50,50);
TH1F* h_yDiffCoor = new TH1F("h_yDiffCoor","Y coordinate difference",100,-50,50);
TH1F* h_zDiffCoor = new TH1F("h_zDiffCoor","Z coordinate difference",100,-50,50);
TH1F* h_xFitDiffCoor = new TH1F("h_xFitDiffCoor","X coordinate difference",100,-50,50);
TH1F* h_yFitDiffCoor = new TH1F("h_yFitDiffCoor","Y coordinate difference",100,-50,50);
TH1F* h_zFitDiffCoor = new TH1F("h_zFitDiffCoor","Z coordinate difference",100,-50,50);
TH1F* h_xFitError = new TH1F("h_xFitError","#sigma_{x}",100,0,100);
TH1F* h_yFitError = new TH1F("h_yFitError","#sigma_{y}",100,0,100);
TH1F* h_zFitError = new TH1F("h_zFitError","#sigma_{z}",100,0,100);
TH1F* h_tFitError = new TH1F("h_tFitError","#sigma_{t}",100,0,1000);
TH2F* h_xFitCoorVsError = new TH2F("h_xFitCoorVsError","",100,-50,50,100,0,100);
TH1F* h_timeDiff = new TH1F("h_timeDiff","Time difference",150,-50,100);
TH2F* h_caussNHits = new TH2F("h_caussNHits","N hits vs. nCaussHits",300,0,300,300,0,300);
TH1F* h_caussEffi = new TH1F("h_caussEffi","CaussalityFilter efficiency",200,0,2);
TH1F* h_tPur = new TH1F("h_tPur","TFilter Purity",200,0,2);
TH1F* h_tEffi = new TH1F("h_tEffi","TFilter Efficiency",200,0,2);
TH1F* h_caussPur = new TH1F("h_caussPur","CaussalityPurity",200,0,2);
TH2F* h_caussDtDr = new TH2F("h_caussDdDr","CaussalityDistanceTime",500,0,500,2500,0,2500);
TH1F* h_caussDtMinusDr = new TH1F("h_caussDtMinusDr","CaussalityDistanceTime",200,-100,100);
TH1F* h_chi2 = new TH1F("h_chi2","#chi^{2}/NDF distribution after the first fit;#chi^{2}/NDF; NoE [#]",100000,0,10000);
TH2F* h_chi2vsDist = new TH2F("h_chi2vsDist","Chi^{2}/NDF vs mismatch distance; Chi^{2}/NDF; Mismatch distance [m]",100,0,100,100,0,100);
TH1F* h_chi2_T = new TH1F("h_chi2_T","#chi^{2}/NDF distribution after the second fit;#chi^{2}/NDF; NoE [#]",1000,0,100);
TH1F* h_nNoiseHits = new TH1F("h_nNoiseHits","Noise hits multiplicity; N_{noise} [#]; NoE [#]",20,0,20);
TH1F* h_vertDistHit = new TH1F("h_vertDistHit","Vertical distance of the hit from cascade",600,-300,300);
TH1F* h_vertDistNoiseHit = new TH1F("h_vertDistNoiseHit","Vertical distance of the hit from cascade",600,-300,300);

TH1F* h_nEvents = new TH1F("h_nEvents","Number of Events",100,0,10);
TH1F* h_nPassEvents = new TH1F("h_nPassEvents","Event-wise Efficiency;log_{10}(E_{sh}[GeV]);Efficiency [%]",100,0,10);
TH1F* h_nPureEvents = new TH1F("h_nPureEvents","Event-wise Purity;log_{10}(E_{sh}[GeV]);Purity [%]",100,0,10);
TH2F* h_efficiencyVsE = new TH2F("h_efficiencyVsE","Hit-wise Efficiency;log_{10}(E_{sh}[GeV]);Efficiency [%]",100,0,10,120,0,1.2);
TH2F* h_purityVsE = new TH2F("h_purityVsE","Hit-wise Purity;log_{10}(E_{sh}[GeV]);Purity [%]",100,0,10,120,0,1.2);
TProfile* efficiencyProfile;
TProfile* purityProfile;

TH1F* h_noiseChargeDistribution = new TH1F("h_noiseChargeDistribution","Noise charge distribution;Q [p.e.]; NoE [#]",500,0,50);
TH1F* h_noiseTimeDistribution = new TH1F("h_noiseTimeDistribution","Noise time distribution; T [ns]; NoE [#]",6000,-3000,3000);

TH1F* h_nHitsAfterCauss = new TH1F("h_nHitsAfterCauss","Number of hits after Causality filter",200,0,200);
TH1F* h_nHitsAfterTFilter = new TH1F("h_nHitsAfterTFilter","Number of hits after TFilter",200,0,200);
TH1F* h_nHitsChange = new TH1F("h_nHitsChange","Change in the number of hits between Causality filter and TFilter; N_{hits} change [#]; NoE [#]",100,-50,50);

TH1F* h_estimateTheta = new TH1F("h_estimateTheta","Theta estimate;#theta_{true}-#theta_{est};NoE [#]",360,-180,180);
TH1F* h_estimatePhi = new TH1F("h_estimatePhi","Phi estimate;#phi_{true}-#phi_{est};NoE [#]",720,-360,360);
TH1F* h_estimateEnergy = new TH1F("h_estimateEnergy","Energy estimate; E_{est}/E_{true} [#];NoE [#]",500,0,5);
TH1F* h_estimateEnergyWrong = new TH1F("h_estimateEnergyWrong","Energy estimate; E_{est}/E_{true} [#];NoE [#]",500,0,5);

TH1F* h_energyDist = new TH1F("h_energyDist","Energy Distribution; log_{10}(E_{sh}[GeV]);NoE [#]",100,0,10);
TH1F* h_energyDistWrong = new TH1F("h_energyDistWrong","Energy Distribution; log_{10}(E_{sh}[GeV]);NoE [#]",100,0,10);
TH1F* h_thetaDist = new TH1F("h_thetaDist","Theta Distribution; #theta [degree]",72,0,360);
TH1F* h_thetaDistWrong = new TH1F("h_thetaDistWrong","Theta Distribution; #theta [degree]",72,0,360);
TH1F* h_phiDist = new TH1F("h_phiDist","Phi Distribution",72,0,360);
TH1F* h_phiDistWrong = new TH1F("h_phiDistWrong","Phi Distribution",72,0,360);

TH2F* h_mcVsEstTheta = new TH2F("h_mcVsEstTheta","Estimated vs. True #theta;#theta_{true} [degree];#theta_{est} [degree]",36,0,180,36,0,180);
TH2F* h_mcVsEstPhi = new TH2F("h_mcVsEstPhi","Estimated vs. True #phi;#phi_{true} [degree]; #phi_{est} [degree]",72,0,360,72,0,360);
TH2F* h_mcVsEstEnergy = new TH2F("h_mcVsEstEnergy","Estimated vs. True Energy; log_{10}(E^{true}_{sh}[GeV];log_{10}(E^{est}_{sh}[GeV]",100,0,10,100,0,10);

TH1F* h_likelihood = new TH1F("h_likelihood","Likelihood minimum value;Likelihood [#];NoE [#]",1000,0,1000);
TH1F* h_mismatchTheta = new TH1F("h_mismatchTheta","Theta mismatch;#theta_{true}-#theta_{rec};NoE [#]",360,-180,180);
TH1F* h_mismatchPhi = new TH1F("h_mismatchPhi","Phi mismatch;#phi_{true}-#phi_{rec};NoE [#]",720,-360,360);
TH1F* h_mismatchEnergy = new TH1F("h_mismatchEnergy","Mismatch energy; E_{rec}/E_{true} [#];NoE [#]",500,-1,5);
TH1F* h_mismatchLogEnergy = new TH1F("h_mismatchLogEnergy","Mismatch log energy; log10(E_{rec}/E_{true}) [#];NoE [#]",500,-1,5);
TH2F* h_mismatchEnergyvsEnergy = new TH2F("h_mismatchEnergyvsEnergy","Mismatch energy vs. energy;log_{10}(E_{sh}[GeV]);E_{rec}/E_{true} [#]",100,0,10,500,0,5);
TH1F* h_mismatchEnergyvsEnergyMedian = new TH1F("h_mismatchEnergyvsEnergyMedian","Median mismatch energy vs. energy;log_{10}(E_{sh}[GeV]);median E_{rec}/E_{true} [#]",100,0,10);
TH1F* h_mismatchPosLong = new TH1F("h_mismatchPosLong","Mismatch Position Longitudinal Direction;Distance [m];NoE [#]",100,0,100);
TH1F* h_mismatchPosPerp = new TH1F("h_mismatchPosPerp","Mismatch Position Perpendicular Direction;Distance [m];NoE [#]",100,0,100);
TH1F* h_mismatchAngle = new TH1F("h_mismatchAngle","Mismatch angle;Mismatch angle [deg];NoE [#]",360,0,360);
TH1F* h_mismatchAngle0Noise = new TH1F("h_mismatchAngle0Noise","Mismatch angle;Mismatch angle [deg];NoE [#]",360,0,360);
TH1F* h_mismatchAngle1Noise = new TH1F("h_mismatchAngle1Noise","Mismatch angle;Mismatch angle [deg];NoE [#]",360,0,360);
TH1F* h_mismatchAngle2Noise = new TH1F("h_mismatchAngle2Noise","Mismatch angle;Mismatch angle [deg];NoE [#]",360,0,360);
TH1F* h_mismatchAngleWell = new TH1F("h_mismatchAngleWell","Mismatch angle of the well reconstructed events;Mismatch angle [deg];NoE [#]",360,0,360);
TH2F* h_mismatchAnglevsEnergy = new TH2F("h_mismatchAnglevsEnergy","Mismatch angle vs. energy;log_{10}(E_{sh}[GeV]);Mismatch angle [deg]",100,0,10,360,0,360);
TH1F* h_mismatchAnglevsEnergyMedian = new TH1F("h_mismatchAnglevsEnergyMedian","Median mismatch angle vs. energy;log_{10}(E_{sh}[GeV]);median Mismatch angle [deg]",100,0,10);
TH2F* h_mismatchPosLongvsEnergy = new TH2F("h_mismatchPosLongvsEnergy","Mismatch Position Longitudinal vs. energy;log_{10}(E_{sh}[GeV]);Distance [m]",100,0,10,100,0,100);
TH2F* h_mismatchPosPerpvsEnergy = new TH2F("h_mismatchPosPerpvsEnergy","Mismatch Position Perpendicular vs. energy;log_{10}(E_{sh}[GeV]);Distance [m]",100,0,10,100,0,100);
TH2F* h_mcVsRecoTheta = new TH2F("h_mcVsRecoTheta","Reconstructed vs. True #theta;#theta_{true} [degree];#theta_{reco} [degree]",180,0,180,180,0,180);
TH2F* h_mcVsRecoPhi = new TH2F("h_mcVsRecoPhi","Reconstructed vs. True #phi;#phi_{true} [degree]; #phi_{reco} [degree]",360,0,360,360,0,360);
TH2F* h_mcVsRecoEnergy = new TH2F("h_mcVsRecoEnergy","Reconstructed vs. True Energy; log_{10}(E^{true}_{sh}[GeV];log_{10}(E^{reco}_{sh}[GeV]",100,0,10,100,0,10);

TH1F* h_thetaError = new TH1F("h_thetaError","Error of the zenith angle obtained from MIGRAD;#sigma_{#theta} [degree]; NoE [#]",100,0,10);
TH1F* h_phiError = new TH1F("h_phiError","Error of the azimuth angle obtained from MIGRAD;#sigma_{#phi} [degree]; NoE [#]",100,0,10);

TH1F* h_energyEstNHits = new TH1F("h_energyEstNHits","Energy estimate based on the number of hits; E_{est}/E_{true} [#];NoE [#]",500,0,5);
TH1F* h_energyEstTotalQ = new TH1F("h_energyEstTotalQ","Energy estimate based on the overall charge; E_{est}/E_{true} [#];NoE [#]",500,0,5);
TH2F* h_energyEstNHitsVsEnergy = new TH2F("h_energyEstNHitsVsEnergy","Energy estimated based on the number of hits vs. energy; log_{10}(E^{true}_{sh}[GeV]E_{est}/E_{true} [#]",100,0,10,500,0,5);
TH2F* h_energyEstTotalQVsEnergy = new TH2F("h_energyEstTotalQVsEnergy","Energy estimated based on the overall charge vs. energy; log_{10}(E^{true}_{sh}[GeV]E_{est}/E_{true} [#]",100,0,10,500,0,5);

TH1F* h_notRecDist = new TH1F("h_notRecDist","Not reconstructed Events (Distance); Distance [m]; NoE [#]",1000,0,1000);
TH1F* h_notRecDistAll = new TH1F("h_notRecDistAll","Not reconstructed Events (Distance); Distance [m]; NoE [#]",1000,0,1000);
TH1F* h_notRecEner = new TH1F("h_notRecEner","Not reconstructed Events (Energy); log_{10}(E_{sh}[GeV])",100,0,10);
TH1F* h_notRecEnerAll = new TH1F("h_notRecEnerAll","Not reconstructed Events (Energy); log_{10}(E_{sh}[GeV])",100,0,10);
TH1F* h_notRecTheta = new TH1F("h_notRecTheta","Not reconstructed Events (Theta); Theta [degree]; NoE [#]",180,0,180);
TH1F* h_notRecThetaAll = new TH1F("h_notRecThetaAll","Not reconstructed Events (Theta); Theta [degree]; NoE [#]",180,0,180);

TH1F* h_vertQRatioUpGoing = new TH1F("vertQRatioUpGoing","Vertical Q-ratio for up-going events; r_{Q} [#]; NoE [#]",1000,0,100);
TH1F* h_vertQRatioDownGoing = new TH1F("vertQRatioDownGoing","Vertical Q-ratio for down-going events; r_{Q} [#]; NoE [#]",1000,0,100);

TH1F* h_zFilter = new TH1F("h_zFilter","Vertical distribution of the reconstructed cascades;Z [m]; NoE [#]",1000,-500,500);
TH1F* h_zFilterAll = new TH1F("h_zFilterAll","Vertical distribution of the reconstructed cascades;Z [m]; NoE [#]",1000,-500,500);
TH1F* h_qFilter = new TH1F("h_qFilter","Distribution of the QRatios; R [%]; NoE [#]",110,0,110);
TH1F* h_branchFilter = new TH1F("h_branchFilter","Distribution of the Branch differences; Up-Down [# hits]; NoE [#]",150,-50,100);
TH1F* h_closeHitsFilter = new TH1F("h_closeHitsFilter","Distribution of number of hit close OMs",20,0,20);

TH1F* h_zFilterCum;
TH1F* h_zFilterCumInv;
TH1F* h_qFilterCum;
TH1F* h_qFilterCumInv;
TH1F* h_branchFilterCum;
TH1F* h_branchFilterCumInv;
TH1F* h_closeHitsFilterCum;
TH1F* h_closeHitsFilterCumInv;

double ExpectedTime(double matrixX, double matrixY, double matrixZ, double matrixTime,int OMID)
{
  	double distanceFromLEDmatrix = TMath::Sqrt(TMath::Power(gOMpositions[OMID].X()-matrixX,2)+TMath::Power(gOMpositions[OMID].Y()-matrixY,2)+TMath::Power(gOMpositions[OMID].Z()-matrixZ,2));
  	// double scattering_correction = (gOMpositions[OMID].Z < LEDmatrixZ) ? (scatteringCorrection/15.0)*(LEDmatrixZ - myOMs[n].Z) : 0;

  	double expectedTime = matrixTime + distanceFromLEDmatrix*ReciprocalSpeedOfLightinWater;// + scattering_correction; //v nanosekundach

  	return  expectedTime;
}

// minimization function
void chi2(Int_t &npar, Double_t* gin, Double_t &f, Double_t* par, Int_t iflag) //keep all these nonsense parameters
{
	double theChi2 = 0;
	double constant = 1.0/(g_pulses.size() - 4);
  	double error = 1.0/4; //error is 3 ns - photomultiplier
  	double theChi;

	for (unsigned int i = 0; i < g_pulses.size(); ++i)
	{
		// chi calculation (measured - theoretical)
		// theChi = constant*error*(g_pulses[i].time - ExpectedTime(par[0], par[1], par[2], par[3], g_pulses[i].OMID))*TMath::Log10(g_pulses[i].charge+1);
		theChi = error*(g_pulses[i].time - ExpectedTime(par[0], par[1], par[2], par[3], g_pulses[i].OMID));
	  	// theChi = constant*error*(g_pulses[i].Times - ExpectedTime(myOMs,g_pulses[i].OMID, par[0], par[1], par[2], par[3]))*log(g_pulses[i].Charge);
		theChi2 += theChi*theChi;
	}
	f = constant*theChi2; // function returns calculated chi2, it is global function which is used in SetFCN()
}

void MEstimator(Int_t &npar, Double_t* gin, Double_t &f, Double_t* par, Int_t iflag)
{
	double sum = 0;
	for (unsigned int i = 0; i < g_pulses.size(); ++i)
	{
		double Tres = (g_pulses[i].time - ExpectedTime(par[0], par[1], par[2], par[3], g_pulses[i].OMID));
		double element = g_pulses[i].charge*TMath::Sqrt(1+TMath::Power(Tres,2)/2);
		sum += element;
	}
	f = sum; // function returns calculated chi2, it is global function which is used in SetFCN()

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
	for (unsigned int i = 0; i < g_pulses.size(); ++i)
	{
		if (g_pulses[i].OMID == OMID)
		{
			inGPulses = true;
			break;
		}
	}
	return !inGPulses;
}

void logLikelihood(Int_t &npar, Double_t* gin, Double_t &f, Double_t* par, Int_t iflag)
{
	// cout << "In log" << endl;
	double logLike = 0;
	double tableParameters[4]{0};
	// cout << "Calculating logLike" << endl;
	// cout << par[0] << " " << par[1] << " " << par[2] << " " << par[3] << " " << par[4] << " " << par[5] << " " << par[6] << endl;
	for (unsigned int i = 0; i < g_pulses.size(); ++i)
	{
		// cout << "GPulses: " << g_pulses[i].OMID << endl;
		GetParameters(par,g_pulses[i].OMID,tableParameters);
		// gOMpositions[g_pulses[i].OMID].Print();
		// cout << i << " " << tableParameters[0] << " " << tableParameters[1] << " " << tableParameters[2] << " " << tableParameters[3] << endl;
		double expectedNPE = GetInterpolatedValue(tableParameters);
		// cout << expectedNPE << " " << expectedNPE*100000000*par[4] << " " << g_pulses[i].charge << " " << TMath::PoissonI(g_pulses[i].charge,expectedNPE*100000000*par[4]) << endl;

		if (TMath::Poisson(g_pulses[i].charge,expectedNPE*110000000*par[4]) > 10e-320)
		{
			logLike -= TMath::Log10(TMath::Poisson(g_pulses[i].charge,expectedNPE*110000000*par[4]));
			// cout << TMath::Log10(TMath::PoissonI(g_pulses[i].charge,expectedNPE*100000000*par[4])) << endl;
		}
		else
		{
			logLike -= -320;
			// cout << -350 << endl;
		}
	}
	// cout << "After g pulses" << endl;
	// for (int i = 0; i < gNOMs; ++i)
	// {
	// 	if (NotInGPulses(i))
	// 	{
	// 		GetParameters(par,i,tableParameters);
	// 		double expectedNPE = GetInterpolatedValue(tableParameters);
	// 		// cout << i << " " << TMath::Poisson(0,expectedNPE*100000000*par[4]) << endl;

	// 		if (TMath::Poisson(0,expectedNPE*110000000*par[4]) > 10e-320)
	// 		{
	// 			logLike -= TMath::Log10(TMath::Poisson(0,expectedNPE*110000000*par[4]));
	// 			// cout << TMath::Log10(TMath::PoissonI(g_pulses[i].charge,expectedNPE*100000000*par[4])) << endl;
	// 			// cout << i << " " << TMath::Log10(TMath::PoissonI(0,expectedNPE*100000000*par[4])) << endl;
	// 		}
	// 		else
	// 		{
	// 			logLike -= -320;
	// 			// cout << i << " " << -350 << endl;
	// 		}
	// 	}
	// }
	// cout << "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF: " << logLike << endl;
	f = logLike;
}

double FitMatrixPosition(TVector3 &R2, double &T2, mcCascade* cascade, bool firstFit = true)
{
	fMinuit->ReleaseParameter(0);
	fMinuit->ReleaseParameter(1);
	fMinuit->ReleaseParameter(2);
	fMinuit->ReleaseParameter(3);

	//cout<<"fitEVent "<<endl;
	fMinuit->SetParameter(0,"LED X",R2.X(),1,-750,750); //estimation of start point
	fMinuit->SetParameter(1,"LED Y",R2.Y(),1,-750,750);
	fMinuit->SetParameter(2,"LED Z",R2.Z(),1,-750,750);
	fMinuit->SetParameter(3,"Time",T2,1,-10000,50000);

	fMinuit->FixParameter(4);
	fMinuit->FixParameter(5);
	fMinuit->FixParameter(6);

	double arglist[10] = {500,0.01}; //500 means 500 iterations and then MIGRAD stops
	fMinuit->ExecuteCommand("MIGRAD",arglist,2); //this line performs minimalisation

	double chi2, edm, errdef;
	int nvpar, nparx;

	fMinuit->GetStats(chi2,edm,errdef,nvpar,nparx);

	R2.SetX(fMinuit->GetParameter(0));
	R2.SetY(fMinuit->GetParameter(1));
	R2.SetZ(fMinuit->GetParameter(2));
	T2 = fMinuit->GetParameter(3);

	if (firstFit)
	{
		h_xFitCoorVsError->Fill(cascade->position[0]-R2[0],fMinuit->GetParError(0));

		h_xFitDiffCoor->Fill(cascade->position[0]-R2[0]);
		h_yFitDiffCoor->Fill(cascade->position[1]-R2[1]);
		h_zFitDiffCoor->Fill(cascade->position[2]-R2[2]);

		h_xFitError->Fill(fMinuit->GetParError(0));
		h_yFitError->Fill(fMinuit->GetParError(1));
		h_zFitError->Fill(fMinuit->GetParError(2));
		h_tFitError->Fill(fMinuit->GetParError(3));
	}
	 // cout<<"Chi2  "<<chi2/(g_pulses.size()-4)<<endl;
	/*cout<<"Pozicia kaskady po Fite "<<R2.X()<<endl;
	cout<<"Pozicia kaskady po Fite "<<R2.Y()<<endl;
	cout<<"Pozicia kaskady po Fite "<<R2.Z()<<endl;
	cout<<"Cas                     "<<T2<<endl;
	cout<<" "<<endl;*/

	return chi2;
	// cout << chi2 << " " << g_pulses.size()-4 << " " << TMath::Prob(chi2,g_pulses.size()-4) << endl;
	// return TMath::Prob(chi2,g_pulses.size()-4);
	// return chi2/(g_pulses.size()-4);
}

double FitMatrixDirection(TVector3 &R2, double &T2, double &energy, double &theta, double &phi, bool &goodError)
{
	// cout << "Fit Matrix direction" << endl;
	fMinuit->ReleaseParameter(4);
	fMinuit->ReleaseParameter(5);
	fMinuit->ReleaseParameter(6);

	// cout << "StartFitting" << endl;
	fMinuit->SetFCN(logLikelihood);
	fMinuit->SetParameter(0,"LED X",R2.X(),1,-750,750); //estimation of start point
	fMinuit->SetParameter(1,"LED Y",R2.Y(),1,-750,750);
	fMinuit->SetParameter(2,"LED Z",R2.Z(),1,-750,750);
	fMinuit->SetParameter(3,"Time",T2,1,-10000,50000);
	fMinuit->SetParameter(4,"Energy",energy,1,0,10000);
	fMinuit->SetParameter(5,"Theta",theta,0.1,0,TMath::Pi());
	fMinuit->SetParameter(6,"Phi",phi,0.1,0,2*TMath::Pi());

	fMinuit->FixParameter(0);
	fMinuit->FixParameter(1);
	fMinuit->FixParameter(2);
	fMinuit->FixParameter(3);

	// cout << "FIT" << endl;
	double arglist[10] = {500,0.01}; //500 means 500 iterations and then MIGRAD stops
	fMinuit->ExecuteCommand("MIGRAD",arglist,2); //this line performs minimalisation

	double chi2, edm, errdef;
	int nvpar, nparx;

	fMinuit->GetStats(chi2,edm,errdef,nvpar,nparx);
	energy = fMinuit->GetParameter(4);
	theta = fMinuit->GetParameter(5);
	phi = fMinuit->GetParameter(6);
	h_thetaError->Fill(fMinuit->GetParError(5));
	h_phiError->Fill(fMinuit->GetParError(6));

	if (fMinuit->GetParError(5) < 0.3 && fMinuit->GetParError(6) < 0.3)
		goodError = true;
	else
		goodError = false;

	// cout << chi2 << " " << g_pulses.size() << " " <<  chi2/g_pulses.size() << " " << fMinuit->GetParameter(5) << " " << fMinuit->GetParameter(6) << endl;
	return chi2/g_pulses.size();
	// return chi2/gNOMs;
}

// Fitter setting
void SetFitter(void)
{
	fMinuit = TVirtualFitter::Fitter(0,7); // the second number is number of parameters
	double arg = -1;
	fMinuit->ExecuteCommand("SET PRINTOUT",&arg,1); // these two lines means that it wont be able to print results on the screen
	fMinuit->ExecuteCommand("SET NOW", &arg ,1);
	fMinuit->SetFCN(chi2);
	// fMinuit->SetFCN(MEstimator);
}

bool IsContained(mcCascade* cascade, int distance = 60)
{
	if (TMath::Sqrt(TMath::Power(cascade->position[0],2)+TMath::Power(cascade->position[1],2)) < distance && TMath::Abs(cascade->position[2]) < 263+(distance-60))
		return true;
	else
		return false;
}

bool IsUncontained(mcCascade* cascade, int near, int far)
{
	double horizontalDist = TMath::Sqrt(TMath::Power(cascade->position[0],2)+TMath::Power(cascade->position[1],2));
	double verticalDist = TMath::Abs(cascade->position[2]);
	if ((horizontalDist < far && horizontalDist > near && verticalDist < 263) || (horizontalDist < far && verticalDist < 263+(far-60) && verticalDist > 263+(near-60)))
		return true;
	else
		return false;
}

double CountCharge(mcCascade*cascade,double &noiseCharge)
{
	double charge = 0;
	for (int i = 0; i < cascade->nHits; ++i)
	{
		charge += cascade->charge[cascade->chID[i]-1];
	}

	noiseCharge = 0;
	for (int i = 0; i < cascade->nNoiseHits; ++i)
	{
		noiseCharge += cascade->noiseCharge[i];
	}
	return charge;
}

double CountCharge(BEvent* event)
{
	double charge = 0;
	for (int i = 0; i < event->NHits(); ++i)
	{
		charge += event->Q(i);
	}

	return charge;
}

int GetNNoiseHits()
{
	int nNoiseHits = 0;
	for (unsigned int i = 0; i < g_pulses.size(); ++i)
	{
		if (g_pulses[i].noise)
			nNoiseHits++;
	}
	return nNoiseHits;
}

double CountCascadeCharge()
{
	double charge = 0;
	for (unsigned int i = 0; i < g_pulses.size(); ++i)
	{
		charge += g_pulses[i].charge;
	}
	return charge;
}

void PrintCascade(mcCascade* cascade)
{
	for (int i = 0; i < cascade->nHits; ++i)
	{
		cout << cascade->chID[i]-1 << " ";
	}
	cout << endl;
}

int CountStrings(mcCascade* cascade)
{
	int nHitStrings = 0;
	int hitStrings[8]{0};

	for (int i = 0; i < cascade->nHits; ++i)
	{
		hitStrings[(cascade->chID[i]-1)/36] = 1;
	}

	for (int i = 0; i < 8; ++i)
	{
		nHitStrings += hitStrings[i];
	}
	// cout << nHitStrings << endl;
	// PrintCascade(cascade);

	return nHitStrings;
}

int CountStrings()
{
	int nHitStrings = 0;
	int hitStrings[8]{0};

	for (unsigned int i = 0; i < g_pulses.size(); ++i)
	{
		hitStrings[(g_pulses[i].OMID)/36] = 1;
	}

	for (int i = 0; i < 8; ++i)
	{
		nHitStrings += hitStrings[i];
	}
	// cout << nHitStrings << endl;
	// PrintCascade(cascade);

	return nHitStrings;
}

bool CheckPurity()
{
	bool pure = true;
	for (unsigned int i = 0; i < g_pulses.size(); ++i)
	{
		if (g_pulses[i].noise)
		{
			pure = false;
			break;
		}
	}
	return pure;
}

int CountHitsWithQHigherThan(mcCascade* cascade)
{
	int nHits = 0;
	for (int i = 0; i < cascade->nHits; ++i)
	{
		if (cascade->charge[cascade->chID[i]-1] > 6.5)
		{
			nHits++;
		}
	}
	return nHits;
}

double EstimateInitialMatrixPositionHQ(mcCascade* cascade, TVector3& position, double& matrixTime)
{
	// cout << "High charge based initial estimation" << endl;
	double maxCharge = numeric_limits<double>::min();
	int maxChargePulseID = -1;
	for (unsigned int i = 0; i < g_pulses.size(); ++i)
	{
		if (g_pulses[i].charge > maxCharge)
		{
			maxCharge = g_pulses[i].charge;
			maxChargePulseID = i;
		}
	}
	position = gOMpositions[g_pulses[maxChargePulseID].OMID];
	matrixTime = g_pulses[maxChargePulseID].time;

	h_xDiffCoor->Fill(cascade->position[0]-position[0]);
	h_yDiffCoor->Fill(cascade->position[1]-position[1]);
	h_zDiffCoor->Fill(cascade->position[2]-position[2]);
	h_timeDiff->Fill(matrixTime);

	return TMath::Sqrt(TMath::Power(cascade->position[0]-position.X(),2)+TMath::Power(cascade->position[1]-position.Y(),2)+TMath::Power(cascade->position[2]-position.Z(),2));
	// position.Print();
	// cout << matrixTime << endl;
}

// minimization function
double chi2(Double_t* par) //keep all these nonsense parameters
{
	double theChi2 = 0;
	double constant = 1.0/(g_pulses.size());
  	double error = 1.0/3; //error is 3 ns - photomultiplier
  	double theChi;

	for (unsigned int i = 0; i < g_pulses.size(); ++i)
	{
		// chi calculation (measured - theoretical)
		// theChi = constant*error*(g_pulses[i].time - ExpectedTime(par[0], par[1], par[2], par[3], g_pulses[i].OMID))*TMath::Log10(g_pulses[i].charge+1);
		// theChi = constant*error*(cascade->time[cascade->chID[i]-1] - ExpectedTime(par[0], par[1], par[2], par[3], cascade->chID[i]-1));
		theChi = error*(g_pulses[i].time - ExpectedTime(par[0], par[1], par[2], par[3], g_pulses[i].OMID));
	  	// theChi = constant*error*(g_pulses[i].Times - ExpectedTime(myOMs,g_pulses[i].OMID, par[0], par[1], par[2], par[3]))*log(g_pulses[i].Charge);
		theChi2 += theChi*theChi;
	}
	return constant*theChi2; // function returns calculated chi2, it is global function which is used in SetFCN()
}

double EstimateInitialMatrixPositionChi(mcCascade* cascade, TVector3& position, double& matrixTime)//, mcCascade* cascade)
{
	// cout << "Grid based initial estimation" << endl;

	double likelihoodValue = 0;
	double cascadeParameters[7];

	double lowestLog = numeric_limits<double>::max();

	double maxCharge = numeric_limits<double>::min();
	int maxChargePulseID = -1;
	for (unsigned int i = 0; i < g_pulses.size(); ++i)
	{
		if (g_pulses[i].charge > maxCharge)
		{
			maxCharge = g_pulses[i].charge;
			maxChargePulseID = i;
		}
	}

	for (int i = 0; i < 7; ++i)
	{
		for (int j = 0; j < 7; ++j)
		{
			for (int k = 0; k < 27; ++k)
			{
				// Int_t &npar, Double_t* gin, Double_t &f, Double_t* par, Int_t iflag
				cascadeParameters[0] = -60+20*i+gOMpositions[269].X();
				cascadeParameters[1] = -60+20*j+gOMpositions[269].Y();
				cascadeParameters[2] = -260+20*k+gOMpositions[269].Z();
				double distance = TMath::Sqrt(TMath::Power(cascadeParameters[0]-gOMpositions[cascade->chID[maxChargePulseID]].X(),2)+TMath::Power(cascadeParameters[1]-gOMpositions[cascade->chID[maxChargePulseID]].Y(),2)+TMath::Power(cascadeParameters[2]-gOMpositions[cascade->chID[maxChargePulseID]].Z(),2));
				cascadeParameters[3] = cascade->time[cascade->chID[maxChargePulseID]-1] - distance*ReciprocalSpeedOfLightinWater;
				likelihoodValue = chi2(cascadeParameters);
				// cout << likelihoodValue << " " << cascadeParameters[0] << " " << cascadeParameters[1] << " " << cascadeParameters[2] << endl;
				if (likelihoodValue < lowestLog)
				{
					lowestLog = likelihoodValue;
					position[0] = cascadeParameters[0];
					position[1] = cascadeParameters[1];
					position[2] = cascadeParameters[2];
					matrixTime = cascadeParameters[3];
				}
			}
		}
	}

	h_xDiffCoor->Fill(cascade->position[0]-position[0]);
	h_yDiffCoor->Fill(cascade->position[1]-position[1]);
	h_zDiffCoor->Fill(cascade->position[2]-position[2]);
	h_timeDiff->Fill(matrixTime);

	// cout << lowestLog << " " << position[0] << " " << position[1] << " " << position[2] << " " << matrixTime << endl;
	return TMath::Sqrt(TMath::Power(cascade->position[0]-position.X(),2)+TMath::Power(cascade->position[1]-position.Y(),2)+TMath::Power(cascade->position[2]-position.Z(),2));
}

double EstimateInitialMatrixPositionMatrix(mcCascade* cascade, TVector3& position, double& matrixTime)
{
	TMatrixD A(g_pulses.size()-1,4);
	TVectorD b(g_pulses.size()-1);

	for (unsigned int i = 0; i < g_pulses.size()-1; ++i)
	{
		double temp = TMath::Power(gOMpositions[g_pulses[i+1].OMID].X(),2) - TMath::Power(gOMpositions[g_pulses[i].OMID].X(),2);
		temp += TMath::Power(gOMpositions[g_pulses[i+1].OMID].Y(),2) - TMath::Power(gOMpositions[g_pulses[i].OMID].Y(),2);
		temp += TMath::Power(gOMpositions[g_pulses[i+1].OMID].Z(),2) - TMath::Power(gOMpositions[g_pulses[i].OMID].Z(),2);
		temp -= (TMath::Power(g_pulses[i+1].time,2) - TMath::Power(g_pulses[i].time,2))/TMath::Power(ReciprocalSpeedOfLightinWater,2);
		b[i] = temp;

		// cout << cascade->chID[i+1]-1 << " " << cascade->chID[i]-1 << endl;
		A[i][0] = 2*(gOMpositions[g_pulses[i+1].OMID].X() - gOMpositions[g_pulses[i].OMID].X());
		A[i][1] = 2*(gOMpositions[g_pulses[i+1].OMID].Y() - gOMpositions[g_pulses[i].OMID].Y());
		A[i][2] = 2*(gOMpositions[g_pulses[i+1].OMID].Z() - gOMpositions[g_pulses[i].OMID].Z());
		A[i][3] = -2*(g_pulses[i+1].time - g_pulses[i].time)/TMath::Power(ReciprocalSpeedOfLightinWater,2);
	}

	TMatrixD B(4,g_pulses.size()-1);
	B.Transpose(A);
	TMatrixD C = (B*A).Invert();
	TVectorD D = A.T()*b;
	TVectorD X = C*D;

	// cout << X[0] << " " << X[1] << " " << X[2] << " " << X[3] << endl;

	h_xDiffCoor->Fill(cascade->position[0]-X[0]);
	h_yDiffCoor->Fill(cascade->position[1]-X[1]);
	h_zDiffCoor->Fill(cascade->position[2]-X[2]);
	h_timeDiff->Fill(X[3]);

	position[0] = X[0];
	position[1] = X[1];
	position[2] = X[2];
	matrixTime = X[3];

	return TMath::Sqrt(TMath::Power(cascade->position[0]-X[0],2)+TMath::Power(cascade->position[1]-X[1],2)+TMath::Power(cascade->position[2]-X[2],2));
}

double EstimateInitialMatrixPositionMatrix(TVector3& position, double& matrixTime)
{
	TMatrixD A(g_pulses.size()-1,4);
	TVectorD b(g_pulses.size()-1);

	for (unsigned int i = 0; i < g_pulses.size()-1; ++i)
	{
		double temp = TMath::Power(gOMpositions[g_pulses[i+1].OMID].X(),2) - TMath::Power(gOMpositions[g_pulses[i].OMID].X(),2);
		temp += TMath::Power(gOMpositions[g_pulses[i+1].OMID].Y(),2) - TMath::Power(gOMpositions[g_pulses[i].OMID].Y(),2);
		temp += TMath::Power(gOMpositions[g_pulses[i+1].OMID].Z(),2) - TMath::Power(gOMpositions[g_pulses[i].OMID].Z(),2);
		temp -= (TMath::Power(g_pulses[i+1].time,2) - TMath::Power(g_pulses[i].time,2))/TMath::Power(ReciprocalSpeedOfLightinWater,2);
		b[i] = temp;

		// cout << cascade->chID[i+1]-1 << " " << cascade->chID[i]-1 << endl;
		A[i][0] = 2*(gOMpositions[g_pulses[i+1].OMID].X() - gOMpositions[g_pulses[i].OMID].X());
		A[i][1] = 2*(gOMpositions[g_pulses[i+1].OMID].Y() - gOMpositions[g_pulses[i].OMID].Y());
		A[i][2] = 2*(gOMpositions[g_pulses[i+1].OMID].Z() - gOMpositions[g_pulses[i].OMID].Z());
		A[i][3] = -2*(g_pulses[i+1].time - g_pulses[i].time)/TMath::Power(ReciprocalSpeedOfLightinWater,2);
	}

	TMatrixD B(4,g_pulses.size()-1);
	B.Transpose(A);
	TMatrixD C = (B*A).Invert();
	TVectorD D = A.T()*b;
	TVectorD X = C*D;

	// cout << X[0] << " " << X[1] << " " << X[2] << " " << X[3] << endl;


	h_timeDiff->Fill(X[3]);

	position[0] = X[0];
	position[1] = X[1];
	position[2] = X[2];
	matrixTime = X[3];

	return 0;
}

int GenerateNoise(mcCascade* cascade, int noiseRateInkHz)
{

	cascade->nNoiseHits = 0;
	for (int i = 0; i < gNOMs; ++i)
	{
		int noisePulses = ranGen.Poisson(noiseRateInkHz/200.0);
		// cout << "OMID: " << i << " nNoisePulses: " << noisePulses << endl;
		for (int j = 0; j < noisePulses; ++j)
		{
			double noiseCharge = 0;//ranGen.Landau(1,2);
			double ranNumber = ranGen.Uniform(1);
			for (unsigned int j = 0; j < gNoiseTable.size(); ++j)
			{
				if (gNoiseTable[j] > ranNumber)
				{
					noiseCharge = j*0.1-0.05;
					break;
				}
			}
			double noiseTime = ranGen.Uniform(5120)-2560;
			cascade->noiseCharge[cascade->nNoiseHits] = noiseCharge;
			h_noiseChargeDistribution->Fill(noiseCharge);
			cascade->noiseTime[cascade->nNoiseHits] = noiseTime;
			h_noiseTimeDistribution->Fill(noiseTime);
			cascade->noiseChID[cascade->nNoiseHits] = i;
			cascade->nNoiseHits++;
			// cout << "\t" << j << " " << noiseTime << " " << noiseCharge << endl;
		}
	}
	return 0;
}

bool bigger_charge(const PulsesVariables& x, const PulsesVariables& y) { return x.charge > y.charge; }

int CaussalityFilter(mcCascade* cascade)
{
	double chargeFilt = 1.5;
	vector<PulsesVariables> vecCharges;

	for (int i = 0; i < cascade->nHits; ++i)
	{
		if (cascade->charge[cascade->chID[i]-1] > chargeFilt)
			vecCharges.push_back(PulsesVariables{cascade->chID[i]-1,cascade->time[cascade->chID[i]-1],cascade->charge[cascade->chID[i]-1],false});
	}
	for (int i = 0; i < cascade->nNoiseHits; ++i)
	{
		if (cascade->noiseCharge[i] > chargeFilt)
			vecCharges.push_back(PulsesVariables{cascade->noiseChID[i],cascade->noiseTime[i],cascade->noiseCharge[i],true});
	}

	sort(vecCharges.begin(), vecCharges.end(), bigger_charge);

	g_pulses.clear();

	for (unsigned int i = 0; i < vecCharges.size(); ++i)
	{
		bool addPulse = true;
		for (unsigned int j = 0; j < g_pulses.size(); ++j)
		{
			double distance = TMath::Sqrt(TMath::Power(gOMpositions[vecCharges[i].OMID].X()-gOMpositions[g_pulses[j].OMID].X(),2)+TMath::Power(gOMpositions[vecCharges[i].OMID].Y()-gOMpositions[g_pulses[j].OMID].Y(),2)+TMath::Power(gOMpositions[vecCharges[i].OMID].Z()-gOMpositions[g_pulses[j].OMID].Z(),2));

			if (j == 0)
			{
				h_caussDtMinusDr->Fill(distance-abs(vecCharges[i].time - g_pulses[j].time)*0.2188);
				h_caussDtDr->Fill(distance,abs(vecCharges[i].time - g_pulses[j].time));
			}
			if(abs(g_pulses[j].time - vecCharges[i].time) >= (distance*ReciprocalSpeedOfLightinWater + 50))
			{
				addPulse = false;
				break;
			}
		}
		bool friendFound = false;
		for (unsigned int j = 0; j < vecCharges.size(); ++j)
		{
			if (i != j && abs(vecCharges[i].OMID-vecCharges[j].OMID) <= 2 && abs(vecCharges[i].time-vecCharges[j].time) < abs(vecCharges[i].OMID-vecCharges[j].OMID)*70+50)
			{
				friendFound = true;
				break;
			}
		}
		if (addPulse && friendFound)
		// if (addPulse)
			g_pulses.push_back(PulsesVariables{vecCharges[i].OMID,vecCharges[i].time,vecCharges[i].charge,vecCharges[i].noise});
	}

	// cout << vecCharges.size() << " " << g_pulses.size() << endl;
	return 0;
}

int CaussalityFilter(BEvent* event)
{
	double chargeFilt = 1.5;
	vector<PulsesVariables> vecCharges;

	// cout << event->GetTotImpulses() << endl;

	for (int i = 0; i < event->GetTotImpulses(); ++i)
	{
		if (event->Q(i) > chargeFilt)
			vecCharges.push_back(PulsesVariables{event->HitChannel(i),event->T(i),event->Q(i),false});
	}

	sort(vecCharges.begin(), vecCharges.end(), bigger_charge);

	g_pulses.clear();

	for (unsigned int i = 0; i < vecCharges.size(); ++i)
	{
		bool addPulse = true;
		for (unsigned int j = 0; j < g_pulses.size(); ++j)
		{
			double distance = TMath::Sqrt(TMath::Power(gOMpositions[vecCharges[i].OMID].X()-gOMpositions[g_pulses[j].OMID].X(),2)+TMath::Power(gOMpositions[vecCharges[i].OMID].Y()-gOMpositions[g_pulses[j].OMID].Y(),2)+TMath::Power(gOMpositions[vecCharges[i].OMID].Z()-gOMpositions[g_pulses[j].OMID].Z(),2));

			if (j == 0)
			{
				h_caussDtMinusDr->Fill(distance-abs(vecCharges[i].time - g_pulses[j].time)*0.2188);
				h_caussDtDr->Fill(distance,abs(vecCharges[i].time - g_pulses[j].time));
			}
			if(abs(g_pulses[j].time - vecCharges[i].time)/ReciprocalSpeedOfLightinWater > (distance + 50))
			{
				addPulse = false;
				break;
			}
		}
		bool friendFound = false;
		for (unsigned int j = 0; j < vecCharges.size(); ++j)
		{
			if (i != j && abs(vecCharges[i].OMID-vecCharges[j].OMID) <= 2 && abs(vecCharges[i].time-vecCharges[j].time) < abs(vecCharges[i].OMID-vecCharges[j].OMID)*70+50)
			{
				friendFound = true;
				break;
			}
		}
		if (addPulse && friendFound)
		// if (addPulse)
			g_pulses.push_back(PulsesVariables{vecCharges[i].OMID,vecCharges[i].time,vecCharges[i].charge,vecCharges[i].noise});
	}

	// cout << vecCharges.size() << " " << g_pulses.size() << endl;
	return 0;
}

void FillPurityEfficiency(mcCascade* cascade, bool afterCaussFilter)
{
	int nNoisePassed = 0;
	int nSignalPassed = 0;
	// cout << "************************" << endl;
	for (unsigned int i = 0; i < g_pulses.size(); ++i)
	{
		// cout << g_pulses[i].OMID << " " << g_pulses[i].charge << " " << g_pulses[i].time << " " << g_pulses[i].noise << endl;
		if (g_pulses[i].noise)
			nNoisePassed++;
		else
			nSignalPassed++;
	}

	if (afterCaussFilter)
	{
		h_efficiencyVsE->Fill(TMath::Log10(cascade->showerEnergy*1000),nSignalPassed/(double)cascade->nHits);
		h_purityVsE->Fill(TMath::Log10(cascade->showerEnergy*1000),1-(nNoisePassed/(double)g_pulses.size()));

		h_caussEffi->Fill(nSignalPassed/(double)cascade->nHits);
		h_caussPur->Fill(1-(nNoisePassed/(double)g_pulses.size()));

		h_caussNHits->Fill(cascade->nHits,nSignalPassed);
	}else{
		h_tEffi->Fill(nSignalPassed/(double)cascade->nHits);
		h_tPur->Fill(1-(nNoisePassed/(double)g_pulses.size()));
	}
}

int TFilter(mcCascade* cascade, TVector3& matrixPosition, double& matrixTime)
{
	int nPulses = 0;
	g_pulses.clear();

	for(int i = 0; i < cascade->nHits; i++)
	{
		if(cascade->charge[cascade->chID[i]-1] > 0)
		{
			double distanceFromLEDmatrix = (matrixPosition - gOMpositions[cascade->chID[i]-1]).Mag();
			// double scattering_correction = (gOMpositions[event->HitChannel(i)].Z() < matrixPosition.Z()) ? (gScatteringCorrection/15.0)*(matrixPosition.Z() - gOMpositions[event->HitChannel(i)].Z()) : 0;
	        double expectedTime = matrixTime + distanceFromLEDmatrix*ReciprocalSpeedOfLightinWater;// + scattering_correction;

	        // h_timeRes->Fill(cascade->time[cascade->chID[i]-1] - expectedTime);
	        if(cascade->time[cascade->chID[i]-1] >= expectedTime-50 && cascade->time[cascade->chID[i]-1] <= expectedTime + 50)
	        {
	        	// if (OMIDAlreadyInGPulses(event->GetImpulse(i)))
	        		// continue;
	        	if (gOMpositions[cascade->chID[i]-1].Z()-matrixPosition.Z() < -140 || gOMpositions[cascade->chID[i]-1].Z()-matrixPosition.Z() > 180)
	        		continue;
	        	g_pulses.push_back(PulsesVariables{cascade->chID[i]-1,cascade->time[cascade->chID[i]-1],cascade->charge[cascade->chID[i]-1],false});
	        	nPulses++;
	        }
	    }
	}
	for(int i = 0; i < cascade->nNoiseHits; i++)
	{
		if(cascade->noiseCharge[i] > 0)
		{
			double distanceFromLEDmatrix = (matrixPosition - gOMpositions[cascade->noiseChID[i]]).Mag();
			// double scattering_correction = (gOMpositions[event->HitChannel(i)].Z() < matrixPosition.Z()) ? (gScatteringCorrection/15.0)*(matrixPosition.Z() - gOMpositions[event->HitChannel(i)].Z()) : 0;
	        double expectedTime = matrixTime + distanceFromLEDmatrix*ReciprocalSpeedOfLightinWater;// + scattering_correction;

	        // h_timeRes->Fill(cascade->noiseTime[i] - expectedTime);
	        if(cascade->noiseTime[i] >= expectedTime-50 && cascade->noiseTime[i] <= expectedTime + 50)
	        {
	        	// if (OMIDAlreadyInGPulses(event->GetImpulse(i)))
	        		// continue;
	        	if (gOMpositions[cascade->noiseChID[i]].Z()-matrixPosition.Z() < -140 || gOMpositions[cascade->noiseChID[i]].Z()-matrixPosition.Z() > 180)
	        		continue;
	        	g_pulses.push_back(PulsesVariables{cascade->noiseChID[i],cascade->noiseTime[i],cascade->noiseCharge[i],true});
	        	nPulses++;
	        }
	    }
	}

	// h_nHitsAfterT->Fill(nPulses);
	return nPulses;
}

double LikelihoodFilterPassed(TVector3 &matrixPosition, double &matrixTime, double &energy, double &theta, double &phi, double &likelihood, bool &passedGoodError)
{
	// cout << "In likelihood" << endl;
	double lowestLog = 10000;
	for (int k = 0; k < 1; ++k)
	{
		for (int l = 0; l < 1; ++l)
		{
			for (int i = 0; i < 1; ++i)
			{
				double cascadeEnergy = TMath::Power(10,1);
				double cascadeTheta = TMath::Pi()/7*k;
				double cascadePhi = TMath::Pi()*2/11*l;
				// double recentLog = FitMatrixDirection(matrixPosition,matrixTime,cascadeEnergy,cascadeTheta,cascadePhi);
				double recentLog = FitMatrixDirection(matrixPosition,matrixTime,energy,theta,phi,passedGoodError);
				// double recentLog = 100000;
				lowestLog = recentLog;
				// if (recentLog < lowestLog)
				// {
				// 	lowestLog = recentLog;
				// 	energy = cascadeEnergy;
				// 	theta = cascadeTheta;
				// 	phi = cascadePhi;
				// }
				// cout << k << " " << l << " " << recentLog  << " " << cascadeEnergy << " " << cascadeTheta << " " << cascadePhi << endl;
			}
			// cout << k << " " << l << endl;
		}
	}
	likelihood = lowestLog;
	// cout << "End likelihood" << endl;
	// cout << lowestLog << " " << energy << " " << theta << " " << phi << endl;
	return lowestLog;
}

double LikelihoodFilterPassedGrid(TVector3 &matrixPosition, double &matrixTime, double &energy, double &theta, double &phi, double &likelihood, bool &passedGoodError)
{
	int nThetaSteps = 4;
	int nPhiSteps = 6;
	int nEnergySteps = 4;
	// cout << "In likelihood" << endl;
	double lowestLog = 10000;
	for (int k = 0; k < nThetaSteps-1; ++k)
	{
		for (int l = 0; l < nPhiSteps-1; ++l)
		{
			for (int i = 0; i < nEnergySteps; ++i)
			{
				double cascadeEnergy = TMath::Power(10,i);
				// double cascadeEnergy = 100*i;
				// double cascadeEnergy = 50;
				double cascadeTheta = TMath::Pi()/(nThetaSteps)*(k+1);
				double cascadePhi = TMath::Pi()*2/nPhiSteps*(l+1);
				double recentLog = FitMatrixDirection(matrixPosition,matrixTime,cascadeEnergy,cascadeTheta,cascadePhi,passedGoodError);
				// double recentLog = FitMatrixDirection(matrixPosition,matrixTime,energy,theta,phi,passedGoodError);
				// double recentLog = 100000;
				// lowestLog = recentLog;
				if (recentLog < lowestLog)
				{
					lowestLog = recentLog;
					energy = cascadeEnergy;
					theta = cascadeTheta;
					phi = cascadePhi;
				}
				// cout << k << " " << l << " " << recentLog  << " " << cascadeEnergy << " " << cascadeTheta << " " << cascadePhi << endl;
			}
			// cout << k << " " << l << endl;
		}
	}
	likelihood = lowestLog;
	// cout << "End likelihood" << endl;
	// cout << lowestLog << " " << energy << " " << theta << " " << phi << endl;
	return lowestLog;
}

double EstimateInitialDirection(mcCascade* cascade, TVector3& position, double& matrixTime, double& energy, double &theta, double &phi)
{
	int nPar = 0;
	double* gin = new double(0);
	int iflag = 0;
	double likelihoodValue = 0;
	double cascadeParameters[7];

	double lowestLog = 1000000;

	cascadeParameters[0] = position[0];
	cascadeParameters[1] = position[1];
	cascadeParameters[2] = position[2];
	cascadeParameters[3] = matrixTime;


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

	double cascadeEnergyTrue = cascade->showerEnergy;
	double cascadeThetaTrue = TMath::Pi()-TMath::ACos(cascade->cosTheta);
	double cascadePhiTrue = cascade->phi+TMath::Pi();
	if (cascadePhiTrue > 2*TMath::Pi())
		cascadePhiTrue -= 2*TMath::Pi();

	h_estimateTheta->Fill((cascadeThetaTrue-theta)/TMath::Pi()*180);
	double phiDiff = cascadePhiTrue-phi;
	if (phiDiff > TMath::Pi())
		phiDiff = 2*TMath::Pi()-phiDiff;
	if (phiDiff < -TMath::Pi())
		phiDiff = -2*TMath::Pi()-phiDiff;
	h_estimatePhi->Fill((phiDiff)/TMath::Pi()*180);
	h_estimateEnergy->Fill(energy/cascadeEnergyTrue);
	h_energyDist->Fill(TMath::Log10(cascadeEnergyTrue*1000));
	h_thetaDist->Fill(cascadeThetaTrue/TMath::Pi()*180);
	h_phiDist->Fill(cascadePhiTrue/TMath::Pi()*180);

	h_mcVsEstTheta->Fill(cascadeThetaTrue/TMath::Pi()*180,theta/TMath::Pi()*180);
	h_mcVsEstPhi->Fill(cascadePhiTrue/TMath::Pi()*180,phi/TMath::Pi()*180);
	h_mcVsEstEnergy->Fill(TMath::Log10(cascadeEnergyTrue*1000),TMath::Log10(energy*1000));


	// if (abs((cascadeThetaTrue-theta)/TMath::Pi()*180) > 20 || abs((phiDiff)/TMath::Pi()*180) > 20)
	if (abs((cascadeThetaTrue-theta)/TMath::Pi()*180) > 20)
	{
		h_energyDistWrong->Fill(TMath::Log10(cascadeEnergyTrue*1000));
		h_thetaDistWrong->Fill(cascadeThetaTrue/TMath::Pi()*180);
		h_phiDistWrong->Fill(cascadePhiTrue/TMath::Pi()*180);
		h_estimateEnergyWrong->Fill(energy/cascadeEnergyTrue);
		return 1;
	}

	return 0;
}

double EstimateInitialDirectionMC(TVector3& position, double& matrixTime, double& energy, double &theta, double &phi)
{
	int nPar = 0;
	double* gin = new double(0);
	int iflag = 0;
	double likelihoodValue = 0;
	double cascadeParameters[7];

	double lowestLog = 1000000;

	cascadeParameters[0] = position[0];
	cascadeParameters[1] = position[1];
	cascadeParameters[2] = position[2];
	cascadeParameters[3] = matrixTime;


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

void FillVertDist(mcCascade* cascade)
{
	for (unsigned int i = 0; i < g_pulses.size(); ++i)
	{
		if (!g_pulses[i].noise)
			h_vertDistHit->Fill(gOMpositions[g_pulses[i].OMID].Z()-cascade->position[2]);
		else
			h_vertDistNoiseHit->Fill(gOMpositions[g_pulses[i].OMID].Z()-cascade->position[2]);
	}
}

double GetVertQRatio(TVector3 matrixPosition)
{
	double upperCharge = 0;
	double lowerCharge = 0;

	double maxDistance = gOMpositions[35].Z()-matrixPosition.Z() < matrixPosition.Z()-gOMpositions[0].Z()? gOMpositions[35].Z()-matrixPosition.Z() : matrixPosition.Z()-gOMpositions[0].Z();

	for (unsigned int i = 0; i < g_pulses.size(); ++i)
	{
		if (TMath::Abs(gOMpositions[g_pulses[i].OMID].Z() - matrixPosition.Z()) > maxDistance)
			continue;
		if (gOMpositions[g_pulses[i].OMID].Z() - matrixPosition.Z() > 0)
		// if (gOMpositions[g_pulses[i].OMID].Z() - cascade->position[2] > 0)
			upperCharge += g_pulses[i].charge;
		else
			lowerCharge += g_pulses[i].charge;
	}
	return upperCharge/lowerCharge;
}

int CountBranchFilter(TVector3& position)
{
	int upper = 0;
	int lower = 0;
	for (unsigned int i = 0; i < g_pulses.size(); ++i)
	{
		if (gOMpositions[g_pulses[i].OMID].Z() > position.Z())
			upper++;
		else
			lower++;
	}
	// cout << upper << " " << lower << endl;
	return upper-lower;
}

struct OMIDdistance
{
	int OMID;
	double distance;
};

bool compareEntries(OMIDdistance first, OMIDdistance second)
{
	return (first.distance <= second.distance);
}

int CountCloseHitsFilter(TVector3& matrixPosition)
{
	int closeHits = 0;

	vector<OMIDdistance> vOMIDdistance;
	// matrixPosition.Print();

	for (int i = 0; i < gNOMs; ++i)
	{
		if (gOMQCal[i] != -1 && gOMpositions[i].Z() > matrixPosition.Z()-20)
		{
			// cout << i << endl;
			// gOMpositions[i].Print();
			double OMdistance = (matrixPosition - gOMpositions[i]).Mag();
			vOMIDdistance.push_back(OMIDdistance{i,OMdistance});
			// cout << OMdistance <<endl;
		}
	}
	sort(vOMIDdistance.begin(),vOMIDdistance.end(),compareEntries);

	int nStudiedHits = vOMIDdistance.size() < 15 ? vOMIDdistance.size() : 15;

	for (int i = 0; i < nStudiedHits; ++i)
	{
		// cout << i << " " << vOMIDdistance[i].OMID << " " << vOMIDdistance[i].distance << endl;
		for (unsigned int j = 0; j < g_pulses.size(); ++j)
		{
			if (g_pulses[j].OMID == vOMIDdistance[i].OMID)
			{
				closeHits++;
				break;
			}
		}
	}

	// PrintGPulses();
	// cout << "closeHits: " << closeHits << endl;

	return closeHits;
}

int DoTheMagicMC(TChain* tree, BEvent* event, BEventMaskMC* eventMask, BSourceEAS* mcEvent)
{
	int nEntries = tree->GetEntries();
	int nEvents = 0;
	int nEnoughStringsEvents = 0;
	int nFirstFitEvents = 0;
	int nGoodErrorEvents = 0;
	int nPassedEvents = 0;
	int nCutPassingEvents = 0;
	int nLogEvents = 0;

	mcCascade* cascade;

	cout << "Number of processed entries: " << nEntries << endl;

	for (int i = 0; i < nEntries; ++i)
	{
		if (i%(nEntries/10) == 0)
		{
			cout << round((double)(i)/nEntries*100) << "% ";
			cout << std::flush;
		}
		tree->GetEntry(i);
		nEvents++;

		CaussalityFilter(event);
		if (CountStrings() > 2 && g_pulses.size() > 5)
		{
			nEnoughStringsEvents++;
			// h_nPassEvents->Fill(TMath::Log10(cascade->showerEnergy*1000));
		}else{
			// h_notRecEner->Fill(TMath::Log10(cascade->showerEnergy*1000));
			// h_notRecDist->Fill(TMath::Sqrt(TMath::Power(cascade->position[0],2)+TMath::Power(cascade->position[1],2)));
			continue;
		}

		double overallCharge = CountCharge(event);

		TVector3 matrixPosition(0,0,0);
		double matrixTime = 0;

		double initDist = EstimateInitialMatrixPositionMatrix(matrixPosition,matrixTime);
		fMinuit->SetFCN(chi2);

		double chi2QResult = FitMatrixPosition(matrixPosition,matrixTime,cascade);
		h_chi2->Fill(chi2QResult);
		// h_chi2vsDist->Fill(chi2QResult,TMath::Sqrt(TMath::Power(cascade->position[0]-matrixPosition.X(),2)+TMath::Power(cascade->position[1]-matrixPosition.Y(),2)+TMath::Power(cascade->position[2]-matrixPosition.Z(),2)));
		if (chi2QResult > 10.83)
		{
			// h_fitDist->Fill(TMath::Sqrt(TMath::Power(cascade->position[0]-matrixPosition.X(),2)+TMath::Power(cascade->position[1]-matrixPosition.Y(),2)+TMath::Power(cascade->position[2]-matrixPosition.Z(),2)));
			// h_notRecEner->Fill(TMath::Log10(cascade->showerEnergy*1000));
			// h_notRecDist->Fill(TMath::Sqrt(TMath::Power(cascade->position[0],2)+TMath::Power(cascade->position[1],2)));
			continue;
		}
		nFirstFitEvents++;

		double vertQRatio = GetVertQRatio(matrixPosition);
		h_vertQRatioUpGoing->Fill(vertQRatio);

		if (matrixPosition.Z() > 252 || CountCascadeCharge()/(overallCharge)*100 < 80 || CountBranchFilter(matrixPosition) < 10 || CountCloseHitsFilter(matrixPosition) < 10)
			continue;
		if (overallCharge < 500 || g_pulses.size() < 70)
			continue;
		nCutPassingEvents++;

		double cascadeEnergy = TMath::Power(10,3.19272+0.0443001*g_pulses.size()-0.0001266*g_pulses.size()*g_pulses.size())/1000;
		double cascadeTheta = 0;
		double cascadePhi = 0;
		double likelihood = 0;


		bool passedGoodError = false;

		if (EstimateInitialDirectionMC(matrixPosition,matrixTime,cascadeEnergy,cascadeTheta,cascadePhi) == 0)
			nGoodErrorEvents++;
		LikelihoodFilterPassed(matrixPosition,matrixTime,cascadeEnergy,cascadeTheta,cascadePhi,likelihood,passedGoodError);

		h_likelihood->Fill(likelihood);

		if (likelihood > 5)
				continue;
			nLogEvents++;



		h_zFilter->Fill(matrixPosition.Z());
		h_qFilter->Fill(CountCascadeCharge()/(overallCharge)*100);
		h_branchFilter->Fill(CountBranchFilter(matrixPosition));
		h_closeHitsFilter->Fill(CountCloseHitsFilter(matrixPosition));
	}

	cout << endl << "Number of entries in the files: " << nEntries << endl;
	cout << "Number of processed events: " << nEvents << endl;
	cout << "Number of enough string events: " << nEnoughStringsEvents << endl;
	cout << "Number of first fit events: " << nFirstFitEvents << endl;
	cout << "Number of good error events: " << nGoodErrorEvents << endl;
	cout << "Number of events passing cuts: " << nCutPassingEvents << endl;
	cout << "Number of events with good likelihood: " << nLogEvents << endl;

	return 0;
}

int DoTheMagicMCCascades(TChain* tree, mcCascade* cascade, int noiseRateInkHz, int initEstTechnique, bool useNoise, bool useChi2, bool likelihoodMultiFit)
{
	int nEntries = tree->GetEntries();
	// int nEntries = 150000;
	int nEvents = 0;
	int nEnoughStringsEvents = 0;
	int nFirstFitEvents = 0;
	int nSecondFitEvents = 0;
	int nGoodErrorEvents = 0;
	int nPassedEvents = 0;
	int nCutPassingEvents = 0;
	int nLogEvents = 0;

	int containmentDistance = 140;
	for (int i = 0; i < nEntries; ++i)
	{
		if (i%(nEntries/10) == 0)
		{
			cout << round((double)(i)/nEntries*100) << "% ";
			cout << std::flush;
		}
		tree->GetEntry(i);

		// if (!IsContained(cascade,containmentDistance) || cascade->showerEnergy < 0 || cascade->showerEnergy > 1000 || i % 250 != 0)
		if (!IsUncontained(cascade,60,200) || cascade->showerEnergy < 0 || cascade->showerEnergy > 10000 || i % 1000 != 0)
		// if (!IsContained(cascade,containmentDistance) || cascade->showerEnergy < 0 || cascade->showerEnergy > 10000)
		{
			continue;
		}
		h_notRecEnerAll->Fill(TMath::Log10(cascade->showerEnergy*1000));
		h_notRecDistAll->Fill(TMath::Sqrt(TMath::Power(cascade->position[0],2)+TMath::Power(cascade->position[1],2)));
		h_notRecThetaAll->Fill((TMath::Pi()-TMath::ACos(cascade->cosTheta))/TMath::Pi()*180);
		nEvents++;
		h_zFilterAll->Fill(cascade->position[2]);
		h_nEvents->Fill(TMath::Log10(cascade->showerEnergy*1000));
		if (useNoise)
			GenerateNoise(cascade,noiseRateInkHz);
		CaussalityFilter(cascade);
		h_xyPosition->Fill(cascade->position[0],cascade->position[1]);
		double overallChargeNoise = 0;
		double overallCharge = CountCharge(cascade,overallChargeNoise);
		h_overallCharge->Fill(overallCharge);
		h_overallChargeWithNoise->Fill(overallCharge+overallChargeNoise);
		h_overallChargeVsEnergy->Fill(TMath::Log10(overallCharge),TMath::Log10(cascade->showerEnergy*1000));
		h_energyEstTotalQ->Fill((TMath::Power(10,1.96907+0.828567*TMath::Log10(overallCharge)))/(cascade->showerEnergy*1000));
		h_energyEstTotalQVsEnergy->Fill(TMath::Log10(cascade->showerEnergy*1000),(TMath::Power(10,1.96907+0.828567*TMath::Log10(overallCharge)))/(cascade->showerEnergy*1000));
		h_nHitStrings->Fill(CountStrings());
		if (CountStrings() > 2 && g_pulses.size() > 5)
		{
			nEnoughStringsEvents++;
			h_nPassEvents->Fill(TMath::Log10(cascade->showerEnergy*1000));
			if (CheckPurity())
				h_nPureEvents->Fill(TMath::Log10(cascade->showerEnergy*1000));
		}else{
			h_notRecEner->Fill(TMath::Log10(cascade->showerEnergy*1000));
			h_notRecDist->Fill(TMath::Sqrt(TMath::Power(cascade->position[0],2)+TMath::Power(cascade->position[1],2)));
			h_notRecTheta->Fill((TMath::Pi()-TMath::ACos(cascade->cosTheta))/TMath::Pi()*180);
			continue;
		}

		int nHitsAfterCauss = g_pulses.size();
		h_nHitsAfterCauss->Fill(nHitsAfterCauss);

		FillPurityEfficiency(cascade,true);

		h_nHits->Fill(cascade->nHits);
		h_nHitsWithNoise->Fill(cascade->nHits+cascade->nNoiseHits);
		h_nHitsQover3p5->Fill(CountHitsWithQHigherThan(cascade));

		TVector3 matrixPosition(0,0,0);
		double matrixTime = 0;

		double initDist = 0;
		double cascadeEnergy = TMath::Power(10,3.30123+0.0447574*g_pulses.size()-0.000135729*g_pulses.size()*g_pulses.size())/1000;
		double cascadeTheta = 0;
		double cascadePhi = 0;
		double likelihood = 0;

		bool passedGoodError = false;
		// cout << cascade->position[0] << " " << cascade->position[1] << " " << cascade->position[2] << endl;

		double cascadeEnergyTrue = cascade->showerEnergy;
		double cascadeThetaTrue = TMath::Pi()-TMath::ACos(cascade->cosTheta);
		double cascadePhiTrue = cascade->phi+TMath::Pi();
		if (cascadePhiTrue > 2*TMath::Pi())
			cascadePhiTrue -= 2*TMath::Pi();

		if (initEstTechnique == 0)
			initDist = EstimateInitialMatrixPositionHQ(cascade,matrixPosition,matrixTime);
		if (initEstTechnique == 1)
			initDist = EstimateInitialMatrixPositionChi(cascade,matrixPosition,matrixTime);
		if (initEstTechnique == 2)
			initDist = EstimateInitialMatrixPositionMatrix(cascade,matrixPosition,matrixTime);

		h_initDist->Fill(initDist);
		h_initDistvsEnergy->Fill(TMath::Log10(cascade->showerEnergy*1000),initDist);

		if (useChi2)
			fMinuit->SetFCN(chi2);
		else
			fMinuit->SetFCN(MEstimator);

		double chi2QResult = FitMatrixPosition(matrixPosition,matrixTime,cascade);
		h_chi2->Fill(chi2QResult);
		h_chi2vsDist->Fill(chi2QResult,TMath::Sqrt(TMath::Power(cascade->position[0]-matrixPosition.X(),2)+TMath::Power(cascade->position[1]-matrixPosition.Y(),2)+TMath::Power(cascade->position[2]-matrixPosition.Z(),2)));

		double chi2Cut = 0;
		if (useChi2)
			chi2Cut = 10.83;
		else
			chi2Cut = 1000000;

		if (chi2QResult > chi2Cut)
		{
			h_fitDist->Fill(TMath::Sqrt(TMath::Power(cascade->position[0]-matrixPosition.X(),2)+TMath::Power(cascade->position[1]-matrixPosition.Y(),2)+TMath::Power(cascade->position[2]-matrixPosition.Z(),2)));
			h_notRecEner->Fill(TMath::Log10(cascade->showerEnergy*1000));
			h_notRecDist->Fill(TMath::Sqrt(TMath::Power(cascade->position[0],2)+TMath::Power(cascade->position[1],2)));
			h_notRecTheta->Fill((TMath::Pi()-TMath::ACos(cascade->cosTheta))/TMath::Pi()*180);
			continue;
		}
		nFirstFitEvents++;

		h_fitDistWell->Fill(TMath::Sqrt(TMath::Power(cascade->position[0]-matrixPosition.X(),2)+TMath::Power(cascade->position[1]-matrixPosition.Y(),2)+TMath::Power(cascade->position[2]-matrixPosition.Z(),2)));
		h_fitDist->Fill(TMath::Sqrt(TMath::Power(cascade->position[0]-matrixPosition.X(),2)+TMath::Power(cascade->position[1]-matrixPosition.Y(),2)+TMath::Power(cascade->position[2]-matrixPosition.Z(),2)));

		int nHitsAfterTFilter = TFilter(cascade,matrixPosition,matrixTime);
		h_nHitsAfterTFilter->Fill(nHitsAfterTFilter);
		h_nHitsChange->Fill(nHitsAfterTFilter-nHitsAfterCauss);
		h_nNoiseHits->Fill(GetNNoiseHits());
		if (nHitsAfterTFilter-nHitsAfterCauss < 0)
		{
			h_notRecEner->Fill(TMath::Log10(cascade->showerEnergy*1000));
			h_notRecDist->Fill(TMath::Sqrt(TMath::Power(cascade->position[0],2)+TMath::Power(cascade->position[1],2)));
			h_notRecTheta->Fill((TMath::Pi()-TMath::ACos(cascade->cosTheta))/TMath::Pi()*180);
			continue;
		}

		fMinuit->SetFCN(chi2);
		double chi2TResult = FitMatrixPosition(matrixPosition,matrixTime,cascade,false);
		h_chi2_T->Fill(chi2TResult);


		if (chi2TResult > 10.83)
		{
			h_notRecEner->Fill(TMath::Log10(cascade->showerEnergy*1000));
			h_notRecDist->Fill(TMath::Sqrt(TMath::Power(cascade->position[0],2)+TMath::Power(cascade->position[1],2)));
			h_notRecTheta->Fill((TMath::Pi()-TMath::ACos(cascade->cosTheta))/TMath::Pi()*180);
			continue;
		}
		nSecondFitEvents++;

		FillPurityEfficiency(cascade,false);
		FillVertDist(cascade);
		h_nHitsVsEnergy->Fill(g_pulses.size(),TMath::Log10(cascade->showerEnergy*1000));
		h_energyEstNHits->Fill((TMath::Power(10,3.30123+0.0447574*g_pulses.size()-0.000135729*g_pulses.size()*g_pulses.size()))/(cascade->showerEnergy*1000));
		h_energyEstNHitsVsEnergy->Fill(TMath::Log10(cascade->showerEnergy*1000),(TMath::Power(10,3.30123+0.0447574*g_pulses.size()-0.000135729*g_pulses.size()*g_pulses.size()))/(cascade->showerEnergy*1000));
		h_fitTDist->Fill(TMath::Sqrt(TMath::Power(cascade->position[0]-matrixPosition.X(),2)+TMath::Power(cascade->position[1]-matrixPosition.Y(),2)+TMath::Power(cascade->position[2]-matrixPosition.Z(),2)));

		TVector3 posTrue(cascade->position[0],cascade->position[1],cascade->position[2]);
		TVector3 posRec(matrixPosition.X(),matrixPosition.Y(),matrixPosition.Z());
		TVector3 posDiff = posRec-posTrue;
		TVector3 cascDirTrue(0,0,1);
		cascDirTrue.SetTheta(cascadeThetaTrue);
		cascDirTrue.SetPhi(cascadePhiTrue);

		h_mismatchPosPerpvsEnergy->Fill(TMath::Log10(cascadeEnergyTrue*1000),posDiff.Perp(cascDirTrue));
		h_mismatchPosPerp->Fill(posDiff.Perp(cascDirTrue));
		h_mismatchPosLongvsEnergy->Fill(TMath::Log10(cascadeEnergyTrue*1000),posDiff*cascDirTrue);
		h_mismatchPosLong->Fill(posDiff*cascDirTrue);

		double vertQRatio = GetVertQRatio(matrixPosition);
		if (TMath::Pi()-TMath::ACos(cascade->cosTheta) < TMath::Pi()/2)
			h_vertQRatioUpGoing->Fill(vertQRatio);
		else
			h_vertQRatioDownGoing->Fill(vertQRatio);

		// if (matrixPosition.Z() > 252 || CountCascadeCharge()/(overallChargeNoise+overallCharge)*100 < 80 || CountBranchFilter(matrixPosition) < 10 || CountCloseHitsFilter(matrixPosition) < 10)
		// 	continue;
		// if (overallChargeNoise+overallCharge < 500 || g_pulses.size() < 70)
		// 	continue;
		nCutPassingEvents++;

		if (1)
		{
			// cout << i << " "  << cascadeThetaTrue << " " << cascadePhiTrue << endl;
			if (!likelihoodMultiFit)
			{
				if (EstimateInitialDirection(cascade,matrixPosition,matrixTime,cascadeEnergy,cascadeTheta,cascadePhi) == 0)
					nGoodErrorEvents++;
				LikelihoodFilterPassed(matrixPosition,matrixTime,cascadeEnergy,cascadeTheta,cascadePhi,likelihood,passedGoodError);
			}else
			{
				LikelihoodFilterPassedGrid(matrixPosition,matrixTime,cascadeEnergy,cascadeTheta,cascadePhi,likelihood,passedGoodError);
			}


			// if (!passedGoodError)
			// 	continue;
			// nGoodErrorEvents++;


			TVector3 cascDirRec(0,0,1);
			cascDirRec.SetTheta(cascadeTheta);
			cascDirRec.SetPhi(cascadePhi);

			if (likelihood > 3)
			{
				h_mismatchAngleWell->Fill(cascDirRec.Angle(cascDirTrue)/TMath::Pi()*180);
				// continue;
			}
			nLogEvents++;


			h_mismatchAngle->Fill(cascDirRec.Angle(cascDirTrue)/TMath::Pi()*180);
			switch (GetNNoiseHits()) {
				case 0:
					h_mismatchAngle0Noise->Fill(cascDirRec.Angle(cascDirTrue)/TMath::Pi()*180);
					break;
				case 1:
					h_mismatchAngle1Noise->Fill(cascDirRec.Angle(cascDirTrue)/TMath::Pi()*180);
					break;
				default:
					h_mismatchAngle2Noise->Fill(cascDirRec.Angle(cascDirTrue)/TMath::Pi()*180);
					break;
			}
			h_mismatchAnglevsEnergy->Fill(TMath::Log10(cascadeEnergyTrue*1000),cascDirRec.Angle(cascDirTrue)/TMath::Pi()*180);
			h_mismatchTheta->Fill((cascadeThetaTrue-cascadeTheta)/TMath::Pi()*180);
			double phiDiff = cascadePhiTrue-cascadePhi;
			if (phiDiff > TMath::Pi())
				phiDiff = 2*TMath::Pi()-phiDiff;
			if (phiDiff < -TMath::Pi())
				phiDiff = -2*TMath::Pi()-phiDiff;
			h_mismatchPhi->Fill((phiDiff)/TMath::Pi()*180);
			h_likelihood->Fill(likelihood);

			h_mismatchEnergy->Fill(cascadeEnergy/cascadeEnergyTrue);
			h_mismatchLogEnergy->Fill(TMath::Log10(cascadeEnergy/cascadeEnergyTrue));
			h_mismatchEnergyvsEnergy->Fill(TMath::Log10(cascadeEnergyTrue*1000),cascadeEnergy/cascadeEnergyTrue);

			h_mcVsRecoTheta->Fill(cascadeThetaTrue/TMath::Pi()*180,cascadeTheta/TMath::Pi()*180);
			h_mcVsRecoPhi->Fill(cascadePhiTrue/TMath::Pi()*180,cascadePhi/TMath::Pi()*180);
			h_mcVsRecoEnergy->Fill(TMath::Log10(cascadeEnergyTrue*1000),TMath::Log10(cascadeEnergy*1000));

		}

		h_zFilter->Fill(matrixPosition.Z());
		h_qFilter->Fill(CountCascadeCharge()/(overallChargeNoise+overallCharge)*100);
		h_branchFilter->Fill(CountBranchFilter(matrixPosition));
		h_closeHitsFilter->Fill(CountCloseHitsFilter(matrixPosition));

	}

	cout << endl << "Number of entries in the files: " << nEntries << endl;
	cout << "Number of processed events: " << nEvents << endl;
	cout << "Number of enough string events: " << nEnoughStringsEvents << endl;
	cout << "Number of first fit events: " << nFirstFitEvents << endl;
	cout << "Number of second fit events: " << nSecondFitEvents << endl;
	cout << "Number of good error events: " << nGoodErrorEvents << endl;
	cout << "Number of events passing cuts: " << nCutPassingEvents << endl;
	cout << "Number of events with good likelihood: " << nLogEvents << endl;

	return 0;
}

int ReadGeometryMCCascades()
{
	TString fileName = "/Data/BaikalData/mc/cascades/array2016_phys_kebkal.dat";

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

int ReadNoiseChargeTable()
{
	ifstream inputFile("/Data/BaikalData/mc/cascades/mc-noise-charge-2016.dat", ios::in);

   	if (!inputFile)
    {
    	cerr << "Calibration file: " << "/Data/BaikalData/showerTable/mc-noise-charge-2016.dat" << " was NOT found. Program termination!" << endl;
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

// gLogTable4D structure: Rho, Z, Phi, cosTheta
int ReadLogTable()
{
	// cout << "4D LogTable reading starts" << endl;
	ifstream fTab;
	fTab.open("/Data/BaikalData/mc/cascades/hq001200_double.dqho2011", ios::in|ios::binary|ios::ate);

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

void DrawHistograms()
{
	TCanvas* c_overallCharge = new TCanvas("c_overallCharge","Results",800,600);
	h_overallCharge->Draw();
	h_overallChargeWithNoise->SetLineColor(kGreen);
	h_overallChargeWithNoise->Draw("same");

	TCanvas* c_nHitStrings = new TCanvas("c_nHitStrings","Results",800,600);
	h_nHitStrings->Draw();

	TCanvas* c_nHits = new TCanvas("c_nHits","Results",800,600);
	h_nHits->Draw();
	h_nHitsWithNoise->SetLineColor(kGreen);
	h_nHitsWithNoise->Draw("same");

	TCanvas* c_nHitsQOver3p5 = new TCanvas("c_nHitsQOver3p5","Results",800,600);
	h_nHitsQover3p5->Draw();

	TCanvas* c_initDist= new TCanvas("c_initDist","Results",800,600);
	h_initDist->Draw();

	TCanvas* c_initDistvsEnergy = new TCanvas("c_initDistvsEnergy","Results",800,600);
	h_initDistvsEnergy->ProfileX()->Draw();

	TCanvas* c_coordinaDiff = new TCanvas("c_coordinaDiff","Results",800,600);
	c_coordinaDiff->Divide(2,2);
	c_coordinaDiff->cd(1);
	h_xDiffCoor->Draw();
	c_coordinaDiff->cd(2);
	h_yDiffCoor->Draw();
	c_coordinaDiff->cd(3);
	h_zDiffCoor->Draw();
	c_coordinaDiff->cd(4);
	h_timeDiff->Draw();

	TCanvas* c_caussNHits = new TCanvas("c_caussNHits","Results",800,600);
	h_caussNHits->Draw("colz");

	TCanvas* c_caussEffi = new TCanvas("c_caussEffi","Results",800,600);
	h_caussEffi->Draw();

	TCanvas* c_caussPur = new TCanvas("c_caussPur","Results",800,600);
	h_caussPur->Draw();

	TCanvas* c_tEffi = new TCanvas("c_tEffi","Results",800,600);
	h_tEffi->Draw();

	TCanvas* c_tPur = new TCanvas("c_tPur","Results",800,600);
	h_tPur->Draw();

	TCanvas* c_caussDist = new TCanvas("c_cuassDist","Results",800,600);
	h_caussDtDr->Draw("colz");

	TCanvas* c_caussDist2 = new TCanvas("c_caussDist2","Results",800,600);
	h_caussDtMinusDr->Draw();

	TCanvas* c_fitDist = new TCanvas("c_fitDist","Results",800,600);
	h_fitDist->Draw();
	h_fitDistWell->SetLineColor(kRed);
	h_fitDistWell->SetLineStyle(9);
	h_fitDistWell->Draw("same");

	TCanvas* c_fitCoordinaDiff = new TCanvas("c_fitCoordinaDiff","Mismatch coordinates after fit",800,600);
	c_fitCoordinaDiff->Divide(2,2);
	c_fitCoordinaDiff->cd(1);
	h_xFitDiffCoor->Draw();
	c_fitCoordinaDiff->cd(2);
	h_yFitDiffCoor->Draw();
	c_fitCoordinaDiff->cd(3);
	h_zFitDiffCoor->Draw();
	c_fitCoordinaDiff->cd(4);

	TCanvas* c_fitCoordinateError = new TCanvas("c_fitCoordinateError","Coordinate Errors after fit",800,600);
	c_fitCoordinateError->Divide(2,2);
	c_fitCoordinateError->cd(1);
	h_xFitError->Draw();
	c_fitCoordinateError->cd(2);
	h_yFitError->Draw();
	c_fitCoordinateError->cd(3);
	h_zFitError->Draw();
	c_fitCoordinateError->cd(4);
	// h_tFitError->Draw();
	h_xFitCoorVsError->Draw("colz");

	TCanvas* c_fitTDist = new TCanvas("c_fitTDist","Results",800,600);
	h_fitTDist->Draw();

	TCanvas* c_purity = new TCanvas("c_purity","Results",800,600);
	h_nPureEvents->Divide(h_nPassEvents);
	h_nPureEvents->Draw();

	TCanvas* c_efficiency = new TCanvas("c_efficiency","Results",800,600);
	h_nPassEvents->Divide(h_nEvents);
	h_nPassEvents->Draw();
	// h_nPassEvents->Draw("Same");

	TCanvas* c_efficiencyVsE = new TCanvas("c_efficiencyVsE","Results",800,600);
	h_efficiencyVsE->Draw("colz");
	efficiencyProfile = h_efficiencyVsE->ProfileX();
	efficiencyProfile->SetLineColor(kRed);
	efficiencyProfile->Draw("same");

	TCanvas* c_purityVsE = new TCanvas("c_purityVsE","Results",800,600);
	h_purityVsE->Draw("colz");
	purityProfile = h_purityVsE->ProfileX();
	purityProfile->SetLineColor(kRed);
	purityProfile->Draw("same");

	TCanvas* c_chi2 = new TCanvas("c_chi2","Results",800,600);
	h_chi2->Draw();

	TCanvas* c_chi2vsDist = new TCanvas("c_chi2vsDist","chi2vsDist",800,600);
	h_chi2vsDist->Draw("colz");

	TCanvas* c_chi2_T = new TCanvas("c_chi2_T","Results",800,600);
	h_chi2_T->Draw();

	TCanvas* c_nHitsAfter = new TCanvas("c_nHitsAfter","Results",800,600);
	h_nHitsAfterTFilter->Draw();
	h_nHitsAfterCauss->SetLineColor(kRed);
	h_nHitsAfterCauss->Draw("Same");

	TCanvas* c_nHitsChange = new TCanvas("c_nHitsChange","Results",800,600);
	h_nHitsChange->Draw();

	TCanvas* c_vertDist = new TCanvas("c_vertDist","Results",800,600);
	h_vertDistHit->Draw();
	h_vertDistNoiseHit->SetLineColor(kRed);
	h_vertDistNoiseHit->Draw("same");

	TCanvas* c_notRecDist = new TCanvas("c_notRecDist","Results",800,600);
	h_notRecDistAll->Sumw2();
	h_notRecDist->Sumw2();
	h_notRecDistAll->SetLineColor(kRed);
	TRatioPlot* rp_notRecDist = new TRatioPlot(h_notRecDist,h_notRecDistAll,"");
	rp_notRecDist->SetH2DrawOpt("Hist");
	rp_notRecDist->Draw("nogrid");
	rp_notRecDist->GetUpperRefYaxis()->SetRangeUser(0,h_notRecDistAll->GetMaximum()*1.1);
	rp_notRecDist->GetUpperRefYaxis()->SetTitle("NoE [#]");
	rp_notRecDist->GetLowerRefYaxis()->SetTitle("Fraction");

	// h_notRecEnerAll->Draw();
	// h_notRecEner->SetLineColor(kRed);
	// h_notRecEner->Draw("same");
	TCanvas* c_notRecEner = new TCanvas("c_notRecEner","Results",800,600);
	h_notRecEnerAll->Sumw2();
	h_notRecEner->Sumw2();
	h_notRecEnerAll->SetLineColor(kRed);
	TRatioPlot* rp_notRecEner = new TRatioPlot(h_notRecEner,h_notRecEnerAll,"");
	rp_notRecEner->SetH2DrawOpt("Hist");
	rp_notRecEner->Draw("nogrid");
	rp_notRecEner->GetUpperRefYaxis()->SetRangeUser(0,h_notRecEnerAll->GetMaximum()*1.1);
	rp_notRecEner->GetUpperRefYaxis()->SetTitle("NoE [#]");
	rp_notRecEner->GetLowerRefYaxis()->SetTitle("Fraction");

	TCanvas* c_notRecTheta = new TCanvas("c_notRecTheta","Results",800,600);
	h_notRecThetaAll->Sumw2();
	h_notRecTheta->Sumw2();
	h_notRecThetaAll->SetLineColor(kRed);
	TRatioPlot* rp_notRecTheta = new TRatioPlot(h_notRecTheta,h_notRecThetaAll,"");
	rp_notRecTheta->SetH2DrawOpt("Hist");
	rp_notRecTheta->Draw("nogrid");
	rp_notRecTheta->GetUpperRefYaxis()->SetRangeUser(0,h_notRecThetaAll->GetMaximum()*1.1);
	rp_notRecTheta->GetUpperRefYaxis()->SetTitle("NoE [#]");
	rp_notRecTheta->GetLowerRefYaxis()->SetTitle("Fraction");
	gPad->BuildLegend(0.75,0.75,0.95,0.95,"");

	TCanvas* c_vertQRatio = new TCanvas("c_vertQRatio","vertQRatio",800,600);
	c_vertQRatio->Divide(1,2);
	c_vertQRatio->cd(1);
	h_vertQRatioDownGoing->DrawNormalized();
	h_vertQRatioUpGoing->SetLineColor(kRed);
	h_vertQRatioUpGoing->DrawNormalized("same");

	c_vertQRatio->cd(2);
	// h_vertQRatioDownGoing->DrawNormalized();
	h_vertQRatioDownGoing->Scale(1/h_vertQRatioDownGoing->Integral());
	h_vertQRatioDownGoing->GetCumulative(false)->Draw();
	h_vertQRatioUpGoing->SetLineColor(kRed);
	h_vertQRatioUpGoing->Scale(1/h_vertQRatioUpGoing->Integral());
	// h_vertQRatioUpGoing->DrawNormalized("Same");
	h_vertQRatioUpGoing->GetCumulative()->Draw("Same");

	TCanvas* c_mcVsEstParameters = new TCanvas("c_mcVsEstParameters","MCvsEstimatedParameters",800,600);
	c_mcVsEstParameters->Divide(2,2);
	c_mcVsEstParameters->cd(1);
	h_mcVsEstTheta->Draw("colz");
	c_mcVsEstParameters->cd(2);
	h_mcVsEstPhi->Draw("colz");
	c_mcVsEstParameters->cd(3);
	h_mcVsEstEnergy->Draw("colz");
	c_mcVsEstParameters->cd(4);
	// h_energyDist->Draw();
	// h_energyDistWrong->SetLineColor(kRed);
	// h_energyDistWrong->Draw("same");
	h_energyDistWrong->Divide(h_energyDist);
	h_energyDistWrong->Draw();

	TCanvas* c_MCminusEstimatedParameters = new TCanvas("c_MCminusEstimatedParameters","MCminusEstimatedParameters",800,600);
	c_MCminusEstimatedParameters->Divide(2,2);
	c_MCminusEstimatedParameters->cd(1);
	h_estimateTheta->Draw();
	c_MCminusEstimatedParameters->cd(2);
	h_estimatePhi->Draw();
	c_MCminusEstimatedParameters->cd(3);
	h_estimateEnergy->Draw();
	h_estimateEnergyWrong->SetLineColor(kRed);
	h_estimateEnergyWrong->Draw("same");
	c_MCminusEstimatedParameters->cd(4);
	// h_thetaDist->Draw();
	// h_thetaDistWrong->SetLineColor(kRed);
	// h_thetaDistWrong->Draw("same");
	// h_phiDist->SetLineColor(kGreen);
	// h_phiDist->Draw("same");
	// h_phiDistWrong->SetLineColor(kYellow);
	// h_phiDistWrong->Draw("same");
	h_thetaDistWrong->Divide(h_thetaDist);
	h_thetaDistWrong->Draw();
	h_phiDistWrong->Divide(h_phiDist);
	h_phiDistWrong->SetLineColor(kRed);
	h_phiDistWrong->Draw("same");

	Double_t x, q;
	q = 0.5;

	for (int i = 0; i < 99; ++i)
	{
		// h_mismatchEnergyvsEnergyMedian->Fill((i)/10.0,h_mismatchEnergyvsEnergy->ProjectionY("_py",i+1,i+1)->GetMean());
		TH1D* tempHist = h_mismatchAnglevsEnergy->ProjectionY("_py",i+1,i+1);

		if (tempHist->GetEntries() > 0)
		{
			tempHist->ComputeIntegral();
			// cout << tempHist->GetQuantiles(1, &x, &q) << endl;
			tempHist->GetQuantiles(1, &x, &q);
			// h_mismatchAnglevsEnergyMedian->Fill(h_mismatchAnglevsEnergyMedian->GetBinCenter(i+1),x);
			h_mismatchAnglevsEnergyMedian->SetBinContent(i,x);
		}
		delete tempHist;
	}
	h_mismatchAnglevsEnergyMedian->Sumw2(1);

	TCanvas* c_mismatchAngleVsEnergy = new TCanvas("c_mismatchAngleVsEnergy","MismatchAngleVsEnergy",800,600);
	h_mismatchAnglevsEnergy->ProfileX()->Draw();
	h_mismatchAnglevsEnergyMedian->SetLineColor(kRed);
	h_mismatchAnglevsEnergyMedian->Draw("same");

	TCanvas* c_MCminusRecoParameters = new TCanvas("c_MCminusRecoParameters","c_MCminusRecoParameters",800,600);
	c_MCminusRecoParameters->Divide(2,2);
	c_MCminusRecoParameters->cd(1);
	h_mismatchTheta->Draw();
	c_MCminusRecoParameters->cd(2);
	h_mismatchPhi->Draw();
	c_MCminusRecoParameters->cd(3);
	h_mismatchEnergy->Draw();
	h_mismatchLogEnergy->SetLineColor(kRed);
	h_mismatchLogEnergy->Draw("Same");
	c_MCminusRecoParameters->cd(4);
	h_mismatchAngle->Draw();
	h_mismatchAngleWell->SetLineColor(kRed);
	h_mismatchAngleWell->Draw("same");
	h_mismatchAngle->ComputeIntegral();
	h_mismatchAngle->GetQuantiles(1, &x,&q);

	cout << "Median mismatch angle: " << x << endl;


	// for (int i = 0; i < 99; ++i)
	// {
	// 	// h_mismatchEnergyvsEnergyMedian->Fill((i)/10.0,h_mismatchEnergyvsEnergy->ProjectionY("_py",i+1,i+1)->GetMean());
	// 	TH1D* tempHist = h_mismatchEnergyvsEnergy->ProjectionY("_py",i+1,i+1);
	// 	if (tempHist->GetEntries())
	// 	{
	// 		// cout << tempHist->GetQuantiles(1, &x, &q) << endl;
	// 		tempHist->ComputeIntegral();
	// 		tempHist->GetQuantiles(1, &x, &q);
	// 		h_mismatchEnergyvsEnergyMedian->Fill((i)/10.0,x);
	// 	}
	// 	delete tempHist;
	// }

	TCanvas* c_mismatchEnergyVsEnergy = new TCanvas("c_mismatchEnergyVsEnergy","MismatchEnergyVsEnergy",800,600);
	h_mismatchEnergyvsEnergy->ProfileX()->Draw();
	// h_mismatchEnergyvsEnergyMedian->SetLineColor(kRed);
	// h_mismatchEnergyvsEnergyMedian->Draw("same");

	TCanvas* c_mismatchPositions = new TCanvas("c_mismatchPositions","MismatchPositions",800,600);
	c_mismatchPositions->Divide(2,2);
	c_mismatchPositions->cd(1);
	h_mismatchPosLong->Draw();
	c_mismatchPositions->cd(2);
	h_mismatchPosPerp->Draw();
	c_mismatchPositions->cd(3);
	h_mismatchPosLongvsEnergy->ProfileX()->Draw();
	c_mismatchPositions->cd(4);
	h_mismatchPosPerpvsEnergy->ProfileX()->Draw();

	TCanvas* c_mcVsRecoParameters = new TCanvas("c_mcVsRecoParameters","MCvsReconstructedParameters",800,600);
	c_mcVsRecoParameters->Divide(2,2);
	c_mcVsRecoParameters->cd(1);
	h_mcVsRecoTheta->Draw("colz");
	c_mcVsRecoParameters->cd(2);
	h_mcVsRecoPhi->Draw("colz");
	c_mcVsRecoParameters->cd(3);
	h_mcVsRecoEnergy->Draw("colz");
	c_mcVsRecoParameters->cd(4);

	TCanvas* c_likelihoodFitErrors = new TCanvas("c_likelihoodFitErrors","LikelihoodFitErrors",800,600);
	c_likelihoodFitErrors->Divide(2,2);
	c_likelihoodFitErrors->cd(1);
	h_thetaError->Draw();
	c_likelihoodFitErrors->cd(2);
	h_phiError->Draw();
	c_likelihoodFitErrors->cd(3);
	h_likelihood->Draw();
	c_likelihoodFitErrors->cd(4);

	TCanvas* c_FilterDistributions = new TCanvas("c_FilterDistributions","FilterDistributions",800,600);
	c_FilterDistributions->Divide(2,2);
	c_FilterDistributions->cd(1);
	h_zFilter->DrawClone();
	h_zFilterAll->SetLineColor(kRed);
	h_zFilterAll->DrawClone("same");
	c_FilterDistributions->cd(2);
	h_qFilter->DrawClone();
	c_FilterDistributions->cd(3);
	h_branchFilter->DrawClone();
	c_FilterDistributions->cd(4);
	h_closeHitsFilter->DrawClone();

	TCanvas* c_FilterDistributionsCumulative = new TCanvas("c_FilterDistributionsCumulative","FilterDistributionsCumulative",800,600);
	c_FilterDistributionsCumulative->Divide(2,2);
	c_FilterDistributionsCumulative->cd(1);
	h_zFilter->Scale(1/h_zFilter->Integral());
	h_zFilterCum = (TH1F*)h_zFilter->GetCumulative();
	h_zFilterCum->Draw();
	c_FilterDistributionsCumulative->cd(2);
	h_qFilter->Scale(1/h_qFilter->Integral());
	h_qFilterCum = (TH1F*)h_qFilter->GetCumulative();
	h_qFilterCum->Draw();
	c_FilterDistributionsCumulative->cd(3);
	h_branchFilter->Scale(1/h_branchFilter->Integral());
	h_branchFilterCum = (TH1F*)h_branchFilter->GetCumulative();
	h_branchFilterCum->Draw();
	c_FilterDistributionsCumulative->cd(4);
	h_closeHitsFilter->Scale(1/h_closeHitsFilter->Integral());
	h_closeHitsFilterCum = (TH1F*)h_closeHitsFilter->GetCumulative();
	h_closeHitsFilterCum->Draw();

	TCanvas* c_FilterDistributionsCumulativeInverse = new TCanvas("c_FilterDistributionsCumulativeInverse","FilterDistributionsCumulativeInverse",800,600);
	c_FilterDistributionsCumulativeInverse->Divide(2,2);
	c_FilterDistributionsCumulativeInverse->cd(1);
	h_zFilterCumInv = (TH1F*)h_zFilter->GetCumulative(false);
	h_zFilterCumInv->Draw();
	c_FilterDistributionsCumulativeInverse->cd(2);
	h_qFilterCumInv = (TH1F*)h_qFilter->GetCumulative(false);
	h_qFilterCumInv->Draw();
	c_FilterDistributionsCumulativeInverse->cd(3);
	h_branchFilterCumInv = (TH1F*)h_branchFilter->GetCumulative(false);
	h_branchFilterCumInv->Draw();
	c_FilterDistributionsCumulativeInverse->cd(4);
	h_closeHitsFilterCumInv = (TH1F*)h_closeHitsFilter->GetCumulative(false);
	h_closeHitsFilterCumInv->Draw();

	TCanvas* c_nNoiseHits = new TCanvas("c_nNoiseHits","NumberNoiseHits",800,600);
	c_nNoiseHits->Divide(2,2);
	c_nNoiseHits->cd(1);
	h_nNoiseHits->Draw();
	c_nNoiseHits->cd(2);
	h_mismatchAngle0Noise->Draw();
	c_nNoiseHits->cd(3);
	h_mismatchAngle1Noise->Draw();
	c_nNoiseHits->cd(4);
	h_mismatchAngle2Noise->Draw();
}

void SaveHistograms(int noiseRateInkHz, int initEstTechnique, bool useNoise, bool useChi2, bool mcData = false, bool likeMultiFit = false)
{
	TFile* outputFile;

	// string filler{"_directionFull"};
	// string filler{"_direction"};
	string filler{""};

	string fitMethod;
	if (useChi2)
		fitMethod = "_chi2";
	else
		fitMethod = "_MEst";
	string noise;
	if (useNoise)
		noise = "";
	else
		noise = "_noNoise";
	string likelihood;
	if (likeMultiFit)
		likelihood = "_multiFit";
	else
		likelihood = "_singleFit";


	if (!mcData)
	{
		outputFile = new TFile(Form("/Data/BaikalData/mc/cascades/studyMCCascades_noiseRate%d_initEst%d%s%s%s%s.root",noiseRateInkHz,initEstTechnique,fitMethod.c_str(),noise.c_str(),likelihood.c_str(),filler.c_str()),"RECREATE");
		// if (useNoise)
		// {
		// 	if (useChi2)
		// 		outputFile = new TFile(Form("/Data/BaikalData/mc/DZH_cascades/studyMCCascades_noiseRate%d_initEst%d_chi2%s.root",noiseRateInkHz,initEstTechnique,filler.c_str()),"RECREATE");
		// 	else
		// 		outputFile = new TFile(Form("/Data/BaikalData/mc/DZH_cascades/studyMCCascades_noiseRate%d_initEst%d_MEst%s.root",noiseRateInkHz,initEstTechnique,filler.c_str()),"RECREATE");
		// }
		// else
		// {
		// 	if (useChi2)
		// 		outputFile = new TFile(Form("/Data/BaikalData/mc/DZH_cascades/studyMCCascades_noiseRate%d_initEst%d_chi2_noNoise%s.root",noiseRateInkHz,initEstTechnique,filler.c_str()),"RECREATE");
		// 	else
		// 		outputFile = new TFile(Form("/Data/BaikalData/mc/DZH_cascades/studyMCCascades_noiseRate%d_initEst%d_MEst_noNoise%s.root",noiseRateInkHz,initEstTechnique,filler.c_str()),"RECREATE");
		// }
	}else
	{
		outputFile = new TFile(Form("/Data/BaikalData/mc/2018may/studyMCCascades.root"),"RECREATE");
	}



	h_overallCharge->Write();
	h_overallChargeWithNoise->Write();
	h_overallChargeVsEnergy->Write();
	h_overallChargeVsEnergy->ProfileX()->Write();
	h_nHitsVsEnergy->Write();
	h_nHitsVsEnergy->ProfileX()->Write();
	h_nHitStrings->Write();
	h_nHits->Write();
	h_nHitsWithNoise->Write();
	h_nHitsQover3p5->Write();
	h_initDist->Write();
	h_initDistvsEnergy->ProfileX()->Write();
	h_xDiffCoor->Write();
	h_yDiffCoor->Write();
	h_zDiffCoor->Write();
	h_xFitDiffCoor->Write();
	h_yFitDiffCoor->Write();
	h_zFitDiffCoor->Write();
	h_timeDiff->Write();
	h_caussNHits->Write();
	h_caussEffi->Write();
	h_caussPur->Write();
	h_tEffi->Write();
	h_tPur->Write();
	h_caussDtDr->Write();
	h_caussDtMinusDr->Write();
	h_fitDist->Write();
	h_fitDistWell->Write();
	h_fitTDist->Write();
	h_nPureEvents->Write();
	h_nPassEvents->Write();
	h_efficiencyVsE->Write();
	efficiencyProfile->Write();
	h_purityVsE->Write();
	purityProfile->Write();
	h_noiseChargeDistribution->Write();
	h_noiseTimeDistribution->Write();
	h_chi2->Write();
	h_chi2vsDist->Write();
	h_chi2_T->Write();
	h_nHitsAfterCauss->Write();
	h_nHitsAfterTFilter->Write();
	h_nHitsChange->Write();
	h_vertDistHit->Write();
	h_vertDistNoiseHit->Write();

	h_estimateEnergy->Write();
	h_estimateTheta->Write();
	h_estimatePhi->Write();
	h_mcVsEstTheta->Write();
	h_mcVsEstPhi->Write();
	h_mcVsEstEnergy->Write();

	h_mcVsRecoTheta->Write();
	h_mcVsRecoPhi->Write();
	h_mcVsRecoEnergy->Write();

	h_likelihood->Write();
	h_mismatchTheta->Write();
	h_mismatchPhi->Write();
	h_mismatchEnergy->Write();
	h_mismatchLogEnergy->Write();
	h_mismatchEnergyvsEnergy->ProfileX()->Write();
	h_mismatchPosLong->Write();
	h_mismatchPosPerp->Write();
	h_mismatchPosPerpvsEnergy->Write();
	h_mismatchPosPerpvsEnergy->ProfileX()->Write();
	h_mismatchPosLongvsEnergy->Write();
	h_mismatchPosLongvsEnergy->ProfileX()->Write();;
	h_mismatchAngle->Write();
	h_mismatchAngleWell ->Write();
	h_mismatchAnglevsEnergy->Write();
	h_mismatchAnglevsEnergy->ProfileX()->Write();
	h_mismatchAnglevsEnergyMedian->Write();

	h_energyEstNHitsVsEnergy->Write();
	h_energyEstNHitsVsEnergy->ProfileX()->Write();
	h_energyEstNHitsVsEnergy->ProfileY()->Write();
	h_energyEstNHits->Write();
	h_energyEstTotalQVsEnergy->Write();
	h_energyEstTotalQVsEnergy->ProfileX()->Write();
	h_energyEstTotalQVsEnergy->ProfileY()->Write();
	h_energyEstTotalQ->Write();

	h_notRecDistAll->Write();
	h_notRecDist->Write();
	h_notRecEnerAll->Write();
	h_notRecEner->Write();
	h_notRecThetaAll->Write();
	h_notRecTheta->Write();

	h_vertQRatioUpGoing->Write();
	h_vertQRatioDownGoing->Write();

	h_zFilter->Write();
	h_zFilterCum->Write();
	h_zFilterCumInv->Write();
	h_qFilter->Write();
	h_qFilterCum->Write();
	h_qFilterCumInv->Write();
	h_branchFilter->Write();
	h_branchFilterCum->Write();
	h_branchFilterCumInv->Write();
	h_closeHitsFilter->Write();
	h_closeHitsFilterCum->Write();
	h_closeHitsFilterCumInv->Write();

	h_nNoiseHits->Write();
	delete outputFile;
}

int studyMCCascade(int noiseRateInkHz = 50, int initEstTechnique = 2, bool useNoise = true, bool chi2 = true, bool likelihoodMultiFit = false)
{
	clock_t begin = clock();
	TChain* mcFiles = new TChain("h11");
	mcFiles->Add("/Data/BaikalData/mc/cascades/ne16_tin_c*_00*.root");
	// mcFiles->Add("/Data/BaikalData/mc/DZH_tauCascades/ntau16_tin_c1_00*.root");

	mcCascade* cascade = new mcCascade;

	ReadGeometryMCCascades();
	ReadNoiseChargeTable();
	ReadLogTable();
	SetFitter();

	mcFiles->SetBranchAddress("L0",&cascade->eventID);
	mcFiles->SetBranchAddress("Esh",&cascade->showerEnergy);
	mcFiles->SetBranchAddress("cost",&cascade->cosTheta);
	mcFiles->SetBranchAddress("fj",&cascade->phi);
	mcFiles->SetBranchAddress("jch",&cascade->nHits);
	mcFiles->SetBranchAddress("xtr",cascade->position);
	mcFiles->SetBranchAddress("Npmt",cascade->chID);
	mcFiles->SetBranchAddress("tre",cascade->time);
	mcFiles->SetBranchAddress("are",cascade->charge);

	DoTheMagicMCCascades(mcFiles,cascade,noiseRateInkHz,initEstTechnique,useNoise,chi2,likelihoodMultiFit);

	DrawHistograms();
	SaveHistograms(noiseRateInkHz,initEstTechnique,useNoise,chi2,false,likelihoodMultiFit);

	clock_t end = clock();
	cout << endl << "Elapsed time : " << double(end - begin)/CLOCKS_PER_SEC << endl;

	return 0;
}

int ReadGeometryMC(TFile* file)
{
	TTree* tree = (TTree*)file->Get("ArrayConfig");

	BGeomTelMC* geomMC = NULL;
	tree->SetBranchAddress("BGeomTel.",&geomMC);

	int nOKOMs = 0;

	tree->GetEntry(0);
	for (int i = 0; i < geomMC->GetNumOMs(); ++i)
	{
		gOMpositions[i] = geomMC->At(i)->GetPosition() - (geomMC->At(270)->GetPosition() + TVector3(0,0,7.6));
		nOKOMs++;
		if ((i > 35 && i < 60) || (i > 71 && i < 84) || i == 130 || i == 131 || i == 245 || i == 246 || i == 247 || i == 256)
			gOMQCal[i] = -1;
	}
	return nOKOMs;
}

int studyMCData(int noiseRateInkHz = 50, int initEstTechnique = 2, bool useNoise = true, bool chi2 = true)
{
	clock_t begin = clock();

	TChain* mcFiles = new TChain("Events");
	const char* filePath = "/Data/BaikalData/mc/2018may/n_cors_n2m_cl2016_x1001.root";

	mcFiles->Add("/Data/BaikalData/mc/2018may/n_cors_n2m_cl2016_x*.root");
	// mcFiles->Add("/Data/BaikalData/mc/2018may/n_cors_n2m_cl2016_x10*.root");
	TFile* file = new TFile(filePath,"READ");

	if (ReadGeometryMC(file) == -1)
	{
		std::cout << "Problem with geometry file!" << std::endl;
		return -1;
	}
	if (ReadLogTable() == -1)
	{
		std::cout << "Problem with 4D LogLikelihood file!" << std::endl;
		return -1;
	}

	delete file;

	SetFitter();

    BEvent* event = NULL;
    mcFiles->SetBranchAddress("BEvent.",&event);
    BEventMaskMC* eventMask = NULL;
    mcFiles->SetBranchAddress("MCEventMask.",&eventMask);
    BSourceEAS* sourceEAS = NULL;
    mcFiles->SetBranchAddress("MCEventSource.",&sourceEAS);
    BMCEvent* mcEvent = NULL;
    mcFiles->SetBranchAddress("BMCEvent.",&mcEvent);

    mcFiles->GetEntry(0);

	DoTheMagicMC(mcFiles,event,eventMask,sourceEAS);

	DrawHistograms();
	SaveHistograms(noiseRateInkHz,initEstTechnique,useNoise,chi2,true);

	clock_t end = clock();
	cout << endl << "Elapsed time : " << double(end - begin)/CLOCKS_PER_SEC << endl;

	return 0;
}
#pragma once
#ifndef PARSEOPTS
#define PARSEOPTS

#include "string"
#include <vector>
#include "TVector3.h"
#include "TVirtualFitter.h"
#include "TRandom2.h"
#include "TTimeStamp.h"

struct EventStats
{
	int nEntries = 0;
	int nSixThrees = 0;
	int nNFilter = 0;
	// int nQFilter = 0;
	int nQFilterChi2 = 0;
	int nTFilter = 0;
	int nTFilterChi2 = 0;
	// int nZFilter = 0;
	// int nTDelayFilter = 0;
	// int nQRatioFilter = 0;
	// int nBranchFilter = 0;
	// int nCloseHitsFilter = 0;
	int nLikelihoodFilter = 0;
};

// structure necessary for the reading of the Zhan-Arys cascade simulations
struct mcCascade
{
	int eventID, nHits;
	float showerEnergy, cosTheta, phi, position[3], charge[288], time[288];
	unsigned short chID[288];
};

// structure holding unified hit
struct UnifiedHit
{
  int OMID;
  double time;
  double charge;
  double expectedCharge;
  bool noise;
  int MCflag;
};

// structure holding unified event format
struct UnifiedEvent
{
	int clusterID = -1;
	int runID = -1;
	int eventID = -1;
	std::vector<UnifiedHit> hits;
	double mcEnergy = -1;
	double mcTheta = -1;
	double mcPhi = -1;
	int mcNTrackHitsAfterTFilter = -1;
	TVector3 mcPosition;
	int nHits = -1;
	int nSignalHits = -1;
	int nNoiseHits = -1;
	double energy = -1;
	double energySigma = -1;
	double correctedEnergy = -1;
	double theta = -1;
	double thetaSigma = -1;
	double phi = -1;
	double phiSigma = -1;
	double directionSigma = -1;
	double declination = -1;
	double rightAscension = -1;
	TTimeStamp eventTime;
	TVector3 position;
	double time = -1;
	double mcFlagID = -1;
	double likelihood = -1;
	int nHitsAfterCaus = -1;
	int nHitsAfterTFilter = -1;
	int nStringsAfterCaus = -1;
	int nStringsAfterTFilter = -1;
	double chi2AfterCaus = -1;
	double chi2AfterTFilter = -1;
	double qTotal = -1;
};

extern int gNEventsProcessed;
extern int gEventID;
extern int gInputType;
extern std::string gProductionID;
extern std::string gFileInputFolder;

extern double gLogTable4D[200][351][21][21];
extern std::vector<double> gNoiseTable;
extern std::vector<double> gNoiseProbability;
extern int gNOMs;
extern const int gNStrings;
extern std::vector<TVector3> gOMpositions;
extern std::vector<double> gOMtimeCal;
extern std::vector<double> gOMqCal;
extern std::vector<UnifiedHit> gPulses; //global vector of PulsesVariables
extern int gVerbose;
extern TVector3 gReferencePosition;
extern const double gRecCinWater;
extern const double gRecC;
extern TVirtualFitter* gMinuit;
extern TRandom2 gRanGen;
extern double gMCMuTimeConstant;
extern double gMCNuTimeConstant;
extern double gMCNoiseRateInkHz;
extern const double gEnergyCorrectionArray[20];

extern int gNCut;
extern double gQTotalCut;
extern double gQCut;
extern double gCausQCut;
extern int gCausFriendCut;
extern double gQCutChi2;
extern double gTCutTimeWindowNs;
extern int gNCutT;
extern double gTCutChi2;
extern double gLikelihoodCut;
extern bool   gUseMultiDirFit;
extern bool gUseEOSRead;
extern bool gUseNewFolderStructure;
extern bool gUseNonHitLikelihoodTerm;
extern bool gUseNoiseHitLikelihoodTerm;
extern bool gUseChargeSatCorrection;

extern const double gLatDet;
extern const double gLonDet;

void parseOpts(int argc, char** argv);
void readRC(const char* rcpath);
void checkParams();

#endif /*PARSEOPTS*/

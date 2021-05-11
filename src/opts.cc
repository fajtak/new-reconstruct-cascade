#include "opts.h"

#include <iostream>
#include <getopt.h>
#include <cstdlib>
#include <stdlib.h>

#include <TEnv.h>

#include "BARS.h"

int gNEventsProcessed = -1;
int gEventID = -1;
std::string gProductionID = "";
int gInputType = 0;
std::string gFileInputFolder = "";

BSectionStatus* gSectionMask = nullptr;
TTree* gSectionMaskTree = nullptr;

double gLogTable4D[200][351][21][21]{0};
std::vector<double> gNoiseTable;
std::vector<double> gNoiseProbability;
const int gNOMs = 288;
const int gNStrings = 8;
std::vector<TVector3> gOMpositions(gNOMs);
std::vector<double> gOMtimeCal(gNOMs,0);
std::vector<double> gOMqCal(gNOMs,0);
std::vector<UnifiedHit> gPulses;
int gVerbose = 0;
TVector3 gReferencePosition;
const double gRecCinWater = 4.57;
const double gRecC = 3.336;
TVirtualFitter* gMinuit;
TRandom2 gRanGen;
double gMCMuTimeConstant = 3900; // seconds per MC file
double gMCNuTimeConstant = 4.4e9;
double gMCNoiseRateInkHz = 50; // in kHz
int gLikelihoodThetaSteps = 4; //
int gLikelihoodPhiSteps = 6;
int gLikelihoodEnergySteps = 4;
const double gEnergyCorrectionArray[20] = { 1.67773, 1.56738, 1.53488, 1.42292, 1.34681, 1.33311, 1.30409, 1.31283, 1.3094, 1.28182, 1.26567, 1.27433, 1.24962, 1.24205, 1.26317, 1.22895, 1.23192, 1.21746, 1.22675, 1.13126 };
double gSigmaMCTimeCal = -1;

int gNCut = -1;
double gQTotalCut = -1;
double gQCut = -1;
double gCausQCut = -1;
int gCausFriendCut = -1;
double gQCutChi2 = -1;
double gTCutQThreshold = -1;
double gTCutDistance = -1;
double gTCutTimeWindowNs = -1;
int gNCutT = -1;
int gNCutDiff = -1;
double gTCutChi2 = -1;
double gLikelihoodCut = -1;
double gZCut = -1;
int gDirFitType = -1;
bool gUseEOSRead = false;
bool gUseNewFolderStructure = false;
bool gUseNonHitLikelihoodTerm = false;
bool gUseNoiseHitLikelihoodTerm = false;
bool gUseChargeSatCorrection = false;
bool gSaveServiceInfo = false;
bool gLaserRun = false;
bool gExcludeTrackHits = false;


//coordinates of the detector converted to radians (assuming east longitude)
// taken from https://baikalforum.jinr.ru/t/detector-geographic-coordinates/726/3
const double gLatDet = TMath::Pi() * 51.764 / 180.0;  //about 3.6 km south from 107 km station
const double gLonDet = TMath::Pi() * 104.414 / 180.0; //about 3.6 km south from 107 km station

using namespace BARS;

// At the end of list you should add 'App::opt_NULL' element
static const struct App::ProgramOption_t options_list[]{
  // --help and --config options will be added automaticly
	{App::opt_Verbose, NOT_REQUIRED},
	{App::opt_Input,   NOT_REQUIRED},
	{App::opt_Output,  NOT_REQUIRED},
	{App::opt_Cluster, NOT_REQUIRED},
	{App::opt_Season,  NOT_REQUIRED},
	{App::opt_Run,     NOT_REQUIRED},
	{
		{
			"chargeSaturation", 'q',
			no_argument,
			"use chargeSaturationCorrection",
			[](char* argv) {gUseChargeSatCorrection = true;},
			[]() {;}
		},

		NOT_REQUIRED
	},
	{
		{
			"noiseHitLog", 'p',
			no_argument,
			"use noiseHitLikelihood term",
			[](char* argv) {gUseNoiseHitLikelihoodTerm = true;},
			[]() {;}
		},

		NOT_REQUIRED
	},
	{
		{
			"nonHitLog", 'g',
			no_argument,
			"use nonHitLikelihood term",
			[](char* argv) {gUseNonHitLikelihoodTerm = true;},
			[]() {;}
		},

		NOT_REQUIRED
	},
	{
		{
			"number", 'n',
			required_argument,
			"set number of events you want to process",
			[](char* argv) {gNEventsProcessed = atoi(argv);},
			[]() {;}
		},

		NOT_REQUIRED
	},
	{
		{
			"tag", 't',
			required_argument,
			"set production ID of the data that will be processed",
			[](char* argv) {gProductionID = argv;},
			[]() {;}
		},
		NOT_REQUIRED
	},
	{
		{
			"mcMu", 'm',
			no_argument,
			"use mc atmospheric muon bundles data as the input for the program",
			[](char* argv) {gInputType = 3;},
			[]() {;}
		},
		NOT_REQUIRED
	},
	{
		{
			"mcNu", 'u',
			no_argument,
			"use mc up-going single muons data as the input for the program",
			[](char* argv) {gInputType = 2;},
			[]() {;}
		},
		NOT_REQUIRED
	},
	{
		{
			"mcCas", 'a',
			no_argument,
			"use mc cascade data from Zhan-Arys as the input for the program",
			[](char* argv) {gInputType = 1;},
			[]() {;}
		},
		NOT_REQUIRED
	},
	{
		{
			"mcInputFolder", 'f',
			required_argument,
			"set input folder with MC data",
			[](char* argv) {gFileInputFolder = argv;},
			[]() {;}
		},
		NOT_REQUIRED
	},
	{
		{
			"eventID", 'e',
			required_argument,
			"set eventID that is going to be processed",
			[](char* argv) {gEventID = atoi(argv);},
			[]() {;}
		},

		NOT_REQUIRED
	},
	{
		{
			"eosRead", 'w',
			no_argument,
			"Uses root://eos.jinr.ru// read style",
			[](char* argv) {gUseEOSRead = true;},
			[]() {;}
		},

		NOT_REQUIRED
	},
	{
		{
			"newFolderStructure", 'l',
			no_argument,
			"Uses new folder structure from the middle of 2019 with clusterX/exp/joint/j01/",
			[](char* argv) {gUseNewFolderStructure = true;},
			[]() {;}
		},

		NOT_REQUIRED
	},
	{
		{
			"saveServiceInfo", 'b',
			no_argument,
			"Saves time residuals and charge saturation service txt files",
			[](char* argv) {gSaveServiceInfo = true;},
			[]() {;}
		},

		NOT_REQUIRED
	},
	{
		{
			"processLaserRun", 'z',
			no_argument,
			"Switches off production of visualizations and some cuts to speed up laser run processing",
			[](char* argv) {gLaserRun = true;},
			[]() {;}
		},

		NOT_REQUIRED
	},
	{
		{
			"excludeTracHitsMC", 'd',
			no_argument,
			"Switches off the reading of the track hits for MC in the TransformToUnified()",
			[](char* argv) {gExcludeTrackHits = true;},
			[]() {;}
		},

		NOT_REQUIRED
	},
	{
		{
			"generateMCTimeCal", 'j',
			required_argument,
			"set sigma for MC Time Calibrations",
			[](char* argv) {gSigmaMCTimeCal = atoi(argv);},
			[]() {;}
		},

		NOT_REQUIRED
	},
  // !!! must be at the end
	{App::opt_NULL, NOT_REQUIRED}
};

// Read resource file, fill parameters
//_____________________________________________________________________
void readRC(const char* rcpath)
{
	TEnv env;
	if (-1 == env.ReadFile(App::RCPath, kEnvAll)) {
		std::cerr << "ERROR: failed to read configuration file at '" << App::RCPath << "'" << std::endl;
		exit(1);
	}

	gNCut = env.GetValue("NCut", 50);
	gQTotalCut = env.GetValue("QTotalCut", 200.0);
	gQCut = env.GetValue("QCut", 100.0);
	gCausQCut = env.GetValue("CausQCut", 1.5);
	gCausFriendCut = env.GetValue("CausFriendCut", 2);
	gTCutTimeWindowNs = env.GetValue("TCutTimeWindow", 50.0);
	gTCutQThreshold = env.GetValue("TCutQThreshold",0.0);
	gTCutDistance = env.GetValue("TCutDistance",1000.0);
	gQCutChi2 = env.GetValue("QCutChi2", 100.0);
	gTCutChi2 = env.GetValue("TCutChi2", 20.0);
	gNCutT = env.GetValue("NCutT", 20);
	gNCutDiff = env.GetValue("NCutDiff",0);
	gZCut = env.GetValue("ZCut",250.0);
	gLikelihoodCut = env.GetValue("LikelihoodCut",3.0); //value 1 is used for the Non-hit and Noise-hit likelihood term
	gDirFitType = env.GetValue("DirFitType",1);
	gUseNonHitLikelihoodTerm = env.GetValue("NonHitLikelihoodTerm",false);
	gUseNoiseHitLikelihoodTerm = env.GetValue("NoiseHitLikelihoodTerm",false);
	gUseChargeSatCorrection = env.GetValue("ChargeSatCorrection",false);
	gMCMuTimeConstant = env.GetValue("MCMuTimeConstant",3900);
	gMCNuTimeConstant = env.GetValue("MCNuTimeConstant",4.4e9);
	gLikelihoodThetaSteps = env.GetValue("LikelihoodThetaSteps",4);
	gLikelihoodPhiSteps = env.GetValue("LikelihoodPhiSteps",6);
	gLikelihoodEnergySteps = env.GetValue("LikelihoodEnergySteps",4);
	gProductionID = env.GetValue("ProductionID","");
}

// Parse options passed to the application.
//_____________________________________________________________________
void parseOpts(int argc, char **argv)
{
	App::ParseProgramOptions(argc, argv, options_list);
}

//_____________________________________________________________________
void checkParams()
{
	App::CheckProgramOptions(options_list);
}

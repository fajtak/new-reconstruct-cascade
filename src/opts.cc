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

double gLogTable4D[200][351][21][21]{0};
std::vector<double> gNoiseTable;
int gNOMs = 288;
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
double gMCTimeConstant = 3900; // seconds per MC file
double gMCNoiseRateInkHz = 50; // in kHz

int gNCut = -1;
double gQTotalCut = -1;
double gQCut = -1;
double gCausQCut = -1;
int gCausFriendCut = -1;
double gQCutChi2 = -1;
double gTCutTimeWindowNs = -1;
int gNCutT = -1;
double gTCutChi2 = -1;
double gLikelihoodCut = -1;
bool gUseMultiDirFit = false;


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
			"mcCas", 'f',
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

	gNCut = env.GetValue("NCut", 70);
	gQTotalCut = env.GetValue("QTotalCut", 200.0);
	gQCut = env.GetValue("QCut", 100.0);
	gCausQCut = env.GetValue("CausQCut", 1.5);
	gCausFriendCut = env.GetValue("CausFriendCut", 2);
	gTCutTimeWindowNs = env.GetValue("TCutTimeWindow", 50.0);
	gQCutChi2 = env.GetValue("QCutChi2", 100.0);
	gTCutChi2 = env.GetValue("TCutChi2", 20.0);
	gNCutT = env.GetValue("NCutT", 20);
	gLikelihoodCut = env.GetValue("LikelihoodCut",3.0);
	gUseMultiDirFit = env.GetValue("MultiDirFit",true);
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

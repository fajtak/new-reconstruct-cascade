#include <iostream>
#include <fstream>
#include <string>

#include "TCanvas.h"
#include "TGraph.h"
// TH1F* h_procSpeed = new TH1F("h_procSpeed","Processing speed in kEvents/s per event; Run [#];kEvents/s [kHz]",600,0,600);
TGraph* g_nEventsPerRun = new TGraph();
TGraph* g_nEventsPerSecond = new TGraph();
TGraph* g_nNFilterEventsPerSecond = new TGraph();
TGraph* g_durationPerRun = new TGraph();
TGraph* g_procSpeed = new TGraph();
TGraph* g_realProcSpeed = new TGraph();
TGraph* g_realProcTime = new TGraph();
TGraph* g_realVsUserTime = new TGraph();
TGraph* g_nRecCasc = new TGraph();
TGraph* g_ratioNFilterNEvents = new TGraph();

int DrawResults()
{
	TCanvas* c_nEventsPerRun = new TCanvas("c_nEventsPerRun","EventsPerRun",800,600);
	g_nEventsPerRun->Draw("AP");
	g_nEventsPerRun->SetTitle("Number of Events per Run;RunID [#]; Events [#]");
	g_nEventsPerRun->SetMarkerStyle(5);
	g_nEventsPerRun->SetLineWidth(4);
	g_nEventsPerRun->SetMarkerColor(4);

	TCanvas* c_durationPerRun = new TCanvas("c_durationPerRun","DurationPerRun",800,600);
	g_durationPerRun->Draw("AP");
	g_durationPerRun->SetTitle("Run duration;RunID [#]; Duration [days]");
	g_durationPerRun->SetMarkerStyle(5);
	g_durationPerRun->SetLineWidth(4);
	g_durationPerRun->SetMarkerColor(3);

	TCanvas* c_nEventsPerSecond = new TCanvas("c_nEventsPerSecond","EventsPerSecond",800,600);
	g_nEventsPerSecond->Draw("AP");
	g_nEventsPerSecond->SetTitle("Number of Detected Events per Second;RunID [#]; Events/Second [# /s]");
	g_nEventsPerSecond->SetMarkerStyle(5);
	g_nEventsPerSecond->SetLineWidth(4);
	g_nEventsPerSecond->SetMarkerColor(2);

	TCanvas* c_nNFilterEventsPerSecond = new TCanvas("c_nNFilterEventsPerSecond","NFilterEventsPerSecond",800,600);
	g_nNFilterEventsPerSecond->Draw("AP");
	g_nNFilterEventsPerSecond->SetTitle("Number of NFilter Events per Second;RunID [#]; Events/Second [# /s]");
	g_nNFilterEventsPerSecond->SetMarkerStyle(5);
	g_nNFilterEventsPerSecond->SetLineWidth(4);
	g_nNFilterEventsPerSecond->SetMarkerColor(3);

	TCanvas* c_procSpeed = new TCanvas("c_procSpeed","ProcessingSpeed",800,600);
	g_procSpeed->Draw("AP");
	g_procSpeed->SetTitle("Processing Speed (Number of Processed Events per Second);RunID [#]; Events/Second [# /s]");
	g_procSpeed->SetMarkerStyle(5);
	g_procSpeed->SetLineWidth(4);
	g_procSpeed->SetMarkerColor(1);

	TCanvas* c_realProcSpeed = new TCanvas("c_realProcSpeed","RealProcessingSpeed",800,600);
	g_realProcSpeed->Draw("AP");
	g_realProcSpeed->SetTitle("Real Processing Speed (Number of Processed Events per Second);RunID [#]; Events/Second [# /s]");
	g_realProcSpeed->SetMarkerStyle(5);
	g_realProcSpeed->SetLineWidth(4);
	g_realProcSpeed->SetMarkerColor(1);

	TCanvas* c_realVsUserTime = new TCanvas("c_realVsUserTime","RealVsUserTime",800,600);
	g_realVsUserTime->Draw("AP");
	g_realVsUserTime->SetTitle("Real/User Processing Time;RunID [#]; Real/User [#]");
	g_realVsUserTime->SetMarkerStyle(5);
	g_realVsUserTime->SetLineWidth(4);
	g_realVsUserTime->SetMarkerColor(6);

	TCanvas* c_realProcTime = new TCanvas("c_realProcTime","RealProcessingTime",800,600);
	g_realProcTime->Draw("AP");
	g_realProcTime->SetTitle("Real Processing Time;RunID [#]; Time [minutes]");
	g_realProcTime->SetMarkerStyle(5);
	g_realProcTime->SetLineWidth(4);
	g_realProcTime->SetMarkerColor(1);

	TCanvas* c_nRecCasc = new TCanvas("c_nRecCasc","NumberReconstructedCascades",800,600);
	g_nRecCasc->Draw("AP");
	g_nRecCasc->SetTitle("Number of Reconstructed Cascades ;RunID [#]; N_{casc} [#]");
	g_nRecCasc->SetMarkerStyle(5);
	g_nRecCasc->SetLineWidth(4);
	g_nRecCasc->SetMarkerColor(4);

	TCanvas* c_ratioNFilterNEvents = new TCanvas("c_ratioNFilterNEvents","RatioNFilterNEvents",800,600);
	g_ratioNFilterNEvents->Draw("AP");
	g_ratioNFilterNEvents->SetTitle("Ratio of NFilter and NEvents ;RunID [#]; Ratio [%]");
	g_ratioNFilterNEvents->SetMarkerStyle(5);
	g_ratioNFilterNEvents->SetLineWidth(4);
	g_ratioNFilterNEvents->SetMarkerColor(6);

	return 0;
}

int programOutputTime(int year, int cluster, TString folderPath)
{
	// ifstream fileIn(Form("/Data/BaikalData/dataVal/exp%d/cluster%d/programOutput_%d_%d.log",year,cluster,year,cluster));
	ifstream fileIn(Form("%s/programOutput_%d_%d.log",folderPath.Data(),year,cluster));

	string oneLine;
	char oneChar;
	int nProcessedRuns = 0;
	int nShorterRuns = 0;
	int nNotJointRuns = 0;
	int nNotBranchesRuns = 0;
	int recentRun = 0;
	int nEvents = 0;
	int nEventsTotal = 0;
	double measTimeHours = 0;
	double measTimeDays = 0;
	double measTimeDaysTotal = 0;

	int nFilter = 0;
	int SixThreeFilter = 0;
	int QFilterChi2 = 0;
	int TFilter = 0;
	int TFilterChi2 = 0;
	int likelihoodFiter = 0;

	double elapsedTime = 0;
	double realTimeMin = 0;
	double realTimeSec = 0;
	double userTimeMin = 0;
	double userTimeSec = 0;
	double realTimeMinSum = 0;
	double userTimeMinSum = 0;

	while(!fileIn.eof())
	{
		fileIn >> oneLine;
		if (oneLine == "shorter")
			nShorterRuns++;
		if (oneLine == "was")
			nNotJointRuns++;
		if (oneLine == "branches")
			nNotBranchesRuns++;
		if (oneLine == "Run:")
		{
			nProcessedRuns++;
			fileIn >> recentRun;
			// cout << recentRun << endl;
			fileIn >> oneLine >> nEvents >> measTimeHours >> measTimeDays;
			// cout << nEvents << " " << measTimeDays << endl;

			if (measTimeDays > 0 && measTimeDays < 5)
			{
				;
			}else
			{
				measTimeDays = nEvents*(measTimeDaysTotal/nEventsTotal);
				measTimeHours = measTimeDays*24;
				// cout << recentRun << " " << measTimeHours << " "  << measTimeDays << endl;
			}

			if (nEvents == 0)
				continue;

			nEventsTotal += nEvents;
			measTimeDaysTotal += measTimeDays;

			g_nEventsPerRun->SetPoint(nProcessedRuns,recentRun,nEvents);
			g_nEventsPerSecond->SetPoint(nProcessedRuns,recentRun,nEvents/measTimeHours/3600);
			g_durationPerRun->SetPoint(nProcessedRuns,recentRun,measTimeDays);

			bool runEnded = false;
			while(!fileIn.eof() && !runEnded)
			{
				fileIn >> oneLine;
				if (oneLine == "Welcome")
				{
					cout << "Unfinished Run: " << recentRun << endl;
					break;
				}
				if (oneLine == "NFilter:")
				{
					fileIn >> nFilter;
				}
				if (oneLine == "SixThreeFilter:")
				{
					fileIn >> SixThreeFilter;
				}
				if (oneLine == "QFilterChi2:")
				{
					fileIn >> QFilterChi2;
				}
				if (oneLine == "TFilter:")
				{
					fileIn >> TFilter;
				}
				if (oneLine == "TFilterChi2:")
				{
					fileIn >> TFilterChi2;
				}
				if (oneLine == "LikelihoodFiter:")
				{
					fileIn >> likelihoodFiter;
				}
				if (oneLine == "time:")
				{
					fileIn >> elapsedTime;
				}
				if (oneLine == "real")
				{
					fileIn >> realTimeMin >> oneChar >> realTimeSec;
					realTimeMinSum += realTimeMin;
				}
				if (oneLine == "user")
				{
					fileIn >> userTimeMin >> oneChar >> userTimeSec;
					userTimeMinSum += userTimeMin;
					runEnded = true;
				}
			}
			g_nNFilterEventsPerSecond->SetPoint(nProcessedRuns,recentRun,nFilter/measTimeHours/3600);

			g_procSpeed->SetPoint(nProcessedRuns,recentRun,nEvents/elapsedTime);
			g_realProcSpeed->SetPoint(nProcessedRuns,recentRun,nEvents/(realTimeMin*60+realTimeSec));
			g_realVsUserTime->SetPoint(nProcessedRuns,recentRun,(realTimeMin*60+realTimeSec)/(userTimeMin*60+userTimeSec));
			g_realProcTime->SetPoint(nProcessedRuns,recentRun,(realTimeMin*60+realTimeSec)/60);
			g_nRecCasc->SetPoint(nProcessedRuns,recentRun,likelihoodFiter);
			g_ratioNFilterNEvents->SetPoint(nProcessedRuns,recentRun,(double)nFilter/nEvents*100);
		}
	}

	cout << "Number of Processed Runs: " << nProcessedRuns << " Number of Events [M#]: " << nEventsTotal/1000000 << " Measurement Time [days] : " << measTimeDaysTotal << endl;
	cout << "Real time processing [hours]: " << realTimeMinSum/60 << " User time processing [hours]: " << userTimeMinSum/60 << endl;
	cout << "Number of all runs: " << nProcessedRuns + nShorterRuns + nNotJointRuns + nNotBranchesRuns << endl;
	cout << "\tNumber of Processed runs: " << nProcessedRuns << endl;
	cout << "\tNumber of runs shorter than 2 hours: " << nShorterRuns << endl;
	cout << "\tNumber of runs without joint.events.root: " << nNotJointRuns << endl;
	cout << "\tNumber of runs without branches: " << nNotBranchesRuns << endl;
	// cout << "Number of Reconstructed Cascades: " << nRecCascTotal << endl;

	// h_procSpeed->Draw();
	DrawResults();

	return 0;
}

			// fileIn >> nEvents >> measTimeHours >> measTimeDays;
			// // cout << nEvents << " " << measTimeDays << endl;
			// if (measTimeDays < 0 || measTimeDays > 10)
			// {
			// 	measTimeDaysTotal += nEvents*(measTimeDaysTotal/nEventsTotal);
			// 	nEventsTotal += nEvents;

			// }else{
			// 	nEventsTotal += nEvents;
			// 	measTimeDaysTotal += measTimeDays;
			// }
			// bool elapsedTimeFound = false;
			// while(!fileIn.eof() && !elapsedTimeFound)
			// {
			// 	fileIn >> oneLine;
			// 	if (oneLine == "LikelihoodFiter:")
			// 	{
			// 		fileIn >> nRecCasc;
			// 		nRecCascTotal += nRecCasc;
			// 	}
			// 	if (oneLine == "time:")
			// 	{
			// 		fileIn >> elapsedTime;
			// 		cout << nEvents << " " << elapsedTime << " " << nEvents/elapsedTime/1000 << endl;
			// 		elapsedTimeFound = true;
			// 		h_procSpeed->SetBinContent(nProcessedRuns,nEvents/elapsedTime/1000);
			// 	}
			// }

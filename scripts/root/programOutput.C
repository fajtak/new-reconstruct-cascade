#include <iostream>
#include <fstream>
#include <string>

TH1F* h_procSpeed = new TH1F("h_procSpeed","Processing speed in kEvents/s per event; Run [#];kEvents/s [kHz]",600,0,600);

int programOutput(int year, int cluster)
{
	ifstream fileIn(Form("/Data/BaikalData/dataVal/exp%d/cluster%d/programOutput_%d_%d.log",year,cluster,year,cluster));
	
	string oneLine;
	int nEventsTotal = 0;
	int nEvents = 0;
	double measTimeHoursTotal = 0;
	double measTimeHours = 0;
	double measTimeDaysTotal = 0;
	double measTimeDays = 0;
	double elapsedTime = 0;
	int nProcessedRuns = 0;
	int nRecCascTotal = 0;
	int nRecCasc = 0;

	while(!fileIn.eof())
	{
		fileIn >> oneLine;
		if (oneLine == "!")
		{
			nProcessedRuns++;
			fileIn >> nEvents >> measTimeHours >> measTimeDays;
			// cout << nEvents << " " << measTimeDays << endl;
			if (measTimeDays < 0 || measTimeDays > 10)
			{
				measTimeDaysTotal += nEvents*(measTimeDaysTotal/nEventsTotal);
				nEventsTotal += nEvents;

			}else{
				nEventsTotal += nEvents;
				measTimeDaysTotal += measTimeDays;				
			}
			bool elapsedTimeFound = false;
			while(!fileIn.eof() && !elapsedTimeFound)
			{
				fileIn >> oneLine;
				if (oneLine == "LikelihoodFiter:")
				{
					fileIn >> nRecCasc;
					nRecCascTotal += nRecCasc;
				}
				if (oneLine == "time:")
				{
					fileIn >> elapsedTime;
					cout << nEvents << " " << elapsedTime << " " << nEvents/elapsedTime/1000 << endl;
					elapsedTimeFound = true;
					h_procSpeed->SetBinContent(nProcessedRuns,nEvents/elapsedTime/1000);
				}
			}
		}
	}

	cout << "Number of Processed Runs: " << nProcessedRuns << " Number of Events: " << nEventsTotal << " Measurement Time [days] : " << measTimeDaysTotal << endl;
	cout << "Number of Reconstructed Cascades: " << nRecCascTotal << endl;

	h_procSpeed->Draw();

	return 0;
}

// #include "BJointMark.h"

const int c_nPed = 5;

const int gNOMs = 288;
const int gNSecs = 25;
double gOMtimeCal[gNOMs] = {0};
double gOMchargeCal[gNOMs] = {0};
TVector3 gOMpositions[gNOMs];
int gSecEntryIDs[gNSecs] = {0};

struct Hit
{
	int OMID;
	double time;
	double expTime;
	double expTimeTrack;
	double expDistTrack;
	double charge;
	double expCharge;
	int eventID;
	TMultiGraph* waveform;
};

vector<Hit> hits;

TString GetSeasonID(int year)
{
	switch (year){
		case 2016:
			return "f";
			break;
		case 2017:
			return "g";
			break;
		case 2018:
			return "h";
			break;
		case 2019:
			return "i";
			break;
		case 2020:
			return "j";
			break;
		default:
			cerr << "Wrong year set!" << endl;
			return "0";
			break;
	}
}

BGeomTel* GetGeomTel(int year, int cluster)
{
	BGeomTel* geomTel = NULL;
	switch (year)
	{
		case 2016:
			switch (cluster)
			{
				case 1:
					geomTel = new BGeomTel2016();
					break;
			}
			break;
		case 2017:
			switch (cluster)
			{
				case 1:
					geomTel = new BGeomTel2017_1();
					break;
				case 2:
					geomTel = new BGeomTel2017_2();
					break;
			}
			break;
		case 2018:
			switch (cluster)
			{
				case 1:
					geomTel = new BGeomTel2018_1();
					break;
				case 2:
					geomTel = new BGeomTel2018_2();
					break;
				case 3:
					geomTel = new BGeomTel2018_3();
					break;
			}
			break;
		case 2019:
			switch (cluster)
			{
				case 1:
					geomTel = new BGeomTel2019_1();
					break;
				case 2:
					geomTel = new BGeomTel2019_2();
					break;
				case 3:
					geomTel = new BGeomTel2019_3();
					break;
				case 4:
					geomTel = new BGeomTel2019_4();
					break;
				case 5:
					geomTel = new BGeomTel2019_5();
					break;
			}
			break;
		case 2020:
			switch (cluster)
			{
				case 1:
					geomTel = new BGeomTel2020_1();
					break;
				case 2:
					geomTel = new BGeomTel2020_2();
					break;
				case 3:
					geomTel = new BGeomTel2020_3();
					break;
				case 4:
					geomTel = new BGeomTel2020_4();
					break;
				case 5:
					geomTel = new BGeomTel2020_5();
					break;
				case 6:
					geomTel = new BGeomTel2020_6();
					break;
				case 7:
					geomTel = new BGeomTel2020_7();
					break;
			}
			break;
	}
	return geomTel;
}

int WaveformVis(TString dataPath, int year, int cluster, int run, int OMID, int pulseID)
{
	BGeomTel* geomConfig = GetGeomTel(year,cluster);

	// BGeomTel2019_2* geomConfig = new BGeomTel2019_2();
	int fileSdC = geomConfig->GetSdc(OMID/12);
	int sectionOMID = OMID%12;

	cout << pulseID << "\t" << OMID/12 << "\t" << OMID << "\t" << fileSdC << "\t" << gSecEntryIDs[OMID/12] << endl;

	TString yearID = GetSeasonID(year);
	if (yearID == "0")
		return -1;

	TString inputFile = Form("%s/%s%04d.raw.events.%d.root",dataPath.Data(),yearID.Data(),run,fileSdC);
	// TString inputFile = Form("/Data/BaikalData/exp16_barsv051/cluster0/0070/f0070.raw.events.%d.sorted.root",fileSdC);

	TFile* file = new TFile(inputFile.Data(),"READ");
	TTree* tree = (TTree*)file->Get("Events");

	if (tree->GetEntries() == 0)
	{
		std::cout << "Files: " << inputFile << " were not found!" << endl;
    	return -1;
	}

	BRawMasterData* rawMasterData = NULL;
    tree->SetBranchAddress("BRawMasterData.", &rawMasterData);
    BRawMasterHeader* rawMasterHeader = NULL;
    tree->SetBranchAddress("BRawMasterHeader.", &rawMasterHeader);

    // TCanvas* myCan = new TCanvas("myCan","Results",800,600);
    // TMultiGraph* mg_waveforms = new TMultiGraph(Form("OM_%d",OMID),"Waveforms;Time [ns];Amplitude [FADC channels]");
    // TLegend* myLegend = new TLegend(0.6,0.7,0.9,0.9);

    int eventID = gSecEntryIDs[OMID/12];

	tree->GetEntry(eventID);

	// cout << rawMasterData->GetNumSamples() << endl;
	TGraph* graphs[rawMasterData->GetNumSamples()];

	for (int i = 0; i < rawMasterData->GetNumSamples(); ++i)
	{
		BRawFADCSample* sample = rawMasterData->GetFADCSample(i);

		if (sample->GetNch() != sectionOMID)
			continue;

		graphs[i] = new TGraph();
	    graphs[i]->SetName(Form("graphs_%d",i));
	    graphs[i]->SetTitle(Form("OMID_%d, Q = %.1f, T = %.1f",OMID,hits[pulseID].charge,hits[pulseID].time));
		// graphs[i]->SetLineColor(i+1);
		graphs[i]->SetLineColor(kRed);
		// graphs[i]->SetMarkerColor(i+1);
		graphs[i]->SetMarkerColor(kRed);
		// graphs[i]->SetMarkerStyle(i+21);
		graphs[i]->SetMarkerStyle(21);
		graphs[i]->SetMarkerSize(0.9);

		// cout << i << " Nbins: " << sample->GetNbins() << " " << sample->GetNch()<< " " << sample->GetOffset() << endl;

		Int_t nbins = sample->GetNbins();
		Short_t *data = sample->GetData();
		Int_t offset = sample->GetOffset();
		double pedestal = 0;
		int amplitude = 0;
		int charge = 0;
		int amplitudeBin = 0;
		for (int k = 0; k < c_nPed; ++k)
		{
			pedestal += data[k];
		}
		pedestal /= c_nPed;
		for (int k = 0; k < nbins; ++k)
		{
			if (data[k]-pedestal > amplitude)
			{
				amplitude = data[k]-pedestal;
				amplitudeBin = k;
			}
			charge += data[k]-pedestal;
		}

		bool saveGraph = false;

		graphs[i]->Set(nbins);
		for(Int_t n = 0; n < nbins; n++) {
			if (TMath::Abs(5*(n+offset-512)-gOMtimeCal[OMID] - hits[pulseID].time) < 200)
				saveGraph = true;
			graphs[i]->SetPoint(n, 5*(n+offset-512)-gOMtimeCal[OMID], data[n]-pedestal);
		}

		if (saveGraph)
			hits[pulseID].waveform->Add(graphs[i]);
		// myLegend->AddEntry(graphs[i],Form("graphs_%d",i),"l");
	}

	// mg_waveforms->Draw("APL");
	// myLegend->Draw();

	// graph->Draw();
	// graph->GetXaxis()->SetTitle("Time [FADC channels]");
	// graph->GetYaxis()->SetTitle("Amplitude [FADC channels]");
	// graph->SetTitle("Comparison of different fitting functions");
	// myCan->Update();

	return 0;
}

int GumbelVis(int pulseID, bool showTrackHits)
{
	// cout << hits[0].OMID << " " << hits[0].charge << " " << hits[0].time << endl;
	TF1* gump = new TF1("gump","[0]*exp(-((x-[1]-9.85)/5/[2] + exp(-(x-[1]-9.85)/5/[2])))",hits[pulseID].expTime-50,hits[pulseID].expTime+150);

	gump->SetParameters(hits[pulseID].expCharge*gOMchargeCal[hits[pulseID].OMID]/6/0.368,hits[pulseID].expTime,2.0);
	TGraph* myGraph = new TGraph(gump);
	myGraph->SetTitle(Form("Gumbel, Q = %.1f, T = %.1f",hits[pulseID].expCharge,hits[pulseID].expTime));
	hits[pulseID].waveform->Add(myGraph);

	if (showTrackHits)
	{
		TF1* gumpTrack = new TF1("gumpTrack","[0]*exp(-((x-[1]-9.85)/5/[2] + exp(-(x-[1]-9.85)/5/[2])))",hits[pulseID].expTimeTrack-50,hits[pulseID].expTimeTrack+150);
		gumpTrack->SetParameters(100*gOMchargeCal[hits[pulseID].OMID]/7/0.368*TMath::Exp(-hits[pulseID].expDistTrack/23),hits[pulseID].expTimeTrack,2.0);
		TGraph* myGraphTrack = new TGraph(gumpTrack);
		myGraphTrack->SetTitle(Form("GumbelTrack, Q = %.1f, T = %.1f, D = %.1f",100*TMath::Exp(-hits[pulseID].expDistTrack/23),hits[pulseID].expTimeTrack,hits[pulseID].expDistTrack));
		myGraphTrack->SetLineColor(kBlue);
		myGraphTrack->SetMarkerColor(kBlue);
		hits[pulseID].waveform->Add(myGraphTrack);
	}

	return 0;
}

int ReadCalibFile(TString dataPath)
{
	TString fileName = Form("%s/calib.txt",dataPath.Data());

	ifstream inputCalibFile;
    inputCalibFile.open(fileName);

    if (!inputCalibFile)
    {
    	cerr << "Calib file: " << fileName << " was NOT found. Program termination!" << endl;
    	return -1;
  	}

  	int OMID = 0;
  	double x,y,z = 0;

  	for(int i = 0; i < gNOMs; i++)
  	{
  		if (inputCalibFile.eof())
  		{
  			cerr << "Calib file: " << fileName << " do not contain information for all 288 OMs. Program termination!" << endl;
  			return -1;
  		}
  		inputCalibFile >> OMID >> gOMchargeCal[OMID] >> gOMtimeCal[OMID] >> x >> y >> z;
  		gOMpositions[OMID].SetXYZ(x,y,z);
	  	// cout << OMID << "\t" << gOMchargeCal[OMID] << "\t" << gOMtimeCal[OMID] << endl;
	  	OMID++;
  	}

	inputCalibFile.close();
	return 0;
}

int ReadTimeResFile(TString dataPath)
{
	TString fileName = Form("%s/timeRes.txt",dataPath.Data());

	ifstream inputTimeResFile;
    inputTimeResFile.open(fileName);

    if (!inputTimeResFile)
    {
    	cerr << "TimeRes file: " << fileName << " was NOT found. Program termination!" << endl;
    	return -1;
  	}

  	int eventID,OMID;
  	double time, expectedTime, expectedTimeTrack, expectedDistanceTrack, charge, expectedCharge;

  	while(!inputTimeResFile.eof())
  	{
  		inputTimeResFile >> eventID >> OMID >> time >> expectedTime >> expectedTimeTrack >> expectedDistanceTrack >> charge >> expectedCharge;
  		if (inputTimeResFile.eof())
  			break;
  		hits.push_back(Hit{OMID,time,expectedTime,expectedTimeTrack,expectedDistanceTrack,charge,expectedCharge});
  		hits.back().waveform = new TMultiGraph(Form("OM_%d",OMID),"Waveforms;Time [ns];Amplitude [FADC channels]");
  	}

  	if (hits.size() == 0)
  	{
  		cerr << "TimeRes file: " << fileName << " do NOT contain any hits. Program termination!" << endl;
    	return -1;
  	}

  	cout << "Nhits = " << hits.size() << " has been read and is going to be processed!" << endl;

	inputTimeResFile.close();
	return 0;
}

int ExtractEventIDs(TString dataPath, int year, int runID, int eventID)
{
	TString yearID = GetSeasonID(year);
	if (yearID == "0")
		return -1;
	TString jointTableFilePath = Form("%s/%s%04d.joint_table.root",dataPath.Data(),yearID.Data(),runID);
	cout << jointTableFilePath << endl;

	TFile* file = new TFile(jointTableFilePath,"READ");
	TTree* tree = (TTree*)file->Get("Marks");

	if (tree->GetEntries() == 0)
	{
		std::cout << "Files: " << jointTableFilePath << " were not found!" << endl;
    	return -1;
	}

	BJointMarkCC* jointMarkCC = NULL;
    tree->SetBranchAddress("BJointMarkCC", &jointMarkCC);

    tree->GetEntry(eventID);

    // cout << jointMarkCC->GetMastersNum() << endl;

    if (jointMarkCC->GetMastersNum() < gNSecs)
    {
    	cout << "WARNING: Joint table: " << jointTableFilePath << " do NOT contain eventID info for all sections but only: " << jointMarkCC->GetMastersNum() << " entries !" << endl;
    	// return -1;
    }

    for (int i = 0; i < jointMarkCC->GetMastersNum(); ++i)
    {
    	if (jointMarkCC->GetIndexes(i) != eventID)
    	{
    		cout << "WARNING: For section " << jointMarkCC->GetMasIndexes(i) - 1 << " master index do not equal to eventID!" << endl;
    	}
    	gSecEntryIDs[jointMarkCC->GetMasIndexes(i) - 1] = jointMarkCC->GetIndexes(i);
    	// cout << jointMarkCC->GetIndexes(i) << "\t" << jointMarkCC->GetMasIndexes(i) << endl;
    }

	return 0;
}

int CreateWaveforms(TString dataPath,int year, int cluster, int run)
{
	// cout << "Waveform Vis" << endl;
	for (int i = 0; i < hits.size(); ++i)
	{
		if (WaveformVis(dataPath,year,cluster,run,hits[i].OMID,i) != 0)
			return -1;
	}
	return 0;
}

int CreatePulses(bool showTrackHits)
{
	// cout << "Gumbel Vis" << endl;
	for (int i = 0; i < hits.size(); ++i)
	{
		GumbelVis(i,showTrackHits);
	}
	return 0;
}

int SaveWaveforms(TString outputFilePath)
{
	TFile* outputFile = new TFile(Form("%s/waveforms.root",outputFilePath.Data()),"RECREATE");

	if (!outputFile->IsOpen())
	{
		std::cout << "Output File: " << outputFilePath << " can not be recreated!" << endl;
    	return -1;
	}
	// TLegend* myLegend = new TLegend(0.6,0.7,0.9,0.9);

	for (int i = 0; i < hits.size(); ++i)
	{
		TCanvas* myCan = new TCanvas(Form("c_OMID_%d",hits[i].OMID),Form("OMID_%d",hits[i].OMID),800,600);
		hits[i].waveform->Draw("APL");

		myCan->BuildLegend();
		myCan->Write();
		myCan->SaveAs(Form("%s/OMID_%d.png",outputFilePath.Data(),hits[i].OMID));
		// myLegend->Draw();
		// hits[i].waveform->Write();
	}
	return 0;
}

int waveformComparison(int year, int cluster, int run, int event, TString dataPath, bool showTrackHits = false)
{
	TString outputFilePath = Form("%s/waveforms/",dataPath.Data());

	if (ReadCalibFile(dataPath) != 0)
		return -1;

	if (ReadTimeResFile(dataPath) != 0)
		return -1;

	if (ExtractEventIDs(dataPath,year,run,event) != 0)
		return -1;

	if (CreateWaveforms(dataPath,year,cluster,run) != 0)
		return -1;

	CreatePulses(showTrackHits);

	SaveWaveforms(outputFilePath);

	return 0;
}
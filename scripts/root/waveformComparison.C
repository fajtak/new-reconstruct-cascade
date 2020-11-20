const int c_nPed = 5;


int waveformVis(int sdc, int eventID, int sampleID, int OMIDshift = 0)
{
	// TString inputFile = Form("/Data/BaikalData/exp16_barsv051/cluster0/0070/f0070.raw.events.%d.sorted.root",sdc);
	TString inputFile = Form("/Data/BaikalData/exp19_barsv051/cluster1/0112/i0112.raw.events.%d.root",sdc);


	TFile* file = new TFile(inputFile,"READ");

	TTree* tree = (TTree*)file->Get("Events");
	BRawMasterData* rawMasterData = NULL;
    tree->SetBranchAddress("BRawMasterData.", &rawMasterData);
    BRawMasterHeader* rawMasterHeader = NULL;
    tree->SetBranchAddress("BRawMasterHeader.", &rawMasterHeader);

    TCanvas* myCan = new TCanvas("myCan","Results",800,600);

    TGraph* graph = new TGraph();
    graph->SetName("graph");
	graph->SetLineColor(4);
	graph->SetMarkerColor(4);
	graph->SetMarkerStyle(21);
	graph->SetMarkerSize(0.9);

	tree->GetEntry(eventID);

	cout << rawMasterData->GetNumSamples() << endl;

	BRawFADCSample* sample = rawMasterData->GetFADCSample(sampleID);

	cout << "Nbins: " << sample->GetNbins() << " " << sample->GetNch()<< " " << sample->GetOffset() << endl;

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

	graph->Set(nbins);
	for(Int_t n = 0; n < nbins; n++) {
		graph->SetPoint(n, n+offset, data[n]-pedestal);
	}

	graph->Draw();
	graph->GetXaxis()->SetTitle("Time [FADC channels]");
	graph->GetYaxis()->SetTitle("Amplitude [FADC channels]");
	graph->SetTitle(Form("OMID_%d",sample->GetNch()+OMIDshift));
	myCan->Update();

	return 0;
}

int waveformVis(int sdc, int eventID)
{
	// TString inputFile = Form("/Data/BaikalData/exp16_barsv051/cluster0/0070/f0070.raw.events.%d.sorted.root",sdc);
	TString inputFile = Form("/Data/BaikalData/exp19_barsv051/cluster1/0112/i0112.raw.events.%d.root",sdc);

	TFile* file = new TFile(inputFile,"READ");

	TTree* tree = (TTree*)file->Get("Events");
	BRawMasterData* rawMasterData = NULL;
    tree->SetBranchAddress("BRawMasterData.", &rawMasterData);
    BRawMasterHeader* rawMasterHeader = NULL;
    tree->SetBranchAddress("BRawMasterHeader.", &rawMasterHeader);

    TCanvas* myCan = new TCanvas("myCan","Results",800,600);
    TMultiGraph* mg_waveforms = new TMultiGraph("mg_waveforms","Waveforms;Time [FADC channels];Amplitude [FADC channels]");
    TLegend* myLegend = new TLegend(0.6,0.7,0.9,0.9);

	tree->GetEntry(eventID);

	cout << rawMasterData->GetNumSamples() << endl;
	TGraph* graphs[rawMasterData->GetNumSamples()];

	for (int i = 0; i < rawMasterData->GetNumSamples(); ++i)
	{
		BRawFADCSample* sample = rawMasterData->GetFADCSample(i);
		graphs[i] = new TGraph();
	    graphs[i]->SetName(Form("graphs_%d",i));
		graphs[i]->SetLineColor(i+1);
		graphs[i]->SetMarkerColor(i+1);
		graphs[i]->SetMarkerStyle(i+21);
		graphs[i]->SetMarkerSize(0.9);

		cout << i << " Nbins: " << sample->GetNbins() << " " << sample->GetNch()<< " " << sample->GetOffset() << endl;

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

		graphs[i]->Set(nbins);
		for(Int_t n = 0; n < nbins; n++) {
			graphs[i]->SetPoint(n, n+offset, data[n]-pedestal);
		}

		mg_waveforms->Add(graphs[i]);
		myLegend->AddEntry(graphs[i],Form("graphs_%d",i),"l");
	}

	mg_waveforms->Draw("APL");
	myLegend->Draw();

	// graph->Draw();
	// graph->GetXaxis()->SetTitle("Time [FADC channels]");
	// graph->GetYaxis()->SetTitle("Amplitude [FADC channels]");
	// graph->SetTitle("Comparison of different fitting functions");
	// myCan->Update();

	return 0;
}

const int gNOMs = 288;
double gOMtimeCal[gNOMs] = {0};
double gOMchargeCal[gNOMs] = {0};
TVector3 gOMpositions[gNOMs];

struct Hit
{
	int OMID;
	double time;
	double expTime;
	double charge;
	double expCharge;
	TMultiGraph* waveform;
};

vector<Hit> hits;

TString calibFilePath = "/Data/BaikalData/exp19_barsv051/cluster1/0112/calib.txt";
TString timeResFilePath = "/Data/BaikalData/exp19_barsv051/cluster1/0112/timeRes.txt";
TString outputFilePath = "/Data/BaikalData/exp19_barsv051/cluster1/0112/waveforms/waveforms.root";
TString pngOutputFilePath = "/Data/BaikalData/exp19_barsv051/cluster1/0112/waveforms/";

TMultiGraph* WaveformVis(int OMID, int pulseID)
{
	BGeomTel2019_2* geomConfig = new BGeomTel2019_2();
	int fileSdC = geomConfig->GetSdc(OMID/12);
	int sectionOMID = OMID%12;
	TString inputFile = Form("/Data/BaikalData/exp19_barsv051/cluster1/0112/i0112.raw.events.%d.root",fileSdC);
	// TString inputFile = Form("/Data/BaikalData/exp16_barsv051/cluster0/0070/f0070.raw.events.%d.sorted.root",fileSdC);

	TFile* file = new TFile(inputFile,"READ");

	TTree* tree = (TTree*)file->Get("Events");
	BRawMasterData* rawMasterData = NULL;
    tree->SetBranchAddress("BRawMasterData.", &rawMasterData);
    BRawMasterHeader* rawMasterHeader = NULL;
    tree->SetBranchAddress("BRawMasterHeader.", &rawMasterHeader);

    // TCanvas* myCan = new TCanvas("myCan","Results",800,600);
    TMultiGraph* mg_waveforms = new TMultiGraph(Form("OM_%d",OMID),"Waveforms;Time [FADC channels];Amplitude [FADC channels]");
    // TLegend* myLegend = new TLegend(0.6,0.7,0.9,0.9);

	tree->GetEntry(364888);

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
			if (TMath::Abs(5*(n+offset-512)-gOMtimeCal[OMID] - hits[pulseID].time) < 100)
				saveGraph = true;
			graphs[i]->SetPoint(n, 5*(n+offset-512)-gOMtimeCal[OMID], data[n]-pedestal);
		}

		if (saveGraph)
			mg_waveforms->Add(graphs[i]);
		// myLegend->AddEntry(graphs[i],Form("graphs_%d",i),"l");
	}

	// mg_waveforms->Draw("APL");
	// myLegend->Draw();

	// graph->Draw();
	// graph->GetXaxis()->SetTitle("Time [FADC channels]");
	// graph->GetYaxis()->SetTitle("Amplitude [FADC channels]");
	// graph->SetTitle("Comparison of different fitting functions");
	// myCan->Update();

	return mg_waveforms;
}

int GumbelVis(int pulseID)
{
	// cout << hits[0].OMID << " " << hits[0].charge << " " << hits[0].time << endl;
	TF1* gump = new TF1("gump","[0]*exp(-((x-[1]-9.85)/5/[2] + exp(-(x-[1]-9.85)/5/[2])))",hits[pulseID].expTime-50,hits[pulseID].expTime+150);

	gump->SetParameters(hits[pulseID].expCharge*gOMchargeCal[hits[pulseID].OMID]/7/0.368,hits[pulseID].expTime,2.0);
	TGraph* myGraph = new TGraph(gump);
	myGraph->SetTitle(Form("Gumbel, Q = %.1f, T = %.1f",hits[pulseID].expCharge,hits[pulseID].expTime));
	hits[pulseID].waveform->Add(myGraph);
	return 0;
}

int ReadCalibFile(TString fileName)
{
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
  		inputCalibFile >> OMID >> gOMchargeCal[OMID] >> gOMtimeCal[OMID] >> x >> y >> z;
  		gOMpositions[OMID].SetXYZ(x,y,z);
	  	// cout << OMID << "\t" << gOMchargeCal[OMID] << "\t" << gOMtimeCal[OMID] << endl;
	  	OMID++;
  	}

	inputCalibFile.close();
	return 0;
}

int ReadTimeResFile(TString fileName)
{
	ifstream inputTimeResFile;
    inputTimeResFile.open(fileName);

    if (!inputTimeResFile)
    {
    	cerr << "TimeRes file: " << fileName << " was NOT found. Program termination!" << endl;
    	return -1;
  	}

  	int eventID,OMID;
  	double timeRes, expectedTime, charge, expectedCharge;

  	while(!inputTimeResFile.eof())
  	{
  		inputTimeResFile >> eventID >> OMID >> timeRes >> expectedTime >> charge >> expectedCharge;
  		if (inputTimeResFile.eof())
  			break;
  		hits.push_back(Hit{OMID,expectedTime+timeRes,expectedTime,charge,expectedCharge});
  	}

	inputTimeResFile.close();

	return 0;
}

int CreateWaveforms()
{
	// cout << "Waveform Vis" << endl;
	for (int i = 0; i < hits.size(); ++i)
	{
		hits[i].waveform = WaveformVis(hits[i].OMID,i);
	}
	return 0;
}

int CreatePulses()
{
	// cout << "Gumbel Vis" << endl;
	for (int i = 0; i < hits.size(); ++i)
	{
		// cout << i << endl;
		GumbelVis(i);
	}
	return 0;
}

int SaveWaveforms()
{
	TFile* outputFile = new TFile(outputFilePath,"RECREATE");

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
		myCan->SaveAs(Form("%s/OMID_%d.png",pngOutputFilePath.Data(),hits[i].OMID));
		// myLegend->Draw();
		// hits[i].waveform->Write();
	}
	return 0;
}

int waveformComparison()
{
	ReadCalibFile(calibFilePath);
	ReadTimeResFile(timeResFilePath);

	CreateWaveforms();
	CreatePulses();

	SaveWaveforms();

	return 0;
}
const int c_nPed = 5;


int waveformVis(int sdc, int eventID, int sampleID, int OMIDshift = 0)
{
	TString inputFile = Form("/Data/BaikalData/exp16_barsv051/cluster0/0070/f0070.raw.events.%d.sorted.root",sdc);

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

	cout << sample->GetNbins() << endl;

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
	TString inputFile = Form("/Data/BaikalData/exp16_barsv051/cluster0/0070/f0070.raw.events.%d.sorted.root",sdc);

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
	

	// cout << rawMasterData->GetNumSamples() << endl;

	// BRawFADCSample* sample = rawMasterData->GetFADCSample(0);

	// cout << sample->GetNbins() << endl;

	// Int_t nbins = sample->GetNbins();
	// Short_t *data = sample->GetData();
	// double pedestal = 0;
	// int amplitude = 0;
	// int charge = 0;
	// int amplitudeBin = 0;
	// for (int k = 0; k < c_nPed; ++k)
	// {
	// 	pedestal += data[k];
	// }
	// pedestal /= c_nPed;
	// for (int k = 0; k < nbins; ++k)
	// {
	// 	if (data[k]-pedestal > amplitude)
	// 	{
	// 		amplitude = data[k]-pedestal;
	// 		amplitudeBin = k;
	// 	}
	// 	charge += data[k]-pedestal;
	// }
	// // if (amplitude-pedestal > 70	|| amplitude-pedestal < 65)
	// if (charge < 0)
	// 	cout << "NEGATIVE CHARGE!" << endl;

	// graph->Set(nbins);
	// for(Int_t n = 0; n < nbins; n++) {
	// 	graph->SetPoint(n, n, data[n]-pedestal);
	// }	
	// fitFunc->SetParameters(0,amplitude/0.37,amplitudeBin,2.2);
	// graph->Fit(fitFunc,"W");
	// graph->Fit(fitFunc2,"W+");
	// graph->Fit(fitFunc3,"W+");
	// // cout << "amplitude = " << amplitude << " charge = " << charge << endl;
	// graph->Draw();
	// graph->GetXaxis()->SetTitle("Time [FADC channels]");	
	// graph->GetYaxis()->SetTitle("Amplitude [FADC channels]");	
	// graph->SetTitle("Comparison of different fitting functions");
	// myCan->Update();

	// TLegend* myLegend = new TLegend(0.6,0.7,0.9,0.9);
	// myLegend->AddEntry(fitFunc,"Gumpbel function, #chi^{2}/NDF = 0.32","l");
	// myLegend->AddEntry(fitFunc2,"Gaussian function, #chi^{2}/NDF = 2.86","l");
	// myLegend->AddEntry(fitFunc3,"Landau function, #chi^{2}/NDF = 2.19","l");
	// myLegend->Draw();

 //    return 0;
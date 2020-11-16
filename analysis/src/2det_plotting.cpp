#include "../inc/2det_plotting.hpp"

/*----------------------------------------------------------------------
 * This file contains a bunch of functions for plotting multi-detector
 * events (as well as some single-detector events too
 *--------------------------------------------------------------------*/

// Conversion of reducedEvent vector (i.e. the sorted events from our detector) 
// into ROOT histograms.
TH2D* plot2DTimeEnergy(std::vector<reducedEvent> results, double minT, double maxT, double minE, double maxE) {
	// 2 dimensional fitting of histogram with energy and hit time.
	// Can use this to later probe slice results.
	
	TH2D* hist2DRoot = new TH2D("detHist", "hit time vs. energy", 30,minT,maxT, 30, minE, minE+maxE);
	for(size_t ii = 0; ii < results.size(); ii++) { // loop through results and fill histogram
		hist2DRoot->Fill(results.at(ii).time, results.at(ii).energy, results.at(ii).weight);
	}
	return hist2DRoot;
}

TH1D* plotTimeHist(std::vector<reducedEvent> results, double minT, double maxT) {
	// Convert a results vector into a time-projection ROOT hist.
	// Probably another function that does this.
	
	TH1D* timeHist = new TH1D("timeHist", "hit time", (int)(maxT-minT), minT, maxT);
	for(size_t ii = 0; ii < results.size(); ii++) { // loop through results and fill histogram
		timeHist->Fill(results.at(ii).time, results.at(ii).weight);
	}
	return timeHist;
}
	
TH1D* plotEneHist(std::vector<reducedEvent> results, double minE, double maxE) {
	// Convert a results vector into an energy-projection ROOT hist.
	
	TH1D* eneHist = new TH1D("eneHist", "energy", 300, minE, maxE);
	for(size_t ii = 0; ii < results.size(); ii++) { // loop through results and fill histogram
		eneHist->Fill(results.at(ii).energy,results.at(ii).weight);
	}
	return eneHist;
}

TH1D* plotThetaHist(std::vector<reducedEvent> results, double minTh, double maxTh) {
	
	TH1D* thetaHist = new TH1D("thetaHist", "theta", 260, minTh, maxTh);
	for (size_t ii = 0; ii < results.size(); ii++) { // loop through results and fill histogram
		thetaHist->Fill(results.at(ii).theta, results.at(ii).weight);
	}
	return thetaHist;
}


TCanvas* plotMapTot(std::vector<reducedEvent> results, double minT, double maxT, double minE, double maxE) {
	// This just calls the three previous functions and sets them onto a single canvas
	// I'm going to ultimately want to copy this structure for later things
		 
	// Generation of histograms on canvas
	TCanvas* c1 = new TCanvas("c1","c1",750,500);
	c1->SetFillColor(0); // set color to white
	
	TPad* pad1 = new TPad("p1","Energy vs Time",0.03,0.03,0.645,0.97,21); // make pad
	pad1->Draw();
	pad1->cd();
	
	TH2D* hist2DRoot = plot2DTimeEnergy(results, minT, maxT, minE, maxE); // Load histogram and draw it
	hist2DRoot->Draw("CONT4Z");
	
	TPad* pad2 = new TPad("p2","Arrival Time",0.675,0.03,0.97,0.485,21); // make second pad
	pad2->Draw();
	pad2->cd();
	
	TH1D* timeHist = plotTimeHist(results, minT, maxT); // Load time histogram
	timeHist->Scale(1/timeHist->Integral());
	timeHist->Draw();
	
	TPad* pad3 = new TPad("p3","Neutron Energy",0.675,0.515,0.97,0.97,21);
	pad3->Draw();
	pad3->cd();
	
	TH1D* eneHist = plotEneHist(results, minE, maxE);
	eneHist->Scale(1/eneHist->Integral());
	eneHist->Draw();

	c1->Update(); // Update canvas for plotting
	
	return c1;
}
	
//TCanvas* plotBlockMapSlice(std::vector<double> wts, std::vector<fixedResultBlock> block, const size_t num, double thr, double time) {
TCanvas* plotMapSlice(std::vector<reducedEvent> results, int nSlices, double minE, double maxE) {
	// This makes a canvas by slicing through and compressing energy levels
	// Slice results data into different energy slices
	
	TCanvas* c2 = new TCanvas("c2","c2",750,1000);
	c2->cd();
	
	TPad* padc2[nSlices];
	TH1D* sliceHist[nSlices];
	double stepsize = (maxE - minE) / ((double)nSlices);
	for (int i = 0; i < nSlices; i++) { // loop through names
		char* tPadName = new char[24]; // initialize names -- bigger than it needs to be to avoid compiler issues
		sprintf(tPadName,"p%d",i);
		char* tPadTitle = new char[48];
		sprintf(tPadTitle,"Arrival Time %f < E < %f m", (double)(i)*stepsize + minE, (double)(i+1)*stepsize + minE);
		
		double padMin = (i+0.0)/nSlices + 0.01; // Initialize TPad		
		double padMax = ((i+1.0)/nSlices) - 0.01;
		
		padc2[i] = new TPad(tPadName,tPadTitle,0.03,padMin,0.97,padMax,21);
		padc2[i]->SetLogy();
		padc2[i]->Draw();
		padc2[i]->cd();
		
		sliceHist[i] = plotEneHist(results, (double)(i)*stepsize + minE, (double)(i+1)*stepsize + minE); // Load slice histogram
		sliceHist[i]->Draw();
	}
	c2->Update();
	return c2;
}

TH1D* dataMCComp(std::vector<double>& dataVector, std::vector<weightedBin>& mcVector) {
	// Compare the weighted bins of a monte carlo simulation with the raw data
	// This is kinda obsolete since we've got chi2 algorithms already
	// This will generate the residual
	
	TH1D* resHist; // Initialization only works with the same sized vectors
	if(dataVector.size() != mcVector.size()) {
        return resHist;
    }
	resHist = new TH1D("residual", "MC - Data Residual", dataVector.size(), 0, dataVector.size());
    
    // Calculate sums for normalizaton
    double sumMC  = std::accumulate(mcVector.begin(), mcVector.end(), 0, [](double sum, weightedBin bin)->double{return sum + bin.wgt;});
    double sumDat = std::accumulate(dataVector.begin(), dataVector.end(), 0);
    
    for (int i = 0; i < mcVector.size(); i++) {
		resHist->SetBinContent(i,mcVector.at(i).wgt / sumMC - dataVector.at(i) / sumDat);
	}
	return resHist;

}
/*	
	TCanvas* c3 = new TCanvas("c3","c3",750,500);
	TPad* pad1c3 = new TPad("p1","Arrival Time (Dagger)",0.03,0.03,0.97,0.97,21);
	pad1c3->Draw();
	pad1c3->cd();
	timeDagHist->Scale(1/timeDagHist->Integral());
	timeDagHist->Draw();
	if (chiFit) {
		if (dataHist->GetSize() == timeDagHist->GetSize()){
			dataHist->Scale(1/dataHist->Integral(0.0,140.0));
			dataHist->SetLineColor(kRed);
			dataHist->Draw("SAME");
		}
	}	
	if (saveMCFit) {
		timeDagHist->SaveAs("./timeDagHist.root");
		printf("\n !!!SAVING DATA AS timeDagHist.root!!! \n");
	}
}*/
std::vector<measurement> calculateDrainTimes(std::vector<double> timeVector,int nBins) {
	// Calculate draining times of dagger from dips
	
	std::vector<measurement> drainTimes;
	std::vector<double> dipEnds = DIPTIMES; // Load dipTimes
	if (dipEnds.size() < 2) { // Don't calculate if there's no dipTimes we can run
		return drainTimes;
	}
	
	TH1D* timeDagHist = new TH1D("timeDag", "hit time (dagger)", nBins, (double)(*dipEnds.begin()), (double)(*dipEnds.end()));
	for (auto tIt = timeVector.begin(); tIt < timeVector.end(); tIt++) { // fill ROOT histogram to use their algorithms
		timeDagHist->Fill(*tIt);
	}
	
	timeDagHist->Scale(1/timeDagHist->Integral()); // Normalize histogram
		
	TF1* expoInv = new TF1("expoInv","[0]*exp(-x/[1])"); // Set parameters
	expoInv->SetParameters(1.0,880.0); // Initialize parameters
	expoInv->SetParLimits(0,0.0,1.0); // Limit parameters to be realistic
	expoInv->SetParLimits(1,0.0,10000.0);
	
	for (size_t i = 0; i < dipEnds.size() - 2; i++) {
		double moveT = 5.0; // time to move dagger
		timeDagHist->Fit("expoInv","Q","",dipEnds[i]+moveT,dipEnds[i+1]-1.0);
		drainTimes.push_back({timeDagHist->GetFunction("expoInv")->GetParameter(1),
							  timeDagHist->GetFunction("expoInv")->GetParError(1)}); // Push back measurement
	}
	return drainTimes;
}


/*
void plotBFieldSlice(std::vector<double> wts, std::vector<fixedResultBlock> block, const size_t num, double thr, double maxF) {
		
	if (bFieldSlice) {
		TCanvas* c4 = new TCanvas("c4","c4",750,1000);
		TPad* padc4[numHistos];
		
		sprintf(fPadName,"p%d",0);
		sprintf(fPadTitle,"B Field E < %f m", threshold);
		padMin = 0.01;
		padMax = 0.99;
		//padMax = (1.0/numHistos) - 0.01;
		padc4[0] = new TPad(fPadName,fPadTitle,0.03,padMin,0.97,padMax,21);
		padc4[0]->Draw();
		//for (int i=1;i<numHistos;i++) { 
		//	sprintf(fPadName,"p%d",i);
		//	sprintf(fPadTitle,"B Field %f < E < %f m", (i-1)*spacing + threshold, i*spacing + threshold);
		//	padMax = ((i+1.0)/numHistos) - 0.01;
		//	padMin = (i+0.0)/numHistos + 0.01;
		//	padc4[i] = new TPad(fPadName,fPadTitle,0.03,padMin,0.97,padMax,21);
		//	padc4[i]->Draw();
		//}
		double normB = fieldBHist[3]->Integral();
		for (int i=0; i<numHistos;i++) {
			padc4[0]->cd();
			fieldBHist[i]->Scale(1/fieldBHist[i]->Integral());
			//fieldBHist[i]->Scale(1/normB);
			sprintf(fPadTitle,"B Field %f < E < %f m", (i-1)*spacing + threshold, i*spacing + threshold);
			fieldBHist[i]->SetTitle(fPadTitle);
			fieldBHist[i]->Draw("Same");		
		}
		c4->Update();
	}
	// check that wts is the same size as block:
	if (wts.size() < block.size()) {
		for (size_t jj = 0; jj < (block.size() - wts.size()); jj++) {
			wts.push_back(1.0);
		}
		for (size_t jj = 0; jj < wts.size(); jj++) { wts.at(jj) = 1.0; } // If we miscalculated set everything to wts=1
	}
	
	double dE = (EMAX - thr) / ((double)num - 1); 
	
	// Need to fill fieldBHist
	char* bName = new char[10];
	char* bTitle = new char[48];
	TH1D* fieldBHist[num];
	
	sprintf(bName,"fBlock0");
	sprintf(bTitle,"hit field E > %f m", thr);
	fieldBHist[0] =  new TH1D(fBName, fBTitle, 300,0,maxF);
	// Loop across the rest of the generated histograms 
	for (size_t ii=1; ii<num;ii++){
		sprintf(bName,"fBlock%d",i);
		sprintf(bTitle,"hit field %f < E < %f m", (ii-1)*dE + thr, ii*dE + thr);
		fieldBHist[i] =  new TH1D(bName, bTitle, 300,0,maxF);
	}
	
	for (size_t ii = 0; ii < block.size(); ii++) {
		for (size_t jj=1; jj < num; jj++) {
			double minSl = ((double)jj - 1.0)* dE + thr;
			double maxSl = ((double)jj)*dE + thr;
			if (((block.at(ii).energy)/(GRAV*MASS_N) > minSl) && ((block.at(ii).energy)/(GRAV*MASS_N) <= maxSl)) {
				fieldBHist[jj]->Fill(fieldstrength(block.at(ii).position),wts.at(ii));
			}
		}
	}
	
	TCanvas* c4 = new TCanvas("c4","c4",750,1000);
	TPad* padc4[num];
	
	char* tPadName = new char[4];
	char* tPadTitle = new char[48];
	sprintf(fPadName,"p%d",0);
	sprintf(fPadTitle,"B Field E < %f m", thr);
	
	double padMin = 0.01;
	double padMax = 0.99;
	//padMax = (1.0/numHistos) - 0.01;
	padc4[0] = new TPad(fPadName,fPadTitle,0.03,padMin,0.97,padMax,21);
	padc4[0]->Draw();
	for (size_t i=1; i<num; i++) { 
		sprintf(fPadName,"p%d",i);
		sprintf(fPadTitle,"B Field %f < E < %f m", (i-1)*spacing + threshold, i*spacing + threshold);
		padMax = ((i+1.0)/numHistos) - 0.01;
		padMin = (i+0.0)/numHistos + 0.01;
		padc4[i] = new TPad(fPadName,fPadTitle,0.03,padMin,0.97,padMax,21);
		padc4[i]->Draw();
	}
	
	for (size_t i=0; i<num; i++) {
		padc4[i]->cd();
		fieldBHist[i]->Scale(1/fieldBHist[i]->Integral());
		sprintf(fPadTitle,"B Field %f < E < %f m", (i-1)*spacing + threshold, i*spacing + threshold);
		fieldBHist[i]->SetTitle(fPadTitle);
		fieldBHist[i]->Draw("Same");		
	}
	c4->Update();
	return;
}
*/


TCanvas* plotChiMap(std::vector<std::vector<double>> variables, std::vector<double> chi2) {
	// This might not work or need to be done in python rather than ROOT
	// Plot the chi2 multi-value curves
	// Should be something like matplotlib's corner
	
	TCanvas* c5 = new TCanvas("c5","c5",750,750);
	int chi2Dim = variables.size(); // dimensions of chi2 scan
	double chi2Min = *std::min_element(chi2.begin(),chi2.end());
	//double chi2Min = chi2.at(mIn);
	
	return c5;
		
	// Find how many chi2 pads we require -- should be like n^2
	int nPads = 0;
	for (int ii = 0; ii < chi2Dim; ii++) {
		for (int jj = 0; jj < ii; jj++) { // More pads for each scan
			nPads += 1;
		}
	}
	// Now generate individual pads	
	TPad* padc5[nPads];
	int mapC = 0; // counter for pad numbers
	for (int ii = 0; ii < chi2Dim; ii++) { // Loop through and find bounds
		for (int jj = 0; jj < ii; jj++) {
			char* tPadName = new char[4]; // initialize names
			sprintf(tPadName,"p%d%d",ii,jj);
			char* tPadTitle = new char[48];
			sprintf(tPadTitle,"chi2 fit, %d vs %d", ii,jj);
				
			// Here we do the pad-sizing based on the number of pads
			double padMini = (((double)ii - 1.0)/(double)(chi2Dim-1)) + 0.01;
			if (padMini > 1.0) { padMini -= 1.0; } // Error catching
			double padMaxi = (((double)ii + 0.0)/(double)(chi2Dim-1)) - 0.01;
			if (padMaxi > 1.0) { padMaxi -= 1.0; } // Error catching
			
			// Same thing but in other dimension
			double padMinj = (((double)jj + 0.0)/(double)(chi2Dim-1)) + 0.01;
			if (padMinj > 1.0) { padMini -= 1.0; }			
			double padMaxj = (((double)jj + 1.0)/(double)(chi2Dim-1)) - 0.01;
			if (padMaxj > 1.0) { padMaxj -= 1.0; }
			
				
			//printf("%f, %f, %f, %f, %d\n", padMini, padMaxi, padMinj, padMaxj, mapC);
			// Generate pad	
			padc5[mapC] = new TPad(tPadName,tPadTitle,padMini,padMinj,padMaxi,padMaxj,21);
			padc5[mapC]->Draw();
			mapC += 1;// Increment counter
		}
	}
		
	// Load our data from file
	TGraph2D* chi2Map[nPads];
	mapC = 0;
	for (int ii = 0; ii < chi2Dim; ii++) { 
		for (int jj = 0; jj < ii; jj++) {
			// TODO: Should attempt to generate this from actual input vectors
			char* c2Filename = new char[24];		
			sprintf(c2Filename,"chi2Fit%d%d.csv", ii, jj); // load 2d object from .csv
			chi2Map[mapC] = new TGraph2D(c2Filename);
				chi2Map[mapC]->SetName(c2Filename);
				chi2Map[mapC]->SetTitle("Chi^2 Fit");
				
			int nBins = chi2Map[mapC]->GetHistogram()->GetSize(); // Rescale bins to only map in a given range
			for (int kk = 0; kk < nBins; kk++) {
				if (chi2Map[mapC]->GetHistogram()->GetBinContent(kk) > (chi2Min+3.0)){ // only map out 3 stdev
					chi2Map[mapC]->GetHistogram()->SetBinContent(kk,0.0);
				}
			}
			
			chi2Map[mapC]->SetMaximum(chi2Min+3.0);
			chi2Map[mapC]->SetMinimum(chi2Min);
			
			padc5[mapC]->cd(); // and draw the map
			chi2Map[mapC]->Draw("COLZ1"); // convert to contour
			
			mapC += 1;
		}
		
	}
	c5->Update();
	return c5;
}

TCanvas* plotSanityCheck(std::vector<reducedEvent> results) {
	// Straightforwards sanity check plotter -- we want to look at 
	// the energy and angular distribution of all detected UCN.
	//
	// This should let us confirm things are loaded properly.
		
	TCanvas* c6 = new TCanvas("c6","c6",750,750); // Load canvas
	
	TPad* pad1 = new TPad("p1","Energy (detected neutrons)",0.03,0.03,0.97,0.485,21); // Running 2 pads for this 
	pad1->Draw();
	pad1->cd();
	
	TH1D* eneHist = plotEneHist(results,0.0,50.0); // Hard coding full possible energy range
	eneHist->Scale(1/eneHist->Integral()); // Normalize histogram
	eneHist->Draw();
		
	TPad* pad2 = new TPad("p2","Theta (detected neutrons)",0.03,0.515,0.97,0.97,21); // And other pad
	pad2->Draw();
	pad2->cd();
	
	TH1D* thetaHist = plotThetaHist(results,0.0,2.7); // Hard coding full possible angle range
	thetaHist->Scale(1/thetaHist->Integral()); // normalize histogram
	thetaHist->Draw();
	
	c6->Update(); // Update canvas
	
	eneHist->SaveAs("./eneHist.root"); // Save histograms for later
	thetaHist->SaveAs("./thetaHist.root");
		
	return c6;
}

TCanvas* plotDetectorPos_2D(std::vector<double> xList, std::vector<double> yList, std::vector<double> zList) {
	// This function records the three-dimensional hit positions of neutrons on 
	// an MC detector
	
	TCanvas* c7 = new TCanvas("c7","c7",500,500); // generate canvas
	c7->SetFillColor(0);
	
	if ((xList.size() != yList.size()) || (xList.size() != zList.size())) {
		return c7;
	}
	TGraph2D* hitPos = new TGraph2D(xList.size()); // Generate a tgraph scatter
	hitPos->SetName("hitPosition");
	hitPos->SetTitle("Hit Position");
	
	for (size_t ii = 0; ii < xList.size(); ii++) { // fill
		hitPos->SetPoint(ii,xList.at(ii),yList.at(ii), zList.at(ii));
	}

	// And save canvases	
	TPad* pad1 = new TPad("p1","Detector Hit Position",0.03,0.03,0.97,0.97,21);	
	pad1->Draw();
	pad1->cd();
	
	hitPos->GetHistogram()->SetMinimum(-1.46); // hardcoded for block, have to change these for other detectors
	hitPos->GetHistogram()->SetMaximum(-1.42); // Histogram is really z_pos
	hitPos->GetXaxis()->SetRangeUser(0.20,0.25); // x and y locations
	hitPos->GetYaxis()->SetRangeUser(0.10,0.15);
	hitPos->Draw("PCOL"); // Scatter points
	
	c7->Update();	
	return c7;
}
	
TH1D* plotPercentDag(std::vector<weightedBin> energyD1, std::vector<weightedBin> energyD2) { 
	// Plot the percentage of neutrons at a given energy hitting detectors 1 and 2
	// WeightedBin here is just a proxy -- need to actually fit energy bands
	
	
	TCanvas* c8 = new TCanvas("c8","c8",750,750);
	double totalHits, percentage;
	TH1D* percentDag = new TH1D("percentDag", "% of neutrons hitting detector 2",300, EMIN, EMAX);
	
	if (energyD1.size() != energyD2.size()) {
		printf("Error! Detectors not same size!\n");
		return percentDag;
	}	
		
	for (int ii = 0; ii < energyD2.size(); ii++)	{
		totalHits = energyD1.at(ii).wgt+energyD2.at(ii).wgt;
		if(totalHits != 0.0){
			percentage = energyD1.at(ii).wgt / totalHits;
		} else{ 
			percentage = 1.0;
		}
		percentDag->SetBinContent(ii, 1.0-percentage);
	}
	percentDag->Draw();
		
	return percentDag;
	
}

	//timeBlockHist->Fit("expo2","WL M","",200.0,240.0+hold_t);
		/*timeBlockHist->Fit("expo2","QL","",200.0,240.0+hold_t);
		 
		blockPop[0]    = timeBlockHist->GetFunction("expo2")->GetParameter(0);
		blockPopErr[0] = timeBlockHist->GetFunction("expo2")->GetParError(0);
		blockSlope[0]  = timeBlockHist->GetFunction("expo2")->GetParameter(1);
		blockErr[0]    = timeBlockHist->GetFunction("expo2")->GetParError(1);
		blockPop[1]    = timeBlockHist->GetFunction("expo2")->GetParameter(2);
		blockPopErr[1] = timeBlockHist->GetFunction("expo2")->GetParError(2);
		blockSlope[1]  = timeBlockHist->GetFunction("expo2")->GetParameter(3);
		blockErr[1]    = timeBlockHist->GetFunction("expo2")->GetParError(3);*/
		
		
		// Block Map Scan Output data 
		/*if (blockSlope[0] < blockSlope[1]) {
			fastIndex = 0;
			slowIndex = 1;
		} else {
			fastIndex = 1;
			slowIndex = 0;
		}
		printf("For all neutrons:\n");
		blockPopErr[fastIndex] = blockPopErr[fastIndex]/blockPop[fastIndex]; // Convert to % error
		blockPop[fastIndex] = blockPop[fastIndex]/blockSlope[fastIndex]*(exp(-(90.0+hold_t)/blockSlope[fastIndex])-exp(-50.0*blockSlope[fastIndex]));
		blockPopErr[fastIndex] = blockPopErr[fastIndex]*blockPop[fastIndex];
		printf("   Fast Decay Constant: %e (%e)\n", blockSlope[fastIndex], blockErr[fastIndex]);
		printf("   Fast Population:     %e (%e)\n", blockPop[fastIndex], blockPopErr[fastIndex]);
		blockPopErr[slowIndex] = blockPopErr[slowIndex]/blockPop[slowIndex]; // Convert to % error
		blockPop[slowIndex] = blockPop[slowIndex]/blockSlope[slowIndex]*(exp(-(90.0+hold_t)/blockSlope[slowIndex])-exp(-50.0/blockSlope[slowIndex]));
		blockPopErr[slowIndex] = blockPopErr[slowIndex]*blockPop[slowIndex];
		printf("   Slow Decay Constant: %e (%e)\n", blockSlope[slowIndex], blockErr[slowIndex]);
		printf("   Slow Population:     %e (%e)\n", blockPop[slowIndex], blockPopErr[slowIndex]);
		printf("\n---------------------------------------------------------------------\n");
	}*/
	//------------------------------------------------------------------
	// Plot weighted information from neutrons hitting blocks
	
/*
	{
	// Start with zeroth pad
	sprintf(tPadName,"p%d",0);
	sprintf(tPadTitle,"Arrival Time E < %f m", threshold); // Name is 
	padMin = 0.01;
	padMax = (1.0/numHistos) - 0.01;
	padc2[0] = new TPad(tPadName,tPadTitle,0.03,padMin,0.97,padMax,21);
	padc2[0]->Draw();
	// Create the remaining pads
	for (int i=1;i<nSlices;i++) { 
		sprintf(tPadName,"p%d",i);
		sprintf(tPadTitle,"Arrival Time %f < E < %f m", (i-1)*spacing + threshold, i*spacing + threshold);
		padMax = ((i+1.0)/numHistos) - 0.01;
		padMin = (i+0.0)/numHistos + 0.01;
	
		padc2[i] = new TPad(tPadName,tPadTitle,0.03,padMin,0.97,padMax,21);
		padc2[i]->SetLogY();
		padc2[i]->Draw();
	} 
	// Fit our remaining histograms
	for (int i=0; i<numHistos;i++) {
		padc2[i]->SetLogy();
		padc2[i]->cd();
		tBlockHist[i]->Draw();
		tBlockHist[i]->Fit("expo2","QL","",50.0,90.0+hold_t);
		slicePop[i][0]    = tBlockHist[i]->GetFunction("expo2")->GetParameter(0);
		slicePopErr[i][0] = tBlockHist[i]->GetFunction("expo2")->GetParError(0);
		sliceSlope[i][0]  = tBlockHist[i]->GetFunction("expo2")->GetParameter(1);
		sliceErr[i][0]    = tBlockHist[i]->GetFunction("expo2")->GetParError(1);
		slicePop[i][1]    = tBlockHist[i]->GetFunction("expo2")->GetParameter(2);
		slicePopErr[i][1] = tBlockHist[i]->GetFunction("expo2")->GetParError(2);
		sliceSlope[i][1]  = tBlockHist[i]->GetFunction("expo2")->GetParameter(3);
		sliceErr[i][1]    = tBlockHist[i]->GetFunction("expo2")->GetParError(3);
		if (sliceSlope[i][0] < sliceSlope[i][1]) {
			fastIndex = 0;
			slowIndex = 1;
		} else {
			fastIndex = 1;
			slowIndex = 0;
		}
		printf("For energy slice between %f and %f cm:\n",(i-1)*spacing + threshold, i*spacing + threshold);
				
		slicePopErr[i][fastIndex] = slicePopErr[i][fastIndex]/slicePop[i][fastIndex]; // Convert to % error
		slicePop[i][fastIndex] = slicePop[i][fastIndex]/sliceSlope[i][fastIndex]*(exp(sliceSlope[i][fastIndex]*(40.0+hold_t))-exp(sliceSlope[i][fastIndex]*50.0));
		slicePopErr[i][fastIndex] = slicePopErr[i][fastIndex]*slicePop[i][fastIndex];
		printf("   Fast Decay Constant: %e (%e)\n", sliceSlope[i][fastIndex], sliceErr[i][fastIndex]);
		printf("   Fast Population:     %e (%e)\n", slicePop[i][fastIndex], slicePopErr[i][fastIndex]);
					
		slicePopErr[i][slowIndex] = slicePopErr[i][slowIndex]/slicePop[i][slowIndex]; // Convert to % error
		slicePop[i][slowIndex] = slicePop[i][slowIndex]/sliceSlope[i][slowIndex]*(exp(sliceSlope[i][slowIndex]*(90.0+hold_t))-exp(sliceSlope[i][slowIndex]*50.0));
		slicePopErr[i][slowIndex] = slicePopErr[i][slowIndex]*slicePop[i][slowIndex];
		printf("   Slow Decay Constant: %e (%e)\n", sliceSlope[i][slowIndex], sliceErr[i][slowIndex]);
		printf("   Slow Population:     %e (%e)\n\n", slicePop[i][slowIndex], slicePopErr[i][slowIndex]);

		gStyle->SetOptFit(1);
		}
		c2->Update();
	}
	//-----------------------------------------------------------------
	
	
	// time should be in terms of maxTime + hold_t
	// Plot the block data by energy slices
		
	double dE = (EMAX - thr) / ((double)num - 1); // energy spacing
	
	char* bName = new char[10];
	char* bTitle = new char[48];
	TH1D* tBlockHist[num];
	
	sprintf(bName,"tBlock0");
	sprintf(bTitle,"hit time E > %f m", thr);
	
	tBlockHist[0] =  new TH1D(bName, bTitle, 300,0,time);
	for (size_t ii=1;ii<num;ii++){
		sprintf(bName,"tBlock%d",ii);
		sprintf(bTitle,"hit time %f < E < %f m", (ii-1)*dE + thr, ii*dE + thr);
		tBlockHist[ii] =  new TH1D(bName, bTitle, 300,0,time);
	}
		
	for (size_t ii = 0; ii < block.size(); ii++) {
		for (size_t jj=1; jj < num; jj++) {
			double minSl = ((double)jj - 1.0)* dE + thr;
			double maxSl = ((double)jj)*dE + thr;
			if (((block.at(ii).energy)/(GRAV*MASS_N) > minSl) && ((block.at(ii).energy)/(GRAV*MASS_N) <= maxSl)) {
				tBlockHist[jj]->Fill((double)block.at(ii).time,wts.at(ii));
			}
		}
	}
	
	TCanvas* c2 = new TCanvas("c2","c2",750,1000);
	c2->cd();
		
	TPad* padc2[num]; // Might have to switch this to vectors
	TH1D* tBlockHist[num];

	TF1* expo2 = new TF1("expo2","[0]*exp(-x/[1]) + [2]*exp(-x/[3])");
	expo2->SetParameters(0.5, 800, 100, 100);
	// Start with zeroth pad
	char* tPadName = new char[4];
	char* tPadTitle = new char[48];
	sprintf(tPadName,"p%d",0);
	sprintf(tPadTitle,"Arrival Time E < %f m", thr);
	
	double padMin = 0.01;
	double padMax = (1.0/num) - 0.01;
	padc2[0] = new TPad(tPadName,tPadTitle,0.03,padMin,0.97,padMax,21);
	padc2[0]->Draw();
	// Create the remaining pads
	for (size_t i=1;i<num;i++) { 
		sprintf(tPadName,"p%d",i);
		sprintf(tPadTitle,"Arrival Time %f < E < %f m", (i-1)*dE + thr, i*dE + thr);
		padMax = ((i+1.0)/num) - 0.01;
		padMin = (i+0.0)/num + 0.01;
	
		padc2[i] = new TPad(tPadName,tPadTitle,0.03,padMin,0.97,padMax,21);
		padc2[i]->Draw();
	} 
	// Fit our histograms
	measurement slicePop[num][2];
	measurement sliceSlope[num][2];
	
	for (int i=0; i<num;i++) {
		padc2[i]->SetLogy();
		padc2[i]->cd();
		tBlockHist[i]->Draw(); // Drawing
		
		tBlockHist[i]->Fit("expo2","QL","",50.0,90.0+holdT); // Fitting
		slicePop.val[i][0]    = tBlockHist[i]->GetFunction("expo2")->GetParameter(0);
		slicePop.err[i][0] = tBlockHist[i]->GetFunction("expo2")->GetParError(0);
		sliceSlope.val[i][0]  = tBlockHist[i]->GetFunction("expo2")->GetParameter(1);
		sliceSlope.err[i][0]    = tBlockHist[i]->GetFunction("expo2")->GetParError(1);
		slicePop.val[i][1]    = tBlockHist[i]->GetFunction("expo2")->GetParameter(2);
		slicePop.err[i][1] = tBlockHist[i]->GetFunction("expo2")->GetParError(2);
		sliceSlope.val[i][1]  = tBlockHist[i]->GetFunction("expo2")->GetParameter(3);
		sliceSlope.err[i][1]    = tBlockHist[i]->GetFunction("expo2")->GetParError(3);
		
		size_t fastIndex = sliceSlope.val[i][0] < sliceSlope.val[i][1]  ? 0 : 1;
		size_t slowIndex = sliceSlope.val[i][0] >= sliceSlope.val[i][1] ? 0 : 1;
		
		printf("For energy slice between %f and %f cm:\n",(i-1)*dE + thr, i*dE + thr);
				
		//slicePop.err[i][fastIndex] = slicePop.err[i][fastIndex]/slicePop.val[i][fastIndex]; // Convert to % error
		double slicePctEF = slicePop.err[i][fastIndex]/slicePop.val[i][fastIndex];
		slicePop.val[i][fastIndex] = slicePop.val[i][fastIndex]/sliceSlope.err[i][fastIndex]*(exp(sliceSlope.val[i][fastIndex]*(40.0+holdT))-exp(sliceSlope.val[i][fastIndex]*50.0));
		//slicePop.err[i][fastIndex] = slicePop.err[i][fastIndex]*slicePop.val[i][fastIndex]; // Convert back to real
		printf("   Fast Decay Constant: %e (%e)\n", sliceSlope.val[i][fastIndex], sliceSlope.err[i][fastIndex]);
		printf("   Fast Population:     %e (%e)\n", slicePop.val[i][fastIndex], slicePop.err[i][fastIndex]);
					
		//slicePop.err[i][slowIndex] = slicePop.err[i][slowIndex]/slicePop.val[i][slowIndex]; // Convert to % error
		double slicePctES = slicePop.err[i][slowIndex]/slicePop.val[i][slowIndex];
		slicePop.val[i][slowIndex] = slicePop.val[i][slowIndex]/sliceSlope.val[i][slowIndex]*(exp(sliceSlope.val[i][slowIndex]*(90.0+holdT))-exp(sliceSlope.val[i][slowIndex]*50.0));
		//slicePop.err[i][slowIndex] = slicePop.err[i][slowIndex]*slicePop.val[i][slowIndex];
		printf("   Slow Decay Constant: %e (%e)\n", sliceSlope.val[i][slowIndex], sliceSlope.err[i][slowIndex]);
		printf("   Slow Population:     %e (%e)\n\n", slicePop.val[i][slowIndex], slicePop.err[i][slowIndex]);

		gStyle->SetOptFit(1);
	}
	c2->Update();
	return;
}*/

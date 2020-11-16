	
	TH1D* timeDagHist = new TH1D("timeDag", "hit time (dagger)", nBins1,-150.0,maxTime); 
	TH1D* eneDagHist = new TH1D("eneDag", "energy (dagger)",300, 0.0, maxEnergy+minEnergy);
	TH1D* thetaDagRoot = new TH1D("thetaDag", "theta (dagger)",300, 0, 1.7); 
	
	
	
	
	// Initialize slice histograms
	// Start with zeroth 
	
	
// CANVASSES: 
// 2: block slice data
// 3: comparison between loaded data and MC simulation (on the dagger)
// 4: B field slice data
// 5: Chi^2 fitting region plot
// 6: energy/angle fitting sanity check
// 7: position on the block

void plotBlockMapTot(std::vector<double> blockWt, std::vector<fixedResultBlock> blockRes, double maxT) {
	
	/*----------------------------------------------------------------*/
	/* !!!PLOTTING AND ANALYSIS!!! */
	/*----------------------------------------------------------------*/
	/* CANVASSES: 
	 * 1: block mapping (block vs. hit time, block hits, etc.
	 * 2: block slice data
	 * 3: comparison between loaded data and MC simulation (on the dagger)
	 * 4: B field slice data
	 * 5: Chi^2 fitting region plot
	 * 6: energy/angle fitting sanity check
	 * 7: position on the block */
	 
	if (blockMapTot) {

		// Generation of histograms on canvas
		TCanvas* c1 = new TCanvas("c1","c1",750,500);
		TPad* pad1c1 = new TPad("p1","Energy vs Time",0.03,0.03,0.645,0.97,21);
		TPad* pad2c1 = new TPad("p2","Arrival Time",0.675,0.03,0.97,0.485,21);
		TPad* pad3c1 = new TPad("p3","Neutron Energy",0.675,0.515,0.97,0.97,21);
		pad1c1->Draw();
		pad2c1->Draw();
		pad3c1->Draw();	
		pad1c1->cd();
		blockHistRoot->Draw("CONT4Z");
		pad2c1->cd();
		timeBlockHist->Scale(1/timeBlockHist->Integral());
		timeBlockHist->Draw();
		timeBlockHist->Fit("expo2","WL M","",200.0,240.0+hold_t);
		/*timeBlockHist->Fit("expo2","QL","",200.0,240.0+hold_t);
		 
		blockPop[0]    = timeBlockHist->GetFunction("expo2")->GetParameter(0);
		blockPopErr[0] = timeBlockHist->GetFunction("expo2")->GetParError(0);
		blockSlope[0]  = timeBlockHist->GetFunction("expo2")->GetParameter(1);
		blockErr[0]    = timeBlockHist->GetFunction("expo2")->GetParError(1);
		blockPop[1]    = timeBlockHist->GetFunction("expo2")->GetParameter(2);
		blockPopErr[1] = timeBlockHist->GetFunction("expo2")->GetParError(2);
		blockSlope[1]  = timeBlockHist->GetFunction("expo2")->GetParameter(3);
		blockErr[1]    = timeBlockHist->GetFunction("expo2")->GetParError(3);*/
		pad3c1->cd();
		eneBlockHist->Scale(1/eneBlockHist->Integral());
		eneBlockHist->Draw();
		c1->Update();
		
		// Block Map Scan Output data 
		if (blockSlope[0] < blockSlope[1]) {
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
	}
	//------------------------------------------------------------------
	// Plot weighted information from neutrons hitting blocks
	
	// Convert vectors to ROOT (weighted) histograms
	TH2D* blockHistRoot = new TH2D("blockHist", "hit time vs. energy", 30, 0, maxT, 30, EMIN, EMIN+EMAX);
	TH1D* timeBlockHist = new TH1D("timeBlock", "hit time (block)", 500,0,maxT);
	TH1D* eneBlockHist = new TH1D("eneBlock", "energy (block)",300, EMIN, EMIN+EMAX);
	for(size_t ii = 0; ii < blockRes.size(); ii++) {
		blockHistRoot->Fill(blockRes.at(ii).time, (blockRes.at(ii).energy)/(GRAV*MASS_N),blockWt.at(ii));
		timeBlockHist->Fill((double)blockRes.at(ii).time,blockWe.at(ii));
		eneBlockHist->Fill((blockRes.at(ii).energy)/(GRAV*MASS_N),blockWt.at(ii));
	}

	// Generation of histogram canvas
	TCanvas* c1 = new TCanvas("c1","c1",750,500); // 2D energy and time
	TPad* pad1c1 = new TPad("p1","Energy vs Time",0.03,0.03,0.645,0.97,21);
	pad1c1->Draw();
	pad1c1->cd();
	blockHistRoot->Draw("CONT4Z");
	
	TPad* pad2c1 = new TPad("p2","Arrival Time",0.675,0.03,0.97,0.485,21); // 1D Time
	pad2c1->Draw();
	pad2c1->cd();
	timeBlockHist->Scale(1/timeBlockHist->Integral());
	timeBlockHist->Draw();
	//timeBlockHist->Fit("expo2","WL M","",200.0,240.0+hold_t);
	
	TPad* pad3c1 = new TPad("p3","Neutron Energy",0.675,0.515,0.97,0.97,21); // 1D energy
	pad3c1->Draw();	
	pad3c1->cd();
	eneBlockHist->Scale(1/eneBlockHist->Integral());
	eneBlockHist->Draw();
		
	c1->Update();
}

void plotBlockMapSlice(std::vector<double> wts, std::vector<fixedResultBlock> block, const size_t num, double thr, double time) {
	
	char* tPadName = new char[4];
	char* tPadTitle = new char[48];
	
	// Slice data on the block
	if (blockMapSlice) {
		TCanvas* c2 = new TCanvas("c2","c2",750,1000);
		c2->cd();
		TPad* padc2[numHistos];
					
		// Start with zeroth pad
		sprintf(tPadName,"p%d",0);
		sprintf(tPadTitle,"Arrival Time E < %f m", threshold);
		padMin = 0.01;
		padMax = (1.0/numHistos) - 0.01;
		padc2[0] = new TPad(tPadName,tPadTitle,0.03,padMin,0.97,padMax,21);
		padc2[0]->Draw();
		// Create the remaining pads
		for (int i=1;i<numHistos;i++) { 
			sprintf(tPadName,"p%d",i);
			sprintf(tPadTitle,"Arrival Time %f < E < %f m", (i-1)*spacing + threshold, i*spacing + threshold);
			padMax = ((i+1.0)/numHistos) - 0.01;
			padMin = (i+0.0)/numHistos + 0.01;
		
			padc2[i] = new TPad(tPadName,tPadTitle,0.03,padMin,0.97,padMax,21);
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
}


//}

double dataMCComp(double blah) {
	// Data vs. MC simulations
	/*if (dataMCComp) {
		
		TCanvas* c3 = new TCanvas("c3","c3",750,500);
		TPad* pad1c3 = new TPad("p1","Arrival Time (Dagger)",0.03,0.03,0.97,0.97,21);
		pad1c3->Draw();
		pad1c3->cd();
		timeDagHist->Scale(1/timeDagHist->Integral());
		timeDagHist->Draw();
		if (chiFit) {
			if (dataHist->GetSize() == timeDagHist->GetSize()){
				dataHist->Scale(1/dataHist->Integral());
				dataHist->SetLineColor(kRed);
				dataHist->Draw("SAME");
			}
		}	
		if (saveMCFit) {
			timeDagHist->SaveAs("./timeDagHist.root");
			printf("\n !!!SAVING DATA AS timeDagHist.root!!! \n");
		}

		if (drainingTimeOn) {timeDagHist->Fit("expo","Q","",45.0,59.0);
			drainFit[0] = timeDagHist->GetFunction("expo")->GetParameter(1);
			drainErr[0] = timeDagHist->GetFunction("expo")->GetParError(1);
			timeDagHist->Fit("expo","Q","",65.0,90.0);
			drainFit[1] = timeDagHist->GetFunction("expo")->GetParameter(1);
			drainErr[1] = timeDagHist->GetFunction("expo")->GetParError(1);
			printf("\n Draining Time: \n%e (Group 2)\n%e (Group 3)\n\n",drainFit[0],drainFit[1]);
		}		
		c3->Update();
	}
	return 3.0;*/
	if (dataMCComp) {
		
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

		if (drainingTimeOn) {timeDagHist->Fit("expo","Q","",45.0,59.0);
			drainFit[0] = timeDagHist->GetFunction("expo")->GetParameter(1);
			drainErr[0] = timeDagHist->GetFunction("expo")->GetParError(1);
			timeDagHist->Fit("expo","Q","",65.0,90.0);
			drainFit[1] = timeDagHist->GetFunction("expo")->GetParameter(1);
			drainErr[1] = timeDagHist->GetFunction("expo")->GetParError(1);
			printf("\n Draining Time: \n%e (Group 2)\n%e (Group 3)\n\n",drainFit[0],drainFit[1]);
		}		
		c3->Update();
	}
}


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
		/*for (int i=1;i<numHistos;i++) { 
			sprintf(fPadName,"p%d",i);
			sprintf(fPadTitle,"B Field %f < E < %f m", (i-1)*spacing + threshold, i*spacing + threshold);
			padMax = ((i+1.0)/numHistos) - 0.01;
			padMin = (i+0.0)/numHistos + 0.01;
			padc4[i] = new TPad(fPadName,fPadTitle,0.03,padMin,0.97,padMax,21);
			padc4[i]->Draw();
		} */
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

//}

void plotChiMap(double blah) {
	if (chiMap) {
		TCanvas* c5 = new TCanvas("c5","c5",750,750);
		
		// Find how many chi2 pads we require
		int nPads=0;
		int nPadsi,nPadsj;
		nPadsi = chi2Dim;
		nPadsj = chi2Dim;
		for (int i=0;i<chi2Dim;i++) { 
			for (int j=0;j<i;j++) {
				if (c2Step[i]!=0 && c2Step[j]!=0){
					nPads=nPads+1;
				}
			} 
			if (c2Step[i]==0) {
					nPadsi=nPadsi-1;
			} 
		}
		printf("%d %d\n",nPads, nPadsi);
		int mapC=0;
		// Create our pads
		TPad* padc5[nPads];
		double padMaxj, padMaxi, padMinj, padMini;
		for (int i=0;i<chi2Dim;i++) { 
			for (int j=0;j<i;j++) {
				if (c2Step[i]!=0.0 && c2Step[j]!=0.0) {
					sprintf(tPadName,"p%d%d",i,j);
					sprintf(tPadTitle,"chi2 fit, variation %d%d", i,j);
					
					padMaxi = (((double)i+0.0)/(double)(nPadsi-1)) - 0.01;
					if (padMaxi > 1.0) { padMaxi = padMaxi - 1.0; }
					padMaxj = (((double)j+1.0)/(double)(nPadsi-1)) - 0.01;
					if (padMaxj > 1.0) { padMaxj = padMaxj - 1.0; }
					
					padMini = (((double)i-1.0)/(double)(nPadsi-1)) + 0.01;
					if (padMini > 1.0) { padMini = padMini - 1.0; }
					padMinj = (((double)j+0.0)/(double)(nPadsi-1)) + 0.01;
					if (padMinj > 1.0) { padMini = padMinj - 1.0; }
					
					printf("%f, %f, %f, %f, %d\n", padMini, padMaxi, padMinj, padMaxj, mapC);
					
					padc5[mapC] = new TPad(tPadName,tPadTitle,padMini,padMinj,padMaxi,padMaxj,21);
					padc5[mapC]->Draw();
					mapC=mapC+1;
				}
			}
		}
		// Load our data
		TGraph2D* chi2Map[nPads];
		mapC=0;
		int nBins;
		char* chi2MinName = new char[48];
		//sprintf(chi2MinName, "%f",
		for (int i=0;i<chi2Dim;i++) { 
			for (int j=0;j<i;j++) {
				if (c2Step[i]!=0.0 && c2Step[j]!=0.0) {
					sprintf(c2Filename,"chi2Fit%d%d.csv",i,j);
					chi2Map[mapC] = new TGraph2D(c2Filename);
						chi2Map[mapC]->SetName(c2Filename);
						chi2Map[mapC]->SetTitle("Chi^2 Fit");
					nBins = chi2Map[mapC]->GetHistogram()->GetSize();
					for (int k=0;k<nBins;k++) {
						if (chi2Map[mapC]->GetHistogram()->GetBinContent(k) > (chi2Min+3.0)){
							chi2Map[mapC]->GetHistogram()->SetBinContent(k,0.0);
						}
					}					
					padc5[mapC]->cd();
					printf("Drawing %d%d...\n", i,j);
					chi2Map[mapC]->SetMaximum(chi2Min+3.0);
					chi2Map[mapC]->SetMinimum(chi2Min);
					chi2Map[mapC]->Draw("COLZ1");
					mapC=mapC+1;
					c5->Update();
				}
			}
		}
		printf("Updating Canvas\n");
		//c5->Update();
		printf("Done Updating!\n");
	}
	
	if (chiMap) {
		TCanvas* c5 = new TCanvas("c5","c5",750,750);
		
		/* Find how many chi2 pads we require */ 
		int nPads=0;
		for (int i=0;i<chi2Dim;i++) { 
			for (int j=0;j<i;j++) {
				if (c2Step[i]!=0 && c2Step[j]!=0) {
					nPads=nPads+1;
				}
			}
		}
		/* Create our pads */
		TPad* padc5[nPads];
		for (int i=0;i<nPads;i++) { 
			
				sprintf(tPadName,"p%d",i);
				sprintf(tPadTitle,"chi2 fit, variation %d", i);
				padMax = ((i+1.0)/nPads) - 0.01;
				padMin = (i+0.0)/nPads + 0.01;
				padc5[i] = new TPad(tPadName,tPadTitle,0.03,padMin,0.97,padMax,21);
				padc5[i]->Draw();
		} 
		/* Load our data */
		TGraph2D* chi2Map[nPads];
		int mapC=0;
		for (int i=0;i<chi2Dim;i++) { 
			for (int j=0;j<i;j++) {
				if (c2Step[i]!=0.0 && c2Step[j]!=0.0) {
					sprintf(c2Filename,"chi2Fit%d%d.csv",i,j);
					chi2Map[mapC] = new TGraph2D(c2Filename);
						chi2Map[mapC]->SetName(c2Filename);
						chi2Map[mapC]->SetTitle("Chi^2 Fit");
					padc5[mapC]->cd();
					chi2Map[mapC]->Draw("CONT4Z");
					mapC=mapC+1;
				}
			}
		}
		c5->Update();
	}
	return;
}

void plotSanityCheck(double blah) {
	if (sanityCheck) {
		TCanvas* c6 = new TCanvas("c6","c6",750,750);
		TPad* pad1c6 = new TPad("p1","Energy (all neutrons)",0.03,0.03,0.97,0.485,21);	
		TPad* pad2c6 = new TPad("p2","Theta (dagger)",0.03,0.515,0.97,0.97,21);
		pad1c6->Draw();
		pad2c6->Draw();
		pad1c6->cd();	
		eneDagHist->Scale(1/eneDagHist->Integral());
		eneDagHist->SaveAs("./eneDagHist.root");
		eneDagHist->Draw();
		pad2c6->cd();
		thetaDagRoot->Scale(1/thetaDagRoot->Integral());
		thetaDagRoot->SaveAs("./thetaDagHist.root");
		thetaDagRoot->Draw();
		c6->Update();
	}

	if (sanityCheck) {
		TCanvas* c6 = new TCanvas("c6","c6",750,750);
		TPad* pad1c6 = new TPad("p1","Energy (all neutrons)",0.03,0.03,0.97,0.485,21);	
		TPad* pad2c6 = new TPad("p2","Theta (dagger)",0.03,0.515,0.97,0.97,21);
		pad1c6->Draw();
		pad2c6->Draw();
		pad1c6->cd();	
		eneDagHist->Scale(1/eneDagHist->Integral());
		eneDagHist->SaveAs("./eneDagHist.root");
		eneDagHist->Draw();
		pad2c6->cd();
		thetaDagRoot->Scale(1/thetaDagRoot->Integral());
		thetaDagRoot->SaveAs("./thetaDagHist.root");
		thetaDagRoot->Draw();
		c6->Update();
	}	
	return;
}

void plotBlockMap(double blah) {
	if (blockMap) {
		TCanvas* c7 = new TCanvas("c7","c7",500,500);
		TPad* pad1c7 = new TPad("p1","Block Hit Position",0.03,0.03,0.97,0.97,21);
		pad1c7->Draw();
		pad1c7->cd();
		hitPos->GetHistogram()->SetMinimum(-1.46);
		hitPos->GetHistogram()->SetMaximum(-1.42);
		hitPos->GetYaxis()->SetRangeUser(0.10,0.15);
		hitPos->GetXaxis()->SetRangeUser(0.20,0.25);
	 	hitPos->Draw("PCOL");
	 	c7->Update();
	 	printf("MAX: (%f, %f, %f);\nMIN: (%f, %f, %f).\n",maxCoords[0],maxCoords[1],maxCoords[2],minCoords[0],minCoords[1],minCoords[2]);
	}
	return;
	
		if (blockMap) {
		TCanvas* c7 = new TCanvas("c7","c7",500,500);
		TPad* pad1c7 = new TPad("p1","Block Hit Position",0.03,0.03,0.97,0.97,21);
		pad1c7->Draw();
		pad1c7->cd();
		hitPos->GetHistogram()->SetMinimum(-1.46);
		hitPos->GetHistogram()->SetMaximum(-1.42);
		hitPos->GetYaxis()->SetRangeUser(0.10,0.15);
		hitPos->GetXaxis()->SetRangeUser(0.20,0.25);
	 	hitPos->Draw("PCOL");
	 	c7->Update();
	 	printf("MAX: (%f, %f, %f);\nMIN: (%f, %f, %f).\n",maxCoords[0],maxCoords[1],maxCoords[2],minCoords[0],minCoords[1],minCoords[2]);
	}
}
	
void plotPercentDag(double blah) { 
	TCanvas* c8 = new TCanvas("c8","c8",750,750);
	double totalHits, percentage;
	TH1D* percentDag = new TH1D("percentDag", "% of neutrons hitting block",300, minEnergy, maxEnergy + minEnergy);
	
	for (int ii = 0; ii < length_dag; ii++)	{
		totalHits = eneDagHist->GetBinContent(ii)+eneDagOnBlock->GetBinContent(ii);
		if(totalHits != 0.0){
			percentage = eneDagHist->GetBinContent(ii)/totalHits;
		} else{ 
			percentage = 1.0;
		}
		percentDag->SetBinContent(ii, 1.0-percentage);
	}
	percentDag->Draw();
	//eneDagHist->Add(eneBlockHist);
	
}

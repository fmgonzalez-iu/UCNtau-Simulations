#include "../inc/blockReader.hpp"

/*----------------------------------------------------------------------
 * 
 * This code features functions for analyzing the effect of the Aluminum
 * block. 
 * 
 * The actual plotting functions are in a separate file.
 * 
 * This is a mess.
 * 
 *--------------------------------------------------------------------*/

double fieldstrength(double x, double y, double z) {  
    // Calculate the strength of the magnetic field at a given position
        
	double R,r;
	double A = sqrt(8.0) * B_REM / M_PI; //parameter related to B -- shows up in expansion

	if (x > 0.0) {
			R = 1.0;
			r = 0.5;
	} else {
			R = 0.5;
			r = 1.0;
	}

	double rho = sqrt(y*y+z*z);
	double r_zeta = sqrt((rho-R)*(rho-R)+x*x);
	double B_halbach = 0.0;

	if (z < -1.0 && r_zeta < r) {
			double zeta = r-r_zeta;
			double eta = r*atan(x/(rho-R));
			double Bsum = 0.0;
			double k_m,k_n,m,n;

			for (m = 1.0; m <= N_TERMS; m += 1.0) {
					k_m = 2*M_PI * (4.0*m - 3.0) / MAG_SPACE;
					for (n = 1.0; n <= N_TERMS; n += 1.0) {
							k_n = 2*M_PI * (4.0*n-3.0) / MAG_SPACE;
							Bsum += pow(-1.0,m)*pow(-1.0,n)/(4.0*m-3.0)/(4.0*n-3.0)*(1-exp(-k_m*MAG_THICK))*(1-exp(-k_n*MAG_THICK))*exp(-(k_n+k_m)*zeta)*cos((k_n-k_m)*eta);
					}
			}
			B_halbach = sqrt(B_HOLD*B_HOLD*(r+R)*(r+R)/(y*y+z*z)+A*A*Bsum);
	}

	return B_halbach;
}
double getNorm(std::vector<measurement> groups) {
	// Get a normalized sum from a single group

	printf("Normalizing group weights!\n");
	double norm = 0.0; // Need to scale all events by some amount
	for (auto it = groups.begin(); it < groups.end(); it++) { // Loop through all measurements, only care about value
		norm += *it.val;
	}
	return norm;
}
// Alumunum (and Oxide) block-------------------------------------------
double alBlockLossProb(double perp_energy) {
	// Hardcoded solid aluminum block loss probability
		
	return 1.0e-3;
}
double calcBlockScatter(double perp_energy) {
	// This is a surface calculation for aluminum and aluminum oxide
	
	const double vAlOxide = (2*M_PI*(HBAR*HBAR/MASS_N))*(AALUMINUM*NALUMINUMAL3O2 + AOXYGEN*NOXYGENAL3O2);
    const double wAlOxide = (HBAR/2)*(NALUMINUMAL3O2*SIGMAALUMINUM + NOXYGENAL3O2*SIGMAOXYGEN);
    const double vAl = (2*M_PI*(HBAR*HBAR/MASS_N))*AALUMINUM*NALUMINUM;
    const double wAl = (HBAR/2)*NALUMINUM*SIGMAALUMINUM;
    
    std::vector<std::complex<double>> potentials = {std::complex<double>(0.0,0.0), 
													std::complex<double>(vAlOxide,-wAlOxide), 
													std::complex<double>(vAl, -wAl)};
    std::vector<double> thickness = {0.0,10.0e-9,10000.0e-9}; // Hardcoded AlOxide thickness
    std::vector<std::complex<double>> mbar = {std::complex<double>(1.0,0.0), std::complex<double>(0.0,0.0), std::complex<double>(0.0,0.0), std::complex<double>(1.0,0.0)};
    
    for (int i = potentials.size(); i > 0; i--) { 
		mbar = matmul(mbar, m(k(perp_energy, potentials[i]), k(perp_energy, potentials[i-1]), thickness[i-1]));
	}
	
	return 1.0 - (std::conj(-mbar[2]/mbar[3])*-mbar[2]/mbar[3]).real();
	
}
double calcZeta(double x, double y, double z) {

	double R, r, rho, r_zeta;
	
	if (x > 0.0) { 
		R = 1.0;
		r = 0.5;
	} else {
		R = 0.5;
		r = 1.0;
	}
	
	rho = sqrt(y*y+z*z);
	r_zeta = sqrt((rho-R)*(rho-R)+x*x);
		
	return r-r_zeta;
}
double calcEta(double x, double y, double z) {
	
	double R, r, rho, r_zeta;
	
	if (x > 0.0) { 
		R = 1.0;
		r = 0.5;
	} else {
		R = 0.5;
		r = 1.0;
	}

	rho = sqrt(y*y+z*z);
	return r*atan(x/(rho-R));
}
//----------------------------------------------------------------------
std::vector<fixedResultDag> parseDagEvents(std::vector<fixedResultDag> dag_events, std::vector<fixedResultBlock> block_events) {
	// This removes events on the dagger that would upscatter off the block.
	
	std::vector<fixedResultDag> surviving_evts;
	// Begin filling events
	size_t bCtr = 0;
	for (auto it = 0; it < dag_events.end(); i++) {
								
		bool upscatter = false; // block scattering
		if ((*it).nHitBlock > 0) { // Check if the dagger has block hits. If so we have to deal.
			// Match first event. The next N will sequentially follow from first to last hit.
			// If the first guess doesn't work we have to scan
			for (size_t j = bCtr; bCtr < events_block.size(); j++) {
				if ((block_events[j].theta == (*it).theta) || (block_events[j].phi == (*it).phi)) {
					bCtr = j;
					break;
				}				
			}
			// Now find the scattering probability for each subsequent hit
			for (size_t j = bCtr; j < bCtr + (*it).nHitBlock; j++ ) { 
					
				if (j >= block_events[j].size()) { // Error handling
					printf(" WARNING: Have too many block hits associated with this dagger hit ");
					break;
				}
				
				if (randAbsProb[j] < calcBlockScatter(block_events[j].perp_energy)) {// Realistic aluminum block
				//if (randAbsProb[j] < alBlockLossProb(block_events[j].perp_energy)) { Single efficiency block
					upscatter = true;
					//block_positions.append(j); // Testing
					break;
				}
			}
		}
		if (!upscatter) {
			surviving_evts.append((*it));
		}
	}
	return surviving_evts;
}
std::vector<fixedResultBlock> parseBlockEvents(std::vector<fixedResultDag> dag_events, std::vector<fixedResultBlock> block_events) {
	// This removes events on the dagger that would upscatter off the block.
	
	std::vector<fixedResultBlock> reduced_events;
	// Begin filling events
	size_t bCtr = 0;
	for (auto it = 0; it < dag_events.end(); i++) {

		if ((*it).nHitBlock > 0) { // Check if the dagger has block hits. If so we have to deal.
			// Match first event. The next N will sequentially follow from first to last hit.
			// If the first guess doesn't work we have to scan
			for (size_t j = bCtr; bCtr < events_block.size(); j++) {
				if ((block_events[j].theta == (*it).theta) || (block_events[j].phi == (*it).phi)) {
					bCtr = j;
					break;
				}				
			}
			// Now find the scattering probability for each subsequent hit
			for (size_t j = bCtr; j < bCtr + (*it).nHitBlock; j++ ) { 
					
				if (j >= block_events[j].size()) { // Error handling
					printf(" WARNING: Have too many block hits associated with this dagger hit ");
					break;
				}
				
				if (randAbsProb[j] < calcBlockScatter(block_events[j].perp_energy)) {// Realistic aluminum block
				//if (randAbsProb[j] < alBlockLossProb(block_events[j].perp_energy)) { Single efficiency block
					reduced_events.append(block_events[j]); // Testing
					break;
				}
			}
		}
	}
	return reduced_events;
}

std::vector<std::vector<double>> chi2_results = block_chi2_sweep(std::vector<evt> reduced_dag, std::vector<double> ref_hist, double* randU01s, double* randDeathTimes) {
		
	// Calculate the chi^2 for a bunch of numbers across a parameter sweep
	std::vector<std::vector<double>> chi2Result;
	int nBins = 10; 
	#pragma omp parallel for collapse(4) // parallelize a square loop to divide into 4 sections
	for(int i = 0; i < nBins+1; i++) {
		for(int j = 0; j < nBins+1; j++) {
			for(int k = 0; k < nBins+1; k++) {
				for(int l = 0; l < nBins+1; l++) {
					std::vector<double> sample;
										
					// Define fitting parameters, relating to parts of the loop -- can load more if needed
					double thresh     = 0.0  + 14.0* i / (double)nBins; 
					double thickBoron = 4.0  + 2.0 * k / (double)nBins;
					double power      = 0.75 + 1.75* l / (double)nBins;
					double cosPower   = 0.0  + 0.5 * j / (double)nBins;
										
//                    std::vector<weightedBin> hist1 = createHistQuantNoOxEPowdEThetaSpline(thickBoron, thresh, power, cosPower, events, randU01s, randDeathTimes);
					// Create a weighted histogram
					std::vector<weightedBin> hist_tmp = createHistQuantSpline(thickBoron, thresh, power, cosPower, reduced_dag, randU01s, randDeathTimes);
					double chisq = calcChisqGagunashvili(refHist, hist_tmp)/(hist_tmp.size()-1); // Calculate Chi^2
										
					printf("%f %f %f %f %f\n", cosPower, thickBoron, thresh, power, chisq); // Print
					fflush(stdout);
					
					// Load vectors
					sample.push_back(chisq);
					sample.push_back(power);
					sample.push_back(thresh);
					sample.push_back(cosPower);
					sample.push_back(thickBoron);
					
					chi2Result.push_back(sample);
				}
			}
		}
	}
	
	return chi2Result;
}
// check if the block is upscattering a neutron, if so find which hit it upscatters on
int blockScatterNum(dag_evt dag_buffer, std::vector<block_evt> events_blok, int counter, double* randUps) {
	
	// load necessary events from input
	int nHits = dag_buffer.nHits[3];
	int hit0 = counter;
	int event = counter;
	int length_blok = events_blok.size();
	double thD  = dag_buffer.theta;
	double phiD = dag_buffer.phi;
	
	// other variables
	//int event;
	double abs = 1.0;
	//printf("%d\n,",length_blok);
	block_evt block_buffer = events_blok[hit0];
	
	// match events, trying to progress normally through the block buffer.	
	if ((block_buffer.theta != thD) || (block_buffer.phi != phiD)) {
		// if given incorrect initial counter, loop through
		//printf("Initial block counter is incorrect, looping through output! \n");
		//printf("%d, %d\n", hit0, event);
		for (event=hit0; event < length_blok; event++) {
			block_buffer = events_blok[event];
			if ((block_buffer.theta == thD) && (block_buffer.phi == phiD)) {
				hit0 = event;
				break;
			}
		}
		//printf("%d, %d\n", hit0-counter, nHits);	
	}
	
	// check on which event we upscatter
	//printf("%d, %d\n", hit0, event);
	while (event < hit0 + nHits) {	
		
		// Error handling 
		if (event >= length_blok || event < 0) {
			printf(" WARNING: Have too many block hits associated with this dagger hit ");
			event = 0;
			break;
		}
		
		// calculate the loss per bounce
		block_buffer = events_blok[event];
		abs = calcBlockScatter(block_buffer.perp_energy);
				
		if (randUps[event] < abs) { 
			printf("%d, %d\n", hit0, event);
			break;
		}
		event = event + 1;
	}
	
	return event-counter;
}

// find if the neutron is actually absorbed or if it's cut off by thresholds
bool checkRealN(dag_evt buffer, double eThres, double theCut, double dagTap, double randDagProb) {
					
	// check if we get absorbed by the bad part of the dagger
	double zeta = calcZeta(buffer.x,buffer.y,buffer.z);
	double boronProb = zeta > dagTap ? 1 : zeta/dagTap;
	
	// filter real N out
	if ((buffer.ecm*JTONEV > eThres) && (buffer.th < theCut) && (randDagProb < boronProb)) {
		return true;
	} else {
		return false;
	}
}


double chi2fitROOT(TH1D* dataHist, TH1D* mcHist, double min1, double min2, double max1, double max2) {
	// chi^2 program because ROOT's arbitrary program is doing something wrong
	// This is fitting both the time constant and the peak height for each dip
	
	
	double chi2Sum = 0.0; 
	if (dataHist->GetSize() != mcHist->GetSize()) { // error catching
		printf("ERROR: chi2fit() has two histograms of incompatable sizes!\n");
		return chi2Sum;
	}
	
	// make sure histograms are properly scaled
	dataHist->Scale(1/dataHist->Integral());
	mcHist->Scale(1/mcHist->Integral());
	
	// TODO turn this into a loop over vector
	TF1* invExp = new TF1("exp","[0]*exp(-x/[1])"); // define logarithmic function for slope
	invExp->SetParameters(0.5, 800);//, 100, 100);
	invExp->SetParLimits(0,0.0,1.0);
	invExp->SetParLimits(1,0.0,2000.0);
	invExp->SetParName(0,"A");
	invExp->SetParName(1,"tau");
	
	int entries = mcHist->GetEntries();	
	
	// find slope of group 2
	dataHist->Fit("invExp","0QL","",min1,min2);
	double s2Dat = dataHist->GetFunction("invExp")->GetParameter(1);
	mcHist->Fit("invExp","0QL","",min1,min2);
	double s2MC = mcHist->GetFunction("invExp")->GetParameter(1);
	
	// find slope of group 3
	dataHist->Fit("expo","0QL","",max1,max2);
	double s3Dat = - 1/dataHist->GetFunction("expo")->GetParameter(1);
	mcHist->Fit("expo","0QL","",max1,max2);
	double s3MC = - 1/mcHist->GetFunction("expo")->GetParameter(1);
	
	// find group yields for the two groups
	g2MC = mcHist->Integral(min1,min2)*entries;
	double g2Dat = dataHist->Integral(min1,min2)*entries;
	g3MC = mcHist->Integral(max1,max2)*entries;
	double g3Dat = dataHist->Integral(max1,max2)*entries;
	
	// sum up chi2 values
	chi2Sum = (s2MC - s2Dat)*(s2MC - s2Dat)/s2Dat + (s3MC - s3Dat)*(s3MC - s3Dat)/s3Dat
			+ (g2MC - g2Dat)*(g2MC - g2Dat)/g2Dat + (g3MC - g3Dat)*(g3MC - g3Dat)/g3Dat;
	
	chi2Sum = chi2Sum / 4.0; // rescale to get chi2/NDF
	return chi2Sum;
}

std::vector<weightedBin> genWeightedDag(std::vector<dag_evt> events_dag, std::vector<block_evt> events_block, double* randDeathTimes){
	// Generate a block weighted dagger list
	
	// For the block we want to parse the dagger and block hits
	std::vector<size_t> block_indices;
	std::vector<dag_evt> real_dag_hits = parseBlockEvents(events_dag, events_block, block_indices, randAbsProb);
	
	std::vector<block_evt> real_block_hits; // Parse the block
	for (auto bInd = block_indices.begin(); bInd < block_indices.end(); bInd++) {
		real_block_hits.append(events_block[block_indices]);
	}
	
	
	
}
std::vector<double> blockChiFit(TH1D* dataHist, std::vector<dag_evt> events_dag, std::vector<block_evt> events_block, double* c2Guess, double* randDeathTimes) {
	// Load TH1D* dataHist somewhere else. 
	// Note that you'll have to GetXaxis()->SetRangeUser(60.0,140.0)
	
	std::vector<double> c2Fitted;

	double c2Guess[6] = {2.296875,12.031250, 1.015625, 1.8, 1.0, 0.017813}; // Hardcoded in here
	double c2Step[6] =  {0.3, 2.0, 0.25, 0.0,0.016076,0.0};
	
	// For the block we want to parse the dagger and block hits
	std::vector<size_t> block_indices;
	std::vector<dag_evt> real_dag_hits = parseBlockEvents(events_dag, events_block, block_indices, randAbsProb);
	
	// IDEA: Step around in an N-D hypercube, and minimize differences in chi2 to find global minima
	for (int step=0; step < 6; step++) {
		
		std::vector<double> chi2Bounds;// Fill the chi2 bounds before doing anything else
		for (size_t dim=0; dim < c2Guess.size(); dim++){
			chi2Bounds.append(c2Guess[dim] - c2Step[dim]);
			chi2Bounds.append(c2Guess[dim] + c2Step[dim]);
		}
		
		// Load our energies from the minimum chi2guesses
		eneMin = c2Guess[0];
		eThMin = c2Guess[1];
		theMin = c2Guess[2];
		theCut = c2Guess[3];
		absMin = c2Guess[4];
		dagCut = c2Guess[5];
		// Start by calculating middle of hypercube
		// Load our relevant dagger events on a TH1D
		TH1D* timeDagFitting = new TH1D("timeDagFit", "hit time (dagger)", nBins,0,maxTime);
	
		for (size_t i = 0; size_t < real_dag_hits.size(); i++) {
			if (real_dag_hits[i].time < randDeathTimes[i] && 
				checkRealN(real_dag_hits[i], c2Guess[1],c2Guess[3],c2Guess[5], randDeathTimes[i])) {
				double histoWeight = calcWeight(real_dag_hits[i].energy, real_dag_hits[i].theta,c2Guess[0],c2Guess[2],c2Guess[1]);
				timeDagFitting->Fill((double)real_dag_hits[i].time+15.59, histoWeight);
			}
		}
		
		chi2Vals[(int)pow(2,chi2Dim)+1] = chi2fitROOT(dataHist, timeDagFitting, 59.0, 75.0, 81.0, 140.0);
		// Flag errors in chi2Vals:
		if (chi2Vals[(int)pow(2,chi2Dim)+1] < 1.0) {
			chi2Vals[(int)pow(2,chi2Dim)+1] = DBL_MAX;
		}
		
		// Save memory and print output
		delete timeDagFitting;
		for (int i = 0; i < chi2Dim; i++ ){
			printf("%f ", c2Guess[i]);
		}
		printf("%f \n",chi2Vals[(int)pow(2,chi2Dim)+1]);
		//printf("%f %f %f %f %f \n",eneMin,eThMin,theMin,absMin,chi2Vals[(int)pow(2,chi2Dim)+1]);
		// start by finding the minumum
		chi2Min = chi2Vals[(int)pow(2,chi2Dim)+1];
		minDex = (int)pow(2,chi2Dim)+1;		
		
		// Look at each endpoint. Potentially skip some.
		for (int b = 0; b < (pow(2,chi2Dim));b++) { 
						
			// Calculate our fit parameters here. Doing binary math
			for (int d = 0; d < chi2Dim; d++ ) {
				c2Vector[d] = chi2Bounds[2*d + ((int)floor(b/pow(2,d)) % 2)];
			}
			/*
			eneMin = chi2Bounds[0 + ((int)floor(b/1) % 2)]; //0 or 1
			eThMin = chi2Bounds[2 + ((int)floor(b/2) % 2)]; //2 or 3
			theMin = chi2Bounds[4 + ((int)floor(b/4) % 2)]; //4 or 5
			theCut = chi2Bounds[6 + ((int)floor(b/8) % 2)]; //6 or 7
			absMin = chi2Bounds[8 + ((int)floor(b/16) % 2)]; //8 or 9
			dagCut = chi2Bounds[10 + ((int)floor(b/32) % 2)]; //10 or 11
			*/
			// Skip if our scan range is 0
			for (int d=0; d < chi2Dim; d++ ) {
				if((chi2Bounds[2*d + ((int)floor(b/pow(2,d)) % 2)] == c2Vector[d]) && (((int)floor(b/pow(2,d)) % 2)==1)) {
					chi2Vals[b] = chi2Vals[(int)pow(2,chi2Dim)+1];
				} else {
					chi2Vals[b] = 0;
				}
			}
			// Check if we haven't found a good chi2 value yet
			if (chi2Vals[b] == 0) {
				// Load our relevant dagger events on a TH1D
				TH1D* timeDagFitting = new TH1D("timeDagFit", "hit time (dagger)", nBins,0,maxTime);
		
				// begin filling events
				blockCounter = 0;
				for (int i=0; i<length_dag;i++) {
					
					// Load the event buffer
					dag_buffer=events_dag[i];
									
					// Calculate event weighting, and check if the neutron isn't removed
					histoWeight = calcWeighting(dag_buffer, c2Vector);
					upscatter = false;
					if (((double)dag_buffer.time + (50.0 + hold_t)) < randDeathTimes[i] && 
							checkRealN(dag_buffer, c2Vector, randDagProb[i])) {
								
						// Check if the dagger has block hits. If so we have to deal with the block.
						if (dag_buffer.nHits[3] > 0) {
							// Match first event. The next N will sequentially follow, from first to last hit.
							block_buffer = events_blok[blockCounter];
							// faster if we find the original hit, but otherwise scan through
							if ((block_buffer.theta != dag_buffer.theta) || (block_buffer.phi != dag_buffer.phi)) {
								for (int j=blockCounter; j < length_blok; j++) {
									block_buffer = events_blok[j];
									if ((block_buffer.theta == dag_buffer.theta) && (block_buffer.phi == dag_buffer.phi)) {
										blockCounter = j;
										break;
									}
								}	
							}			
							// Now find the scattering probability for each
							for (int j = blockCounter; j < blockCounter + dag_buffer.nHits[3]; j++ ) {	
								// Error handling
								if (j >= length_blok) {
									printf(" WARNING: Have too many block hits associated with this dagger hit ");
									break;
								}
								
								block_buffer = events_blok[j];
								absMin = calcBlockScatter(block_buffer.perp_energy);
								//absMin = 1.0;
								
								// Quit loop on first hit that upscatters, and load block data
								if (randAbsProb[j] < absMin) {
									
									// For chi2 fit we just care about filling timeDagFitting()
									upscatter = true;
									//blockCounter = j+1;
									break;
								}
							}
						}
					
						// fill timeDagFitting of non-upscattered  (real) particles
						if (!upscatter) {
							timeDagFitting->Fill((double)dag_buffer.time+15.59, histoWeight);
							//timeDagFitting->AddBinContent((double)dag_buffer.time+15.59, histoWeight);
						}
					}
					// increment block buffer by number of hits that occured
					blockCounter = blockCounter + dag_buffer.nHits[3];
				}
				/*for (int i = 0; i<length_dag; i++) { 
					hitEvt = 0;
					
					// Load the event buffer
					dag_buffer=events_dag[i];
									
					// Calculate event weighting, and check if the neutron isn't removed
					histoWeight = calcWeighting(dag_buffer, c2Vector);
					if (((double)dag_buffer.time + (50.0 + hold_t)) < randDeathTimes[i] && 
							checkRealN(dag_buffer, c2Vector, randDagProb[i])) {
								
						// Check if the dagger has block hits. If so we have to deal with the block.
						if (dag_buffer.nHits[3] > 0) {
							hitEvt = blockScatterNum(dag_buffer, events_blok, blockCounter, randAbsProb);
						}
						
						// If there is no upscattering, hitEvt = nHits + 1, so fill dagger
						if (hitEvt <= dag_buffer.nHits[3]) {
							timeDagFitting->Fill((double)dag_buffer.time+15.0, histoWeight);
						} 
					}
					// increment blockCounter by number of hits
					blockCounter = blockCounter + (dag_buffer.nHits[3]);
				}
				*/ 
				//chi2Vals[b] = chi2fit(dataHist, timeDagFitting, 60.0, 76.0, 80.0, 130.0);
				chi2Vals[b] = chi2fit(dataHist, timeDagFitting, 59.0, 75.0, 81.0, 140.0);
				//chi2Vals[b] = chi2Vals[b]  + chi2fit(dataHist, timeDagFitting, 80.0, 100.0, 100.0, 130.0);
				/*
				dataHist->GetXaxis()->SetRangeUser(60.0,76.0);
				//dataHist->Scale(1/dataHist->Integral(60.0,76.0));
				timeDagFitting->GetXaxis()->SetRangeUser(60.0,76.0);
				//timeDagFitting->Scale(1/timeDagFitting->Integral(60.0,76.0));
				//chi2Vals[b] = dataHist->Chi2Test(timeDagFitting,"UW CHI2/NDF") * (double)timeDagFitting->Integral(60.0,76.0);
				chi2Vals[b] = timeDagFitting->Chi2Test(dataHist,"WW CHI2/NDF");// * (double)dataHist->Integral(60.0,76.0);
				
				
				//dataHist->Scale(dataHist->Integral(60.0,76.0));
				//timeDagFitting->Scale(timeDagFitting->Integral(60.0,76.0));
				
				// Get ready to comment this out when it doesnt work!
				dataHist->GetXaxis()->SetRangeUser(82.0,140.0);
				//dataHist->Scale(1/dataHist->Integral(86.0,120.0));
				timeDagFitting->GetXaxis()->SetRangeUser(82.0,140.0);
				//timeDagFitting->Scale(1/timeDagFitting->Integral(86.0,120.0));
				//chi2Vals[b] = chi2Vals[b] + dataHist->Chi2Test(timeDagFitting,"UW CHI2/NDF") * (double)timeDagFitting->Integral(82.0,130.0);
				chi2Vals[b] = chi2Vals[b] + timeDagFitting->Chi2Test(dataHist,"WW CHI2/NDF");// * (double)dataHist->Integral(82.0,130.0);
				*/
				if (chi2Vals[b] < 1.0) {
				//	chi2Vals[b] = DBL_MAX;
				}
									
				//dataHist->Scale(dataHist->Integral(86.0,120.0));
				//timeDagFitting->Scale(timeDagFitting->Integral(86.0,120.0));
				/* // Do actual chi2 calc
				timeDagFitting->GetXaxis()->SetRangeUser(58.0,130.0);
				chi2Vals[b] = dataHist->Chi2Test(timeDagFitting,"UW CHI2/NDF");*/
				delete timeDagFitting;
				for (int k = 0; k < chi2Dim; k++) {
					printf("%f ",c2Vector[k]);
				}
				printf("%f \n", chi2Vals[b]);
				
				// check where minimum is
				if (chi2Vals[b] < chi2Min) {
					chi2Min = chi2Vals[b];
					minDex = b;
				}
				//printf("%f %f %f %f %f \n",eneMin,eThMin,theMin,absMin,chi2Vals[b]);
			
				/*
					// Check if we upscatter off the block
					upscatter = false;
					while (dag_buffer.nHits[3] > 0) {
						if (randAbsProb[blockCounter] < absMin) {
							upscatter = true;
							blockCounter = blockCounter + dag_buffer.nHits[3];
							dag_buffer.nHits[3] = 0;
						} else {
							upscatter = false;
							blockCounter = blockCounter + 1;
							dag_buffer.nHits[3] = dag_buffer.nHits[3] - 1;
						}
					}
					
					// check if we get absorbed by the bad part of the dagger
				    if(dag_buffer.pos[0] > 0) {
				        zeta = 0.5 - sqrt(dag_buffer.pos[0]*dag_buffer.pos[0] + pow(fabs(dag_buffer.pos[2] - dag_buffer.z_offset) - 1.0, 2));
				    }
				    else {
						zeta = 1.0 - sqrt(dag_buffer.pos[0]*dag_buffer.pos[0] + pow(fabs(dag_buffer.pos[2] - dag_buffer.z_offset) - 0.5, 2));
				    }
				    boronProb = zeta > dagCut ? 1 : zeta/dagCut;
				    // Calculate weighting
					if (!upscatter) {
						angleTheta = dag_buffer.theta;
						energyCm = dag_buffer.energy*JTONEV;
						if (thetaOn) {
							histoWeight = (pow(cos(angleTheta),theMin)*sin(angleTheta))/(cos(angleTheta)*sin(angleTheta)) * pow((energyCm),eneMin)/energyCm;
						} else {
							histoWeight = (pow(cos(angleTheta),theMin)*sin(angleTheta)) * pow((energyCm),eneMin)/energyCm;
						}
						if ((std::isfinite(histoWeight)) && (energyCm > eThMin) && (angleTheta < theCut) && (randDagProb[i] < boronProb)) {
							timeDagFitting->Fill(dag_buffer.time+15.0, histoWeight);
						}
					}
				}*/
				/*timeDagFitting->GetXaxis()->SetRangeUser(42.0,220.0);
				chi2Vals[b] = dataHist->Chi2Test(timeDagFitting,"UW CHI2/NDF");
				delete timeDagFitting;
				for (int k = 0; k < chi2Dim; k++) {
					printf("%f ",c2Guess[k]);
				}
				printf("%f \n", chi2Vals[b]);
				//printf("%f %f %f %f %f \n",eneMin,eThMin,theMin,absMin,chi2Vals[b]);
				// check where minimum is
				if (chi2Vals[b] < chi2Min) {
					chi2Min = chi2Vals[b];
					minDex = b;
				}*/
			}
		}
		// Divide step size by 2
		for (int dim=0; dim<chi2Dim; dim++) {
			c2Step[dim] = c2Step[dim]/2;
			// Move center (if necessary)
			if (minDex != (pow(2,chi2Dim)+1)) {
				c2Guess[dim] = (c2Guess[dim] -(1 - 2*((int)floor(minDex/pow(2,dim)) % 2)) * c2Step[dim]);
			}
		}
	}
	
	printf("Min Chi^2:");
	for (int dim = 0; dim<chi2Dim; dim++) {
		printf("%f ", c2Guess[dim]);
	}
	printf("%f\n\n", chi2Min);
	printf("----------------------------------------------------------\n\n");
	dataHist->GetXaxis()->SetRangeUser(0.0,140.0);
	//printf("Min Chi^2: %f %f %f %f %f %f %f\n",eneMin,eThMin,theMin,theCut,absMin,dagCut,chi2Min);
	/*------------------------------------------------------------*/
	/* ----------------MAP CHI^2 TANGENT SURFACES---------------- */
	/*------------------------------------------------------------*/
	if (chiMap) {
		printf("Generating tangent surface plots! \n");			
		// Increase scan size (increase our smallest step by dimensions)
		for (int dim=0; dim<chi2Dim; dim++) {
			c2Step[dim] = pow(2,numSteps-2)*c2Step[dim];
			//c2Step[dim] = 2.0 * c2Step[dim];
			c2Mapper[dim] = c2Guess[dim];
		}
		
		int stepx, stepy;
		int scanRangex, scanRangey = numSteps;
		double c2Temp, c2Val = 0.0;
		
		
		// Scan across our dimensions to produce our chi2 plots
		for (int dimx=0; dimx<chi2Dim; dimx++) {
			for (int dimy=0; dimy<dimx; dimy++) { // Scanning for only a given quadrant
				if (c2Step[dimx]!=0.0 && c2Step[dimy]!=0.0) { //Ignore if there's no step!
					
					// Load our file for scanning
					sprintf(c2Filename,"chi2Fit%d%d.csv",dimx,dimy);
					printf("Opening File %s\n\n",c2Filename);
					chi2OutFile = fopen(c2Filename,"w");
					
					// Initialize our step range
					scanRangex = numSteps;
					stepx = -scanRangex;
											
					// Loop across our steps
					while (stepx < scanRangex) { 
						c2Temp = chi2Min;
						// Re-initialize in y
						scanRangey = numSteps;	
						stepy = -scanRangey;
						while (stepy < scanRangey) {
						
						//for (int stepx=-scanRange; stepx<=scanRange+1; stepx++) {
						//for (int stepy=-scanRange; stepy<=scanRange+1; stepy++) {
							
							// Figure out what our positions are.
							for (int dim=0;dim<chi2Dim;dim++) {
								c2Vector[dim] = c2Guess[dim];
							}
							c2Vector[dimx] = c2Guess[dimx] + (double)(stepx/((double)numSteps/2.0))*c2Step[dimx];
							c2Vector[dimy] = c2Guess[dimy] + (double)(stepy/((double)numSteps/2.0))*c2Step[dimy];
							// if we get a negative c2Vector, move up.
							while (c2Vector[dimx] < 0.0) {
								stepx = stepx + 1;
								scanRangex = scanRangex +1;
								c2Vector[dimx] = c2Guess[dimx] + (double)(stepx/((double)numSteps/2.0))*c2Step[dimx];
							}
							while (c2Vector[dimy] < 0.0) {
								stepy = stepy + 1;
								scanRangey = scanRangey +1;
								c2Vector[dimy] = c2Guess[dimy] + (double)(stepy/((double)numSteps/2.0))*c2Step[dimy];
							}
							
							// Putting in a fast-scanning thing to get rid of non-visible angle cuts
							/*if (c2Vector[3] > 1.571 && (dimx == 3 || dimy == 3)) {
								stepy = scanRangey;
								stepx = scanRangex;
							}*/
							
							
							// Fill only one histogram at a time
							TH1D* timeDagFitting = new TH1D("timeDagFit", "hit time (dagger)", nBins,0,maxTime);
							
							blockCounter = 0;
							for (int i=0; i<length_dag;i++) {
								
								// Load the event buffer
								dag_buffer=events_dag[i];
												
								// Calculate event weighting, and check if the neutron isn't removed
								histoWeight = calcWeighting(dag_buffer, c2Vector);
								upscatter = false;
								if (((double)dag_buffer.time + (50.0 + hold_t)) < randDeathTimes[i] && 
										checkRealN(dag_buffer, c2Vector, randDagProb[i])) {
											
									// Check if the dagger has block hits. If so we have to deal with the block.
									if (dag_buffer.nHits[3] > 0) {
										// Match first event. The next N will sequentially follow, from first to last hit.
										block_buffer = events_blok[blockCounter];
										// faster if we find the original hit, but otherwise scan through
										if ((block_buffer.theta != dag_buffer.theta) || (block_buffer.phi != dag_buffer.phi)) {
											for (int j=blockCounter; j < length_blok; j++) {
												block_buffer = events_blok[j];
												if ((block_buffer.theta == dag_buffer.theta) && (block_buffer.phi == dag_buffer.phi)) {
													blockCounter = j;
													break;
												}
											}	
										}			
										// Now find the scattering probability for each
										for (int j = blockCounter; j < blockCounter + dag_buffer.nHits[3]; j++ ) {	
											// Error handling
											if (j >= length_blok) {
												printf(" WARNING: Have too many block hits associated with this dagger hit ");
												break;
											}
											
											block_buffer = events_blok[j];
											absMin = calcBlockScatter(block_buffer.perp_energy);
											//absMin = 1.0;
											
											// Quit loop on first hit that upscatters, and load block data
											if (randAbsProb[j] < absMin) {
												
												// For chi2 fit we just care about filling timeDagFitting()
												upscatter = true;
												//blockCounter = j+1;
												break;
											}
										}
									}
								
									// fill timeDagFitting of non-upscattered  (real) particles
									if (!upscatter) {
										timeDagFitting->Fill((double)dag_buffer.time+16.0, histoWeight);
									}
								}
								// increment block buffer by number of hits that occured
								blockCounter = blockCounter + dag_buffer.nHits[3];
							}
							
							// Calculate chi2 values for a region around
							c2Val = chi2fit(dataHist, timeDagFitting, 59.0, 75.0, 81.0, 140.0);
							
							delete timeDagFitting;
							for (int dim = 0; dim<chi2Dim; dim++) {
								printf("%f ", c2Vector[dim]);
							}
							printf("%f \n", c2Val);
						
							fprintf(chi2OutFile, "%lg %lg %lg \n",c2Vector[dimx],c2Vector[dimy],c2Val);
														
							// if we're still descending want to move another step
							if (c2Val < c2Temp) {
								scanRangex = scanRangex + 1;
								scanRangey = scanRangey + 1;
								// do stuff if we've found a new minimum
								if (c2Val < chi2Min ) {
									chi2Min = c2Val;
									printf("Shifting Minimum Chi2 point!\n");
									for (int dim = 0; dim<chi2Dim; dim++) {
										c2Mapper[dim] = c2Vector[dim];
									}
								}
							// Break if we're not changing anything 
							} else if (fabs(c2Val-c2Temp) < 1e-9) {
								printf("No change in chi2 value! Leaving fitting loop!\n");
								//stepx = scanRangex;
								//stepy = scanRangey;
							}
													
							// catch in case we somehow got a infinite loop
							if (stepy > 4*numSteps+2) {
								printf("ERROR: Found too many steps beyond expected!\n");
								printf("       This means we found a local minimum earlier!\n");
								printf("       Ending function!\n\n");
								stepy = scanRangey;
							}
							
							// Increment
							stepy = stepy+1;
							c2Temp = c2Val;
						}
						if (stepx > 4*numSteps+2) {
							printf("ERROR: Found too many steps beyond expected!\n");
							printf("       This means we found a local minimum earlier!\n");
							printf("       Ending function!\n\n");
							stepx = scanRangex;
							
						}
						
						stepx = stepx+1;
					}
					
					fclose(chi2OutFile);	
					
				}
			}
		}
	}
	
	return c2Fitted;
	
}


int LifetimeShift2() {
  //  ROOT::EnableThreadSafety();
    initxorshift();
    //TApplication theApp("App",0,0);
    
    /*----------------------------------------------------------------*/
    /* !!!DECLARATION OF VARIABLES!!! */
    /*----------------------------------------------------------------*/
    /* Buffer declarations */
	const size_t buff_len = 4 + 11*8 + 4*4 + 4; // Change if MC data changes
	const size_t blok_len = (4 + 10*8 + 4*4 + 4);
	
	// Load events structures
    std::vector<dag_evt> events_dag;
    dag_evt event_dag, dag_buffer;
	std::vector<block_evt> events_blok;
    block_evt event_blok, block_buffer;
		
    // IO variables
    int test; // test debug parameter
    int length_dag, length_blok;
    
    if(argc <= 3) {
        printf("Error! Usage: ./block_lifetime_fit folder t1 t2 t... \n");
        return 0;
    }
    int numArgs = argc - 2; // number of arguments -- find now 
    double hold_t[numArgs];
	
    // Analysis and weighting variables
    bool upscatter = false;
    int blockCounter, length_od, length_ob = 0;
    double angleTheta, zeta, boronProb;
    double neutronYield[numArgs] = {0.0};
    double neutronError[numArgs] = {0.0};
    double totalYield[numArgs] = {0.0};
    double histoWeight, eneMin, eThMin, theMin, theCut, absMin, dagCut;
    double groupYield[3][numArgs] = {0.0};
    double groupError[3][numArgs] = {0.0};
    double iniState[3];
    double binBreak[2];
    
    // Yield finding initialization 
	binBreak[0] = 40.0; // Change these numbers to change the position of bin dividers (only 3dip for now)
	binBreak[1] = 60.0;
		
	/* Chi2 fitting dimensions table: 
	 * 1 = Energy (power)
	 * 2 = Energy (threshold)
	 * 3 = Angle  (power) 
	 * 4 = Angle  (cutoff)
	 * 5 = Upscatter (prob)
	 * 6 = B10 Bad Region (m) */
	int chi2Dim = 6;
	//double c2Guess[chi2Dim] = {2.15, 0.088, 1.10, 2.0, 1.0, 0.0}; //"Old," not as good numbers
	//double c2Guess[chi2Dim] = {2.0, 0.0, 1.0, 2.0, 1.0, 0.0}; // No weighting or truncation
	//double c2Guess[chi2Dim] = {2.30448, 7.6071, 1.264079, 2.0, 1.0, 0.0160763}; //"Nathan's suggested" numbers
	
	//double c2Guess[chi2Dim] = {3.504480, 0.048571, 2.564079, 1.0, 1.0, 0.052076}; //I found these numbers!	 -- They're super unphysical!
    char dagFile[48];
    char blokFile[48];
        
    /*----------------------------------------------------------------*/
    /* !!!BOOLEAN LOGIC to quickly change program analysis!!! */
    /*----------------------------------------------------------------*/
    bool debugMode = false; // Turns on debug mode
    bool quickLoad =  true; // Only loads the first 10 dagger and block events
    bool thetaOn   =  true; // Modify thetaOn to change plotting E and theta
    bool realBlock =  true; // Turns on/off a real aluminum block spectrum
    
	//printf("%d\n",numArgs);
	eneMin = c2Guess[0];
	eThMin = c2Guess[1];
	theMin = c2Guess[2];
	theCut = c2Guess[3];
	absMin = c2Guess[4];
	dagCut = c2Guess[5];
	
	// This loop opens files and calculates the yields
	for (int ii=0; ii<numArgs; ii++) {
		
		char* buf;
		char* buf_blok;
		buf = new char[buff_len];
		buf_blok = new char[blok_len];

		sprintf(dagFile,"%sdag.out.%ss", argv[1],argv[ii+2],'c');

		hold_t[ii]= std::stod(argv[ii+2]);
		std::ifstream binfile(dagFile, std::ios::in | std::ios::binary);
		
		if(!binfile.is_open()) {
			printf("Error! Could not open file %s\n", dagFile);
			return 1;
		}

		int test=1;  		
		while(!binfile.eof()) {
			binfile.read(buf, buff_len);
			if(binfile.eof()) {
				break;
			}

			// Load dagger events
       		event_dag.time = *((double *)(&buf[0] + sizeof(unsigned int)));
	        event_dag.energy = *((double *)(&buf[0] + sizeof(unsigned int)+sizeof(double)));
	        event_dag.perp_energy = *((double *)(&buf[0] + sizeof(unsigned int)+2*sizeof(double)));
	        for(int jj = 0; jj < 3; jj++) {
				event_dag.pos[jj] = *((double *)(&buf[0] + sizeof(unsigned int)+(3 + jj)*sizeof(double)));
			}
			event_dag.z_offset = *((double *)(&buf[0] + sizeof(unsigned int) + 6*sizeof(double)));
			for(int jj = 0; jj < 4; jj++) {
				event_dag.nHits[jj] = *((int *)(&buf[0] + sizeof(unsigned int) + 7*sizeof(double) + jj*sizeof(int)));
			}
			
			for(int jj = 0; jj < 3; jj++) {
				iniState[jj] = *((double *)(&buf[0] + sizeof(unsigned int) + (7+jj)*sizeof(double) + 4*sizeof(int))); 
			}
			event_dag.theta = acos(iniState[2]/sqrt(iniState[0]*iniState[0] + iniState[1]*iniState[1] + iniState[2]*iniState[2]));
			event_dag.phi = acos(iniState[5]/(event_dag.theta));
			
			if (debugMode) {
				if (test < 10) { 
					printf( " %e %e %e \n ", event_dag.time, event_dag.energy, event_dag.perp_energy);
					printf( " %e %e %e %e \n ", event_dag.pos[0], event_dag.pos[1], event_dag.pos[2], event_dag.z_offset);
					printf( " %u %u %u %u \n ", event_dag.nHits[0], event_dag.nHits[1], event_dag.nHits[2], event_dag.nHits[3]);
					printf( " %e %e %e \n\n", iniState[0], iniState[1], iniState[2]); 
					test = test+1;
				} else if (quickLoad) {
					break;
				}
	        }
	        events_dag.push_back(event_dag);
	    }
	    binfile.close();
	
	    length_dag = events_dag.size() - length_od;
	    printf("Read %u (dagger) Events in file %s!\n", length_dag, dagFile);
	   
		sprintf(blokFile,"%sblok.out.%ss", argv[1],argv[ii+2],'c');	    
	    std::ifstream binfile2(blokFile, std::ios::in | std::ios::binary);
		if(!binfile2.is_open()) {
			printf("Error! Could not open file %s\n", blokFile);
			return 1;
		}	
	    test = 1;
	       	    
	    while(!binfile2.eof()) {
	        binfile2.read(buf_blok, blok_len);
	        if(binfile2.eof()) {
	            break;
	        }
	        event_blok.time = *((double *)(&buf_blok[0] +sizeof(unsigned int)));
	        event_blok.energy = *((double *)(&buf_blok[0] + +sizeof(unsigned int) + sizeof(double)));
	        for(int jj = 0; jj < 4; jj++) {
				event_blok.nHits[jj] = *((int *)(&buf_blok[0] + sizeof(unsigned int) + 2*sizeof(double) + jj*sizeof(int)));
			}
			for(int jj = 0; jj < 3; jj++) {
				event_blok.position[jj] = *((double *)(&buf_blok[0] + sizeof(unsigned int) + (2 + jj)*sizeof(double) + 4*sizeof(int)));
			}
			event_blok.perp_energy = *((double *)(&buf_blok[0] + sizeof(unsigned int) + 5*sizeof(double) + 4*sizeof(int)));
			for(int jj = 0; jj < 3; jj++) {
				iniState[jj] = *((double *)(&buf_blok[0] + sizeof(unsigned int) + (6+jj)*sizeof(double) + 4*sizeof(int))); 
			}
	        event_blok.theta = acos(iniState[2]/sqrt(iniState[0]*iniState[0] + iniState[1]*iniState[1] + iniState[2]*iniState[2]));
	        event_blok.phi = acos(iniState[5]/(event_blok.theta));
	
			if (debugMode) {
				if (test < 10) { 
					printf( " %e %e \n ", event_blok.time, event_blok.energy);
					printf( " %u %u %u %u \n ", event_blok.nHits[0], event_blok.nHits[1], event_blok.nHits[2], event_blok.nHits[3]);
					printf( " %e %e %e \n ", event_blok.position[0], event_blok.position[1], event_blok.position[2]);
					printf( " %e %e \n\n ", event_blok.perp_energy, event_blok.theta);
					printf( " %e %e %e \n\n", iniState[0], iniState[1], iniState[2]); 
					test = test+1;
				} else if (quickLoad) {
					break;
				}
	        }			
	        events_blok.push_back(event_blok);
	    }
	    binfile2.close();
	   	    
	    length_blok = events_blok.size() - length_ob;
	    printf("Read %u (block) Events in file %s!\n", length_blok, blokFile);
	    
	    // Create random numbers for weighting absorption spectrum (dagger)
	    double* randDeathTimes = new double[length_dag];
		double* randAbsProb = new double[length_blok];
	    double* randDagProb = new double[length_dag];
	    for(unsigned long i = 0; i < length_dag; i++) {
	        
	        randDeathTimes[i] = 99999999999;
	    }
	    for(unsigned long i = 0; i < length_blok; i++) {
			randAbsProb[i] = 1.0 - nextU01();
			if (randAbsProb[i] < 0.0 || randAbsProb[i] > 1.0) {  // check and make sure random number generator doesn't crash
				i = i-1;
			}
		}
		for(unsigned long i = 0; i  < length_dag; i++) {
			randDagProb[i] = nextU01();
		}
		
		// Calculate Neutron weighting
		printf("Calculating Weighting!\n");
		blockCounter = 0;
		for (int jj=0; jj<length_dag; jj++) {
			dag_buffer=events_dag[jj];
			angleTheta = dag_buffer.theta;
			
			// Calculate event weighting. Weighting based on the DAGGER energy, since the block will have diff. hits. Might want to change to more reasonable number (in energy)
			/*if (thetaOn) { 
				histoWeight = (pow(cos(angleTheta),theMin)*sin(angleTheta))/(cos(angleTheta)*sin(angleTheta)) * (pow((dag_buffer.energy),eneMin)/dag_buffer.energy);
			} else {
				histoWeight = (pow(cos(angleTheta),theMin)*sin(angleTheta)) * (pow((dag_buffer.energy),eneMin)/dag_buffer.energy);
			}*/
			histoWeight = 1.0;
			// check if we get absorbed by the bad part of the dagger
		    if(dag_buffer.pos[0] > 0) {
		        zeta = 0.5 - sqrt(dag_buffer.pos[0]*dag_buffer.pos[0] + pow(fabs(dag_buffer.pos[2] - dag_buffer.z_offset) - 1.0, 2));
		    }
		    else {
				zeta = 1.0 - sqrt(dag_buffer.pos[0]*dag_buffer.pos[0] + pow(fabs(dag_buffer.pos[2] - dag_buffer.z_offset) - 0.5, 2));
		    }
		    boronProb = zeta > dagCut ? 1 : zeta/dagCut;
		    
		    // load our dagger things
			if ((((double)dag_buffer.time + (50.0 + hold_t[ii])) < randDeathTimes[jj]) && (dag_buffer.energy*JTONEV > eThMin) && (angleTheta < theCut) && dag_buffer.time > 0.0 && (randDagProb[jj] < boronProb)){
				upscatter = false;
				// Check if the dagger event actually hit the block. Counter is now index of FIRST block hit
				if (dag_buffer.nHits[3] > 0){
					// Match first event. The next N will sequentially follow, from first to last hit.
					block_buffer = events_blok[blockCounter];
					if ((block_buffer.theta != dag_buffer.theta) || (block_buffer.phi != dag_buffer.phi)) {
						for (int kk=0; kk<length_blok; kk++) {
							block_buffer = events_blok[kk];
							if ((block_buffer.theta == dag_buffer.theta) && (block_buffer.phi == dag_buffer.phi)) {
								blockCounter = kk;
								break;
							}
						}	
					}
					for (int kk=blockCounter; kk<(blockCounter+dag_buffer.nHits[3]); kk++) {
						// Error handling
						if (kk >= length_blok) {
							printf(" WARNING: Have too many block hits associated with this dagger hit ");
							break;
						}
						block_buffer = events_blok[kk];
						if (realBlock) {
							absMin = calcBlockScatter(block_buffer.perp_energy);
							printf("absMin = %f \n", absMin);
						} else {
							absMin = 0.0001;
						}
						
						// Quit loop on first hit that upscatters, and load block data
						if (randAbsProb[kk] < absMin) {
							// If we've upscattered, leave the FOR loop, go to calculate yields
							upscatter = true;
							blockCounter = kk+1;
							break;
						}
					}
				}
			    totalYield[ii] = totalYield[ii] + histoWeight;
			    //neutronError[ii] = neutronError[ii] + histoWeight*histoWeight;
				// Fill dagger weights/ neutron yields
				if (!upscatter) {					
					if (std::isfinite(histoWeight)) {
						neutronYield[ii] = neutronYield[ii] + histoWeight;
						neutronError[ii] = neutronError[ii] + histoWeight*histoWeight;
						//totalYield[ii] = totalYield[ii] + histoWeight;
						if (dag_buffer.time < binBreak[0]) { 
							groupYield[0][ii] = groupYield[0][ii] + histoWeight;
							groupError[0][ii] = groupError[0][ii] + histoWeight*histoWeight;
						} else if (dag_buffer.time >= binBreak[0] && dag_buffer.time < binBreak[1]) {
							groupYield[1][ii] = groupYield[1][ii] + histoWeight;
							groupError[1][ii] = groupError[1][ii] + histoWeight*histoWeight;
						} else if (dag_buffer.time >=binBreak[1]) {
							groupYield[2][ii] = groupYield[2][ii] + histoWeight;
							groupError[2][ii] = groupError[2][ii] + histoWeight*histoWeight;
						}				
					}
				} //else {
					//totalYield[ii] = totalYield[ii] + histoWeight;
				//}
			}
		}
				
				
		//* normalize the errors */
		neutronError[ii] = sqrt(neutronError[ii]);
		for (int jj=0; jj<3; jj++) { 
			groupError[jj][ii] = sqrt(groupError[jj][ii]);
		}
		printf("Done Weighting!\n");
		
		// Clear memory cache
		events_dag.erase(events_dag.begin(),events_dag.end());
		events_blok.erase(events_blok.begin(),events_blok.end());
		delete[] randDeathTimes;
		delete[] randAbsProb;
		delete[] randDagProb;
		delete[] buf;
		delete[] buf_blok;

	}
		
	printf("Normalizing group weights!\n");
	double normalization = totalYield[0];
	//double normalization = (totalYield[0] - neutronYield[0]) / totalYield[0];
	for(int ii =0 ;ii<numArgs;ii++){
		neutronYield[ii] = /*(totalYield[ii]-*/(neutronYield[ii])/totalYield[ii]/normalization;
		neutronError[ii] = neutronError[ii]/normalization;
		for(int jj=0; jj<3; jj++) {
			groupYield[jj][ii] = groupYield[jj][ii]/totalYield[ii]/normalization;
			groupError[jj][ii] = groupError[jj][ii]/normalization;
		}
		
		
		printf("\n%e, %e, %e, %e\n", neutronYield[ii], groupYield[0][ii],groupYield[1][ii],groupYield[2][ii]);
		printf("%e, %e, %e, %e\n", neutronError[ii], groupError[0][ii],groupError[1][ii],groupError[2][ii]);
	}
			
	if (numArgs > 2) {
		printf("Saving output file! \n");
		
		FILE* outfile;
		outfile = fopen("block_yields.csv","w");
		for (int ii = 0; ii<numArgs; ii++){
			fprintf(outfile,"%f, ",hold_t[ii]);
			fprintf(outfile, "%f, %f, %f, %f, ",     neutronYield[ii], neutronError[ii], groupYield[0][ii], groupError[0][ii]);
			fprintf(outfile, "%f, %f, %f, %f, \n", groupYield[1][ii], groupError[1][ii], groupYield[2][ii] ,groupError[2][ii]);
		}
		fclose(outfile);
	} else {
		     
	    double totalLifetime;
	    double blockLifetime;
	    double totalError;
	    
	    double groupLife[3];
	    double dT = hold_t[0] - hold_t[1];
	    totalLifetime = -dT/log(1/(exp(dT/877.7)*neutronYield[1]/neutronYield[0]));
		blockLifetime = totalLifetime - 877.7;
		//blockLifetime = 1/(1/totalLifetime - 1/877.7);
		//totalError = -dT / (neutronError[1]*neutronYield[1]/neutronYield[0]) / pow(log(1/(exp(-dT/877.7)*(neutronError[1]*neutronYield[1]/neutronYield[0]))),2);
		totalError = sqrt(((neutronError[1]/neutronYield[1])*(neutronError[1]/neutronYield[1]) + (neutronError[0]/neutronYield[0])*(neutronError[0]/neutronYield[0]))/totalLifetime)/(totalLifetime);//(log(neutronYield[1]/neutronYield[0])*log(neutronYield[1]/neutronYield[0]));
		
		for (int ii=0;ii<3;ii++){
			groupLife[ii] = -dT/log(1/(exp(dT/877.7)*groupYield[ii][1]/groupYield[ii][0]));
		}
		
		printf("The total lifetime is: %f +- %f !\n", totalLifetime, totalError);
		printf("This corresponds to a block shift of %f !\n", blockLifetime);
		printf("Individual group lifetimes:\n   %f (group 1) \n   %f (group 2) \n   %f (group 3)\n", groupLife[0], groupLife[1], groupLife[2]);
	}
	
    return 0;
}




void fillBlockHistos(){
if (realRefl) { 
		printf("Using realistic Al Block reflection data. Filling daggger and block histograms! \n\n");
		
		// Look at each dagger event separately
		for (int i=0; i<length_dag;i++) {
			
			upscatter = false;
			hitEvt = 0;
			// Load the event buffer
			dag_buffer=events_dag[i];
			if (pseOn) {
				dag_buffer.time = dag_buffer.time - (180.0 + hold_t);
			}
			
			// Calculate event weighting, and check if the neutron isn't removed
			histoWeight = calcWeighting(dag_buffer, c2Guess);
			//printf("%f\n", histoWeight);
			if (((double)dag_buffer.time + (50.0 + hold_t)) < randDeathTimes[i] && 
					checkRealN(dag_buffer, c2Guess, randDagProb[i])) {
						
				// Check if the dagger has block hits. If so we have to deal with the block.
				if (dag_buffer.nHits[3] > 0) {
					// Match first event. The next N will sequentially follow, from first to last hit.
					block_buffer = events_blok[blockCounter];
					// faster if we find the original hit, but otherwise scan through
					if ((block_buffer.theta != dag_buffer.theta) || (block_buffer.phi != dag_buffer.phi)) {
						for (int j=blockCounter; j < length_blok; j++) {
							block_buffer = events_blok[j];
							if ((block_buffer.theta == dag_buffer.theta) && (block_buffer.phi == dag_buffer.phi)) {
								blockCounter = j;
								break;
							}
						}	
					}			
					// Now find the scattering probability for each
					for (int j = blockCounter; j < blockCounter + dag_buffer.nHits[3]; j++ ) {	
						// Error handling
						if (j >= length_blok) {
							printf(" WARNING: Have too many block hits associated with this dagger hit ");
							break;
						}
						
						block_buffer = events_blok[j];
						absMin = calcBlockScatter(block_buffer.perp_energy);
						//absMin = 1.0;
						
						// Quit loop on first hit that upscatters, and load block data
						if (randAbsProb[j] < absMin) {							
							upscatter = true;
							
							// now we want to load block hits
							if (blockMapTot) {
								blockHistRoot->Fill(block_buffer.time, (block_buffer.energy)*JTONEV,histoWeight);
								timeBlockHist->Fill((double)block_buffer.time,histoWeight);
								eneBlockHist->Fill((block_buffer.energy)*JTONEV,histoWeight);
							}
							field_block = fieldstrength(block_buffer.position);
							if (blockMapSlice || bFieldSlice) {
								if ((block_buffer.energy)*JTONEV <= threshold) {
										tBlockHist[0]->Fill((double)block_buffer.time,histoWeight);				
										fieldBHist[0]->Fill(field_block,histoWeight);
								} else {
									for (int j=1; j < numHistos; j++) {
										if (((block_buffer.energy)*JTONEV > (((double)j-1.0)*spacing + threshold)) && ((block_buffer.energy)*JTONEV <= (((double)j)*spacing + threshold))) {
											tBlockHist[j]->Fill((double)block_buffer.time,histoWeight);
											fieldBHist[j]->Fill(field_block,histoWeight);
										}
									}
								}			
							}
							if (blockMap) {
								hitPos->SetPoint(i,block_buffer.position[0],block_buffer.position[1],block_buffer.position[2]);
								if (block_buffer.position[2] > maxCoords[2]) {
									for (int q = 0; q<3; q++) {
										maxCoords[q] = block_buffer.position[q];
									}
								}
								if (block_buffer.position[2] < minCoords[2]) {
									for (int q = 0; q<3; q++) {
										minCoords[q] = block_buffer.position[q];
									}
								}
								hitB->SetPoint(i,block_buffer.position[0],block_buffer.position[1],field_block);
							}
							break;
						} else {
							doubleCounts = doubleCounts + 1;
						}
					}
				}
				
				// fill dagger histograms depending if we have upscattered
				if (sanityCheck) {
					eneDagHist->Fill((dag_buffer.energy)*JTONEV,histoWeight);
				}
				
				if (!upscatter) {
					timeDagHist->Fill((double)dag_buffer.time+15.59, histoWeight);
					if (sanityCheck) {
						thetaDagRoot->Fill((double(dag_buffer.theta)),histoWeight);
					}
				}
			} else {
				brokenCounter = brokenCounter + 1;
			}
			
			// Increment block counter by number of hits.
			blockCounter = blockCounter + (dag_buffer.nHits[3]);
		}

		printf( " Calculating yields! \n");
		if (timeDagHist->GetEntries() > 0) {
			neutronYield = timeDagHist->IntegralAndError(0,nBins,yieldError,"");
			lossYield = timeBlockHist->IntegralAndError(0,nBins,lossError,"");
			// Normalize losses to total neutron yields
			lossYield = lossYield / thetaDagRoot->Integral();
			lossError = lossError / thetaDagRoot->Integral();
			for (int i=0; i<sizeof(groupYield)/sizeof(*groupYield); i++) {
				groupYield[i] = timeDagHist->IntegralAndError(timeDagHist->GetXaxis()->FindBin(gBins[2*i]),timeDagHist->GetXaxis()->FindBin(gBins[2*i+1]),groupError[i],"");
			}
		} else {
			printf("No neutrons on dagger! Yields set to 0! \n");
			neutronYield = 0.0;
			for (int i=0; i<sizeof(groupYield)/sizeof(*groupYield); i++) {
				groupYield[i] = 0.0;
			}
		}
		
		
	} else { 
		printf("Filling Dagger histograms! \n");
		for (int i=0; i<length_dag; i++) {
			dag_buffer=events_dag[i];
			
			if (pseOn) {
				dag_buffer.time = dag_buffer.time - (180.0 + hold_t);
			}
			
			// Calculate out the weighting for theta
			angleTheta = dag_buffer.theta;
			if (thetaOn) { 
				histoWeight = (pow(cos(angleTheta),theMin)*sin(angleTheta))/(cos(angleTheta)*sin(angleTheta)) * pow((dag_buffer.energy),eneMin)/dag_buffer.energy;
			} else {
				histoWeight = (pow(cos(angleTheta),theMin)*sin(angleTheta)) * pow((dag_buffer.energy),eneMin)/dag_buffer.energy;
			}
			// check if we get absorbed by the bad part of the dagger
		    if(dag_buffer.pos[0] > 0) {
		        zeta = 0.5 - sqrt(dag_buffer.pos[0]*dag_buffer.pos[0] + pow(fabs(dag_buffer.pos[2] - dag_buffer.z_offset) - 1.0, 2));
		    }
		    else {
				zeta = 1.0 - sqrt(dag_buffer.pos[0]*dag_buffer.pos[0] + pow(fabs(dag_buffer.pos[2] - dag_buffer.z_offset) - 0.5, 2));
		    }
		    boronProb = zeta > dagCut ? 1 : zeta/dagCut;
			// Check if we're (a) loading a living neutron and (b) above threshold
			if (((double)dag_buffer.time + (50.0 + hold_t)) < randDeathTimes[i] && (dag_buffer.energy*JTONEV > eThMin) && (angleTheta < theCut) && dag_buffer.time > 0.0 && (randDagProb[i] < boronProb)){
				upscatter = false;
				// Check if we've been killed by the block
				while (dag_buffer.nHits[3] > 0) {
					if (randAbsProb[blockCounter] < absMin) {
						upscatter = true;
						blockCounter = blockCounter + dag_buffer.nHits[3];
						dag_buffer.nHits[3] = 0;
					} else {
						upscatter = false;
						blockCounter = blockCounter + 1;
						dag_buffer.nHits[3] = dag_buffer.nHits[3] - 1;
					}
				}
			    /* Calculate weighting */
				/* Fill Necessary histograms. */
				if (!upscatter) { 
					timeDagHist->Fill((double)dag_buffer.time, histoWeight);
					if (sanityCheck) {
						thetaDagRoot->Fill((double(dag_buffer.theta)),histoWeight);
					}
					
				}
				if (sanityCheck) {
					eneDagHist->Fill((dag_buffer.energy)*JTONEV,histoWeight);
				}
			}	
		}
		// Calculate Yields based on histogram integrals
		neutronYield = timeDagHist->IntegralAndError(0,nBins,yieldError,"");
		
		for (int i=0; i<sizeof(groupYield)/sizeof(*groupYield); i++) {
			groupYield[i] = timeDagHist->IntegralAndError(timeDagHist->GetXaxis()->FindBin(gBins[2*i]),timeDagHist->GetXaxis()->FindBin(gBins[2*i+1]),groupError[i],"");
		}
		
		printf("Filling Block histograms! \n");
		blockCounter = 0;
		// Fill Necessary Block Data
		for (int i=0; i<length_blok; i++) {
			block_buffer=events_blok[i];
			deathTime = randDeathTimes[i+length_dag];
			angleTheta = block_buffer.theta;
			upscatter = false;
			// Check Block scattering
			while (block_buffer.nHits[3] > 1) {
				if (randAbsProb[blockCounter] < absMin) {
					if (randAbsProb[blockCounter] < 0.0 || randAbsProb[blockCounter] > 1.0) { // fix buffer overflows or something?
						upscatter = true;
						block_buffer.nHits[3] = 1;
						blockCounter = 0;
					//	brokenCounter = brokenCounter + 1;
					} else {
						upscatter = true;
						blockCounter = blockCounter + block_buffer.nHits[3];
						block_buffer.nHits[3] = 1;
					}
					
				} else {
					upscatter = false;
					blockCounter = blockCounter + 1;
					block_buffer.nHits[3] = block_buffer.nHits[3] - 1;
				}
			}
			// Load histograms that don't scatter and don't die
			if ((!upscatter) && ((double)block_buffer.time < deathTime) && (block_buffer.position[2] < -1.0) && (block_buffer.energy*JTONEV> eThMin) && (angleTheta < theCut)) {
				
				if (thetaOn) { 
					histoWeight = (pow(cos(angleTheta),theMin)*sin(angleTheta))/(cos(angleTheta)*sin(angleTheta)) * pow((block_buffer.energy),eneMin)/block_buffer.energy;
				} else {
					histoWeight = (pow(cos(angleTheta),theMin)*sin(angleTheta))*pow((block_buffer.energy),eneMin)/block_buffer.energy;
				}
				
				if (blockMapTot) {
					blockHistRoot->Fill(block_buffer.time, (block_buffer.energy)*JTONEV,histoWeight);
					timeBlockHist->Fill((double)block_buffer.time,histoWeight);
					eneBlockHist->Fill((block_buffer.energy)*JTONEV,histoWeight);
				}
				field_block = fieldstrength(block_buffer.position);
				if (blockMapSlice || bFieldSlice) {
					if ((block_buffer.energy)*JTONEV <= threshold) {
							tBlockHist[0]->Fill((double)block_buffer.time,histoWeight);				
							fieldBHist[0]->Fill(field_block,histoWeight);
					} else {
						for (int j=1; j < numHistos; j++) {
							if (((block_buffer.energy)*JTONEV > (((double)j-1.0)*spacing + threshold)) && ((block_buffer.energy)*JTONEV <= (((double)j)*spacing + threshold))) {
								tBlockHist[j]->Fill((double)block_buffer.time,histoWeight);
								fieldBHist[j]->Fill(field_block,histoWeight);
							}
						}
					}			
				}
				if (blockMap) {
					hitPos->SetPoint(i,block_buffer.position[0],block_buffer.position[1],block_buffer.position[2]);
					if (block_buffer.position[2] > maxCoords[2]) {
						for (int q = 0; q<3; q++) {
							maxCoords[q] = block_buffer.position[q];
						}
					}
					if (block_buffer.position[2] < minCoords[2]) {
						for (int q = 0; q<3; q++) {
							minCoords[q] = block_buffer.position[q];
						}
					}
					hitB->SetPoint(i,block_buffer.position[0],block_buffer.position[1],field_block);
				}
			} else if (!upscatter && (double)block_buffer.time < deathTime && (block_buffer.energy*JTONEV > eThMin) && (angleTheta < theCut)) { 
				// Check if we have a block hit outside of the block
				brokenCounter = brokenCounter + 1;
			} else if (upscatter  && (double)block_buffer.time < deathTime && (block_buffer.energy*JTONEV > eThMin) && (angleTheta < theCut)) {
				doubleCounts = doubleCounts + 1;
			}
		}
	}
	
	printf("Broken Trajectories: %d \nDouble Hits: %d\n", brokenCounter, doubleCounts);
	printf("\n---------------------------------------------------------------------\n");
	printf(" !!! OUTPUT DATA !!! ");
	printf("\n---------------------------------------------------------------------\n");
	
	
	printf("The neutron yield is: %e (+- %e) !\n", neutronYield, yieldError);
	printf("The loss yield is: %e (+- %e) ! \n", lossYield, lossError);
	printf("This includes a yield of: \n");
	printf("%e (+- %e) (%f +- %f percent) (group 1) \n", groupYield[0],groupError[0], groupYield[0]/neutronYield, groupError[0]/neutronYield);
	printf("%e (+- %e) (%f +- %f percent) (group 2) \n", groupYield[1],groupError[1], groupYield[1]/neutronYield, groupError[1]/neutronYield);
	printf("%e (+- %e) (%f +- %f percent) (group 3) \n", groupYield[2],groupError[2], groupYield[2]/neutronYield, groupError[2]/neutronYield);
	
	
	theApp.Run();
    
    delete[] randDeathTimes;
    delete[] randAbsProb;
    delete[] buf;
	
    return 0;
}

//======================================================================
// ---------------------- Here's the actual MAIN program!! ----------------------------------
int blockLifetimeShift(int argc, char** argv) {
  //  ROOT::EnableThreadSafety();
    initxorshift();

    if(argc <= 3) {
        printf("Error! Usage: ./block_lifetime_fit folder t1 t2 t... \n");
        return 0;
    }
    int numArgs = argc - 2; // number of arguments -- find now 
    std::vector<double> hold_t;
	
    // Analysis and weighting variables 
    /*bool upscatter = false;
    int blockCounter, length_od, length_ob = 0;
    double angleTheta, zeta, boronProb;
    double neutronYield[numArgs] = {0.0};
    double neutronError[numArgs] = {0.0};
    double histoWeight, eneMin, eThMin, theMin, theCut, absMin, dagCut;
    double groupYield[3][numArgs] = {0.0};
    double groupError[3][numArgs] = {0.0};
    double binBreak[2];*/
    
    // Yield finding initialization 
	//binBreak[0] = 40.0; // Change these numbers to change the position of bin dividers (only 3dip for now)
	//binBreak[1] = 60.0;
	//std::vector<double> binBreak = ENDTIMES; // Load binBreak from setup.h
	
		
	// Chi2 fitting dimensions table: 
	// 1 = Energy (power)
	// 2 = Energy (threshold)
	// 3 = Angle  (power) 
	// 4 = Angle  (cutoff)
	//5 = Upscatter (prob)
	// 6 = B10 Bad Region (m)
		
	//int chi2Dim = 6;
	std::vector<double> c2Guess;
	
	
	//double c2Guess[chi2Dim] = {2.15, 0.088, 1.10, 2.0, 1.0, 0.0}; //"Old," not as good numbers
	//double c2Guess[chi2Dim] = {2.0, 0.0, 1.0, 2.0, 1.0, 0.0}; // No weighting or truncation
	double c2Guess[chi2Dim] = {2.30448, 0.076071, 1.264079, 2.0, 1.0e-3, 0.0160763}; //"Nathan's suggested" numbers
	//double c2Guess[chi2Dim] = {3.504480, 0.048571, 2.564079, 1.0, 1.0, 0.052076}; //I found these numbers!	 -- They're super unphysical!
    
    std::vector<std::vector<measurement>> corrYield;
	for (size_t ii = 0; ii < numArgs; ii++) {
		
		hold_t[ii]= std::stod(argv[ii+2]);
		
		// Load vectors from file
		char dagFile[48];
		sprintf(dagFile,"%sdag.out.%ss", argv[1],argv[ii+2],'c');	
		std::vector<fixedResultDag> events_dag = readFileDagRes(dagFile);
		if(events_dag.size() <= 0) {
			return 1;
		}
		
		char blockFile[48];
		sprintf(blockFile,"%sblok.out.%ss", argv[1],argv[ii+2],'c');	    
		std::vector<fixedResultBlock> events_block = readFileBlockRes(blockFile);
		if (events_block.size() <= 0) {
			return 1;
		}
		
		corrYield.push_back(blockUpscatterYield(events_dag,events_block,hold_t[ii]));
	} 
	
	double norm = getNorm(groupYield.at(0));
	for (size_t ii = 0; ii < groupYield.size(); ii++) {
		std::vector<measurement> yBuff = groupYield.at(ii);
	}
		
	return 0;
}
 
// ---------------------- Here's the actual MAIN program!! -----------------------------------
int blockPlottingStuff(int argc, char** argv) {
	// This is the block plotting function, rewritten functionally

	initxorshift(); // RNG

    ROOT::EnableThreadSafety();	// Initialize ROOT and our random number generator
    TApplication theApp("App",0,0); // This is so that we can plot ROOT stuff
    
    if((argc != 3) && (argc != 4) && (argc != 5)) {
        printf("Error! Usage: ./block_spectrum_map fname_dag fname_blok (hold_t) (fitting file) \n");
        return 0;
    }
    
    // Load vectors from file
	std::vector<fixedResultDag> events_dag = readFileDagRes(argv[1]);
	if(events_dag.size() <= 0) {
		return 1;
	}
		
	std::vector<fixedResultBlock> events_block = readFileBlockRes(argv[2]);
	if (events_block.size() <= 0) {
		return 1;
	}
    
    // Initialize random numbers
	double* randU01s = new double[events.size()*2*NRECORDS];
	for(unsigned long i = 0; i < events.size()*2*NRECORDS; i++) {
		randU01s[i] = nextU01();
	}
	double* randDeathTimes = new double[events.size()];
	for(unsigned long i = 0; i < events.size(); i++) {
		randDeathTimes[i] = -TAU_N*log(nextU01());
	}
	
	// Chi2 Fitting happens with the dagger -- Need "REDUCED" dagger events
	std::vector<fixedResultDag>   reduced_dag   = parseDagEvents(events_dag, events_block);
	std::vector<fixedResultBlock> reduced_block = parseBlockEvents(events_dag, events_block);
	
	std::vector<std::vector<double>> chi2_results = block_chi2_sweep(reduced_dag,ref_hist,randU01s,randDeathTimes);
	
	std::vector<double> ch2_red; // Find minimum chi^2
	for (auto cIt = chi2_results.begin(); cIt < chi2_results.end(); cIt++){
		c2_red.append((*cIt).at(0));
	}
	int minc2Ind = std::min_element(c2_red.begin(),c2_red.end());
	std::vector<double> results = chi2_results[minc2Ind];
	
	// calculate weighted dagger and block events
	
    std::vector<double> wtsDag = calcHistoWeight(events_dag,events_block); // Reweight events to fit ideal spectrum
    if(wtsDag.size() != events_dag.size()) {
		printf("Error! Weighting vector mis-sized!\n");
		return 1;
	}
   
	std::vector<double> wtsBlock = calcHistoWeightBlock(events_dag,events_blocK);
	if(wtsBlock.size() != events_block.size()) {
		printf("Error! Weighting vector mis-sized!\n");
		return 1;
	}
    
    // Here's where all the actual plotting things go
    
    theApp.Run();
    return 0;
}    
    
    /*
	
	//----------------------------------------------------------------
    // !!!INITIAL PARAMETERS to quickly change analysis!!! 
    //----------------------------------------------------------------
    // ROOT plotting initialization 
    int nBins1=500;
    nBins = 500;
	numHistos = 5; // Change this number if we want to change number of slices!
	threshold = 0.16; // Change this number if we want to change our slice threshold energy!
	maxEnergy=0.345; // Energies from spectrum generation
	minEnergy=0.018; 
	maxTime = 500.0;
	maxField = 0.2; // Few enough neutrons can penetrate past 0.2 that this should be set.
	spacing = (maxEnergy - threshold)/((numHistos-1));
	
	// Yield finding initialization
	binBreak[0] = 40.0; // Change these numbers to change the position of bin dividers (only 3dip for now)
	binBreak[1] = 60.0;                                                     
	binBuff = 0.0; // Change this number to modify the buffer space on either side of the bin
	gBins[0] = 0.0;
	gBins[5] = maxTime;
	for (int i = 0; i < 2; i++) { // this loop will calculate bin positions. Trust me!
		for (int j = 0; j < 2; j++) {
			gBins[(j+1)+2*i] = binBreak[i] + (double)(2*j-1)*binBuff;
		}
	}
	
	// Chi2 Minimization guesses -- variables declared here since we need chi2Dim declared previously
	
		
	//----------------------------------------------------------------
    // !!!LOAD DATA -- actual function begins here!!! 
    //----------------------------------------------------------------
    // Input/Output data 
    if((argc != 3) && (argc != 4) && (argc != 5)) {
        printf("Error! Usage: ./block_spectrum_map fname_dag fname_blok (hold_t) (fitting file) \n");
        return 0;
    }
	if (argc >= 4) {
		hold_t = std::stod(argv[3]);
		//if (pseOn) {
		//	hold_t = hold_t+200.0;
		//}
	} else {
		hold_t = 0.0;
		decayOn = false;
	}
	
	//nBins = 700+hold_t;
	printf("%d\n",nBins);
	if (argc >= 5) {
		chiFit = true; // Turn on chi^2 fitting if we have a fitting file.
	}
   
    // Have to declare these random number strings here (since they're variable length)
    double* randDeathTimes = new double[length_dag+length_blok];
    double* randAbsProb = new double[2*length_blok];
    double* randDagProb = new double[length_dag];
    // Create random numbers for weighting absorption spectrum 
    for(unsigned long i = 0; i < length_dag+length_blok; i++) {
		if (decayOn) {
			randDeathTimes[i] = -877.7*log(nextU01());
		} else {
			randDeathTimes[i] = DBL_MAX;
		}
    }
    for(unsigned long i = 0; i < 2*length_blok; i++) {
		randAbsProb[i] = nextU01();
		if (randAbsProb[i] < 0.0 || randAbsProb[i] > 1.0) {  // check and make sure random number generator doesn't crash
			i = i-1;
		}
	}
	for(unsigned long i = 0; i  < length_dag; i++) {
		randDagProb[i] = nextU01();
		//randDagProb[i] = 1.0;
		//if (i > length_dag) {
		//	printf("%d,%f \n",i,randDagProb[i]);
		}
	}
    //----------------------------------------------------------------
	// !!!INITIALIZE all ROOT histograms!!! 
	//----------------------------------------------------------------
		
	// Create hitpos and B field scatter plots
	TGraph2D* hitPos = new TGraph2D(length_blok);
		hitPos->SetName("hitPosition");
		hitPos->SetTitle("Hit Position");
	TGraph2D* hitB = new TGraph2D(length_blok);
		hitB->SetName("hitBField");
		hitB->SetTitle("Hit B Field");
	TMultiGraph* posMaps = new TMultiGraph();
	// Define two-exponential decay function
	
	// Chi2 ROOT initializations
	TFile* rootFile;
	TH1D* dataHist;
	printf("Base ROOT histograms loaded\n");*/
	

//}

//double reweightHistos(std::vector<double> doubleblah) {
	// Load our energies from the minimum chi2guesses
/*	eneMin = c2Guess[0];
	eThMin = c2Guess[1];
	theMin = c2Guess[2];
	theCut = c2Guess[3];
	absMin = c2Guess[4];
	dagCut = c2Guess[5];
	
	//----------------------------------------------------------------
	// !!!FILL all ROOT dagger histograms!!! 
	//----------------------------------------------------------------
	printf("Filling Dagger histograms! \n");
	blockCounter = 0;
*/



{	
	printf("Filling Block histograms! \n");
	blockCounter = 0;
	// Fill Necessary Block Data 
	for (int i=0; i<length_blok; i++) {
		block_buffer=events_blok[i];
		deathTime = randDeathTimes[i+length_dag];
		angleTheta = block_buffer.theta;
		upscatter = false;
		// Check Block scattering
		while (block_buffer.nHits[3] > 1) {
			if (randAbsProb[blockCounter] < absMin) {
				if (randAbsProb[blockCounter] < 0.0 || randAbsProb[blockCounter] > 1.0) { // fix buffer overflows or something?
					upscatter = true;
					block_buffer.nHits[3] = 1;
					blockCounter = 0;
				//	brokenCounter = brokenCounter + 1;
				} else {
					upscatter = true;
					blockCounter = blockCounter + block_buffer.nHits[3];
					block_buffer.nHits[3] = 1;
				}
				
			} else {
				upscatter = false;
				blockCounter = blockCounter + 1;
				block_buffer.nHits[3] = block_buffer.nHits[3] - 1;
			}
		}
		// Load histograms that don't scatter and don't die
		if ((!upscatter) && ((double)block_buffer.time < deathTime) && (block_buffer.position[2] < -1.0) && (block_buffer.energy/(GRAV*MASS_N) > eThMin) && (angleTheta < theCut)) {
			
			if (thetaOn) { 
				histoWeight = (pow(cos(angleTheta),theMin)*sin(angleTheta))/(cos(angleTheta)*sin(angleTheta)) * pow((block_buffer.energy),eneMin)/block_buffer.energy;
			} else {
				histoWeight = (pow(cos(angleTheta),theMin)*sin(angleTheta))*pow((block_buffer.energy),eneMin)/block_buffer.energy;
			}
			
			if (blockMapTot) {
				blockHistRoot->Fill(block_buffer.time, (block_buffer.energy)/(GRAV*MASS_N),histoWeight);
				timeBlockHist->Fill((double)block_buffer.time,histoWeight);
				eneBlockHist->Fill((block_buffer.energy)/(GRAV*MASS_N),histoWeight);
			}
			field_block = fieldstrength(block_buffer.position);
			if (blockMapSlice || bFieldSlice) {
				if ((block_buffer.energy)/(GRAV*MASS_N) <= threshold) {
						tBlockHist[0]->Fill((double)block_buffer.time,histoWeight);				
						fieldBHist[0]->Fill(field_block,histoWeight);
				} else {
					for (int j=1; j < numHistos; j++) {
						if (((block_buffer.energy)/(GRAV*MASS_N) > (((double)j-1.0)*spacing + threshold)) && ((block_buffer.energy)/(GRAV*MASS_N) <= (((double)j)*spacing + threshold))) {
							tBlockHist[j]->Fill((double)block_buffer.time,histoWeight);
							fieldBHist[j]->Fill(field_block,histoWeight);
						}
					}
				}			
			}
			if (blockMap) {
				hitPos->SetPoint(i,block_buffer.position[0],block_buffer.position[1],block_buffer.position[2]);
				if (block_buffer.position[2] > maxCoords[2]) {
					for (int q = 0; q<3; q++) {
						maxCoords[q] = block_buffer.position[q];
					}
				}
				if (block_buffer.position[2] < minCoords[2]) {
					for (int q = 0; q<3; q++) {
						minCoords[q] = block_buffer.position[q];
					}
				}
				hitB->SetPoint(i,block_buffer.position[0],block_buffer.position[1],field_block);
			}
		} else if (!upscatter && (double)block_buffer.time < deathTime && (block_buffer.energy/(GRAV*MASS_N) > eThMin) && (angleTheta < theCut)) { 
			// Check if we have a block hit outside of the block 
			brokenCounter = brokenCounter + 1;
		} else if (upscatter  && (double)block_buffer.time < deathTime && (block_buffer.energy/(GRAV*MASS_N) > eThMin) && (angleTheta < theCut)) {
			doubleCounts = doubleCounts + 1;
		}
	}
	
	printf("Broken Trajectories: %d \nDouble Hits: %d\n", brokenCounter, doubleCounts);
	printf("\n---------------------------------------------------------------------\n");
	printf(" !!! OUTPUT DATA !!! ");
	printf("\n---------------------------------------------------------------------\n");
	}

//}

//-------------------------------------------------------------------------------------------------------------------------
//	 !!!PLOTTING AND ANALYSIS!!!
//-------------------------------------------------------------------------------------------------------------------------

	
	
//measurement getYield(double group, blah) {
	
			
	/*printf("The neutron yield is: %e (+- %e) !\n", neutronYield, yieldError);
	printf("This includes a yield of: \n");
	printf("%e (+- %e) (%f +- %f percent) (group 1) \n", groupYield[0],groupError[0], groupYield[0]/neutronYield, groupError[0]/neutronYield);
	printf("%e (+- %e) (%f +- %f percent) (group 2) \n", groupYield[1],groupError[1], groupYield[1]/neutronYield, groupError[1]/neutronYield);
	printf("%e (+- %e) (%f +- %f percent) (group 3) \n", groupYield[2],groupError[2], groupYield[2]/neutronYield, groupError[2]/neutronYield);
	
	
	theApp.Run();
    
    delete[] randDeathTimes;
    delete[] randAbsProb;
    delete[] buf;
	
    return neutronYield;*/
//}

	//sprintf(rootName,"%s",argv[4]);// Load file from ROOT, and force it to fit inside our ROI
	//printf("Loading %s\n",rootName);
	//rootFile = TFile::Open(rootName);
	//dataHist = (TH1D*)rootFile->Get("summedDipS");
	//dataHist->GetXaxis()->SetRangeUser(60.0,140.0);
	/*for (int i = 0; i<length_dag; i++) { 
			
			dag_buffer = events_dag[i];
			upscatter = false;
			
			// Calculate weighting 
			angleTheta = dag_buffer.theta;
			energyCm = dag_buffer.energy*JTONEV;
			if (thetaOn) {
				histoWeight = (pow(cos(angleTheta),theMin)*sin(angleTheta))/(cos(angleTheta)*sin(angleTheta)) * pow((energyCm),eneMin)/energyCm;
			} else {
				histoWeight = (pow(cos(angleTheta),theMin)*sin(angleTheta)) * pow((energyCm),eneMin)/energyCm;
			}
			
			// check if we get absorbed by the bad part of the dagger
		    if(dag_buffer.pos[0] > 0) {
		        zeta = 0.5 - sqrt(dag_buffer.pos[0]*dag_buffer.pos[0] + pow(fabs(dag_buffer.pos[2] - dag_buffer.z_offset) - 1.0, 2));
		    } else {
				zeta = 1.0 - sqrt(dag_buffer.pos[0]*dag_buffer.pos[0] + pow(fabs(dag_buffer.pos[2] - dag_buffer.z_offset) - 0.5, 2));
		    }
		    boronProb = zeta > dagCut ? 1 : zeta/dagCut;
			
			// check if we survived the thresholds
			if ((std::isfinite(histoWeight)) && (energyCm > eThMin) &&
				(angleTheta < theCut) && (randDagProb[i] < boronProb)) {
					
				// check if we upscattered
				while (dag_buffer.nHits[3] > 0) {
					if (randAbsProb[blockCounter] < absMin) {
						upscatter = true;
						blockCounter = blockCounter + dag_buffer.nHits[3];
						dag_buffer.nHits[3] = 0;
					} else {
						upscatter = false;
						blockCounter = blockCounter + 1;
						dag_buffer.nHits[3] = dag_buffer.nHits[3] - 1;
					}
				}
			}
		    
			if (!upscatter) {
				
				
				timeDagFitting->Fill(dag_buffer.time+15.0, histoWeight);
			}
			
		}*/
		
		/*
		for (int i=0; i<length_dag;i++) {
			dag_buffer=events_dag[i];
			if (pseOn) {
				dag_buffer.time = dag_buffer.time - (180.0 + hold_t);
			}
			
			// Calculate event weighting. Weighting based on the DAGGER energy, since the block will have diff. hits. Might want to change to more reasonable number (in energy) 
			angleTheta = dag_buffer.theta;
			if (thetaOn) { 
				histoWeight = (pow(cos(angleTheta),theMin)*sin(angleTheta))/(cos(angleTheta)*sin(angleTheta)) * (pow((dag_buffer.energy),eneMin)/dag_buffer.energy)*JTONEV;
			} else {
				histoWeight = (pow(cos(angleTheta),theMin)*sin(angleTheta)) * (pow((dag_buffer.energy),eneMin)/dag_buffer.energy)*JTONEV;
			}
			// check if we get absorbed by the bad part of the dagger 
		    if(dag_buffer.pos[0] > 0) {
		        zeta = 0.5 - sqrt(dag_buffer.pos[0]*dag_buffer.pos[0] + pow(fabs(dag_buffer.pos[2] - dag_buffer.z_offset) - 1.0, 2));
		    }
		    else {
				zeta = 1.0 - sqrt(dag_buffer.pos[0]*dag_buffer.pos[0] + pow(fabs(dag_buffer.pos[2] - dag_buffer.z_offset) - 0.5, 2));
		    }
			boronProb = zeta > dagCut ? 1 : zeta/dagCut;
			// Check that (a) the neutron is alive, (b) that it's above threshold, and (c) not in bad dagger part. 
			// These are all the neutrons we actually want to count and/or will weight.
			// Fix the randDeathTimes when you get a chance... 
			if (((double)dag_buffer.time + (50.0 + hold_t)) < randDeathTimes[i] && (dag_buffer.energy*JTONEV > eThMin) && (angleTheta < theCut) && dag_buffer.time > 0.0 && (randDagProb[i] < boronProb)){
				
				//* Separate hits into block and dagger events
				upscatter = false;
				//blockCounter = 0;
				// Check if the dagger event actually hit the block. Counter is now index of FIRST block hit 
				if (dag_buffer.nHits[3] > 0){
					// Match first event. The next N will sequentially follow, from first to last hit.
					block_buffer = events_blok[blockCounter];
					if ((block_buffer.theta != dag_buffer.theta) || (block_buffer.phi != dag_buffer.phi)) {
						for (int j=0; j < length_blok; j++) {
							block_buffer = events_blok[j];
							if ((block_buffer.theta == dag_buffer.theta) && (block_buffer.phi == dag_buffer.phi)) {
								blockCounter = j;
								break;
							}
						}	
					}
					for (int j = blockCounter; j < blockCounter + dag_buffer.nHits[3]; j++ ) {
						// Error handling
						if (j >= length_blok) {
							printf(" WARNING: Have too many block hits associated with this dagger hit ");
							break;
						}
						
						block_buffer = events_blok[j];
						//absMin = calcBlockScatter(block_buffer.perp_energy);
						absMin = 1.0;
						
						// Quit loop on first hit that upscatters, and load block data
						if (randAbsProb[j] < absMin) {
							// If we're upscattering, load block data, based on boolean logics
							if (blockMapTot) {
								blockHistRoot->Fill(block_buffer.time, (block_buffer.energy)*JTONEV,histoWeight);
								timeBlockHist->Fill((double)block_buffer.time,histoWeight);
								eneBlockHist->Fill((block_buffer.energy)*JTONEV,histoWeight);
							}
							field_block = fieldstrength(block_buffer.position);
							if (blockMapSlice || bFieldSlice) {
								if ((block_buffer.energy)*JTONEV <= threshold) {
										tBlockHist[0]->Fill((double)block_buffer.time,histoWeight);				
										fieldBHist[0]->Fill(field_block,histoWeight);
								} else {
									for (int j=1; j < numHistos; j++) {
										if (((block_buffer.energy)*JTONEV > (((double)j-1.0)*spacing + threshold)) && ((block_buffer.energy)*JTONEV <= (((double)j)*spacing + threshold))) {
											tBlockHist[j]->Fill((double)block_buffer.time,histoWeight);
											fieldBHist[j]->Fill(field_block,histoWeight);
										}
									}
								}			
							}
							if (blockMap) {
								hitPos->SetPoint(i,block_buffer.position[0],block_buffer.position[1],block_buffer.position[2]);
								if (block_buffer.position[2] > maxCoords[2]) {
									for (int q = 0; q<3; q++) {
										maxCoords[q] = block_buffer.position[q];
									}
								}
								if (block_buffer.position[2] < minCoords[2]) {
									for (int q = 0; q<3; q++) {
										minCoords[q] = block_buffer.position[q];
									}
								}
								hitB->SetPoint(i,block_buffer.position[0],block_buffer.position[1],field_block);
							}
							// If we've upscattered, leave the FOR loop, go to calculate yields
							upscatter = true;
							blockCounter = j+1;
							break;
						} else { 
							doubleCounts = doubleCounts + 1;
						}
					}
				}
				// Load dagger histograms
				if (!upscatter) {
					timeDagHist->Fill((double)dag_buffer.time+15.0, histoWeight);
					if (sanityCheck) {
						thetaDagRoot->Fill((double(dag_buffer.theta)),histoWeight);
					}
					
				}
				if (sanityCheck) {
					eneDagHist->Fill((dag_buffer.energy)*JTONEV,histoWeight);
				}
				
			} else {
				brokenCounter = brokenCounter + 1;
			}
		}
		
		*/
/*
			// Check if the dagger has block hits. If so we have to deal with the block.
			if (dag_buffer.nHitBlock > 0) {
				// Match first event. The next N will sequentially follow, from first to last hit.
				// faster if we find the original hit, but otherwise scan through
				if ((events_block[blockCounter].theta != events_dag[i].theta) || (events_block[blockCounter].phi != events_dag[i].phi)) {
					for (size_t j=blockCounter; j < events_block.size(); j++) {
						block_buffer = events_blok[j];
						if ((block_buffer.theta == dag_buffer.theta) && (block_buffer.phi == dag_buffer.phi)) {
							blockCounter = j;
							break;
						}
					}	
				}			
				// Now find the scattering probability for each
				for (int j = blockCounter; j < blockCounter + dag_buffer.nHits[3]; j++ ) {	
					// Error handling
					if (j >= length_blok) {
						printf(" WARNING: Have too many block hits associated with this dagger hit ");
						break;
					}
					
					block_buffer = events_blok[j];
					absMin = calcBlockScatter(block_buffer.perp_energy);
					//absMin = 1.0;
					
					// Quit loop on first hit that upscatters, and load block data
					if (randAbsProb[j] < absMin) {
						
						// For chi2 fit we just care about filling timeDagFitting()
						upscatter = true;
						//blockCounter = j+1;
						break;
					}
				}
			}		//chi2Vals[(int)pow(2,chi2Dim)+1] = chi2Vals[(int)pow(2,chi2Dim)+1]  + chi2fit(dataHist, timeDagFitting, 80.0, 100.0, 100.0, 130.0);
		//chi2Vals[(int)pow(2,chi2Dim)+1] = chi2fit(dataHist, timeDagFitting, 60.0, 80.0);
		/*
		dataHist->GetXaxis()->SetRangeUser(60.0,76.0);
		//dataHist->Scale(1/dataHist->Integral(60.0,76.0));
		timeDagFitting->GetXaxis()->SetRangeUser(60.0,76.0);
		//timeDagFitting->Scale(1/timeDagFitting->Integral(60.0,76.0));
		//chi2Vals[(int)pow(2,chi2Dim)+1] = dataHist->Chi2Test(timeDagFitting,"UW CHI2/NDF") * (double)timeDagFitting->Integral(60.0,76.0);
		chi2Vals[(int)pow(2,chi2Dim)+1] = timeDagFitting->Chi2Test(dataHist,"WW CHI2/NDF");// * (double)timeDagFitting->Integral(60.0,76.0);
		
		
		//dataHist->Scale(dataHist->Integral(60.0,76.0));
		//timeDagFitting->Scale(timeDagFitting->Integral(60.0,76.0));
		
		// Get ready to comment this out when it doesnt work!
		dataHist->GetXaxis()->SetRangeUser(82.0,140.0);
		//dataHist->Scale(1/dataHist->Integral(86.0,120.0));
		timeDagFitting->GetXaxis()->SetRangeUser(82.0,140.0);
		//timeDagFitting->Scale(1/timeDagFitting->Integral(86.0,120.0));
		//chi2Vals[(int)pow(2,chi2Dim)+1] = chi2Vals[(int)pow(2,chi2Dim)+1] + dataHist->Chi2Test(timeDagFitting,"UW CHI2/NDF") * (double)timeDagFitting->Integral(82.0,130.0);
		chi2Vals[(int)pow(2,chi2Dim)+1] = chi2Vals[(int)pow(2,chi2Dim)+1] + timeDagFitting->Chi2Test(dataHist,"WW CHI2/NDF");// * (double)dataHist->Integral(82.0,130.0);
		*/
				
		//dataHist->Scale(dataHist->Integral(86.0,120.0));
		//timeDagFitting->Scale(timeDagFitting->Integral(86.0,120.0));
		/*
		// Do a chi2 fit for each element -- Try something weird
		dataHist->GetXaxis()->SetRangeUser(60.0,76.0);
		timeDagFitting->GetXaxis()->SetRangeUser(60.0,76.0);
		chi2Vals[(int)pow(2,chi2Dim)+1] = dataHist->Chi2Test(timeDagFitting,"UW CHI2/NDF");
		c2peak[0] = chi2Vals[(int)pow(2,chi2Dim)+1];
		
		// Get ready to comment this out when it doesnt work!
		dataHist->GetXaxis()->SetRangeUser(86.0,120.0);
		timeDagFitting->GetXaxis()->SetRangeUser(86.0,120.0);
		chi2Vals[(int)pow(2,chi2Dim)+1] = chi2Vals[(int)pow(2,chi2Dim)+1] + dataHist->Chi2Test(timeDagFitting,"UW CHI2/NDF");
		c2peak[1] = chi2Vals[(int)pow(2,chi2Dim)+1];
		*/
		//chi2Vals[(int)pow(2,chi2Dim)+1] = c2peak[0] + c2peak[1];
		

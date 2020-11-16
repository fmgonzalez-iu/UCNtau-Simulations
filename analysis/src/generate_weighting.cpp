#include "../inc/generate_weighting.hpp"

/*----------------------------------------------------------------------
 * These functions generate our weighted "histogram." 
 * 
 * By binning and reweighting data from our analysis file, we can fit
 * MC histograms to data. Minimizing the chi^2 residual of these helps 
 * find an ideal spectrum, which then lets us do actual physics with the
 * results!
 * 
 * 
 * TODO: I think evt needs to be well-defined? Call it like nResultsEvt"?
 * TODO: Get fitting histogram's size loaded into here
 *--------------------------------------------------------------------*/ 

double weightE(double eventE, double ePower, double threshold, double maxE) {
	// Power law weight calculation debugging E thresholds
	
	
	double w = -1.0;
	if ((eventE * JTONEV > maxE) || (eventE * JTONEV < threshold)) { // Check if UCN is outside fitting range
		return w;
	}
	if (ePower == 0.0) { // Assuming EdE -- no external weighting
		w = (eventE * JTONEV-threshold) / (maxE - threshold); // EdE
		return w;
	}
	
	// Need to double check if we're one or zero indexed on power
	w = pow(eventE * JTONEV / maxE, ePower - 1); // One indexed
	// w = pow(eventE * JTONEV / maxE, ePower); // Zero indexed
	return w;
}

double weightTh(double eventTh, double cosPower, double threshold, double maxTh) {
	// Lambertian forward-directed debugging theta
	
	double w = -1.0;
	if ((eventTh < threshold) || (eventTh > maxTh)) {
		return w;
	}
	
	w = pow(cos(eventTh), cosPower);
	return w;
	
}

double calcWeight(double eventE, double eventTh, double ePower, double cosPower, double threshold) {
	// Power law weight calculation
	// To assume EdE weighting set ePower to zero.
	// If you don't want theta calculation set cosPower to 0.0.
	
	double w = -1.0;
	if ((threshold >= 0.0) && (eventE * JTONEV < threshold)) { // Check if UCN is below threshold
		return w;
	}		
	if (ePower == 0.0) { // Assuming EdE -- no external weighting
		w = (eventE * JTONEV-threshold)/(EMAX-threshold); // weight is EdE
		return w;
	}
	
	if (eventTh > 0.0) { // Need to double check if we're one or zero indexed on power
		w = pow(eventE * JTONEV / EMAX, ePower-1)*pow(cos(eventTh), cosPower);
		//w = pow(eventE / (0.3444*GRAV*MASS_N), ePower)*pow(cos(eventTh),cosPower); // NB the 0.344 factor thingy is close to JTONEV/EMAX
	} else {
		w = pow(eventE * JTONEV / EMAX, ePower-1);
		//w = pow(eventE * JTONEV / EMAX, ePower);
	}
	return w;	
}

std::vector<weightedBin> createHistQuantMultilayer(double thickOxide, double thickBoron, double threshold, double enePower, double cosPower, std::vector<noabsResult>& events,std::vector<double> ref, double* randU01s, double* randDeathTimes) {
	// Weighting the hists by energy and aborption layer
	// For this to be general, want all these parameters to be loaded. 
	// However, thickOxide, enePower, cosPower can be set to zero,
	// which will speed up processing.

	std::vector<weightedBin> hist; // Initialize histogram as same size as ref
	weightedBin zero = {0.0, 0.0, 0.0}; // weightedBin has wgt, wgtSqr, num
	hist.resize(ref.size(), zero);

	int numCount = 0; // initialize counters
	int numMiss  = 0;
	
	double dipEnds[NDIPS] = DIPTIMES; // Dip times
	
	for(size_t i = 0; i < events.size(); i++) { // loop through events
		
		double weight = calcWeight(events[i].energy,events[i].theta,enePower,cosPower,threshold);
	
		if (weight < 0.0) { // minimum energy threshold
			continue;
		}
				
		for(size_t j = 0; j < NRECORDS; j++) { // loop through the dagger hits
			// Check timing is reasonable
			
			if(events[i].times[j] < dipEnds[1]) { // Ignore events with dagger @ pk1
				continue;
			}
			if(events[i].times[j] - dipEnds[1] > randDeathTimes[i]) { // Check if particle survives
				break;
			}
			if(events[i].times[j] >= dipEnds[NDIPS-1]) { // End of unload
				break;
			}
			// Check if we get an absorption
			if(absorbMultilayer(events[i].ePerps[j], randU01s[i*2*NRECORDS + j], thickOxide, thickBoron)) {
				size_t tInd = size_t(events[i].times[j] - dipEnds[1]);
				if(tInd > (hist.size() - 1)) { // Error checker
					printf("Boo!\n");
					break; // This would segfault by calling a bigger index than there are numbers in hist
				}

				hist[tInd].wgt    += weight; // Increment events
				hist[tInd].wgtSqr += weight*weight;
				hist[tInd].num    += 1;
                
                numCount += 1;
				break;
			}
			if(j == NRECORDS - 1) {
				numMiss += 1;
			}
		}
	}
	//printf("%d %d\n", numCount, numMiss); // Print the counts at end of run
	return hist;
}

std::vector<weightedBin> createHistQuantSpline(double thickOxide, double thickBoron, double threshold, double enePower, double cosPower, std::vector<noabsResult>& events,std::vector<double> ref, double* randU01s, double* randDeathTimes) {
	// Weighting the hists by energy and aborption layer
	// For this to be general, want all these parameters to be loaded. 
	// However, thickOxide, enePower, cosPower can be set to zero,
	// which will speed up processing.
	// This uses a spline fit model for boron/oxide thickness.

	gsl_spline *spline; // This is our spline
	gsl_interp_accel *acc;
	createSplineQuantOxide(thickOxide, thickBoron, &spline, &acc);
	
	std::vector<weightedBin> hist; // Initialize histogram as same size as ref
	weightedBin zero = {0.0, 0.0, 0.0}; // weightedBin has wgt, wgtSqr, num
	hist.resize(ref.size(), zero);

	int numCount = 0; // initialize counters
	int numMiss  = 0;
	
	double dipEnds[NDIPS] = ENDTIMES; // Dip times
	
	for(size_t i = 0; i < events.size(); i++) { // loop through events
		
		double weight = calcWeight(events[i].energy,events[i].theta,enePower,cosPower,threshold);
		if (weight < 0.0) { // minimum energy threshold
			continue;
		}
				
		for(size_t j = 0; j < NRECORDS; j++) { // loop through the dagger hits
			// Check timing is reasonable
			if(events[i].times[j] < dipEnds[1]) { // Ignore events with dagger @ pk1
				continue;
			}
			//if(events[i].times[j] - dipEnds[i] > randDeathTimes[i]) { // Check if particle survives
			//	printf("time + DT : %f,%f\n",events[i].times[j], randDeathTimes[i]);
			//	break;
			//}
			if(events[i].times[j] >= (dipEnds[NDIPS-1])) { // End of unload
				break;
			}
			
			if(absorbSpline(events[i].ePerps[j], randU01s[i*2*NRECORDS + j], spline, acc)) {
				size_t tInd = size_t(events[i].times[j] - dipEnds[1]);
				if (tInd > (hist.size() - 1)) { // Error checker
					printf("Boo!\n");
					break; // This would segfault by calling a bigger index than there are numbers in hist
				}
				
				hist[tInd].wgt    += weight; // Increment events
				hist[tInd].wgtSqr += weight*weight;
				hist[tInd].num    += 1;
                
                numCount += 1;
				break;
			}
			if(j == NRECORDS - 1) {
				numMiss += 1;
			}
		}
	}
	//printf("%d %d\n", numCount,numMiss);
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
    return hist;

}
/*
std::vector<double> calcHistoWeight(std::vector<fixedResultDag> dag, std::vector<fixedResultBlock> block) {
	// BLOCK stuff
	// Return a vector with histogram weights
	
	std::vector<double> wts;
	size_t bBuff = 0;
	for (auto dBuff = dag.begin(); dBuff < dag.end(); dBuff++) { // Loop through dagger events
		double wt = -1.0;
		double zeta = dBuff.x > 0 ? (0.5 - sqrt(dBuff.x*dBuff.x + pow(fabs(dBuff.z - dBuff.z_offset) - 1.0, 2))) 
									: (1.0 - sqrt(dBuff.x*dBuff.x + pow(fabs(dBuff.z - dBuff.z_offset) - 0.5,2))); // Calculate zeta for this event
		double boronProb = zeta > ZETACUT ? 1 : zeta / ZETACUT; // Find probability of B10 with weird cuts
		double tDeath = -TAU_N*log(nextU01());
		double bProb  = nextU01(); 
		bool upscatter = true;
		if ((dBuff.time > 0.0) && (dBuff.time + (50.0+holdT) < tDeath)) &&
			((dBuff.energy / (GRAV*MASS_N) > EMIN) && (dBuff.theta < THEMIN) && (bProb < boronProb)) { // Cut out unphysical UCN
				
			upscatter = false;
			while (dBuff.nBlockHits > 0) { // Loop through block hits on dagger
				if (nextU01() < alBlockLossProb(block.at(bBuff).perp_energy)) {
					bBuff += (size_t)dBuff.nBlockHits;
					upscatter = true;
					break;
				} else {
					bBuff += 1;
					if (bBuff == block.size()) { break; }
					dBuff.nBlockHits -= 1;
				}
			}
		
			if (!upscatter) { // If we didn't upscatter let's add to yields
				wt = (pow(cos(dBuff.theta),THEMIN)*sin(dBuff.theta))/(cos(dBuff.theta)*sin(dBuff.theta)) * pow((dBuff.energy),EMIN)/dBuff.energy;
			}
		}
		wts.push_back(wt);
	}
	return wts;
}*/

TH1D* createHistQuantSplineRoot(double thickOxide, double thickBoron, double threshold, double enePower, double cosPower, std::vector<noabsResult>& events,std::vector<double> ref, double* randU01s, double* randDeathTimes) {
    // Weighting the hists by energy and aborption layer
	// For this to be general, want all these parameters to be loaded. 
	// However, thickOxide, enePower, cosPower can be set to zero,
	// which will speed up processing.
	// This uses a spline fit model for boron/oxide thickness.
    
    gsl_spline *spline; // This is our spline
	gsl_interp_accel *acc;
	createSplineQuantOxide(thickOxide, thickBoron, &spline, &acc);
    
    char hName[256]; // This is our ROOT histogram
    sprintf(hName, "MC%f%f%f", (float)thickOxide, (float)thickBoron, (float)threshold);
    TH1D* hist = new TH1D(hName, "Monte Carlo", (int)ref.size(), 0, (int)ref.size());
    printf("Loaded hist\n");
    double dipEnds[NDIPS] = DIPTIMES; // Dip times
    printf("dipEnds loaded\n");
    for(size_t i = 0; i < events.size(); i++) { // loop through events
		
		double weight = calcWeight(events[i].energy,events[i].theta,enePower,cosPower,threshold);
		
		if (weight < 0.0) { // minimum energy threshold
			continue;
		}
		
		for(size_t j = 0; j < NRECORDS; j++) { // loop through the dagger hits
			// Check timing is reasonable
			if(events[i].times[j] < dipEnds[1]) { // Ignore events with dagger @ pk1
				continue;
			}
			if(events[i].times[j] - dipEnds[1] > randDeathTimes[i]) { // Check if particle survives
				break;
			}
			if(events[i].times[j] >= dipEnds[NDIPS-1]) { // End of unload
				break;
			}
			// Check if we get an absorption	
			if(absorbSpline(events[i].ePerps[j], randU01s[i*2*NRECORDS + j], spline, acc)) {
                hist->Fill(events[i].times[j]-dipEnds[1], weight);// Fill ROOT hist
                break;
            }
        }
    }
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
    hist->Draw(hName);
    return hist;
}

/*std::vector<fixedResultBlock> parseBlockEvents(std::vector<fixedResultDag> dag_events, std::vector<fixedResultBlock> block_events) {
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
}*/

std::vector<weightedBin> createHistQuantSpline_2det(
		double thickOxide1, double thickBoron1, // Assume both detectors have oxide + boron 
		double thickOxide2, double thickBoron2,
		double thresh, double power, double cosPower,
		std::vector<noabsCleanDagResult> &eventsD, std::vector<noabsCleanResult> &eventsC, 
		int refSize, double* randU01s, double* randDeathTimes
		) {
	// Weighting the hists by energy and aborption layer
	//
	// For this to be general, want all these parameters to be loaded. 
	// However, thickOxide, enePower, cosPower can be set to zero,
	// which will speed up processing.
	// This uses a spline fit model for boron/oxide thickness.
	//
	// We generate two splines for our two different detectors
	//
	// Turns out splines use too much memory/don't efficiently allocate memory,
	// so we'll have to do the slower multilayer calculation with less overhead.

	gsl_spline *spline1; // This is our spline for the dagger
	gsl_interp_accel *acc1;
	createSplineQuantOxide(thickOxide1, thickBoron1, &spline1, &acc1);
		
	gsl_spline *spline2; // This is our spline for the AC
	gsl_interp_accel *acc2;
	createSplineQuantOxide(thickOxide2,thickBoron2, &spline2, &acc2);

	std::vector<weightedBin> hist; // Initialize histogram as same size as ref
	weightedBin zero = {0.0, 0.0, 0.0}; // weightedBin has wgt, wgtSqr, num
	hist.resize(refSize, zero); // Note that we're only fitting according to dagger (detector 1)

	int numCount = 0; // initialize counters
	int numMiss  = 0;
	int numClean = 0;
	int numGiant = 0;
	int numPeak1 = 0;
	int numACP1  = 0;
	int uncleaned = 0;
	double dipEnds[NDIPS] = ENDTIMES; // Dip times
	size_t evtCtr = 0; // Event counter for detector 2
	
	for(size_t i = 0; i < eventsD.size(); i++) { // loop through events
		
		double weight = calcWeight(eventsD[i].energy,eventsD[i].theta,power,cosPower,thresh);
		if (weight < 0.0) { // minimum energy threshold
			continue;
		}
				
		bool cleaned = false; // pass flag if we're cleaned
		for(size_t j = 0; j < NRECORDS; j++) { // loop through the dagger hits
			if (eventsD[i].nClean[j] > 0) { // If there's cleaner hits now we have to deal
				for (size_t k = evtCtr; k < eventsC.size(); k++) { // Loop through cleaner hits
					if ((fabs(eventsC[k].theta - eventsD[i].theta) < std::numeric_limits<float>::epsilon())
						&& (fabs(eventsC[k].energy - eventsD[i].energy) < std::numeric_limits<float>::epsilon()*JTONEV)) {
					
					//if ((eventsC[k].theta == eventsD[i].theta) && (eventsC[k].energy == eventsD[i].energy)) {
						evtCtr = k; // This keeps track of where we are initialized
						break;
					}
				}
				
				// Now find the scattering probability for each subsequent hit
				for (size_t k = evtCtr; k < evtCtr + eventsD[i].nClean[j]; k++ ) { 
					// Bad counter check
					if (k >= eventsC.size()) {
						//evtCtr = 0;
						printf( "Cleaner/Dagger size mismatch!\n");
						break;
					}
					// Check timing:
					if ((eventsC[k].t < eventsD[i].times[j]) || (eventsD[i].times[j] < 0.0)) {
						
						if (eventsC[k].y > 0) { // Giant cleaner, 100% efficiency
							cleaned = true;
							evtCtr = k; // if we're cleaned, start searching from this event in the future
							numGiant +=1;
							if (eventsC[k].t > -HOLDTIME && eventsC[k].t  < 0 ) {
								uncleaned += 1;
							}
							break;
						}
						// Active cleaner
						
						if (absorbSpline(eventsC[k].ePerp, randU01s[i*2*NRECORDS + k+j], spline2,acc2)) {
							
							//printf("%e\n",eventsC[k].ePerp);
						//if (absorbMultilayer(eventsC[i].ePerp, nextU01(), thickOxide2, thickBoron2)) {
						//if (absorbMultilayer(eventsC[i].ePerp, randU01s[i*2*NRECORDS + k+j], thickOxide2, thickBoron2)) {
						//if (absorbSpline(eventsC[k].ePerp, nextU01(), spline2,acc2)) {
							cleaned = true;
							evtCtr = k;
							if (eventsC[k].t > 0 && eventsC[k].t < (dipEnds[1]-HOLDTIME)) {
								numACP1 += 1;
							}
							if (eventsC[k].t > -HOLDTIME && eventsC[k].t  < 0 ) {
								//printf("%f\n",eventsC[k].t);
								uncleaned += 1;
							}	
							break;
						}
					} else { // If timing is off
						evtCtr = k;
						if (eventsD[i].nClean[j] != 50) {
							printf("bad timing!\n");
						}
						break;
					}
				}
				if (cleaned) { // note that cleaned only resets at the end of dagger hits
					numClean += 1;
					//evtCtr += eventsD[i].nClean[j];
					break;
				}
			}
										
			// Check timing is reasonable
			//if(eventsD[i].times[j] < (dipEnds[1]-HOLDTIME)) { // Ignore events with dagger @ pk1
			//	continue;
			//}
			// if(eventsD[i].times[j] - dipEnds[1] > randDeathTimes[i]) { // Check if particle survives
			// 	printf("time + DT : %f,%f\n",events[i].times[j], randDeathTimes[i]);
			// 	break;
			// }
			if(eventsD[i].times[j] >= (dipEnds[NDIPS-1]-HOLDTIME)) { // End of unload
				break;
			}
			if (eventsD[i].ePerps[j] <= 0.0) {
				// More than 50 cleaner hits here
				break;
			}
			//printf("%e\n", eventsD[i].ePerps[j]);
			//if(absorbMultilayer(eventsD[i].ePerps[j], randU01s[i*2*NRECORDS + evtCtr+ j], thickOxide1, thickBoron1)) {
			//if(absorbMultilayer(eventsD[i].ePerps[j], nextU01(), thickOxide1, thickBoron1)) {
			
			if(absorbSpline(eventsD[i].ePerps[j], randU01s[i*2*NRECORDS + evtCtr+ j],  spline1, acc1)) {
				if(eventsD[i].times[j] < (dipEnds[1]-HOLDTIME)) {
					numPeak1 += 1;
					break;
				}
				
				
				size_t tInd = size_t(eventsD[i].times[j] - (dipEnds[1] -HOLDTIME));//dipEnds[0]);
				if (tInd > (hist.size() - 1)) { // Error checker
					printf("Boo!\n");
					printf("%f %f\n",eventsD[i].times[j],dipEnds[1]);
					break; // This would segfault by calling a bigger index than there are numbers in hist
				}
				
				hist[tInd].wgt    += weight; // Increment events
				hist[tInd].wgtSqr += weight*weight;
				hist[tInd].num    += 1;
                
                numCount += 1;
				break;
			}	
			
			if(j == NRECORDS - 1) { // Check on the last run
				numMiss += 1;
				break; // Completely unnecessary but w/e
			}
		}	
	}
	printf("%d %d %d %d %d %d\n", numCount,numMiss,numClean,numGiant, uncleaned, numACP1);
	
	
    gsl_spline_free(spline1);
    gsl_interp_accel_free(acc1);
    
    gsl_spline_free(spline2);
    gsl_interp_accel_free(acc2);
    return hist;
	
}

/*std::vector<reducedEvent> parseEvents_2det(
		double thickOxide1, double thickBoron1,
		double thickOxide2, double thickBoron2,
		double thresh, double power, double cosPower,
		std::vector<noabsCleanDagResult> &eventsD, std::vector<noabsCleanResult> &eventsC,
		int refSize, double* randU01s, double* randDeathTimes,
*/
std::vector<weightedBin> createHistQuantSpline_2det_debug(
		double thickOxide1, double thickBoron1, // Assume both detectors have oxide + boron 
		double thickOxide2, double thickBoron2,
		double thresh, double power, double cosPower,
		std::vector<noabsCleanDagResult> &eventsD, std::vector<noabsCleanResult> &eventsC, 
		int refSize, double* randU01s, double* randDeathTimes,
		TH1D* energyHist, TH1D* angleHist) {
	// Weighting the hists by energy and aborption layer
	//
	// For this to be general, want all these parameters to be loaded. 
	// However, thickOxide, enePower, cosPower can be set to zero,
	// which will speed up processing.
	// This uses a spline fit model for boron/oxide thickness.
	//
	// We generate two splines for our two different detectors
	
	gsl_spline *spline1; // This is our spline for the dagger
	gsl_interp_accel *acc1;
	createSplineQuantOxide(thickOxide1, thickBoron1, &spline1, &acc1);
		
	gsl_spline *spline2; // This is our spline for the AC
	gsl_interp_accel *acc2;
	createSplineQuantOxide(thickOxide2,thickBoron2, &spline2, &acc2);

	std::vector<weightedBin> hist; // Initialize histogram as same size as ref
	weightedBin zero = {0.0, 0.0, 0.0}; // weightedBin has wgt, wgtSqr, num
	hist.resize(refSize, zero); // Note that we're only fitting according to dagger (detector 1)

	int numCount = 0; // initialize counters
	int numMiss  = 0;
	int numClean = 0;
	int numGiant = 0;
	int numPeak1 = 0;
	int numACP1  = 0;
	int uncleaned = 0;
	double dipEnds[NDIPS] = ENDTIMES; // Dip times
	size_t evtCtr = 0; // Event counter for detector 2
	
	for(size_t i = 0; i < eventsD.size(); i++) { // loop through events
		
		double weight = calcWeight(eventsD[i].energy,eventsD[i].theta,power,cosPower,thresh);
		if (weight < 0.0) { // minimum energy threshold
			continue;
		}
				
		bool cleaned = false; // pass flag if we're cleaned
		for(size_t j = 0; j < NRECORDS; j++) { // loop through the dagger hits
			if (eventsD[i].nClean[j] > 0) { // If there's cleaner hits now we have to deal
				for (size_t k = evtCtr; k < eventsC.size(); k++) { // Loop through cleaner hits
					if ((fabs(eventsC[k].theta - eventsD[i].theta) < std::numeric_limits<float>::epsilon())
						&& (fabs(eventsC[k].energy - eventsD[i].energy) < std::numeric_limits<float>::epsilon()*JTONEV)) {
					
					//if ((eventsC[k].theta == eventsD[i].theta) && (eventsC[k].energy == eventsD[i].energy)) {
						evtCtr = k; // This keeps track of where we are initialized
						break;
					}
				}
				
				// Now find the scattering probability for each subsequent hit
				for (size_t k = evtCtr; k < evtCtr + eventsD[i].nClean[j]; k++ ) { 
					// Bad counter check
					if (k >= eventsC.size()) {
						//evtCtr = 0;
						printf( "Cleaner/Dagger size mismatch!\n");
						break;
					}
					// Check timing:
					if ((eventsC[k].t < eventsD[i].times[j]) || (eventsD[i].times[j] < 0.0)) {
						
						if (eventsC[k].y > 0) { // Giant cleaner, 100% efficiency
							cleaned = true;
							evtCtr = k; // if we're cleaned, start searching from this event in the future
							numGiant +=1;
							if (eventsC[k].t > -HOLDTIME && eventsC[k].t  < 0 ) {
								uncleaned += 1;
							}
							break;
						}
						// Active cleaner
						
						if (absorbSpline(eventsC[k].ePerp, randU01s[i*2*NRECORDS + k+j], spline2,acc2)) {
							
							//printf("%e\n",eventsC[k].ePerp);
						//if (absorbMultilayer(eventsC[i].ePerp, nextU01(), thickOxide2, thickBoron2)) {
						//if (absorbMultilayer(eventsC[i].ePerp, randU01s[i*2*NRECORDS + k+j], thickOxide2, thickBoron2)) {
						//if (absorbSpline(eventsC[k].ePerp, nextU01(), spline2,acc2)) {
							cleaned = true;
							evtCtr = k;
							if (eventsC[k].t > 0 && eventsC[k].t < (dipEnds[1]-HOLDTIME)) {
								numACP1 += 1;
							}
							if (eventsC[k].t > -HOLDTIME && eventsC[k].t  < 0 ) {
								//printf("%f\n",eventsC[k].t);
								uncleaned += 1;
							}	
							break;
						}
					} else { // If timing is off
						evtCtr = k;
						if (eventsD[i].nClean[j] != 50) {
							printf("bad timing!\n");
						}
						break;
					}
				}
				if (cleaned) { // note that cleaned only resets at the end of dagger hits
					numClean += 1;
					//evtCtr += eventsD[i].nClean[j];
					break;
				}
			}
										
			// Check timing is reasonable
			//if(eventsD[i].times[j] < (dipEnds[1]-HOLDTIME)) { // Ignore events with dagger @ pk1
			//	continue;
			//}
			// if(eventsD[i].times[j] - dipEnds[1] > randDeathTimes[i]) { // Check if particle survives
			// 	printf("time + DT : %f,%f\n",events[i].times[j], randDeathTimes[i]);
			// 	break;
			// }
			if(eventsD[i].times[j] >= (dipEnds[NDIPS-1]-HOLDTIME)) { // End of unload
				break;
			}
			if (eventsD[i].ePerps[j] <= 0.0) {
				// More than 50 cleaner hits here
				break;
			}
			//printf("%e\n", eventsD[i].ePerps[j]);
			//if(absorbMultilayer(eventsD[i].ePerps[j], randU01s[i*2*NRECORDS + evtCtr+ j], thickOxide1, thickBoron1)) {
			//if(absorbMultilayer(eventsD[i].ePerps[j], nextU01(), thickOxide1, thickBoron1)) {
			
			if(absorbSpline(eventsD[i].ePerps[j], randU01s[i*2*NRECORDS + evtCtr+ j],  spline1, acc1)) {
				if(eventsD[i].times[j] < (dipEnds[1]-HOLDTIME)) {
					numPeak1 += 1;
					break;
				}
				
				
				size_t tInd = size_t(eventsD[i].times[j] - (dipEnds[1] -HOLDTIME));//dipEnds[0]);
				if (tInd > (hist.size() - 1)) { // Error checker
					printf("Boo!\n");
					printf("%f %f\n",eventsD[i].times[j],dipEnds[1]);
					break; // This would segfault by calling a bigger index than there are numbers in hist
				}
				
				hist[tInd].wgt    += weight; // Increment events
				hist[tInd].wgtSqr += weight*weight;
				hist[tInd].num    += 1;
                
                numCount += 1;
				break;
			}	
			
			if(j == NRECORDS - 1) { // Check on the last run
				numMiss += 1;
				break; // Completely unnecessary but w/e
			}
		}	
	}
	printf("%d %d %d %d %d %d\n", numCount,numMiss,numClean,numGiant, uncleaned, numACP1);
	
	
    gsl_spline_free(spline1);
    gsl_interp_accel_free(acc1);
    
    gsl_spline_free(spline2);
    gsl_interp_accel_free(acc2);
    return hist;
	


}

std::vector<measurement> fitTCRoot(TH1D* hist, std::vector<double> fitOffsets, std::vector<double> fitEnds) {
	// This just uses ROOT histograms to fit time constants. 
	// Could do a non-root based thing here, but whatever.
	// Assumes we've already generated from createHistQuantSplineRoot
	std::vector<measurement> results;
	
	if (fitOffsets.size() != fitEnds.size()) {
		printf("ERROR: Need fitOffsets and fitEnds to be the same size!\n");
		return results;
	}
	
	for (size_t i = 0; i < fitOffsets.size(); i++) { // Fit each peak
		// ROOT fitting 
		TF1* fit = new TF1("mcFit", "[0]*exp(-(x-[2])/[1])", fitOffsets[i], fitEnds[i]); // assume exponential
        fit->SetParameter(2, fitOffsets[i]);
        fit->SetParLimits(2, fitOffsets[i], fitOffsets[i]);
        fit->SetParameter(0, 1000);
        fit->SetParameter(1, 20);
        hist->Fit(fit, "LQ", "", fitOffsets[i], fitEnds[i]);
        
		measurement tcs = {fit->GetParameter(1), fit->GetParError(1)};
		results.push_back(tcs);
	}
	return results;
}

//std::vector<weightedBin> 
/*
std::vector<double> calcHistoWeightBlock(std::vector<fixedResultBlock> block) {
	// This calculates the histogram weights but for the block this time
	
	std::vector<double> wts;
	for (size_t it = 0; it < block.size(); it++) { // Loop block events
		double tDeath = -TAU_N*log(nextU01());
		bool us = false;
		if (nextU01() < alBlockLossProb(block.at(it).perp_energy)) {
			us = true;
			
			size_t jt = it+1; //If we've upscattered, move to the next UCN
			if (jt < block.size()) {
				while (block.at(jt).nHitBlock > 1) {
					jt += 1;
					wts.push_back(-1.0);
					if (jt == block.size()) { break; }
				}
			} 
			it = jt - 1;
		}
		if ((!us) && ((double)block.at(it).time < tDeath) &&
			(block.at(it).energy/(GRAV*MASS_N) > EMIN) && (block.at(it).theta < THEMIN)) { // Make sure we're weighted right

			double wt = (pow(cos(block.at(it).theta),THEMIN)*sin(block.at(it).theta))/(cos(block.at(it).theta)*sin(block.at(it).theta)) * pow((block.at(it).energy),EMIN)/block.at(it).energy;
			wts.push_back(wt);
		} else { 
			wts.push_back(-1.0); 
		}
	}
	return wts;		
}*/




//-----------------------------------------------------------------------------------------------------------------------------
// Here there be dragons
//-----------------------------------------------------------------------------------------------------------------------------


/*double calcWeightEdE(double eventE, double energy, double threshold) {
	// EdE Weight Calculation
	
	double w = -1.0; // Return for an underthreshold neutron
	if (threshold >= 0.0) { // Check if UCN is below threshold
		if (eventE * JTONEV < threshold) {
			return w;
		}
		w = (eventE * JTONEV-threshold)/(EMAX-threshold); // Fitting an energy and threshold
	} else {
		w = eventE * JTONEV / EMAX; // Just an energy
	}

	return w;
}*/
/*
std::vector<weightedBin> createHistQuantMultilayerEdESpline(double thickOxide, double thickBoron, double threshold, std::vector<evt>& events, double* randU01s, double* randDeathTimes) {
	// weighting the hists by energy and aborption layers, assuming EdE
	// Uses a spline fit model for boron/oxide thickness 

	gsl_spline *spline; // This is our spline
	gsl_interp_accel *acc;
	createSplineQuantOxide(thickOxide, thickBoron, &spline, &acc);

	std::vector<weightedBin> hist; // Initialize histogram as same size 
	weightedBin zero = {0.0, 0.0, 0.0}; // weightedBin has wgt, wgtSqr, num
	hist.resize(184, zero);

	int numCount = 0; // initialize counters
	int numMiss  = 0;
	
	double dipEnds[NDIPS] = DIPTIMES; // Dip times
	
	for(size_t i = 0; i < events.size(); i++) {
		double weight = 0.0;
		if(events[i].energy*JTONEV < threshold) {
			continue;
		}
//        weight = events[i].energy*JTONEV/34.5;
		weight = (events[i].energy*JTONEV-threshold)/(34.5-threshold);

// Check timing is reasonable
		for(int j = 0; j < NRECORDS; j++) {        
            if(events[i].times[j] < 41) {
                continue;
            }
            if(events[i].times[j] - 41 > randDeathTimes[i]) {
                break;
            }
            if(events[i].times[j] >= 225) {
                break;
            }
//            if(absorbMultilayer(events[i].ePerp[j], randU01s[i*100 + j], thickOxide, thickBoron)) {
            if(absorbSpline(events[i].ePerp[j], randU01s[i*2*NRECORDS + j], spline, acc)) {
                if(int(events[i].times[j])-41 > 183) {
                    printf("Boo!\n");
                }
                hist[int(events[i].times[j])-41].wgt += weight;
                hist[int(events[i].times[j])-41].wgtSqr += weight*weight;
                hist[int(events[i].times[j])-41].num += 1;
                numCount += 1;
                break;
            }
        }
    }
    printf("%d\n", numCount);
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
    return hist;
}



std::vector<weightedBin> createHistQuantMultilayerEPowdESpline(double thickOxide, double thickBoron, double threshold, double power, std::vector<evt>& events, double* randU01s, double* randDeathTimes) {
	//----------------------------------------------------------------------
// EPowdE
//----------------------------------------------------------------------


    gsl_spline *spline;
    gsl_interp_accel *acc;
    createSplineQuantOxide(thickOxide, thickBoron, &spline, &acc);
    std::vector<weightedBin> hist;
    weightedBin zero = {0.0, 0.0};
    hist.resize(184, zero);
    int numCount = 0;
    for(unsigned long i = 0; i < events.size(); i++) {
        double weight = 0.0;
        if(events[i].energy*JTONEV < threshold) {
            continue;
        }
        weight = pow(events[i].energy*JTONEV/34.5, power);
//        weight = (events[i].energy*JTONEV-threshold)/(34.5-threshold);

        for(int j = 0; j < NRECORDS; j++) {
            if(events[i].ePerp[j] == 0) {
                continue;
            }
            if(events[i].times[j] < 41) {
                continue;
            }
            if(events[i].times[j] - 41 > randDeathTimes[i]) {
                break;
            }
            if(events[i].times[j] >= 225) {
                break;
            }
//            if(absorbMultilayer(events[i].ePerp[j], randU01s[i*100 + j], thickOxide, thickBoron)) {
            if(absorbSpline(events[i].ePerp[j], randU01s[i*2*NRECORDS + j], spline, acc)) {
                if(int(events[i].times[j])-41 > 183) {
                    printf("Boo!\n");
                }
                hist[int(events[i].times[j])-41].wgt += weight;
                hist[int(events[i].times[j])-41].wgtSqr += weight*weight;
                hist[int(events[i].times[j])-41].num += 1;
                numCount += 1;
                break;
            }
        }
    }
//    printf("%d\n", numCount);
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
    return hist;
}
std::vector<weightedBin> createHistQuantNoOxEPowdEThetaSpline(double thickBoron, double threshold, double power, double cosPower, std::vector<evt>& events, double* randU01s, double* randDeathTimes) {
 

//----------------------------------------------------------------------
// No oxide layer
//----------------------------------------------------------------------

    double sum = 0;
    double sumWgt = 0;
    gsl_spline *spline;
    gsl_interp_accel *acc;
    createSplineQuantOxide(0.0, thickBoron, &spline, &acc);
    std::vector<weightedBin> hist;
    weightedBin zero = {0.0, 0.0};
    hist.resize((int)ceil(END-START), zero);
    int numCount = 0;
    for(unsigned long i = 0; i < events.size(); i++) {
        double weight = 0.0;
        if(events[i].energy*JTONEV < threshold) {
            continue;
        }
//        if(events[i].energy > (0.3444*GRAV*MASS_N)) {
//            continue;
//        }
        weight = pow(events[i].energy*JTONEV/34.5, power)*pow(cos(events[i].theta), cosPower);
//        weight = (events[i].energy*JTONEV-threshold)/(34.5-threshold);

        for(int j = 0; j < NRECORDS; j++) {
            if(events[i].ePerp[j] == 0) {
                continue;
            }
            if(events[i].times[j] < START) {
                continue;
            }
            if(events[i].times[j] - START > randDeathTimes[i]) {
                break;
            }
            if(events[i].times[j] >= END) {
                break;
            }
//            if(absorbMultilayer(events[i].ePerp[j], randU01s[i*100 + j], thickOxide, thickBoron)) {
            if(absorbSpline(events[i].ePerp[j], randU01s[i*2*NRECORDS + j], spline, acc)) {
                sum += (j+1)*weight;
                sumWgt += weight;
                if(int(events[i].times[j])-START > (END-START)) {
                    printf("Boo!\n");
                }
                hist[int(events[i].times[j])-START].wgt += weight;
                hist[int(events[i].times[j])-START].wgtSqr += weight*weight;
                hist[int(events[i].times[j])-START].num += 1;
                numCount += 1;
                break;
            }
        }
    }
//    printf("%d\n", numCount);
    printf("%f\n", sum/sumWgt);
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
    return hist;
}

std::vector<weightedBin> createHistQuantNoOxEPowdEThetaSpline_C(double thickBoron, double threshold, double power, double cosPower, std::vector<evt>& events, double* randU01s, double* randDeathTimes) {
    gsl_spline *spline;
    gsl_interp_accel *acc;
    createSplineQuantOxide(0.0, thickBoron, &spline, &acc);
    std::vector<weightedBin> hist;
    weightedBin zero = {0.0, 0.0};
    hist.resize((int)ceil(END-START), zero);
    int numCount = 0;
    for(unsigned long i = 0; i < events.size(); i++) {
        double weight = 0.0;
        if(events[i].energy*JTONEV < threshold) {
            continue;
        }
        weight = pow(events[i].energy/(0.3444*GRAV*MASS_N), power-1)*pow(cos(events[i].theta), cosPower);
//        weight = (events[i].energy*JTONEV-threshold)/(34.5-threshold);

        for(int j = 0; j < NRECORDS; j++) {
            if(events[i].ePerp[j] == 0) {
                continue;
            }
            if(events[i].times[j] < START) {
                continue;
            }
            if(events[i].times[j] - START > randDeathTimes[i]) {
                break;
            }
            if(events[i].times[j] >= END) {
                break;
            }
//            if(absorbMultilayer(events[i].ePerp[j], randU01s[i*100 + j], thickOxide, thickBoron)) {
            if(absorbSpline(events[i].ePerp[j], randU01s[i*2*NRECORDS + j], spline, acc)) {
                if(int(events[i].times[j])-START > (END-START)) {
                    printf("Boo!\n");
                }
                hist[int(events[i].times[j])-START].wgt += weight;
                hist[int(events[i].times[j])-START].wgtSqr += weight*weight;
                hist[int(events[i].times[j])-START].num += 1;
                numCount += 1;
                break;
            }
        }
    }
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
    return hist;
}
* 
* TH1D* createHistQuantMultilayerEdESplineRoot(double thickOxide, double thickBoron, double threshold, std::vector<evt>& events, double* randU01s, double* randDeathTimes) {
    gsl_spline *spline;
    gsl_interp_accel *acc;
    createSplineQuantOxide(thickOxide, thickBoron, &spline, &acc);
    char hName[256];
    sprintf(hName, "MC%f%f%f", thickOxide, thickBoron, threshold);
    TH1D* hist = new TH1D(hName, "Monte Carlo", 184, 0, 184);
    for(unsigned long i = 0; i < events.size(); i++) {
        double weight = 0.0;
        if(events[i].energy*JTONEV < threshold) {
            continue;
        }
//        weight = events[i].energy*JTONEV/34.5;
        weight = (events[i].energy*JTONEV-threshold)/(34.5-threshold);

        for(int j = 0; j < NRECORDS; j++) {
            if(events[i].times[j] < 41) {
                continue;
            }
            if(events[i].times[j] - 41 > randDeathTimes[i]) {
                break;
            }
            if(events[i].times[j] >= 225) {
                break;
            }
//            if(absorbMultilayer(events[i].ePerp[j], randU01s[i*100 + j], thickOxide, thickBoron)) {
            if(absorbSpline(events[i].ePerp[j], randU01s[i*2*NRECORDS + j], spline, acc)) {
                if(int(events[i].times[j])-41 > 183) {
                    printf("Boo!\n");
                }
                hist->Fill(events[i].times[j]-41, weight);
                break;
            }
        }
    }
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
    return hist;
}
* 
* 
void fitTCQuantMultilayerEdESplineRoot(double thickOxide, double thickBoron, double threshold, std::vector<evt>& events, double* randU01s, double* randDeathTimes, double* tcs, double* tcsErr) {
    double fitOffsets[8] = {4, 42, 62, 82, 102, 122, 142, 162};
    double fitEnds[8] = {38, 58, 78, 98, 118, 138, 158, 178};
    gsl_spline *spline;
    gsl_interp_accel *acc;
    createSplineQuantOxide(thickOxide, thickBoron, &spline, &acc);
    char hName[256];
    sprintf(hName, "MC%f%f%f", thickOxide, thickBoron, threshold);
    TH1D* hist = new TH1D(hName, "Monte Carlo", 184, 0, 184);
    for(unsigned long i = 0; i < events.size(); i++) {
        double weight = 0.0;
        if(events[i].energy*JTONEV < threshold) {
            continue;
        }
//        weight = events[i].energy*JTONEV/34.5;
        weight = (events[i].energy*JTONEV-threshold)/(34.5-threshold);

        for(int j = 0; j < NRECORDS; j++) {
            if(events[i].times[j] < 41) {
                continue;
            }
            if(events[i].times[j] - 41 > randDeathTimes[i]) {
                break;
            }
            if(events[i].times[j] >= 225) {
                break;
            }
//            if(absorbMultilayer(events[i].ePerp[j], randU01s[i*100 + j], thickOxide, thickBoron)) {
            if(absorbSpline(events[i].ePerp[j], randU01s[i*2*NRECORDS + j], spline, acc)) {
                if(int(events[i].times[j])-41 > 183) {
                    printf("Boo!\n");
                }
                hist->Fill(events[i].times[j]-41, weight);
                break;
            }
        }
    }
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
    for(int i = 0; i < 8; i++) {
        TF1* fit = new TF1("mcFit", "[0]*exp(-(x-[2])/[1])", fitOffsets[i], fitEnds[i]);
        fit->SetParameter(2, fitOffsets[i]);
        fit->SetParLimits(2, fitOffsets[i], fitOffsets[i]);
        fit->SetParameter(0, 1000);
        fit->SetParameter(1, 20);
        hist->Fit(fit, "LQ", "", fitOffsets[i], fitEnds[i]);
        tcs[i] = fit->GetParameter(1);
        tcsErr[i] = fit->GetParError(1);
        delete fit;
    }
}
*/



//void createHist(double absProb, double absCut, double threshold, double saturation, std::vector<evt>& events, double* randU01s) {

		//if(events[i].energy*JTONEV < threshold) {
		//	continue;
		//}
			
		
		//(events[i].energy * JTONEV-threshold)/(34.5-threshold); // calculate weighting
//        weight = events[i].energy*JTONEV/34.5;
//if(events[i].times[j] < 41) { 
//if(events[i].times[j] - 41 > randDeathTimes[i]) {
//if(events[i].times[j] >= 225) { 

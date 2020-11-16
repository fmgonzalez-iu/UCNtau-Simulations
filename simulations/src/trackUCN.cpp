#include "../../general/setup.h"
#include "../../general/typedefs.h"

#include <limits>
#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>
#include <assert.h>

#include "../inc/data_io.hpp"
#include "../inc/trackUCN.hpp"
#include "../inc/symplectic.hpp"
#include "../inc/geometry.hpp"
#include "../inc/quant_refl.hpp"

//#include "../setup.h"

extern "C" {
    #include "../inc/fields_nate.h"
    #include "../inc/xorshift.h"
}
/*----------------------------------------------------------------------
 * This file contains UCN tracking and result outputs.
 * 
 * See the ./inc/trackUCN.hpp file for structuring of results
 * 
 *--------------------------------------------------------------------*/

noabsResult daggerHitTimes(std::vector<double> state, double dt, trace tr, double holdT) {
	// Track a bunch of dagger hits
	
	// Tangent and norm vectors for reflection off infinitely thin dagger
    std::vector<double> tang = {0.0, 0.0, 1.0};
    std::vector<double> normPlus = {0.0, 1.0, 0.0};
    std::vector<double> normMinus = {0.0, -1.0, 0.0};
    
    noabsResult res; // Initialize result (for post-processed reweighting
    res.theta = acos(state[5]/sqrt(state[3]*state[3] + state[4]*state[4] + state[5]*state[5])); // Initial theta term
    
    double t = 0;
    potential(&state[0], &state[1], &state[2], &(res.energy), &t, &tr);
    res.energy = res.energy - MINU + (state[3]*state[3] + state[4]*state[4] + state[5]*state[5])/(2*MASS_N); // Initial energy
    
    double settlingTime; // Load initial settling time (How long is the neutron in the trap at the beginning?)
    do {
        settlingTime = -FILLINGEXPTIME*log(nextU01());
    } while(settlingTime >= FILLINGTIME);
    
    settlingTime = settlingTime + CLEANINGTIME;
	//----- Do we need this function??? ---
    std::vector<double> prevState(6);  // Start by moving around with cleaner down
    int numSteps = settlingTime/dt;
    double energy;
    for(int i = 0; i < numSteps; i++) {
        prevState = state;
        symplecticStep(state, dt, energy, t, tr);
        if((prevState[2] < -1.5+0.38 && state[2] > -1.5+0.38 && state[1] > 0) || (prevState[2] > -1.5+0.38 && state[2] < -1.5+0.38 && state[1] > 0)) { //cleaned (giant)
            res.energy = -1.0;
            return res;
        }
        t = t + dt;
    }
    
    int nHit = 0;
    int nHitHouseLow = 0;
    int nHitHouseHigh = 0;
    while(true) { // Symplectic stepping until we hit NRECORDS on the dagger
        prevState = state;
        symplecticStep(state, dt, energy, t, tr);
        t = t + dt;
        if(t >= 3000.0) { // Time out = 3000s
            for(int i = nHit; i < NRECORDS; i++) {
                res.times[i] = t - settlingTime;
                res.ePerps[i] = 0.0;
                res.zetas[i] = 0.0;
            }
            return res;
        }
        // Kill UCN if cleaned (giant only)
        if((prevState[2] < -1.5+0.38 && state[2] > -1.5+0.38 && state[1] > 0) || (prevState[2] > -1.5+0.38 && state[2] < -1.5+0.38 && state[1] > 0)) { //cleaned
            for(int i = nHit; i < NRECORDS; i++) {
                res.times[i] = t - settlingTime;
                res.ePerps[i] = 0.0;
                res.zetas[i] = 0.0;
            }
            return res;
        }

        if((prevState[1] < 0 && state[1] > 0) || (prevState[1] > 0 && state[1] < 0)) { // Cross dagger plane
            double fracTravel = fabs(prevState[1])/(fabs(state[1]) + fabs(prevState[1]));
            double predX = prevState[0] + fracTravel * (state[0] - prevState[0]);
            double predZ = prevState[2] + fracTravel * (state[2] - prevState[2]);
            
            double zOff = zOffDipCalc(t - settlingTime);
            
            if(checkDagHit(predX, 0.0, predZ, zOff)) { // Did we hit the dagger?
                res.times[nHit] = t - settlingTime;	// Upload this to our res structure
                res.ePerps[nHit] = state[4]*state[4]/(2*MASS_N);
                res.zetas[nHit] = calcDagZeta(predX, 0.0, predZ, zOff);
                nHit += 1;
                if(nHit >= NRECORDS) {
                    break;
                }
                // We're going to keep reflecting forever -- used for finding boron thickness
                if(prevState[1] > 0 && prevState[4] < 0) {
                    reflect(prevState, normPlus, tang); 
                }
                else {
                    reflect(prevState, normMinus, tang);
                }
                state = prevState;
            }
            else if(checkHouseHitLow(predX, 0.0, predZ, zOff)) { // Or did we hit the dagger casing?
                nHitHouseLow += 1;
                if(prevState[1] > 0 && prevState[4] < 0) {
                    reflect(prevState, normPlus, tang);
                }
                else {
                    reflect(prevState, normMinus, tang);
                }
                state = prevState;
            }
            else if(checkHouseHitHigh(predX, 0.0, predZ, zOff)) {
                nHitHouseHigh += 1;
                if(prevState[1] > 0 && prevState[4] < 0) {
                    reflect(prevState, normPlus, tang);
                }
                else {
                    reflect(prevState, normMinus, tang);
                }
                state = prevState;
            }
        }
    }
    return res;
}

noabsCleanDagResultZ noAbsCleanTrack_multidata(std::vector<double> state, double dt, trace tr, double holdT, std::ofstream &bin2) {//char* fName) {
	// Track a bunch of dagger hits
	//std::ofstream bin2(fName, std::ios::out | std::ios::binary);  // Load a binary file -- I can't put a pointer to the file outside of the function so this is my hack
	// Tangent and norm vectors for reflection off infinitely thin dagger
    std::vector<double> normPlus = {0.0, 1.0, 0.0};
    std::vector<double> normMinus = {0.0, -1.0, 0.0};
    std::vector<double> tang = {0.0, 0.0, 1.0};
    
    noabsCleanDagResultZ res; // Initialize result (for post-processed reweighting (convert to zeta form))
    noabsCleanResult resC;
    res.theta  = acos(state[5]/sqrt(state[3]*state[3] + state[4]*state[4] + state[5]*state[5])); // Initial theta term
    resC.theta = res.theta; // Initial theta term (cleaner)
       
    double t = 0;
    potential(&state[0], &state[1], &state[2], &(res.energy), &t, &tr);
    res.energy  = res.energy - MINU + (state[3]*state[3] + state[4]*state[4] + state[5]*state[5])/(2*MASS_N); // Initial energy
    resC.energy = res.energy; // We'll match these in post
    
    int nClean = 0; // Counter for cleaner
    int nHit = 0; // Counter for hits
    double settlingTime; // Load initial settling time (How long is the neutron in the trap at the beginning?)
    do {
        settlingTime = -FILLINGEXPTIME*log(nextU01());
    } while(settlingTime >= FILLINGTIME);
    
    double tOffset = settlingTime + CLEANINGTIME + holdT; // Amount of time neutron is in trap before end of hold
    settlingTime = settlingTime + CLEANINGTIME;
	
    std::vector<double> prevState(6);  // Start by moving around with cleaner down
    int numSteps = settlingTime/dt;
    double energy;
    for(int i = 0; i < numSteps; i++) {
        prevState = state;
        symplecticStep(state, dt, energy, t, tr);
        
        int code = checkClean(state, prevState, CLEANINGHEIGHT); // Functional cleaning
        if (code != 0) {			
			// Now we save a cleaner hit
			resC.ePerp = state[5]*state[5]/(2*MASS_N); // Perpendicular energy
			resC.t = t - tOffset;
			resC.x = prevState[0];
			resC.y = prevState[1];
			resC.z = prevState[2];
			
			writenoabsCleanRes(bin2,resC); // This saves resC, now can fix in post
								
			std::vector<double> normZPlus = {0.0, 0.0, 1.0};// Need to have UCNs scatter if we hit the cleaner
			std::vector<double> normZMinus = {0.0, 0.0, -1.0};
			std::vector<double> cleanTang = {0.0, 1.0, 0.0};
			if(prevState[2] < -1.5 + CLEANINGHEIGHT && state[2] > -1.5 + CLEANINGHEIGHT) { //moving up
				reflect(prevState, normZMinus, cleanTang);
			} else if(prevState[2] > -1.5 + CLEANINGHEIGHT && state[2] < -1.5 + CLEANINGHEIGHT) {
				reflect(prevState, normZPlus, cleanTang);
			} else {
				printf("Boo!\n");
			}
			state = prevState;
			nClean += 1;
		}
        t = t + dt;
        
        if (nClean + nHit >= NRECORDS) { // A check on nClean -- for a "cleaner" timeout
			for(int i = nHit; i < NRECORDS; i++) { // Fill the remainder of res with zeros
                res.times[i] = t - tOffset;
                res.ePerps[i] = 0.0;
                res.zetas[i]  = 0.0;
                res.nClean[i] = nClean;
            }       
            return res;
		}
    }
        
    // Now look at hold + unload
    while(true) { // Symplectic stepping until we hit NRECORDS on the dagger
        prevState = state;
        symplecticStep(state, dt, energy, t, tr);
        t = t + dt;
        if(t >= 3000.0) { // Time out = 3000s
            for(int i = nHit; i < NRECORDS; i++) {
                res.times[i] = t - tOffset;
                res.ePerps[i] = 0.0; // Should have 0 probability 
                res.zetas[i]  = 0.0;
                res.nClean[i] = nClean;
            }
            printf("timeout nclean = %d\n",nClean);
            return res;
        }
        
        int code = checkClean(state, prevState, RAISEDCLEANINGHEIGHT); // Giant Cleaner
        if (code == 1) { 
			// Now we save a cleaner hit
			resC.ePerp = state[5]*state[5]/(2*MASS_N); // Perpendicular energy
			resC.t = t - tOffset;
			resC.x = prevState[0];
			resC.y = prevState[1];
			resC.z = prevState[2];
			writenoabsCleanRes(bin2,resC); // This saves resC, now can fix in post
							
			std::vector<double> normZPlus = {0.0, 0.0, 1.0};// Need to have UCNs scatter if we hit the cleaner
			std::vector<double> normZMinus = {0.0, 0.0, -1.0};
			std::vector<double> cleanTang = {0.0, 1.0, 0.0};
			if(prevState[2] < -1.5 + RAISEDCLEANINGHEIGHT && state[2] > -1.5 + RAISEDCLEANINGHEIGHT) { // positive z momentum
				reflect(prevState, normZMinus, cleanTang);
			} else if(prevState[2] > -1.5 + RAISEDCLEANINGHEIGHT && state[2] < -1.5 + RAISEDCLEANINGHEIGHT) { // negative z momentum
				reflect(prevState, normZPlus, cleanTang);
			} else {
				printf("Boo!\n");
			}
			state = prevState;
			nClean += 1;
		}
		
		// We're re-inserting the active cleaner here.
		double acHeight = calcCleanHeight(t - (settlingTime + holdT),RAISEDCLEANINGHEIGHT,CLEANINGHEIGHT);
		code = checkClean(state,prevState,acHeight);
        if (code == 2) { // Now we save an active cleaner hit
			resC.ePerp = state[5]*state[5]/(2*MASS_N); // Perpendicular energy
			resC.t = t - tOffset;
			resC.x = prevState[0];
			resC.y = prevState[1];
			resC.z = prevState[2];
			writenoabsCleanRes(bin2,resC); // This saves resC, now can fix in post
						
			std::vector<double> normZPlus = {0.0, 0.0, 1.0};// Need to have UCNs scatter if we hit the cleaner
			std::vector<double> normZMinus = {0.0, 0.0, -1.0};
			std::vector<double> cleanTang = {0.0, 1.0, 0.0};
			if(prevState[2] < -1.5 + acHeight && state[2] > -1.5 + acHeight) { // positive z momentum
				reflect(prevState, normZMinus, cleanTang);
			} else if(prevState[2] > -1.5 + acHeight && state[2] < -1.5 + acHeight) { // negative z momentum
				reflect(prevState, normZPlus, cleanTang);
			} else {
				printf("Boo!\n");
			}
			state = prevState;
			nClean += 1;
		}
        
		// Now that we've checked for cleaning, we can check the dagger hit
		if((prevState[1] < 0 && state[1] > 0) || (prevState[1] > 0 && state[1] < 0)) { // Cross dagger plane
			double fracTravel = fabs(prevState[1])/(fabs(state[1]) + fabs(prevState[1]));
			double predX = prevState[0] + fracTravel * (state[0] - prevState[0]);
			double predZ = prevState[2] + fracTravel * (state[2] - prevState[2]);
            
            //double zOff = zOffDipCalc(t - settlingTime);
            double zOff = zOffDipCalc(t - tOffset);
            if(checkDagHit(predX, 0.0, predZ, zOff)) { // Did we hit the dagger?
                res.times[nHit] = t - tOffset;	// Upload this to our res structure
                res.ePerps[nHit] = state[4]*state[4]/(2*MASS_N);
                res.zetas[nHit]  = calcDagZeta(predX, 0.0, predZ, zOff); 
                res.nClean[nHit] = nClean;
                nHit += 1;
                if(nHit + nClean >= NRECORDS) {
					if (nHit < NRECORDS) {
						for(int i = nHit; i < NRECORDS; i++) { // Fill the remainder of res with zeros
			                res.times[i] = t - tOffset;
			                res.ePerps[i] = 0.0;
			                res.zetas[i]  = 0.0;
			                res.nClean[i] = nClean;
			            } 	
                    } break;
                }
                // We're going to keep reflecting forever -- used for finding boron thickness
                if(prevState[1] > 0 && prevState[4] < 0) { // negative y momentum
                    reflect(prevState, normPlus, tang); 
                } else { // positive y momentum
                    reflect(prevState, normMinus, tang);
                }
                state = prevState;
            } else if(checkHouseHitLow(predX, 0.0, predZ, zOff)) { // Or did we hit the dagger casing?
                if(prevState[1] > 0 && prevState[4] < 0) {
                    reflect(prevState, normPlus, tang);
                } else {
                    reflect(prevState, normMinus, tang);
                }
                state = prevState;
            } else if(checkHouseHitHigh(predX, 0.0, predZ, zOff)) {
                if(prevState[1] > 0 && prevState[4] < 0) {
                    reflect(prevState, normPlus, tang);
                } else {
                    reflect(prevState, normMinus, tang);
                }
                state = prevState;
            }
        }
        
        if (nClean+nHit >= NRECORDS) { // And finally a check on nClean -- for a "cleaner" timeout
			for(int i = nHit; i < NRECORDS; i++) { // Fill the remainder of res with zeros
                res.times[i] = t - tOffset;
                res.ePerps[i] = 0.0;
                res.zetas[i]  = 0.0;
                res.nClean[i] = nClean;
            }       
            return res;
		}
    }
    return res;
}

fixedResult fixedEffDaggerHitTime(std::vector<double> state, double dt, trace tr, double holdT) {
	// For a fixed efficiency, find when a neutron hits the dagger
	
	// Tangent and norm vectors for reflection off infinitely thin dagger
    std::vector<double> tang = {0.0, 0.0, 1.0};
    std::vector<double> normPlus = {0.0, 1.0, 0.0};
    std::vector<double> normMinus = {0.0, -1.0, 0.0};
    
    fixedResult res; // Initialize starting angle and energy
    res.theta = acos(state[5]/sqrt(state[3]*state[3] + state[4]*state[4] + state[5]*state[5]));
    double t = 0;
    potential(&state[0], &state[1], &state[2], &(res.eStart), &t, &tr);
    res.eStart = res.eStart - MINU + (state[3]*state[3] + state[4]*state[4] + state[5]*state[5])/(2*MASS_N);
    
    double deathTime = -TAU_N*log(nextU01()); // Randomly choose when the neutron decays (for speed)

	double settlingTime; // randomly decide when the neutron enters the trap
	do { 
		settlingTime = -FILLINGEXPTIME * log(nextU01());
	} while(settlingTime >= FILLINGTIME);
	
	settlingTime = settlingTime + CLEANINGTIME;

	res.settlingT = settlingTime;
    
#ifdef KILLDECAYEDUCN // If we want to speed up calculation we can kill UCN here by creating this global var
	if(deathTime < FIRSTDIPTIME) {
		res.energy = res.eStart;
		res.t = t;
		res.ePerp = state[4]*state[4]/(2*MASS_N);
		res.x = state[0];
		res.y = state[1];
		res.z = state[2];
		res.zOff = -1;
		res.nHit = 0;
		res.nHitHouseLow = 0;
		res.nHitHouseHigh = 0;
		res.deathTime = deathTime;
        return res;
	}
#endif

	std::vector<double> prevState(6); // Progress through settlingTime
	int numSteps = settlingTime/dt;
	double energy;
	for(int i = 0; i < numSteps; i++) {
		prevState = state;
		symplecticStep(state, dt, energy, t, tr);
		if((prevState[2] < -1.5+CLEANINGHEIGHT && state[2] > -1.5+CLEANINGHEIGHT && state[1] > 0) || (prevState[2] > -1.5+CLEANINGHEIGHT && state[2] < -1.5+CLEANINGHEIGHT && state[1] > 0)) { //cleaned
			res.energy = energy; // Save cleaned neutrons
			res.t = t - settlingTime;
			res.ePerp = state[5]*state[5]/(2*MASS_N);
			res.x = state[0];
			res.y = state[1];
			res.z = state[2];
			res.zOff = -2;
			res.nHit = 0;
			res.nHitHouseLow = 0;
			res.nHitHouseHigh = 0;
			res.deathTime = deathTime;
			return res;
		}
		t = t + dt;
	}

	int nHit = 0;
	int nHitHouseLow = 0;
	int nHitHouseHigh = 0;
	while(true) {
		prevState = state;
		symplecticStep(state, dt, energy, t, tr);
		t = t + dt;
		if(t - settlingTime > deathTime) { // If the neutron decays, save where it dies
			res.energy = energy;
			res.t = t - settlingTime;
			res.ePerp = state[4]*state[4]/(2*MASS_N);
			res.x = state[0];
			res.y = state[1];
			res.z = state[2];
			res.zOff = -1;
			res.nHit = nHit;
			res.nHitHouseLow = nHitHouseLow;
			res.nHitHouseHigh = nHitHouseHigh;
			res.deathTime = deathTime;
			return res;
		}
		if((prevState[2] < -1.5+RAISEDCLEANINGHEIGHT && state[2] > -1.5+RAISEDCLEANINGHEIGHT && state[1] > 0) || (prevState[2] > -1.5+RAISEDCLEANINGHEIGHT && state[2] < -1.5+RAISEDCLEANINGHEIGHT && state[1] > 0)) { //cleaned
			res.energy = energy; // Kill cleaned neutrons too
			res.t = t - settlingTime;
			res.ePerp = state[5]*state[5]/(2*MASS_N);
			res.x = state[0];
			res.y = state[1];
			res.z = state[2];
			res.zOff = -3;
			res.nHit = nHit;
			res.nHitHouseLow = nHitHouseLow;
			res.nHitHouseHigh = nHitHouseHigh;
			res.deathTime = deathTime;
			return res;
		}
		if(std::isnan(energy)) { // If we have a non-real energy (i.e. we got outside the trap) also kill the neutron
			res.energy = energy;
			res.t = t - settlingTime;
			res.ePerp = 0.0;
			res.x = state[0];
			res.y = state[1];
			res.z = state[2];
			res.zOff = -4;
			res.nHit = nHit;
			res.nHitHouseLow = nHitHouseLow;
			res.nHitHouseHigh = nHitHouseHigh;
			res.deathTime = deathTime;
			return res;
        }
		if((prevState[1] < 0 && state[1] > 0) || (prevState[1] > 0 && state[1] < 0)) { // Check if we passed through the dagger
			double fracTravel = fabs(prevState[1])/(fabs(state[1]) + fabs(prevState[1]));
			double predX = prevState[0] + fracTravel * (state[0] - prevState[0]);
			double predZ = prevState[2] + fracTravel * (state[2] - prevState[2]);

			double zOff = zOffDipCalc(t - settlingTime);

			if(checkDagHit(predX, 0.0, predZ, zOff)) { // Check the dagger
				nHit += 1;
				if(absorbMultilayer(state[4]*state[4]/(2*MASS_N), BTHICK, predX, 0.0, predZ, zOff)) { // Calculate B10 abs off dagger
					res.energy = energy;
					res.t = t - settlingTime;
					res.ePerp = state[4]*state[4]/(2*MASS_N);
					res.x = predX;
					res.y = 0.0;
					res.z = predZ;
					res.zOff = zOff;
					res.nHit = nHit;
					res.nHitHouseLow = nHitHouseLow;
					res.nHitHouseHigh = nHitHouseHigh;
					res.deathTime = deathTime;
					return res;
				}
				if(prevState[1] > 0 && prevState[4] < 0) { // If we didn't absorb, reflect
					reflect(prevState, normPlus, tang);
				} else {
					reflect(prevState, normMinus, tang);
				}
				state = prevState;
			} else if(checkHouseHitLow(predX, 0.0, predZ, zOff)) { // Check the dagger casing too
				nHitHouseLow += 1;
				if(prevState[1] > 0 && prevState[4] < 0) {
					reflect(prevState, normPlus, tang);
				} else {
					reflect(prevState, normMinus, tang);
				}
				state = prevState;
			} else if(checkHouseHitHigh(predX, 0.0, predZ, zOff)) { // More dagger casing check
				nHitHouseHigh += 1;
				if(prevState[1] > 0 && prevState[4] < 0) {
					reflect(prevState, normPlus, tang);
				} else {
					reflect(prevState, normMinus, tang);
				}
				state = prevState;
			}
		}
	}
}

fixedResult fixedEffDaggerHitTime_reinsert(std::vector<double> state, double dt, trace tr, double holdT) {
	// Fixed dagger thickness, reinsert small cleaner at end of run during counting
	
	// Initialize reflection vectors
    std::vector<double> tang = {0.0, 0.0, 1.0};
    std::vector<double> normPlus = {0.0, 1.0, 0.0};
    std::vector<double> normMinus = {0.0, -1.0, 0.0};
    
    fixedResult res; // Initialize angle and energy
    res.theta = acos(state[5]/sqrt(state[3]*state[3] + state[4]*state[4] + state[5]*state[5]));
    double t = 0;
    potential(&state[0], &state[1], &state[2], &(res.eStart), &t, &tr);
    res.eStart = res.eStart - MINU + (state[3]*state[3] + state[4]*state[4] + state[5]*state[5])/(2*MASS_N);
    
    double deathTime = -TAU_N*log(nextU01());
    
    double settlingTime;
    do {
        settlingTime = -FILLINGEXPTIME*log(nextU01());
    } while(settlingTime >= FILLINGTIME);
    
    settlingTime = settlingTime + CLEANINGTIME;
    
    res.settlingT = settlingTime;

#ifdef KILLDECAYEDUCN // If we want to speed up calculation we can kill UCN here by creating this global var    
	if(deathTime < FIRSTDIPTIME) {
		res.energy = res.eStart;
		res.t = t - (holdT+settlingTime);
		res.ePerp = state[4]*state[4]/(2*MASS_N);
		res.x = state[0];
		res.y = state[1];
		res.z = state[2];
		res.zOff = -1;
		res.nHit = 0;
		res.nHitHouseLow = 0;
		res.nHitHouseHigh = 0;
		res.deathTime = deathTime;
		return res;
	}
#endif
    
	std::vector<double> prevState(6);
	int numSteps = settlingTime/dt;
	double energy;
	for(int i = 0; i < numSteps; i++) { // Loop through initial cleaning time
		prevState = state;
		symplecticStep(state, dt, energy, t, tr);
		int code = checkClean(state, prevState, RAISEDCLEANINGHEIGHT); // Functional cleaning
		if(code == 1) { // Giant cleaner
#ifdef GC_EFF // Compiler flag, if we have an inefficient cleaner we can reflect
			if(nextU01() < GC_EFF) {
#endif
//        if((prevState[2] < -1.5+CLEANINGHEIGHT && state[2] > -1.5+CLEANINGHEIGHT && state[1] > 0) || (prevState[2] > -1.5+CLEANINGHEIGHT && state[2] < -1.5+CLEANINGHEIGHT && state[1] > 0)) { //cleaned
				res.energy = energy;
	            res.t = t - (holdT+settlingTime);
	            res.ePerp = state[5]*state[5]/(2*MASS_N);
	            res.x = state[0];
	            res.y = state[1];
	            res.z = state[2];
	            res.zOff = -2;
	            res.nHit = 0;
	            res.nHitHouseLow = 0;
	            res.nHitHouseHigh = 0;
	            res.deathTime = deathTime;
	            return res;
#ifdef GC_EFF // Reflect off of cleaner
			} else {
				std::vector<double> normZPlus = {0.0, 0.0, 1.0};
				std::vector<double> normZMinus = {0.0, 0.0, -1.0};
				std::vector<double> cleanTang = {0.0, 1.0, 0.0};
				/*if(state[2] > -1.5 + CLEANINGHEIGHT && prevState[2] < -1.5 + CLEANINGHEIGHT) { //moving up
					reflect(prevState, normZMinus, cleanTang);
				} else if(state[2] < -1.5 + CLEANINGHEIGHT && prevState[2] > -1.5 + CLEANINGHEIGHT) {
					reflect(prevState, normZPlus, cleanTang);
                } else {
					printf("Boo!\n");
				}*/
				if(state[2] > -1.5 + RAISEDCLEANINGHEIGHT && prevState[2] < -1.5 + RAISEDCLEANINGHEIGHT) { //moving up
					reflect(prevState, normZMinus, cleanTang);
				} else if(state[2] < -1.5 + RAISEDCLEANINGHEIGHT && prevState[2] > -1.5 + RAISEDCLEANINGHEIGHT) {
					reflect(prevState, normZPlus, cleanTang);
                } else {
					printf("Boo!\n");
				}
				state = prevState;
			}
#endif
        } else if(code == 2) { // Active cleaner
#ifdef AC_EFF // Compiler flag, if we have an inefficient cleaner we can reflect
			if(nextU01() < AC_EFF) {
#endif
                res.energy = energy;
                res.t = t - (holdT+settlingTime);
                res.ePerp = state[5]*state[5]/(2*MASS_N);
                res.x = state[0];
                res.y = state[1];
                res.z = state[2];
                res.zOff = -2;
                res.nHit = 0;
                res.nHitHouseLow = 0;
                res.nHitHouseHigh = 0;
                res.deathTime = deathTime;
                return res;
#ifdef AC_EFF // If AC doesn't absorb, we need to reflect
			} else {
				std::vector<double> normZPlus = {0.0, 0.0, 1.0};
				std::vector<double> normZMinus = {0.0, 0.0, -1.0};
				std::vector<double> cleanTang = {0.0, 1.0, 0.0};
				/*if(state[2] > -1.5 + CLEANINGHEIGHT && prevState[2] < -1.5 + CLEANINGHEIGHT) { //moving up
					reflect(prevState, normZMinus, cleanTang);
				} else if(state[2] < -1.5 + CLEANINGHEIGHT && prevState[2] > -1.5 + CLEANINGHEIGHT) {
					reflect(prevState, normZPlus, cleanTang);
                } else {
					printf("Boo!\n");
				}*/
				if(state[2] > -1.5 + RAISEDCLEANINGHEIGHT && prevState[2] < -1.5 + RAISEDCLEANINGHEIGHT) { //moving up
					reflect(prevState, normZMinus, cleanTang);
				} else if(state[2] < -1.5 + RAISEDCLEANINGHEIGHT && prevState[2] > -1.5 + RAISEDCLEANINGHEIGHT) {
					reflect(prevState, normZPlus, cleanTang);
                } else {
					printf("Boo!\n");
				}
				state = prevState;
			}
#endif
        }
        t = t + dt;
    }
    
    int nHit = 0; // Now precess through the holding and unload time
    int nHitHouseLow = 0;
    int nHitHouseHigh = 0;
    while(true) {
        prevState = state;
        symplecticStep(state, dt, energy, t, tr);
        t = t + dt;
#ifdef KILLDECAYEDUCN // If we want to speed up calculation we can kill UCN here    
        if(t - settlingTime > deathTime) { // Kill UCN if they would die on this step
            res.energy = energy;
            res.t = t - (holdT+settlingTime);
            res.ePerp = state[4]*state[4]/(2*MASS_N);
            res.x = state[0];
            res.y = state[1];
            res.z = state[2];
            res.zOff = -1;
            res.nHit = nHit;
            res.nHitHouseLow = nHitHouseLow;
            res.nHitHouseHigh = nHitHouseHigh;
            res.deathTime = deathTime;
            return res;
        }
#endif
        int code = checkClean(state, prevState, RAISEDCLEANINGHEIGHT);
        if (code == 1) { // Giant cleaner
        //if((prevState[2] < -1.5+RAISEDCLEANINGHEIGHT && state[2] > -1.5+RAISEDCLEANINGHEIGHT && state[1] > 0) || (prevState[2] > -1.5+RAISEDCLEANINGHEIGHT && state[2] < -1.5+RAISEDCLEANINGHEIGHT && state[1] > 0)) { //cleaned
#ifdef GC_EFF // Check giant cleaner efficiency
			if(nextU01() < GC_EFF) {
#endif
	            res.energy = energy;
	            res.t = t - (holdT + settlingTime);
	            res.ePerp = state[5]*state[5]/(2*MASS_N);
	            res.x = state[0];
	            res.y = state[1];
	            res.z = state[2];
	            res.zOff = -3;
	            res.nHit = nHit;
	            res.nHitHouseLow = nHitHouseLow;
	            res.nHitHouseHigh = nHitHouseHigh;
	            res.deathTime = deathTime;
	            return res;
#ifdef GC_EFF	       
			} else {
                std::vector<double> normZPlus = {0.0, 0.0, 1.0};
                std::vector<double> normZMinus = {0.0, 0.0, -1.0};
                std::vector<double> cleanTang = {0.0, 1.0, 0.0};
                if(state[2] > -1.5 + RAISEDCLEANINGHEIGHT && prevState[2] < -1.5 + RAISEDCLEANINGHEIGHT) { //moving up
                    reflect(prevState, normZMinus, cleanTang);
                }
                else if(state[2] < -1.5 + RAISEDCLEANINGHEIGHT && prevState[2] > -1.5 + RAISEDCLEANINGHEIGHT) {
                    reflect(prevState, normZPlus, cleanTang);
                }
                else {
                    printf("Boo!\n");
                }
                state = prevState;
            }
#endif
        }
#if AC_REINSERT
        double acHeight = calcCleanHeight(t - (settlingTime + holdT),RAISEDCLEANINGHEIGHT,CLEANINGHEIGHT);
#else
		double acHeight = RAISEDCLEANINGHEIGHT;
#endif
        code = checkClean(state,prevState,acHeight);
        // Active cleaner re-insertion height
        /*if(((t - settlingTime <= HOLDTIME) && (code == 2)) || // Raised AC during hold
			((t - settlingTime > HOLDTIME) && (checkClean(state, prevState, CLEANINGHEIGHT) == 2))) { // Lower active cleaner at end
		*/
		if (code == 2) {
#ifdef AC_EFF
            if(nextU01() < AC_EFF) {
#endif
                res.energy = energy;
                res.t = t - (holdT + settlingTime);
                res.ePerp = state[5]*state[5]/(2*MASS_N);
                res.x = state[0];
                res.y = state[1];
                res.z = state[2];
                res.zOff = -5;
                res.nHit = nHit;
                res.nHitHouseLow = nHitHouseLow;
                res.nHitHouseHigh = nHitHouseHigh;
                res.deathTime = deathTime;
                return res;
#ifdef AC_EFF
            } else {
                std::vector<double> normZPlus = {0.0, 0.0, 1.0};
                std::vector<double> normZMinus = {0.0, 0.0, -1.0};
                std::vector<double> cleanTang = {0.0, 1.0, 0.0};
                if(state[2] > -1.5 + acHeight && prevState[2] < -1.5 + acHeight) { //moving up
                    reflect(prevState, normZMinus, cleanTang);
                }
                else if(state[2] < -1.5 + acHeight && prevState[2] > -1.5 + acHeight) {
                    reflect(prevState, normZPlus, cleanTang);
                }
                else {
                    printf("Boo!\n");
                }
                state = prevState;
            }
#endif
        }
        if(std::isnan(energy)) { // flag for NaN in energy (outside of trap!)
            res.energy = energy;
            res.t = t - (holdT + settlingTime);
            res.ePerp = 0.0;
            res.x = state[0];
            res.y = state[1];
            res.z = state[2];
            res.zOff = -4;
            res.nHit = nHit;
            res.nHitHouseLow = nHitHouseLow;
            res.nHitHouseHigh = nHitHouseHigh;
            res.deathTime = deathTime;
            return res;
        }
        if((prevState[1] < 0 && state[1] > 0) || (prevState[1] > 0 && state[1] < 0)) { // Now go to dagger
            double fracTravel = fabs(prevState[1])/(fabs(state[1]) + fabs(prevState[1]));
            double predX = prevState[0] + fracTravel * (state[0] - prevState[0]);
            double predZ = prevState[2] + fracTravel * (state[2] - prevState[2]);
            
            double zOff = zOffDipCalc(t - (holdT+settlingTime));
            
            if(checkDagHit(predX, 0.0, predZ, zOff)) {
                nHit += 1;
                if(absorbMultilayer(state[4]*state[4]/(2*MASS_N), BTHICK, predX, 0.0, predZ, zOff)) { // B10 absorption
                    res.energy = energy;
                    res.t = t - (holdT+settlingTime);
                    res.ePerp = state[4]*state[4]/(2*MASS_N);
                    res.x = predX;
                    res.y = 0.0;
                    res.z = predZ;
                    res.zOff = zOff;
                    res.nHit = nHit;
                    res.nHitHouseLow = nHitHouseLow;
                    res.nHitHouseHigh = nHitHouseHigh;
                    res.deathTime = deathTime;
                    return res;
                }
                if(prevState[1] > 0 && prevState[4] < 0) {
                    reflect(prevState, normPlus, tang);
                }
                else {
                    reflect(prevState, normMinus, tang);
                }
                state = prevState;
            }
            else if(checkHouseHitLow(predX, 0.0, predZ, zOff)) { // Casing checks
                nHitHouseLow += 1;
                if(prevState[1] > 0 && prevState[4] < 0) {
                    reflect(prevState, normPlus, tang);
                }
                else {
                    reflect(prevState, normMinus, tang);
                }
                state = prevState;
            }
            else if(checkHouseHitHigh(predX, 0.0, predZ, zOff)) { // Other casing check
                nHitHouseHigh += 1;
                if(prevState[1] > 0 && prevState[4] < 0) {
                    reflect(prevState, normPlus, tang);
                }
                else {
                    reflect(prevState, normMinus, tang);
                }
                state = prevState;
            }
        }
    }
}

fixedResult fixedEffDaggerHitTime_PSE(std::vector<double> state, double dt, trace tr, double holdT) {
	// Assuming we have a fixed result with dagger lowering for PSE. 
	// TODO: Clean this code up since it's a copy of the other two.
	
    std::vector<double> tang = {0.0, 0.0, 1.0};
    std::vector<double> normPlus = {0.0, 1.0, 0.0};
    std::vector<double> normMinus = {0.0, -1.0, 0.0};
    fixedResult res;
    res.theta = acos(state[5]/sqrt(state[3]*state[3] + state[4]*state[4] + state[5]*state[5]));
    double t = 0;
    potential(&state[0], &state[1], &state[2], &(res.eStart), &t, &tr);
    res.eStart = res.eStart - MINU + (state[3]*state[3] + state[4]*state[4] + state[5]*state[5])/(2*MASS_N);
    
    double deathTime = -TAU_N*log(nextU01());
    
    double settlingTime;
    do {
        settlingTime = -70*log(nextU01());
    } while(settlingTime >= 150);
    
    t = -settlingTime;
    
    res.settlingT = settlingTime;
    
    if(deathTime < FIRSTDIPTIME) {
        res.energy = res.eStart;
        res.t = t;
        res.ePerp = state[4]*state[4]/(2*MASS_N);
        res.x = state[0];
        res.y = state[1];
        res.z = state[2];
        res.zOff = -1;
        res.nHit = 0;
        res.nHitHouseLow = 0;
        res.nHitHouseHigh = 0;
        res.deathTime = deathTime;
        return res;
    }

    std::vector<double> prevState(6);
    int nHit = 0;
    int nHitHouseLow = 0;
    int nHitHouseHigh = 0;
    double energy;
    
    while(true) {
        prevState = state;
        symplecticStep(state, dt, energy, t, tr);
        t = t + dt;
        
        double cleanHeight = t < CLEANINGTIME ? -1.5+0.38 : -1.5+0.38+0.05;
        double zOff = zOffDipCalc((t < 0 ? 0 : t));
        
        if(t > deathTime) {
            res.energy = energy;
            res.t = t;
            res.ePerp = state[4]*state[4]/(2*MASS_N);
            res.x = state[0];
            res.y = state[1];
            res.z = state[2];
            res.zOff = -1;
            res.nHit = nHit;
            res.nHitHouseLow = nHitHouseLow;
            res.nHitHouseHigh = nHitHouseHigh;
            res.deathTime = deathTime;
            return res;
        }
        if((prevState[2] < cleanHeight && state[2] > cleanHeight && state[1] > 0) || (prevState[2] > cleanHeight && state[2] < cleanHeight && state[1] > 0)) { //cleaned
            res.energy = energy;
            res.t = t;
            res.ePerp = state[5]*state[5]/(2*MASS_N);
            res.x = state[0];
            res.y = state[1];
            res.z = state[2];
            res.zOff = -3;
            res.nHit = nHit;
            res.nHitHouseLow = nHitHouseLow;
            res.nHitHouseHigh = nHitHouseHigh;
            res.deathTime = deathTime;
            return res;
        }
        if(std::isnan(energy)) {
            res.energy = energy;
            res.t = t;
            res.ePerp = 0.0;
            res.x = state[0];
            res.y = state[1];
            res.z = state[2];
            res.zOff = -4;
            res.nHit = nHit;
            res.nHitHouseLow = nHitHouseLow;
            res.nHitHouseHigh = nHitHouseHigh;
            res.deathTime = deathTime;
            return res;
        }
        if((prevState[1] < 0 && state[1] > 0) || (prevState[1] > 0 && state[1] < 0)) {
            double fracTravel = fabs(prevState[1])/(fabs(state[1]) + fabs(prevState[1]));
            double predX = prevState[0] + fracTravel * (state[0] - prevState[0]);
            double predZ = prevState[2] + fracTravel * (state[2] - prevState[2]);
            
            if(checkDagHit(predX, 0.0, predZ, zOff)) {
                nHit += 1;
                if(absorbMultilayer(state[4]*state[4]/(2*MASS_N), BTHICK, predX, 0.0, predZ, zOff)) {
                    res.energy = energy;
                    res.t = t;
                    res.ePerp = state[4]*state[4]/(2*MASS_N);
                    res.x = predX;
                    res.y = 0.0;
                    res.z = predZ;
                    res.zOff = zOff;
                    res.nHit = nHit;
                    res.nHitHouseLow = nHitHouseLow;
                    res.nHitHouseHigh = nHitHouseHigh;
                    res.deathTime = deathTime;
                    return res;
                }
                if(prevState[1] > 0 && prevState[4] < 0) {
                    reflect(prevState, normPlus, tang);
                }
                else {
                    reflect(prevState, normMinus, tang);
                }
                state = prevState;
            }
            else if(checkHouseHitLow(predX, 0.0, predZ, zOff)) {
                nHitHouseLow += 1;
                if(prevState[1] > 0 && prevState[4] < 0) {
                    reflect(prevState, normPlus, tang);
                }
                else {
                    reflect(prevState, normMinus, tang);
                }
                state = prevState;
            }
            else if(checkHouseHitHigh(predX, 0.0, predZ, zOff)) {
                nHitHouseHigh += 1;
                if(prevState[1] > 0 && prevState[4] < 0) {
                    reflect(prevState, normPlus, tang);
                }
                else {
                    reflect(prevState, normMinus, tang);
                }
                state = prevState;
            }
        }
    }
}

cleanResult cleanTime(std::vector<double> state, double dt, trace tr, double holdT){
    std::vector<double> tang = {0.0, 0.0, 1.0};
    std::vector<double> normPlus = {0.0, 1.0, 0.0};
    std::vector<double> normMinus = {0.0, -1.0, 0.0};
    std::vector<double> normZPlus = {0.0, 0.0, 1.0};
    std::vector<double> normZMinus = {0.0, 0.0, -1.0};
    std::vector<double> cleanTang = {0.0, 1.0, 0.0};
    cleanResult res;
    res.theta = acos(state[5]/sqrt(state[3]*state[3] + state[4]*state[4] + state[5]*state[5]));
    double t = 0;
    
//    double deathTime = -877.7*log(nextU01());
    
    double settlingTime;
    do {
        settlingTime = -70*log(nextU01());
    } while(settlingTime >= 150);
    
    settlingTime = settlingTime + CLEANINGTIME;

    std::vector<double> prevState(6);
    int numSteps = settlingTime/dt;
    double energy;
    for(int i = 0; i < numSteps; i++) {
        prevState = state;
        symplecticStep(state, dt, energy, t, tr);
        int code = checkClean(state, prevState, CLEANINGHEIGHT);
        if(code == 1) {
            if(nextU01() < 0.75) {
                res.energy = energy;
                res.t = t-settlingTime;
                res.x = state[0];
                res.y = state[1];
                res.z = state[2];
                res.code = -2;
                return res;
            }
            else {
                if(state[2] > -1.5 + CLEANINGHEIGHT && prevState[2] < -1.5 + CLEANINGHEIGHT) { //moving up
                    reflect(prevState, normZMinus, cleanTang);
                }
                else if(state[2] < -1.5 + CLEANINGHEIGHT && prevState[2] > -1.5 + CLEANINGHEIGHT) {
                    reflect(prevState, normZPlus, cleanTang);
                }
                else {
                    printf("Boo!\n");
                }
                state = prevState;
            }
        }
        if(code == 2) {
            if(nextU01() < 0.2) {
                res.energy = energy;
                res.t = t-settlingTime;
                res.x = state[0];
                res.y = state[1];
                res.z = state[2];
                res.code = -5;
                return res;
            }
            else {
                if(state[2] > -1.5 + CLEANINGHEIGHT && prevState[2] < -1.5 + CLEANINGHEIGHT) { //moving up
                    reflect(prevState, normZMinus, cleanTang);
                }
                else if(state[2] < -1.5 + CLEANINGHEIGHT && prevState[2] > -1.5 + CLEANINGHEIGHT) {
                    reflect(prevState, normZPlus, cleanTang);
                }
                else {
                    printf("Boo!\n");
                }
                state = prevState;
            }
        }
        t = t + dt;
    }
    
    res.energy = energy;
    res.t = t-settlingTime;
    res.x = state[0];
    res.y = state[1];
    res.z = state[2];
    res.code = -1;
    return res;
    
    while(true) {
        prevState = state;
        symplecticStep(state, dt, energy, t, tr);
        t = t + dt;
//        if(t - settlingTime > deathTime) {
//            res.energy = energy;
//            res.t = t-settlingTime;
//            res.x = state[0];
//            res.y = state[1];
//            res.z = state[2];
//            res.code = -1;
//            return res;
//        }
        int code = checkClean(state, prevState, RAISEDCLEANINGHEIGHT);
        if(code == 1) {
            if(nextU01() < 0.75) {
                res.energy = energy;
                res.t = t-settlingTime;
                res.x = state[0];
                res.y = state[1];
                res.z = state[2];
                res.code = -3;
                return res;
            }
            else {
                if(state[2] > -1.5 + RAISEDCLEANINGHEIGHT && prevState[2] < -1.5 + RAISEDCLEANINGHEIGHT) { //moving up
                    reflect(prevState, normZMinus, cleanTang);
                }
                else if(state[2] < -1.5 + RAISEDCLEANINGHEIGHT && prevState[2] > -1.5 + RAISEDCLEANINGHEIGHT) {
                    reflect(prevState, normZPlus, cleanTang);
                }
                else {
                    printf("Boo!\n");
                }
                state = prevState;
            }
        }
        if(code == 2) {
            if(nextU01() < 0.2) {
                res.energy = energy;
                res.t = t-settlingTime;
                res.x = state[0];
                res.y = state[1];
                res.z = state[2];
                res.code = -6;
                return res;
            }
            else {
                if(state[2] > -1.5 + RAISEDCLEANINGHEIGHT && prevState[2] < -1.5 + RAISEDCLEANINGHEIGHT) { //moving up
                    reflect(prevState, normZMinus, cleanTang);
                }
                else if(state[2] < -1.5 + RAISEDCLEANINGHEIGHT && prevState[2] > -1.5 + RAISEDCLEANINGHEIGHT) {
                    reflect(prevState, normZPlus, cleanTang);
                }
                else {
                }
                state = prevState;
            }
        }
        if(std::isnan(energy)) {
            res.energy = energy;
            res.t = t-settlingTime;
            res.x = state[0];
            res.y = state[1];
            res.z = state[2];
            res.code = -4;
            return res;
        }
        if((prevState[1] < 0 && state[1] > 0) || (prevState[1] > 0 && state[1] < 0)) {
            double fracTravel = fabs(prevState[1])/(fabs(state[1]) + fabs(prevState[1]));
            double predX = prevState[0] + fracTravel * (state[0] - prevState[0]);
            double predZ = prevState[2] + fracTravel * (state[2] - prevState[2]);
            
            double zOff = zOffDipCalc(t - settlingTime);
            
            if(checkDagHit(predX, 0.0, predZ, zOff)) {
                if(absorbMultilayer(state[4]*state[4]/(2*MASS_N), BTHICK, predX, 0.0, predZ, zOff)) {
                    res.energy = energy;
                    res.t = t-settlingTime;
                    res.x = state[0];
                    res.y = state[1];
                    res.z = state[2];
                    res.code = 0;
                    return res;
                }
                if(prevState[1] > 0 && prevState[4] < 0) {
                    reflect(prevState, normPlus, tang);
                }
                else {
                    reflect(prevState, normMinus, tang);
                }
                state = prevState;
            }
            else if(checkHouseHitLow(predX, 0.0, predZ, zOff)) {
                if(prevState[1] > 0 && prevState[4] < 0) {
                    reflect(prevState, normPlus, tang);
                }
                else {
                    reflect(prevState, normMinus, tang);
                }
                state = prevState;
            }
            else if(checkHouseHitHigh(predX, 0.0, predZ, zOff)) {
                if(prevState[1] > 0 && prevState[4] < 0) {
                    reflect(prevState, normPlus, tang);
                }
                else {
                    reflect(prevState, normMinus, tang);
                }
                state = prevState;
            }
        }
    }
}

noabsResult noabsCleanTime(std::vector<double> state, double dt, trace tr, double holdT){
    std::vector<double> tang = {0.0, 0.0, 1.0};
    std::vector<double> normPlus = {0.0, 1.0, 0.0};
    std::vector<double> normMinus = {0.0, -1.0, 0.0};
    std::vector<double> normZPlus = {0.0, 0.0, 1.0};
    std::vector<double> normZMinus = {0.0, 0.0, -1.0};
    std::vector<double> cleanTang = {0.0, 1.0, 0.0};
    noabsResult res;
    
    double t = 0;
    
    res.theta = acos(state[5]/sqrt(state[3]*state[3] + state[4]*state[4] + state[5]*state[5]));
    
    potential(&state[0], &state[1], &state[2], &(res.energy), &t, &tr);
    res.energy = res.energy - MINU + (state[3]*state[3] + state[4]*state[4] + state[5]*state[5])/(2*MASS_N);
    
    double settlingTime;
    do {
        settlingTime = -70*log(nextU01());
    } while(settlingTime >= 150);
    
    settlingTime = settlingTime + CLEANINGTIME;
    
    int nHit = 0;

    std::vector<double> prevState(6);
    int numSteps = settlingTime/dt;
    double energy;
    for(int i = 0; i < numSteps; i++) {
        if(nHit >= NRECORDS) {
            return res;
        }
        prevState = state;
        symplecticStep(state, dt, energy, t, tr);
        int code = checkClean(state, prevState, CLEANINGHEIGHT);
        if(code == 1) {
            res.times[nHit] = t - settlingTime;
            res.ePerps[nHit] = state[5]*state[5]/(2*MASS_N);
            res.zetas[nHit] = -1;
            if(state[2] > -1.5 + CLEANINGHEIGHT && prevState[2] < -1.5 + CLEANINGHEIGHT) { //moving up
                reflect(prevState, normZMinus, cleanTang);
            }
            else if(state[2] < -1.5 + CLEANINGHEIGHT && prevState[2] > -1.5 + CLEANINGHEIGHT) {
                reflect(prevState, normZPlus, cleanTang);
            }
            else {
                printf("Boo!\n");
            }
            state = prevState;
            nHit += 1;
        }
        if(code == 2) {
            res.times[nHit] = t - settlingTime;
            res.ePerps[nHit] = state[5]*state[5]/(2*MASS_N);
            res.zetas[nHit] = -2;
            if(state[2] > -1.5 + CLEANINGHEIGHT && prevState[2] < -1.5 + CLEANINGHEIGHT) { //moving up
                reflect(prevState, normZMinus, cleanTang);
            }
            else if(state[2] < -1.5 + CLEANINGHEIGHT && prevState[2] > -1.5 + CLEANINGHEIGHT) {
                reflect(prevState, normZPlus, cleanTang);
            }
            else {
                printf("Boo!\n");
            }
            state = prevState;
            nHit += 1;
        }
        t = t + dt;
    }
    
    while(nHit < NRECORDS) {
        prevState = state;
        symplecticStep(state, dt, energy, t, tr);
        t = t + dt;
        int code = checkClean(state, prevState, RAISEDCLEANINGHEIGHT);
        if(code == 1) {
            res.times[nHit] = t - settlingTime;
            res.ePerps[nHit] = state[5]*state[5]/(2*MASS_N);
            res.zetas[nHit] = -1;
            if(state[2] > -1.5 + RAISEDCLEANINGHEIGHT && prevState[2] < -1.5 + RAISEDCLEANINGHEIGHT) { //moving up
                reflect(prevState, normZMinus, cleanTang);
            }
            else if(state[2] < -1.5 + RAISEDCLEANINGHEIGHT && prevState[2] > -1.5 + RAISEDCLEANINGHEIGHT) {
                reflect(prevState, normZPlus, cleanTang);
            }
            else {
                printf("Boo!\n");
            }
            state = prevState;
            nHit += 1;
        }
        if(code == 2) {
            res.times[nHit] = t - settlingTime;
            res.ePerps[nHit] = state[5]*state[5]/(2*MASS_N);
            res.zetas[nHit] = -2;
            if(state[2] > -1.5 + RAISEDCLEANINGHEIGHT && prevState[2] < -1.5 + RAISEDCLEANINGHEIGHT) { //moving up
                reflect(prevState, normZMinus, cleanTang);
            }
            else if(state[2] < -1.5 + RAISEDCLEANINGHEIGHT && prevState[2] > -1.5 + RAISEDCLEANINGHEIGHT) {
                reflect(prevState, normZPlus, cleanTang);
            }
            else {
                printf("Boo!\n");
            }
            state = prevState;
            nHit += 1;
        }
        if(std::isnan(energy)) {
            res.times[nHit] = std::numeric_limits<float>::quiet_NaN();
            res.ePerps[nHit] = std::numeric_limits<float>::quiet_NaN();
            res.zetas[nHit] = std::numeric_limits<float>::quiet_NaN();
            break;
        }
        if((prevState[1] < 0 && state[1] > 0) || (prevState[1] > 0 && state[1] < 0)) {
            double fracTravel = fabs(prevState[1])/(fabs(state[1]) + fabs(prevState[1]));
            double predX = prevState[0] + fracTravel * (state[0] - prevState[0]);
            double predZ = prevState[2] + fracTravel * (state[2] - prevState[2]);
            
            double zOff = zOffDipCalc(t - settlingTime);
            
            if(checkDagHit(predX, 0.0, predZ, zOff)) {
                res.times[nHit] = t - settlingTime;
                res.ePerps[nHit] = state[4]*state[4]/(2*MASS_N);
                res.zetas[nHit] = calcDagZeta(predX, 0.0, predZ, zOff);
                nHit += 1;
                if(prevState[1] > 0 && prevState[4] < 0) {
                    reflect(prevState, normPlus, tang);
                }
                else {
                    reflect(prevState, normMinus, tang);
                }
                state = prevState;
            }
            else if(checkHouseHitLow(predX, 0.0, predZ, zOff)) {
                if(prevState[1] > 0 && prevState[4] < 0) {
                    reflect(prevState, normPlus, tang);
                }
                else {
                    reflect(prevState, normMinus, tang);
                }
                state = prevState;
            }
            else if(checkHouseHitHigh(predX, 0.0, predZ, zOff)) {
                if(prevState[1] > 0 && prevState[4] < 0) {
                    reflect(prevState, normPlus, tang);
                } else {
                    reflect(prevState, normMinus, tang);
                }
                state = prevState;
            }
        }
    }
    
    return res;
}

bounceResult timeToBounce(std::vector<double> state, double dt, trace tr, double holdT) {
	
	bounceResult res; // Initialize result and save state
	double t = 0; // tracking through time
	int nHit = 0; // number of records to count
	res.theta = acos(state[5]/sqrt(state[3]*state[3] + state[4]*state[4] + state[5]*state[5]));
	potential(&state[0], &state[1], &state[2], &(res.energy), &t, &tr); // calculate potential
	
    res.energy = res.energy - MINU + (state[3]*state[3] + state[4]*state[4] + state[5]*state[5])/(2*MASS_N);
    double energy;
    
    for (int ii = 0; ii < NRECORDS; ii++) { // Initialize results (in case we break something)
		res.times[ii] = -1.0;
		res.refl[ii] = 0; 
	}
    
    while (nHit < NRECORDS) {
		std::vector<double> prevState = state;
		if (std::isnan(energy)) { break; }
		symplecticStep(state, dt, energy, t, tr);
		t = t + dt;
				
		if (hitTD(prevState, state)) {
			res.times[nHit] = t;
			res.refl[nHit] = 2;
			nHit += 1;
		}
		/*
		if (checkMagReflection(prevState, state)) { // Count reflections
			res.times[nHit] = t;
			if (hitTD(prevState, state)) {
				res.refl[nHit] = 2; // TD is tagged 2,
			} else {
				res.refl[nHit] = 1; // Other is tagged 1
			}
			nHit += 1;
		}*/
		if (t >= 150.0) { break; } // Only track for 300s filling.
	}
	return res;
}

void trackAndPrint(std::vector<double> state, double dt, trace tr, long traj, std::ofstream &bin){
    double t = 0;
    double energy;
    int counter = 0;
    for(int i = 0; i < (int)(60.0/dt); i++) {
		std::vector<double> prevState = state;
		symplecticStep(state, dt, energy, t, tr);
		t = t + dt;
		if (std::isnan(energy)) { break; }
        //if (checkMagReflection(prevState,state)) {
        //if (hitTD(prevState,state)) {
		//	counter += 1;
        printf("%e %e %e %e %e \n", t, state[0], state[1], state[2],energy);
		//}	
    }
    writePos(bin,traj,t,state);
}

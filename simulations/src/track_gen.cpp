#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>
#include <stdio.h>
//#include "../inc/constants.h"
#include "../inc/track_gen.hpp"
#include "../../general/setup.h"

extern "C" {
    #include "../inc/xorshift.h"
    #include "../inc/fields_nate.h"
}

/*----------------------------------------------------------------------
 * This source file contains the generation function for our UCN.
 * 
 * The two randomPointTrapEdE use a normal EdE generating function, while
 * the randomPointTrapOptimum require us to set external optimum parameters
 *  
 * I've modifited so it should just be one randomPointTrap, but I'm prone
 * to making mistakes.
 *--------------------------------------------------------------------*/

std::vector<double> randomPointTrap(trace tr) {
	// This is a randomPointTrap designed to incorporate whatever
	// Should be the only one needed
	
	std::vector<double> state(6);
	
	double maxEnergy = GRAV*MASS_N*(0.01*EMAX); // Max generated UCN energy
    double maxP = sqrt(2*MASS_N*maxEnergy);
		
	double t = 0.0; // potential time buffer
	double energy; // Energy generating
	
	while(true) {
        energy = maxEnergy * nextU01();
        if(energy < EMIN / JTONEV) {
            continue;
        }
        if (EPOW == 1.0) { // Speeds up calculation so we don't need to do powers
			if (nextU01() < energy/maxEnergy) {
				break;
			}
		} else { // If we're exponentially weighting energy
			if(nextU01() < pow(energy/maxEnergy, EPOW)) {
				break;
			}
		}
	}

	double totalU; // Find the potential of our initial state
	do {
		state[2] = -1.464413669130002; // Generate above the trap door
		state[0] = nextU01()*0.15 - 0.075; // Random XY position
		state[1] = nextU01()*0.15 - 0.075;
		potential(&state[0], &state[1], &state[2], &totalU, &t, &tr);
		totalU = totalU - MINU;
	} while(totalU >= energy);
    
	double targetP = sqrt(2.0*MASS_N*(energy - totalU)); // Find initial momentum
	
	// Optimum
	double theta;
	if (THETAPOW == 0.0)  {
		double u1 = nextU01();
		theta = asin(sqrt(u1));
	} else { // Optimally want forward-directed angles
		while(true) {
			double u1 = nextU01();
			theta = asin(sqrt(u1));
			if(nextU01() < pow(cos(theta), THETAPOW)) { // Check forward-directedness
				break;
			}
		}
	}
	double u2 = nextU01();
    double phi = 2 * M_PI * u2;
    
	state[3] = sin(theta)*cos(phi); // Convert angle to momentum direction
    state[4] = sin(theta)*sin(phi);
    state[5] = cos(theta);
	
	// Rotate Momentum Vector
	double pLen = sqrt(state[3]*state[3] + state[4]*state[4] + state[5]*state[5]);
	state[3] = (targetP/pLen)*state[3];
    state[4] = (targetP/pLen)*state[4];
    state[5] = (targetP/pLen)*state[5];
        
    return state;

}

/*std::vector<double> randomPointTrapEdE(trace tr) {
    std::vector<double> state(6);
    double maxEnergy = GRAV*MASS_N*0.3444;
    double maxP = sqrt(2*MASS_N*maxEnergy);
    
    double t = 0.0;    
    
    double energy;
    while(true) {
        energy = maxEnergy * nextU01();
        if(energy < 0.5/JTONEV) {
            continue;
        }
        if(nextU01() < energy/maxEnergy) {
            break;
        }
    }
    
    double totalU;
    do {
        state[2] = -1.464413669130002;
        state[0] = nextU01()*0.15 - 0.075;
        state[1] = nextU01()*0.15 - 0.075;
        potential(&state[0], &state[1], &state[2], &totalU, &t, &tr);
        totalU = totalU - MINU;
    } while(totalU >= energy);
    
    double targetP = sqrt(2.0*MASS_N*(energy - totalU));
    
    double u1 = nextU01();
    double u2 = nextU01();
    double theta = asin(sqrt(u1));
    double phi = 2 * M_PI * u2;
    
    state[3] = sin(theta)*cos(phi);
    state[4] = sin(theta)*sin(phi);
    state[5] = cos(theta);
    
    double pLen = sqrt(state[3]*state[3] + state[4]*state[4] + state[5]*state[5]);
    
    state[3] = (targetP/pLen)*state[3];
    state[4] = (targetP/pLen)*state[4];
    state[5] = (targetP/pLen)*state[5];
    
    return state;
}

std::vector<double> randomPointTrapOptimum(trace tr) {
    std::vector<double> state(6);
    double maxEnergy = GRAV*MASS_N*0.3444;
    double maxP = sqrt(2*MASS_N*maxEnergy);
    
    double t = 0.0;    
    
    double energy;
    while(true) {
        energy = maxEnergy * nextU01();
        if(energy < EMIN/JTONEV) {
            continue;
        }
        if(nextU01() < pow(energy/maxEnergy, EPOW)) {
            break;
        }
    }
    
    double totalU;
    do {
        state[2] = -1.464413669130002;
        state[0] = nextU01()*0.15 - 0.075;
        state[1] = nextU01()*0.15 - 0.075;
        potential(&state[0], &state[1], &state[2], &totalU, &t, &tr);
        totalU = totalU - MINU;
    } while(totalU >= energy);
    
    double targetP = sqrt(2.0*MASS_N*(energy - totalU));
    
    double theta;
    while(true) {
        double u1 = nextU01();
        theta = asin(sqrt(u1));
        if(nextU01() < pow(cos(theta), THETAPOW)) {
            break;
        }
    }
    
    double u2 = nextU01();
    double phi = 2 * M_PI * u2;
    
    state[3] = sin(theta)*cos(phi);
    state[4] = sin(theta)*sin(phi);
    state[5] = cos(theta);
    
    double pLen = sqrt(state[3]*state[3] + state[4]*state[4] + state[5]*state[5]);
    
    state[3] = (targetP/pLen)*state[3];
    state[4] = (targetP/pLen)*state[4];
    state[5] = (targetP/pLen)*state[5];
    
    return state;
}

std::vector<double> randomPointTrapOptimumCleanable(trace tr) {
    std::vector<double> state(6);
    double maxEnergy = GRAV*MASS_N*0.45;
    double maxP = sqrt(2*MASS_N*maxEnergy);
    
    double t = 0.0;    
    
    double energy;
    while(true) {
        energy = maxEnergy * nextU01();
        if(energy < EMIN/JTONEV) {
            continue;
        }
        if(nextU01() < pow(energy/maxEnergy, EPOW)) {
            break;
        }
    }
    
    double totalU;
    do {
        state[2] = -1.464413669130002;
        state[0] = nextU01()*0.15 - 0.075;
        state[1] = nextU01()*0.15 - 0.075;
        potential(&state[0], &state[1], &state[2], &totalU, &t, &tr);
        totalU = totalU - MINU;
    } while(totalU >= energy);
    
    double targetP = sqrt(2.0*MASS_N*(energy - totalU));
    
    double theta;
    while(true) {
        double u1 = nextU01();
        theta = asin(sqrt(u1));
        if(nextU01() < pow(cos(theta), THETAPOW)) {
            break;
        }
    }
    
    double u2 = nextU01();
    double phi = 2 * M_PI * u2;
    
    state[3] = sin(theta)*cos(phi);
    state[4] = sin(theta)*sin(phi);
    state[5] = cos(theta);
    
    double pLen = sqrt(state[3]*state[3] + state[4]*state[4] + state[5]*state[5]);
    
    state[3] = (targetP/pLen)*state[3];
    state[4] = (targetP/pLen)*state[4];
    state[5] = (targetP/pLen)*state[5];
    
    return state;
}

std::vector<double> randomPointTrapEdECleanable(trace tr) {
    std::vector<double> state(6);
    double maxEnergy = GRAV*MASS_N*0.45;
    double maxP = sqrt(2*MASS_N*maxEnergy);
    
    double t = 0.0;    
    
    double energy;
    while(true) {
        energy = maxEnergy * nextU01();
        if(energy < 0.5/JTONEV) {
            continue;
        }
        if(nextU01() < energy/maxEnergy) {
            break;
        }
    }
    
    double totalU;
    do {
        state[2] = -1.464413669130002;
        state[0] = nextU01()*0.15 - 0.075;
        state[1] = nextU01()*0.15 - 0.075;
        potential(&state[0], &state[1], &state[2], &totalU, &t, &tr);
        totalU = totalU - MINU;
    } while(totalU >= energy);
    
    double targetP = sqrt(2.0*MASS_N*(energy - totalU));
    
    double u1 = nextU01();
    double u2 = nextU01();
    double theta = asin(sqrt(u1));
    double phi = 2 * M_PI * u2;
    
    state[3] = sin(theta)*cos(phi);
    state[4] = sin(theta)*sin(phi);
    state[5] = cos(theta);
    
    double pLen = sqrt(state[3]*state[3] + state[4]*state[4] + state[5]*state[5]);
    
    state[3] = (targetP/pLen)*state[3];
    state[4] = (targetP/pLen)*state[4];
    state[5] = (targetP/pLen)*state[5];
    
    return state;
}

std::vector<double> randomPointTrapOptimumOnlyCleanable(trace tr) {
    std::vector<double> state(6);
    double maxEnergy = GRAV*MASS_N*0.45;
    double maxP = sqrt(2*MASS_N*maxEnergy);
    
    double t = 0.0;    
    
    double energy;
    while(true) {
        energy = maxEnergy * nextU01();
        if(energy < ECLEAN) {
            continue;
        }
        if(nextU01() < pow(energy/maxEnergy, EPOW)) {
            break;
        }
    }
    
    double totalU;
    do {
        state[2] = -1.464413669130002;
        state[0] = nextU01()*0.15 - 0.075;
        state[1] = nextU01()*0.15 - 0.075;
        potential(&state[0], &state[1], &state[2], &totalU, &t, &tr);
        totalU = totalU - MINU;
    } while(totalU >= energy);
    
    double targetP = sqrt(2.0*MASS_N*(energy - totalU));
    
    double theta;
    while(true) {
        double u1 = nextU01();
        theta = asin(sqrt(u1));
        if(nextU01() < pow(cos(theta), THETAPOW)) {
            break;
        }
    }
    
    double u2 = nextU01();
    double phi = 2 * M_PI * u2;
    
    state[3] = sin(theta)*cos(phi);
    state[4] = sin(theta)*sin(phi);
    state[5] = cos(theta);
    
    double pLen = sqrt(state[3]*state[3] + state[4]*state[4] + state[5]*state[5]);
    
    state[3] = (targetP/pLen)*state[3];
    state[4] = (targetP/pLen)*state[4];
    state[5] = (targetP/pLen)*state[5];
    
    return state;
}*/

#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <assert.h>
#include "../../general/setup.h"
//#include "../inc/constants.h"
#include "../inc/geometry.hpp"

extern "C" {
	#include "../inc/xorshift.h"
}

/*----------------------------------------------------------------------
 * This file contains a set of functions that contain geometric 
 * properties of the UCNtau experiment.
 * 
 * TODO: Plug in dimensions of small cleaner, dagger, into setup.h
 * -------------------------------------------------------------------*/

std::vector<double> cross(std::vector<double> a, std::vector<double> b) {
	// Cross product of two vectors (must be size 3)
	assert(a.size()>=3);
	assert(b.size()>=3);
	
	std::vector<double> res(3);
	res[0] = a[1]*b[2]-a[2]*b[1];
	res[1] = a[0]*b[2]-a[2]*b[0];
	res[2] = a[0]*b[1]-a[1]*b[0];
	
	return res;
}

std::vector<double> matmul(std::vector<std::vector<double>> mat, std::vector<double> vec) {
	// Multiply a matrix by a vector
	// Don't actually use anywhere, I think I did it wrong...
	
	std::vector<double> buff(mat.size());
	// Loop through the size of matrix:
	for (size_t ii = 0; ii < mat.size(); ii++) {
		assert(mat.at(ii).size() == vec.size()); // Make sure the size is right
		double counter = 0.0;
		for (size_t jj = 0; jj < vec.size(); jj++) {
			buff.at(ii) += mat.at(ii).at(jj) * vec.at(jj);
		}
	}
	
	return buff;
}

void normalize(std::vector<double> &a) {
	// normalize a 3D vector so that it's magnitude is one.
	
	assert(a.size()==3);
	double len = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
	a[0] = a[0]/len;
	a[1] = a[1]/len;
	a[2] = a[2]/len;
}

double zOffDipCalc(double t) {
	// Calculate the dagger height
		
	// From ../general/setup.h
	double acc = DAGGERACC;
	double vel = DAGGERVEL;
	double spr = DAGGERSPR;
	
	//sdouble holdT = HOLDTIME;
	double dipHeights[NDIPS] = HEIGHTS; // From file
	//double dipEnds[NDIPS] = ENDTIMES; // Hardcoded endtime
	double dipEnds[NDIPS] = DIPTIMES; // Endtime from t
	
	if(t > dipEnds[NDIPS-1]) {
		return 0.01; // height is 1cm at end of unload
	}
    
	int i; 
	for(i = 0; i < NDIPS; i++) {
		if(dipEnds[i] > t) {
			break;
		}
	}
	
	if(i == 0) {
		return dipHeights[0];
	}
	
	double target = dipHeights[i];
	double start = dipHeights[i-1];
	
	double distance = fabs(target - start)*1000000;
	
	double time = (distance - vel*(vel/acc)*spr)/(vel*spr);
	
	double sign = dipHeights[i] - dipHeights[i-1] > 0 ? 1 : -1;
	sign = dipHeights[i] == dipHeights[i-1] ? 0 : sign;
	
	if(time > 0) {
//		printf("Flat\n");
		if((t - dipEnds[i-1]) < (vel/acc)) {
			return dipHeights[i-1]
				+ sign*(
					acc*spr*(t - dipEnds[i-1])*(t - dipEnds[i-1])/2.0
				)/1000000;
		} else if((t - dipEnds[i-1]) >= (vel/acc) && (t - dipEnds[i-1]) < ((vel/acc) + time)) {
			return dipHeights[i-1]
				+ sign*(
					(vel/acc)*vel*spr/2.0 + (t - dipEnds[i-1] - vel/acc)*vel*spr
				)/1000000;
		} else if((t - dipEnds[i-1]) >= ((vel/acc) + time) && (t - dipEnds[i-1]) < (2*(vel/acc) + time)){
			return dipHeights[i-1]
				+ sign*(
					(vel/acc)*vel*spr/2.0
					+ time*vel*spr
					+ vel*spr*(t - dipEnds[i-1] - time - vel/acc)
					- acc*spr*(t - dipEnds[i-1] - time - vel/acc)*(t - dipEnds[i-1] - time - vel/acc)/2.0
				)/1000000;
		} else {
			return dipHeights[i];
		}
	} else {
		double halfTime = sqrt(distance/(acc*spr));
//			printf("UnFlat\n");
		if(t - dipEnds[i-1] < halfTime) {
			return dipHeights[i-1]
				+ sign*(
					acc*spr*(t - dipEnds[i-1])*(t - dipEnds[i-1])/2.0
				)/1000000;
		} else if(t - dipEnds[i-1] >= halfTime && t - dipEnds[i-1] < 2*halfTime){
			return dipHeights[i-1]
				+ (dipHeights[i]-dipHeights[i-1])/2
				+ sign*(
					(acc*spr*halfTime*(t - dipEnds[i-1] - halfTime))
					- acc*spr*(t - dipEnds[i-1] - halfTime)*(t - dipEnds[i-1] - halfTime)/2.0
				)/1000000;
		} else {
			return dipHeights[i];
		}
	}
}

double calcCleanHeight(double t, double h1, double h2) {
	// Here I'm calculating the height of the cleaner assuming it moves 
	
	double vel = CLEANERVEL; // cleaner velocity from setup
	double tMax = fabs(h2 - h1) / vel; // Time spent accelerating + decelerating
	
	if (t < 0.0) { // Before cleaner movement, at height 1
		return h1;
	} else if (t > tMax) { // after cleaner movement at height 2
		return h2;
	} else { // Intermediate time
		return h1 + (h2 - h1) * t / tMax;
	}
	
}

//double zOffDipCalc(double t) {
//	double z;
//	int i;
//
//	double holdT = HOLDTIME;
//	double speed;
//	
//	double dipHeights[NDIPS] = HEIGHTS; //3 dip
//	double dipEnds[NDIPS] = ENDTIMES; //3 dip
//
//	if(t > dipEnds[NDIPS-1]) {
//		return 0.01;
//	}
//	
//	for(i = 0; i < NDIPS; i++) {
//		if(dipEnds[i] > t) {
//			break;
//		}
//	}
//	
//	if(i == 0) {
//		return dipHeights[0];
//	}
//	
//	if(dipHeights[i-1] > dipHeights[i]) {
//		speed = -0.49/13.0;
//		z = dipHeights[i-1] + speed*(t-dipEnds[i-1]);
//		return z > dipHeights[i] ? z : dipHeights[i];
//	}
//	else {
//		speed = 0.49/13.0;
//		z = dipHeights[i-1] + speed*(t-dipEnds[i-1]);
//		return z < dipHeights[i] ? z : dipHeights[i];
//	}
//}

void reflect(std::vector<double> &state, std::vector<double> norm, std::vector<double> tang) {
	// Reflection off a solid surface
	
	double pTarget = sqrt(state[3]*state[3] + state[4]*state[4] + state[5]*state[5]); // Cons. of energy
	std::vector<double> tangPrime = cross(norm, tang); // need norm/tangent and derivative
	
	double u1 = nextU01(); // Specular reflection (use RNG)
	double u2 = nextU01();
	
	double theta = asin(sqrt(u1));
	double phi = 2 * M_PI * u2;
	
	double pN = cos(theta);
	double pT = sin(theta)*cos(phi);
	double pTprime = sin(theta)*sin(phi);
	
	std::vector<double> newPdir(3); // Find momentum direction + accel.
	newPdir[0] = pN*norm[0] + pT*tang[0] + pTprime*tangPrime[0];
	newPdir[1] = pN*norm[1] + pT*tang[1] + pTprime*tangPrime[1];
	newPdir[2] = pN*norm[2] + pT*tang[2] + pTprime*tangPrime[2];
	
	double pLen = sqrt(newPdir[0]*newPdir[0] + newPdir[1]*newPdir[1] + newPdir[2]*newPdir[2]);
	
	state[3] = newPdir[0] * pTarget / pLen; // Update momentum vectors, scaling
	state[4] = newPdir[1] * pTarget / pLen;
	state[5] = newPdir[2] * pTarget / pLen;
}

bool checkDagHit(double x, double y, double z, double zOff) {
	// Check if we hit the dagger
	
	double zeta;
	if(x > 0) {
		zeta = 0.5 - sqrt(x*x + pow(fabs(z - zOff) - 1.0, 2));
	} else {
		zeta = 1.0 - sqrt(x*x + pow(fabs(z - zOff) - 0.5, 2));
	} if((x > -0.3524) && (x < 0.0476) && (zeta > 0.0) && (z < (-1.5 + zOff + 0.2))) {
		return true;
	}
	return false;
}

int checkClean(std::vector<double> state, std::vector<double> prevState, double cleanHeight) {
	// Check if we hit the cleaner
	
	if((prevState[2] < -1.5 + cleanHeight && state[2] > -1.5 + cleanHeight) || 
		(prevState[2] > -1.5 + cleanHeight && state[2] < -1.5 + cleanHeight)) { // Did UCN pass through the cleaning plane?
			
		if(state[1] > 0) { // Giant cleaner is just 1/2 of trap
			return 1;
		}
		if((state[0] > SMALLCLEANERXMIN && state[0] < SMALLCLEANERXMAX)
			&& (state[1] > SMALLCLEANERYMIN && state[1] < SMALLCLEANERYMAX)) { // Small active cleaner
			return 2;
		}
		/*if(state[1] > -(0.218041/2 + 0.335121 + 0.3556) && state[1] < -(0.218041/2 + 0.335121) // Small active cleaner
			&& state[0] > (0.115529/2 + 0.212841 - 0.6604) && state[0] < (0.115529/2 + 0.212841)) {
			return 2;
		}*/
	}
	
	return 0;
}

double calcDagZeta(double x, double y, double z, double zOff) {
	// Curvature of the dagger (it's asymmetric!)
	
	double zeta;
	if(x > 0) {
		zeta = 0.5 - sqrt(x*x + pow(fabs(z - zOff) - 1.0, 2));
	}
	else {
		zeta = 1.0 - sqrt(x*x + pow(fabs(z - zOff) - 0.5, 2));
	}
	return zeta;
}

bool checkBlockHit(std::vector<double> &state, std::vector<double> prevState, double &ePerp) {
	// Re-written from original FORTRAN code like a year and a half later
	// This is a lot of hard-coded in stuff, because FORTRAN was weird
	// Should return a reflected state vector
	
	// Need to not crash pls
	assert(state.size()==6);
	assert(prevState.size()==6);
	
	// Hard-coding in Al-block rotation matrices from a Mathematica notebook.
	double blockPos[3] = {0.2175882546128424, 0.1264454150954313, -1.436798256096196}; // block raised by 0.7 mm (most physical)
	
	double rotation[3][3] = {{0.4981846187446881, -0.7463198129200769, -0.4413823244195138},
							{0.8070370832294234, 0.5852378206028751, -0.07866452547073647},
							{0.3170227597175610, -0.3170224470581768, 0.8938643835182667}}; // Rotation Matrix
							
	double invRotation[3][3] = {{0.4981846187446881, 0.8070370832294234, 0.3170224470581768},
								{-0.7463198129200769, 0.5852378206028750, -0.3170224470581768},
								{-0.4413823244195138, -0.07866452547073647, 0.8938643835182667}}; // Should be inverse of rotation
	
	double blockSize[3] = {0.0127 * BLOCKSCALE, 0.0127 * BLOCKSCALE, 0.00635 * BLOCKSCALE};
	
	//Shift the neutron position to "block coordinates"
	std::vector<double> shiftState(6);
	shiftState = {0.0,0.0,0.0, 0.0,0.0,0.0}; // Initialize at zero
	std::vector<double> shiftPrevState(6);
	shiftPrevState = {0.0,0.0,0.0, 0.0,0.0,0.0};

	// Matrix multiplication
	for (size_t ii = 0; ii < 3; ii++) {
		for (size_t jj = 0; jj < 3; jj++) {
			// Position shifted by block pos
			shiftState[ii] += rotation[ii][jj] * (state[jj] - blockPos[jj]);
			shiftPrevState[ii] += rotation[ii][jj] * (prevState[jj] - blockPos[jj]);
			// Velocity just rotated
			shiftState[ii+3] += rotation[ii][jj] * state[jj+3];
			shiftPrevState[ii+3] += rotation[ii][jj] * prevState[jj+3];
		}
	}

	//Check if the neutron is inside the block
	if ((fabs(shiftState[0]) > blockSize[0]) || (fabs(shiftState[1]) > blockSize[1]) || (fabs(shiftState[2]) > blockSize[2])) {	
		return false;
	} else {
		// We've hit the block, now we need to do reflection calculations
		std::vector<double> fracTravel(3);
		fracTravel = {0.0,0.0,0.0};
		for (size_t ii = 0; ii < 3; ii++) {
			// FracTravel depends on the sign of reflection
			if (shiftState[ii] >= 0) {
				fracTravel[ii] = (blockSize[ii] - shiftPrevState[ii])
									/ fabs(shiftState[ii] - shiftPrevState[ii]);
			} else {
				fracTravel[ii] = - (blockSize[ii] + shiftPrevState[ii])
									/ fabs(shiftState[ii] - shiftPrevState[ii]);
			}
			// Catch to make sure we don't have something reflecting inside the block
			if (fabs(fracTravel[ii]) > 1.0) { fracTravel[ii] = 0.0; }
		}
		
		// If we're not upscattering, reflect off the block
		// Hardcoding in different planes of reflection
		if ((fabs(fracTravel[2]) > fabs(fracTravel[0])) && (fabs(fracTravel[2]) > fabs(fracTravel[1]))) {
			//printf("Reflection in z!");
			ePerp = shiftPrevState[5]*shiftPrevState[5] / (2.0*MASS_N);
			if (fracTravel[2] < 0) {
				reflect(shiftPrevState, {0.0, 0.0, 1.0}, {1.0, 0.0, 0.0});
			} else {
				//printf("Managed to sneak under block!\n");
				reflect(shiftPrevState, {0.0, 0.0, 1.0}, {1.0, 0.0, 0.0});
				
				// Matrix multiplication
				for (size_t ii = 0; ii < 3; ii++) {
					for (size_t jj = 0; jj < 3; jj++) {
						// Position shifted by block pos
						state[ii] += invRotation[ii][jj] * (shiftPrevState[jj] + blockPos[jj]);
						// Velocity just rotated
						state[ii+3] += invRotation[ii][jj] * shiftPrevState[jj+3];
					}
				}
				return false; // Skip this since it's an unphysical case
			}
			
		} else if ((fabs(fracTravel[0]) > fabs(fracTravel[1])) && (fabs(fracTravel[0]) > fabs(fracTravel[2]))) {
			//printf("Reflection in x!\n");
			ePerp = shiftPrevState[3]*shiftPrevState[3]/(2.0*MASS_N);
			
			if (fracTravel[0] >= 0) {
				reflect(shiftPrevState, {-1.0, 0.0, 0.0}, {0.0, 0.0, 1.0}); // - x
			} else {
				reflect(shiftPrevState, {1.0, 0.0, 0.0}, {0.0, 0.0, 1.0}); // + x
			}
		} else if ((fabs(fracTravel[1]) > fabs(fracTravel[0])) && (fabs(fracTravel[1]) > fabs(fracTravel[2]))) {
			//printf("Reflection in y!\n");
			ePerp = shiftPrevState[4]*shiftPrevState[4]/(2.0*MASS_N);
			
			if (fracTravel[1] >= 0) {
				reflect(shiftPrevState, {0.0, -1.0, 0.0}, {0.0, 0.0, 1.0}); // - y		
			} else {
				reflect(shiftPrevState, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}); // + y
			}
		} else { 
			//printf("Boo!\n");
			return false;
		}
		
		// Matrix multiplication to return original state
		for (size_t ii = 0; ii < 3; ii++) {
			for (size_t jj = 0; jj < 3; jj++) {
				// Position shifted by block pos
				state[ii] += invRotation[ii][jj] * (shiftPrevState[jj] + blockPos[jj]);
				// Velocity just rotated
				state[ii+3] += invRotation[ii][jj] * shiftPrevState[jj+3];
			}
		}
		return true;
	}
}

std::vector<double> convertToTrap(std::vector<double> state) {
	// Convert from cartesian to trap coordinates
	
	std::vector<double> trapCoords;
	
	// Copy state vector to easier to read things
	double x  = state[0];
	double y  = state[1];
	double z  = state[2];
	double px = state[3];
	double py = state[4];
	double pz = state[5];
		
	double R, r; // get major and minor axes radii
	if (x > 0) {
		R = 1.0;
		r = 0.5;
	} else {
		R = 0.5;
		r = 1.0;
	}
	double rho = sqrt(y*y+z*z);
	
	// Spacial coordinates
	double zeta = r - sqrt(x*x + pow(rho - R,2)); // normal to trap
	double eta  = r - atan(x / (rho - R)); // tangent along holding field
	double xi   = (r + R) * atan(y / 2); // other tangent
	
	trapCoords.push_back(zeta); // Load vector
	trapCoords.push_back(eta); 
	trapCoords.push_back(xi);
	
	// Momentum coordinates
	double zetaDot = -2*(x*px) - (y*py + z*pz) / rho; // normal to trap
	double etaDot  = (1 / (R*R - 2*R*rho + x*x + rho*rho))*((R-rho)*px + (x/rho)*(y*py + z*pz));
	double xiDot   = 2*(r+R)*py / (y*y + 4);
																
	trapCoords.push_back(zetaDot);
	trapCoords.push_back(etaDot);
	trapCoords.push_back(xiDot);
	
	return trapCoords;
		
}

bool checkMagReflection(std::vector<double> s1, std::vector<double> s2) {
	// Want to see if a UCN reflects off the trap
	
	double zetaDot1 = - 2*(s1[0]*s1[3]) - (s1[1]*s1[4] + s1[2]*s1[5])
										/ sqrt(s1[1]*s1[1] + s1[2]*s1[2]) ; // Momentum in zeta direction (1)
	double zetaDot2 = - 2*(s2[0]*s2[3]) - (s2[1]*s2[4] + s2[2]*s2[5])
										/ sqrt(s2[1]*s2[1] + s2[2]*s2[2]) ; // Momentum in zeta direction (2)
	
	if ((std::signbit(zetaDot1) == 1) && (std::signbit(zetaDot2) == 0)) { // If zd1 is negative and zd2 is positive, we've flipped!
		return true;
	} else {
		return false;
	}
		
}

bool hitTD(std::vector<double> s1, std::vector<double> s2) {
	// Want to regenerate if there's a mag reflection off the trap door.
	
	double tdMin = -0.15;
	double tdMax = 0.15;
	
	if ((( tdMin < s2[0]) && (tdMax > s2[0]))
		&& (( tdMin < s2[1]) && (tdMax > s2[1]))) { // Check if we're in the td region
		
		if (checkMagReflection(s1,s2)) { return true; }
		
	}
	return false;
			
}

bool checkHouseHitLow(double x, double y, double z, double zOff) {
	// Check hits off the dagger housing
	
	if(z >= (-1.5 + zOff + 0.2) && z < (-1.5 + zOff + 0.2 + 0.14478) && 
		fabs(x + 0.1524) < (0.40 + 2.0179*(z + 1.5 - zOff - 0.2))/2.0) {
		return true;
	}
	return false;
}

bool checkHouseHitHigh(double x, double y, double z, double zOff) {
	// Check hits off the dagger actuator
	
	if(z >= (-1.5 + zOff + 0.2 + 0.14478) && z < (-1.5 + zOff + 0.2 + 0.2667) && 
		fabs(x + 0.1524) < 0.69215/2.0) {
		return true;
	}
	return false;
}

std::vector<double> initializeLyapState(std::vector<double> ref) {
	// Initialize lyuponev separation state
	
	double lenP = sqrt(ref[3]*ref[3] + ref[4]*ref[4] + ref[5]*ref[5]); // cons. of energy
	
	std::vector<double> pair(6);
	std::vector<double> randomP(3);
	std::vector<double> paraP(3);
	
	randomP[0] = nextU01(); // Random momentum
	randomP[1] = nextU01();
	randomP[2] = nextU01();
	normalize(randomP);
	
	paraP[0] = ref[3]; // Reference momentum
	paraP[1] = ref[4];
	paraP[2] = ref[5];	
	normalize(paraP);
	
	std::vector<double> perpP = cross(paraP, randomP); // Perpendicular
	normalize(perpP);
	
	double alpha = EPSILON*EPSILON*PSCALE*PSCALE/(lenP*2);
	double beta = sqrt(EPSILON*EPSILON*PSCALE*PSCALE - (EPSILON*EPSILON*PSCALE*PSCALE/(2*lenP))*(EPSILON*EPSILON*PSCALE*PSCALE/(2*lenP)));
	
	pair[0] = ref[0];
	pair[1] = ref[1];
	pair[2] = ref[2];
	
	pair[3] = ref[3] - alpha*paraP[0] + beta*perpP[0];
	pair[4] = ref[4] - alpha*paraP[1] + beta*perpP[1];
	pair[5] = ref[5] - alpha*paraP[2] + beta*perpP[2];
	
	return pair;
}

void resetStates(std::vector<double> ref, std::vector<double> &pair) {
	std::vector<double> dist(6);
	dist[0] = pair[0]-ref[0];
	dist[1] = pair[1]-ref[1];
	dist[2] = pair[2]-ref[2];
	dist[3] = pair[3]-ref[3];
	dist[4] = pair[4]-ref[4];
	dist[5] = pair[5]-ref[5];
	
	double scaling = EPSILON/distance(ref, pair);
	
	pair[0] = ref[0] + scaling*dist[0];
	pair[1] = ref[1] + scaling*dist[1];
	pair[2] = ref[2] + scaling*dist[2];
	pair[3] = ref[3] + scaling*dist[3];
	pair[4] = ref[4] + scaling*dist[4];
	pair[5] = ref[5] + scaling*dist[5];
}

double distance(std::vector<double> ref, std::vector<double> pair) {
	double dist = 0;
	dist += (pair[0]-ref[0])*(pair[0]-ref[0])/(XSCALE*XSCALE);
	dist += (pair[1]-ref[1])*(pair[1]-ref[1])/(YSCALE*YSCALE);
	dist += (pair[2]-ref[2])*(pair[2]-ref[2])/(ZSCALE*ZSCALE);
	dist += (pair[3]-ref[3])*(pair[3]-ref[3])/(PSCALE*PSCALE);
	dist += (pair[4]-ref[4])*(pair[4]-ref[4])/(PSCALE*PSCALE);
	dist += (pair[5]-ref[5])*(pair[5]-ref[5])/(PSCALE*PSCALE);
	
	dist = sqrt(dist);
	return dist;
}

/*void shiftTD(double *x, double *y, double *z, double t) {
	// Shift UCN hitting the trap door up by some amount
	
    int iLow = (int)(t/SAMPDT);
    double frac = (t - iLow*SAMPDT)/(SAMPDT);
    int iHi = iLow + 1;
    iLow = (iLow % tr->num + tr->num) % tr->num;
    iHi = (iHi % tr->num + tr->num) % tr->num;

    *x = *x + HEATMULT*(tr->x[iLow] + frac*(tr->x[iHi] - tr->x[iLow]));
    *y = *y + HEATMULT*(tr->y[iLow] + frac*(tr->y[iHi] - tr->y[iLow]));
    *z = *z + HEATMULT*(tr->z[iLow] + frac*(tr->z[iHi] - tr->z[iLow]));
}*/

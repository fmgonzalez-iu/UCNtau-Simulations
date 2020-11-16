#pragma once
#ifndef GENERATEWEIGHTING_H
#define GENERATEWEIGHTING_H

#define _USE_MATH_DEFINES
//#include <vector>
//#include <gsl/gsl_spline.h>
//#include "TH1D.h"
//#include "TF1.h"
//#include "TROOT.h"

//#include <numeric>
//#include <cmath>
//#include <cstring>
//#include <complex>
//#include <assert.h>
#include "./surface_multilayer.hpp"
extern "C" {
	#include "../../general/setup.h"
	#include "../../general/typedefs.h"
}

/*----------------------------------------------------------------------
 * Include file for generate_weighting.cpp
 *--------------------------------------------------------------------*/ 
double weightE(double eventE, double ePower, double threshold, double maxE); 
double weightTh(double eventTh, double cosPower, double threshold, double maxTh);
double calcWeight(double eventE, double eventTh, double ePower, double cosPower, double threshold);

std::vector<weightedBin> createHistQuantMultilayer(double thickOxide, double thickBoron, double threshold, double enePower, double cosPower, std::vector<noabsResult>& events,std::vector<double> ref, double* randU01s, double* randDeathTimes); 
std::vector<weightedBin> createHistQuantSpline(double thickOxide, double thickBoron, double threshold, double enePower, double cosPower, std::vector<noabsResult>& events,std::vector<double> ref, double* randU01s, double* randDeathTimes);
TH1D* createHistQuantSplineRoot(double thickOxide, double thickBoron, double threshold, double enePower, double cosPower, std::vector<noabsResult>& events,std::vector<double> ref, double* randU01s, double* randDeathTimes);
std::vector<measurement> fitTCRoot(TH1D* hist, std::vector<double> fitOffsets, std::vector<double> fitEnds);
std::vector<weightedBin> createHistQuantSpline_2det(
		double thickOxide1, double thickBoron1, // Assume both detectors have oxide + boron 
		double thickOxide2, double thickBoron2,
		double thresh, double power, double cosPower,
		std::vector<noabsCleanDagResult> &eventsD, std::vector<noabsCleanResult> &eventsC, 
		int refSize, double* randU01s, double* randDeathTimes);
		
std::vector<weightedBin> createHistQuantSpline_2det(
		double thickOxide1, double thickBoron1, // Assume both detectors have oxide + boron 
		double thickOxide2, double thickBoron2,
		double thresh, double power, double cosPower,
		std::vector<noabsCleanDagResult> &eventsD, std::vector<noabsCleanResult> &eventsC, 
		int refSize, double* randU01s, double* randDeathTimes);
		
//std::vector<weightedBin> createHistQuantMultilayerEdE(double thickOxide, double thickBoron, double threshold, std::vector<evt>& events, double* randU01s, double* randDeathTimes);
//std::vector<weightedBin> createHistQuantMultilayerEdESpline(double thickOxide, double thickBoron, double threshold, std::vector<evt>& events, double* randU01s, double* randDeathTimes);
//std::vector<weightedBin> createHistQuantMultilayerEPowdESpline(double thickOxide, double thickBoron, double threshold, double power, std::vector<evt>& events, double* randU01s, double* randDeathTimes);
//std::vector<weightedBin> createHistQuantNoOxEPowdEThetaSpline(double thickBoron, double threshold, double power, double cosPower, std::vector<evt>& events, double* randU01s, double* randDeathTimes);
//std::vector<weightedBin> createHistQuantNoOxEPowdEThetaSpline_C(double thickBoron, double threshold, double power, double cosPower, std::vector<evt>& events, double* randU01s, double* randDeathTimes);

//TH1D* createHistQuantMultilayerEdESplineRoot(double thickOxide, double thickBoron, double threshold, std::vector<evt>& events, double* randU01s, double* randDeathTimes);
//void fitTCQuantMultilayerEdESplineRoot(double thickOxide, double thickBoron, double threshold, std::vector<evt>& events, double* randU01s, double* randDeathTimes, double* tcs, double* tcsErr);

#endif /* GENERATEWEIGHTING_H */

#pragma once
#ifndef DETPLOTTING_H
#define DETPLOTTING_H

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

#include "TH1D.h"
#include "TH2D.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TGraph2D.h"

#include "./surface_multilayer.hpp"
extern "C" {
	#include "../../general/setup.h"
	#include "../../general/typedefs.h"
}

/*----------------------------------------------------------------------
 * Include file for generate_weighting.cpp
 *--------------------------------------------------------------------*/ 
TH2D* plot2DTimeEnergy(std::vector<reducedEvent> results, double minT, double maxT, double minE, double maxE);

TH1D* plotTimeHist(std::vector<reducedEvent> results, double minT, double maxT);
TH1D* plotEneHist(std::vector<reducedEvent> results, double minE, double maxE);
TH1D* plotThetaHist(std::vector<reducedEvent> results, double minTh, double maxTh);
TH1D* plotPercentDag(std::vector<weightedBin> energyD1, std::vector<weightedBin> energyD2);
TH1D* dataMCComp(std::vector<double>& dataVector, std::vector<weightedBin>& mcVector);

std::vector<measurement> calculateDrainTimes(std::vector<double> timeVector,int nBins);

TCanvas* plotChiMap(std::vector<std::vector<double>> variables, std::vector<double> chi2);
TCanvas* plotSanityCheck(std::vector<reducedEvent> results);
TCanvas* plotDetectorPos_2D(std::vector<double> xList, std::vector<double> yList, std::vector<double> zList);
TCanvas* plotMapTot(std::vector<reducedEvent> results, double minT, double maxT, double minE, double maxE);
TCanvas* plotMapSlice(std::vector<reducedEvent> results, int nSlices, double minE, double maxE);



/*(double weightE(double eventE, double ePower, double threshold, double maxE); 
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
*/	
//std::vector<weightedBin> createHistQuantMultilayerEdE(double thickOxide, double thickBoron, double threshold, std::vector<evt>& events, double* randU01s, double* randDeathTimes);
//std::vector<weightedBin> createHistQuantMultilayerEdESpline(double thickOxide, double thickBoron, double threshold, std::vector<evt>& events, double* randU01s, double* randDeathTimes);
//std::vector<weightedBin> createHistQuantMultilayerEPowdESpline(double thickOxide, double thickBoron, double threshold, double power, std::vector<evt>& events, double* randU01s, double* randDeathTimes);
//std::vector<weightedBin> createHistQuantNoOxEPowdEThetaSpline(double thickBoron, double threshold, double power, double cosPower, std::vector<evt>& events, double* randU01s, double* randDeathTimes);
//std::vector<weightedBin> createHistQuantNoOxEPowdEThetaSpline_C(double thickBoron, double threshold, double power, double cosPower, std::vector<evt>& events, double* randU01s, double* randDeathTimes);

//TH1D* createHistQuantMultilayerEdESplineRoot(double thickOxide, double thickBoron, double threshold, std::vector<evt>& events, double* randU01s, double* randDeathTimes);
//void fitTCQuantMultilayerEdESplineRoot(double thickOxide, double thickBoron, double threshold, std::vector<evt>& events, double* randU01s, double* randDeathTimes, double* tcs, double* tcsErr);

#endif /* DETPLOTTING_H */

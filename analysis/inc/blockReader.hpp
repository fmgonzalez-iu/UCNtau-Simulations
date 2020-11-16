#ifndef BLOCKREADER_H
#define BLOCKREADER_H

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <numeric>
#define _USE_MATH_DEFINES

#include <cmath>
#include <cstring>
#include <complex>
#include <gsl/gsl_spline.h>
#include <cstdlib>
#include "TH1D.h"
#include "TF1.h"
#include "TH2D.h"
#include "TAxis.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TApplication.h"
#include "TGraph2D.h"
#include "TMultiGraph.h"
#include "TPolyMarker3D.h"
#include "TObjArray.h"
#include "TFile.h"

extern "C" {
    #include "xorshift.h"
}

std::vector<double> refHist = {276.0,2147.0,6357.0,6724.0,6445.0,6515.0,6303.0,5883.0,5937.0,5750.0,5862.0,5741.0,5499.0,5298.0,5118.0,5336.0,5208.0,4989.0,4914.0,4840.0,4773.0,4761.0,4603.0,4670.0,4585.0,4490.0,4504.0,4344.0,4234.0,4190.0,4215.0,4159.0,4099.0,4054.0,4042.0,3958.0,3875.0,3758.0,3770.0,6777.0,19591.0,26111.0,24937.0,23352.0,21985.0,21214.0,20137.0,19300.0,18771.0,17531.0,16790.0,16027.0,15480.0,15079.0,14571.0,13997.0,13510.0,13170.0,12590.0,18787.0,32668.0,30531.0,29039.0,27367.0,25942.0,24249.0,22440.0,21718.0,20811.0,19545.0,18862.0,17545.0,16773.0,15964.0,15590.0,14640.0,14114.0,13309.0,12663.0,19992.0,29095.0,27203.0,25432.0,24246.0,22596.0,21316.0,20074.0,18979.0,17967.0,16720.0,16213.0,15231.0,14874.0,13783.0,13126.0,12620.0,12152.0,11748.0,11099.0,19469.0,31086.0,29762.0,28316.0,26476.0,24893.0,23408.0,22080.0,21091.0,19686.0,18667.0,17434.0,16717.0,15647.0,15110.0,14155.0,13375.0,12639.0,12072.0,11831.0,18263.0,23110.0,22100.0,21373.0,20302.0,19436.0,18349.0,17475.0,16732.0,15687.0,14945.0,14201.0,13216.0,12934.0,12358.0,11613.0,11052.0,10341.0,9803.0,9612.0,15876.0,20286.0,19403.0,17943.0,16981.0,15999.0,14919.0,13821.0,12959.0,11977.0,11107.0,10590.0,9923.0,9119.0,8469.0,8032.0,7199.0,6634.0,6475.0,5908.0,8792.0,10150.0,8870.0,7618.0,6546.0,5978.0,5135.0,4541.0,3926.0,3486.0,3176.0,2779.0,2399.0,2159.0,1892.0,1670.0,1493.0,1387.0,1236.0,1038.0,959.0,787.0,739.0,663.0,589.0};

typedef struct weightedBin {
    double wgt;
    double wgtSqr;
    double num;
} weightedBin;



typedef struct dag_evt {
	double time;
    double energy;
    double perp_energy;
    double pos[3];
    double z_offset;
    double theta;
    int nHits[4];
} dag_evt;


typedef struct block_evt {
	double time;
	double energy;
	int nHits[4];
	double position[3];
	double theta;
} block_evt;

typedef struct block_evt {
	double hit_time;
	double energy_b;
	double ey;
	double bx;
	double by;
	double bz; 
	
} block_evt;

#endif /* BLOCKREADER_H */

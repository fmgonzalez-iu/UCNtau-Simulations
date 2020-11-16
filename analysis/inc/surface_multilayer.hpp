#pragma once
#ifndef SURFACEMULTILAYER_H
#define SURFACEMULTILAYER_H

#define _USE_MATH_DEFINES
#include <numeric>
#include <cmath>
//#include <cstring>
#include <complex>

//#include "TH1D.h"
//#include "TF1.h"
//#include "TROOT.h"
//#include <assert.h>

#include "./quant_refl.hpp"

extern "C" {
	#include "../../general/setup.h"
	//#include "./xorshift.h"
//	#include "../../general/typedefs.h"
}
//#include "../inc/chisq_spectrum.hpp"
//#include "../inc/quant_refl.hpp"

/*----------------------------------------------------------------------
 * Include file for surfaceMultilayer.cpp
 *--------------------------------------------------------------------*/

double absorbProbQuantNoOxide(double ePerp, double thickBoron);
double absorbProbQuantOxide(double ePerp, double thickOxide, double thickBoron);
bool absorbMultilayer(double ePerp, double u, double thickOxide, double thickBoron);
bool absorbSpline(double ePerp, double u, gsl_spline *spline, gsl_interp_accel *acc);
void createSplineQuantOxide(double thickOxide, double thickBoron, gsl_spline **spline, gsl_interp_accel **acc);

#endif /* SURFACEMULTILAYER_H */

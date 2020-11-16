#pragma once
#ifndef QUANTREFL_H
#define QUANTREFL_H

#include <numeric>
#define _USE_MATH_DEFINES
#include <cmath>
//#include <cstring>
#include <complex>
#include <vector>
#include <gsl/gsl_spline.h>
//#include <gsl/gsl_spline.h>
#include "TH1D.h"
#include "TF1.h"
#include "TROOT.h"
//#include <assert.h>

//#include "../inc/chisq_spectrum.hpp"
extern "C" {
	#include "../../general/setup.h"
	#include "./xorshift.h"
}

/*----------------------------------------------------------------------
 * Include file for quant_refl.cpp
 *--------------------------------------------------------------------*/

std::vector<std::complex<double>> matmul(std::vector<std::complex<double>> a, std::vector<std::complex<double>> b);
std::complex<double> k(double ePerp, std::complex<double> u);
std::complex<double> gamma(std::complex<double> kn, std::complex<double> knm1);
std::vector<std::complex<double>> m(std::complex<double> kn, std::complex<double> knm1, double z);

#endif /* QUANTREFL_H */

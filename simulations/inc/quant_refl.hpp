#ifndef QUANT_REFL_H
#define QUANT_REFL_H
#pragma once
#include <vector>
#include <complex>
#include "../../general/typedefs.h"

std::vector<std::complex<double>> matmul(std::vector<std::complex<double>> a, std::vector<std::complex<double>> b);
std::complex<double> k(double ePerp, std::complex<double> u);
std::complex<double> gamma(std::complex<double> kn, std::complex<double> knm1);
std::vector<std::complex<double>> m(std::complex<double> kn, std::complex<double> knm1, double z);
double absorbProbQuantOxide(double ePerp, double thickBoron);
bool absorbMultilayer(double ePerp, double thickBoron, double x, double y, double z, double zOff);

#endif /* QUANT_REFL_H */

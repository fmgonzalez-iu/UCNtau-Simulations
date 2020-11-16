#include "../inc/quant_refl.hpp"

/*----------------------------------------------------------------------
 * This contains quantum mechanical UCN reflection functions. Assuming 
 * a UCN-matter interaction acts as a particle in a 1D potential, we can
 * solve the transmission and upscattering probabilities.
 * 
 * Additionally there's some matrix algebra required for these 
 * calculations to work -- they're included here too.
 *--------------------------------------------------------------------*/

std::vector<std::complex<double>> matmul(std::vector<std::complex<double>> a, std::vector<std::complex<double>> b) {
	// 2D matrix multiplication
	
	// Initialize a residual matrix (0% probability of transmission
    std::vector<std::complex<double>> res = {std::complex<double>(0,0), std::complex<double>(0,0),
											 std::complex<double>(0,0), std::complex<double>(0,0)};
    if(a.size() != 4 || b.size() != 4) { // Check that both vectors have the right size
        return res;
    }
    
    res[0] = a[0]*b[0] + a[1]*b[2]; // Matrix multiplication (2*2)
    res[1] = a[0]*b[1] + a[1]*b[3];
    res[2] = a[2]*b[0] + a[3]*b[2];
    res[3] = a[2]*b[1] + a[3]*b[3];
    
    return res;
}

std::complex<double> k(double ePerp, std::complex<double> u) {
	// k is complex wavevector of a particle with given ePerp on potential u
    return std::sqrt((2*MASS_N/(HBAR*HBAR))*(ePerp - u));
}

std::complex<double> gamma(std::complex<double> kn, std::complex<double> knm1) {
	// Ratio of two complex numbers
    return knm1/kn;
}

std::vector<std::complex<double>> m(std::complex<double> kn, std::complex<double> knm1, double z) {
	// Transmission matrix M for wavevectors -- see Golub for derivation
	
    std::vector<std::complex<double>> res = {std::complex<double>(0,0), std::complex<double>(0,0), 
											 std::complex<double>(0,0), std::complex<double>(0,0)};
    
    res[0] = (1.0/2.0)*(1.0 + gamma(kn,knm1))*std::exp(std::complex<double>(0,1)*(knm1-kn)*z);
    res[1] = (1.0/2.0)*(1.0 - gamma(kn,knm1))*std::exp(-std::complex<double>(0,1)*(knm1+kn)*z);
    res[2] = (1.0/2.0)*(1.0 - gamma(kn,knm1))*std::exp(std::complex<double>(0,1)*(knm1+kn)*z);
    res[3] = (1.0/2.0)*(1.0 + gamma(kn,knm1))*std::exp(-std::complex<double>(0,1)*(knm1-kn)*z);
    return res;
}

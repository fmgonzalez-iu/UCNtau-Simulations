#include "../inc/surface_multilayer.hpp"

/*----------------------------------------------------------------------
 * These functions describe the absorption probability of a UCN reflecting
 * off a given material with a thickness of Oxide and Boron.
 * 
 * These functions are energy dependent, acting as a 1D particle in a 
 * finite potential.
 *--------------------------------------------------------------------*/ 

//----------------------------------------------------------------------
// Absorption multilayers
//----------------------------------------------------------------------
double absorbProbQuantNoOxide(double ePerp, double thickBoron) {
	// Calculation of absorption prob for boron coated ZnS
	
	const double vboron = (2*M_PI*(HBAR*HBAR)/MASS_N)*ABORON*NBORON; // Define real/imaginary potentials
	const double wboron = (HBAR/2)*NBORON*SIGMABORON;
	const double vzns = (2*M_PI*(HBAR*HBAR)/MASS_N)*AZINC*NZINC + (2*M_PI*(HBAR*HBAR)/MASS_N)*ASULFUR*NSULFUR;
	const double wzns = (HBAR/2)*NZINC*SIGMAZINC + (HBAR/2)*NSULFUR*SIGMASULFUR;
    
	std::vector<std::complex<double>> pots = {std::complex<double>(0, 0), // Create a vector of potentials
											  std::complex<double>(vboron, -wboron),
											  std::complex<double>(vzns, -wzns)};
	std::vector<std::complex<double>> mbar = {std::complex<double>(1,0), std::complex<double>(0,0), std::complex<double>(0,0), std::complex<double>(1,0)};
	std::vector<double> zs = {0.0, thickBoron*1e-9, 10000e-9};
	
	for(int i = pots.size()-1; i > 0; i--) { // Loop through potentials and calculate probabilities
		mbar = matmul(mbar, m(k(ePerp, pots[i]), k(ePerp, pots[i-1]), zs[i-1]));
	}
    
	return 1.0 - (std::conj(-mbar[2]/mbar[3])*-mbar[2]/mbar[3]).real();
}


double absorbProbQuantOxide(double ePerp, double thickOxide, double thickBoron) {
	// Calculation of absorption prob for boron coated ZnS with oxide layer

	const double voxide = (2*M_PI*(HBAR*HBAR)/MASS_N)*ABORON*NBORONB2O3 + (2*M_PI*(HBAR*HBAR)/MASS_N)*AOXYGEN*NOXYGENB2O3;  // Define real/imaginary potentials
	const double woxide = (HBAR/2)*NBORONB2O3*SIGMABORON + (HBAR/2)*NOXYGENB2O3*SIGMAOXYGEN;
	const double vboron = (2*M_PI*(HBAR*HBAR)/MASS_N)*ABORON*NBORON;
	const double wboron = (HBAR/2)*NBORON*SIGMABORON;
	const double vzns = (2*M_PI*(HBAR*HBAR)/MASS_N)*AZINC*NZINC + (2*M_PI*(HBAR*HBAR)/MASS_N)*ASULFUR*NSULFUR;
	const double wzns = (HBAR/2)*NZINC*SIGMAZINC + (HBAR/2)*NSULFUR*SIGMASULFUR;

	std::vector<std::complex<double>> pots = {std::complex<double>(0, 0), // Create a vector of potentials
											  //std::complex<double>(vcarbon, -wcarbon),
											  std::complex<double>(voxide, -woxide),
											  std::complex<double>(vboron, -wboron),
											  std::complex<double>(vzns, -wzns)};

	std::vector<std::complex<double>> mbar = {std::complex<double>(1,0), std::complex<double>(0,0), std::complex<double>(0,0), std::complex<double>(1,0)};
	std::vector<double> zs = {0.0, thickOxide*1e-9, thickOxide*1e-9 + thickBoron*1e-9, 10000e-9};

	for(int i = pots.size()-1; i > 0; i--) { // Loop through potentials and calculate probabilities
		mbar = matmul(mbar, m(k(ePerp, pots[i]), k(ePerp, pots[i-1]), zs[i-1]));
	}
    
	return 1.0 - (std::conj(-mbar[2]/mbar[3])*-mbar[2]/mbar[3]).real();
}

bool absorbMultilayer(double ePerp, double u, double thickOxide, double thickBoron) {
	// Check if particle was absorbed
	#ifdef USEOXIDE
	double abs = absorbProbQuantOxide(ePerp, thickOxide, thickBoron);
	#else
	double abs = absorbProbQuantNoOxide(ePerp, thickBoron);
	#endif
    if(u < abs) { // and check!
        return true;
    }
    return false; 
}

bool absorbSpline(double ePerp, double u, gsl_spline *spline, gsl_interp_accel *acc) {
	// Evaluate spline to check if particle was absorbed
	double abs = gsl_spline_eval(spline, ePerp, acc);
	if(u < abs) { // and check!
        return true;
    }
    return false;
}

void createSplineQuantOxide(double thickOxide, double thickBoron, gsl_spline **spline, gsl_interp_accel **acc) {
	// Spline plot for fitting the oxide and boron layers
	
    *acc = gsl_interp_accel_alloc(); // gnu spline library
    *spline = gsl_spline_alloc(gsl_interp_akima, 1001); // Use Akima spline, 1001 points
    
    double xs[1001]; 
    double ys[1001];
    xs[0] = 0.0;
    ys[0] = 0.0;
    for(int i = 1; i < 1001; i++) { // Generate positions
        xs[i] = 0.0 + EMAX / JTONEV * (i/1000.0) * (i/1000.0); // X positions
        ys[i] = absorbProbQuantOxide(xs[i], thickOxide, thickBoron); // Y positions are probability
    }
    gsl_spline_init(*spline, xs, ys, 1001); // Interpolate spline
}

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
#include "TH1D.h"
#include "TF1.h"
#include "TROOT.h"

extern "C" {
    #include "xorshift.h"
}

#define START (300+20+41)
#define END (300+20+41+184)
//#define START 0
//#define END 1000
//#define START 41
//#define END 225

#define NRECORDS 50

#define MASS_N 1.674927471e-27
#define GRAV 9.80665e0
#define JTONEV 6.2415091e27
#define HBAR 1.054571800e-34

#define NBORONB2O3 8.824e27
#define NBORON 1.37e29
#define ABORON -0.1e-15
#define SIGMABORON 2200*3.835e-25

#define NOXYGENB2O3 1.32e28
#define AOXYGEN 5.803e-15
#define SIGMAOXYGEN 4.232e-28
//#define SIGMAOXYGEN 0.0

#define NCARBON 1.133e29
#define ACARBON 6.6460e-15
#define SIGMACARBON 5.551e-28
//#define SIGMACARBON 0.0

#define NZINC 2.527e28
#define AZINC 5.68e-15
#define SIGMAZINC 5.241e-28
//#define SIGMAZINC 0.0

#define NSULFUR 2.527e28
#define ASULFUR 2.847e-15
#define SIGMASULFUR 1.556e-28
//#define SIGMASULFUR 0.0

typedef struct evt {
    double energy;
    double theta;
    float times[NRECORDS];
    float ePerp[NRECORDS];
} evt;

typedef struct measurement {
    double val;
    double err;
} measurement;

std::vector<std::complex<double>> matmul(std::vector<std::complex<double>> a, std::vector<std::complex<double>> b) {
    std::vector<std::complex<double>> res = {std::complex<double>(0,0), std::complex<double>(0,0), std::complex<double>(0,0), std::complex<double>(0,0)};
    if(a.size() != 4 || b.size() != 4) {
        return res;
    }
    res[0] = a[0]*b[0] + a[1]*b[2];
    res[1] = a[0]*b[1] + a[1]*b[3];
    res[2] = a[2]*b[0] + a[3]*b[2];
    res[3] = a[2]*b[1] + a[3]*b[3];
    return res;
}

std::complex<double> k(double ePerp, std::complex<double> u) {
    return std::sqrt((2*MASS_N/(HBAR*HBAR))*(ePerp - u));
}

std::complex<double> gamma(std::complex<double> kn, std::complex<double> knm1) {
    return knm1/kn;
}

std::vector<std::complex<double>> m(std::complex<double> kn, std::complex<double> knm1, double z) {
    std::vector<std::complex<double>> res = {std::complex<double>(0,0), std::complex<double>(0,0), std::complex<double>(0,0), std::complex<double>(0,0)};
    res[0] = (1.0/2.0)*(1.0 + gamma(kn,knm1))*std::exp(std::complex<double>(0,1)*(knm1-kn)*z);
    res[1] = (1.0/2.0)*(1.0 - gamma(kn,knm1))*std::exp(-std::complex<double>(0,1)*(knm1+kn)*z);
    res[2] = (1.0/2.0)*(1.0 - gamma(kn,knm1))*std::exp(std::complex<double>(0,1)*(knm1+kn)*z);
    res[3] = (1.0/2.0)*(1.0 + gamma(kn,knm1))*std::exp(-std::complex<double>(0,1)*(knm1-kn)*z);
    return res;
}

double absorbProbQuantOxide(double ePerp, double thickOxide, double thickBoron) {
//    const double voxide = (2*M_PI*(HBAR*HBAR)/MASS_N)*ABORON*NBORONB2O3 + (2*M_PI*(HBAR*HBAR)/MASS_N)*AOXYGEN*NOXYGENB2O3;
//    const double woxide = (HBAR/2)*NBORONB2O3*SIGMABORON + (HBAR/2)*NOXYGENB2O3*SIGMAOXYGEN;
    const double vboron = (2*M_PI*(HBAR*HBAR)/MASS_N)*ABORON*NBORON;
    const double wboron = (HBAR/2)*NBORON*SIGMABORON;
    const double vzns = (2*M_PI*(HBAR*HBAR)/MASS_N)*AZINC*NZINC + (2*M_PI*(HBAR*HBAR)/MASS_N)*ASULFUR*NSULFUR;
    const double wzns = (HBAR/2)*NZINC*SIGMAZINC + (HBAR/2)*NSULFUR*SIGMASULFUR;
    
    std::vector<std::complex<double>> pots = {std::complex<double>(0, 0),
//                                              std::complex<double>(vcarbon, -wcarbon),
//                                              std::complex<double>(voxide, -woxide),
                                              std::complex<double>(vboron, -wboron),
                                              std::complex<double>(vzns, -wzns)};
    std::vector<std::complex<double>> mbar = {std::complex<double>(1,0), std::complex<double>(0,0), std::complex<double>(0,0), std::complex<double>(1,0)};
//    std::vector<double> zs = {0.0, thickOxide*1e-9, thickOxide*1e-9 + thickBoron*1e-9, 10000e-9};
    std::vector<double> zs = {0.0, thickBoron*1e-9, 10000e-9};
    
    for(int i = pots.size()-1; i > 0; i--) {
        mbar = matmul(mbar, m(k(ePerp, pots[i]), k(ePerp, pots[i-1]), zs[i-1]));
    }
    
    return 1.0 - (std::conj(-mbar[2]/mbar[3])*-mbar[2]/mbar[3]).real();
}

bool absorbMultilayer(double ePerp, double u, double thickOxide, double thickBoron) {
    if(u < absorbProbQuantOxide(ePerp, thickOxide, thickBoron)) {
        return true;
    }
    return false; 
}

bool absorbSpline(double ePerp, double u, gsl_spline *spline, gsl_interp_accel *acc) {
    double abs = gsl_spline_eval(spline, ePerp, acc);
    if(u < abs) {
        return true;
    }
    return false;
}

void createSplineQuantOxide(double thickOxide, double thickBoron, gsl_spline **spline, gsl_interp_accel **acc) {
    *acc = gsl_interp_accel_alloc();
    *spline = gsl_spline_alloc(gsl_interp_akima, 1001);
    double xs[1001];
    double ys[1001];
    xs[0] = 0.0;
    ys[0] = 0.0;
    for(int i = 1; i < 1001; i++) {
        xs[i] = 0.0 + 50.0/JTONEV*(i/1000.0)*(i/1000.0);
        ys[i] = absorbProbQuantOxide(xs[i], thickOxide, thickBoron);
    }
    gsl_spline_init(*spline, xs, ys, 1001);
}

measurement sumCtsSpline(double start, double end, double threshold, double power, double cosPower, std::vector<evt>& events, double* randU01s, double* randDeathTimes, gsl_spline *spline, gsl_interp_accel *acc) {
    measurement meas = {0.0, 0.0};
    double sum = 0.0;
    double sumErr = 0.0;
    int numCount = 0;
    for(unsigned long i = 0; i < events.size(); i++) {
        double weight = 0.0;
        if(events[i].energy*JTONEV < threshold) {
            continue;
        }
        weight = pow(events[i].energy/(0.3444*GRAV*MASS_N), power-1)*pow(cos(events[i].theta), cosPower);
//        weight = (events[i].energy*JTONEV-threshold)/(34.5-threshold);

        for(int j = 0; j < NRECORDS; j++) {
            if(events[i].times[j] > randDeathTimes[i]) {
                break;
            }
            if(events[i].times[j] >= end) {
                break;
            }
//            if(absorbMultilayer(events[i].ePerp[j], randU01s[i*100 + j], thickOxide, thickBoron)) {
            if(absorbSpline(events[i].ePerp[j], randU01s[i*2*NRECORDS + j], spline, acc)) {
                if(events[i].times[j] > start) {
                    meas.val += weight;
                    meas.err += weight*weight;
                    numCount += 1;
                }
                break;
            }
        }
    }
    
    printf("%d\n", numCount);
    meas.err = sqrt(meas.err);

    return meas;
}

int main(int argc, char** argv) {
//    ROOT::EnableThreadSafety();
    initxorshift();
    //buff_len = 4 + 3*8 + 4
    const size_t buff_len = 4 + 2*8 + 2*NRECORDS*4 + 4;
    char* buf = new char[buff_len];
    if(argc != 2) {
        printf("Error! Usage: ./chisq_spectrum_fit fname\n");
    }

    std::ifstream binfile(argv[1], std::ios::in | std::ios::binary);
    if(!binfile.is_open()) {
        printf("Error! Could not open file %s\n", argv[1]);
        return 1;
    }
    
    std::vector<evt> events;
    evt event;
    while(!binfile.eof()) {
        binfile.read(buf, buff_len);
        if(binfile.eof()) { //Breaks on last read of file (i.e. when 0 bytes are read and EOF bit is set)
            break;
        }
        if(*((unsigned int *)&buf[0]) != 2*8 + 2*NRECORDS*4) {
            fprintf(stderr, "Error! Aliased read on\n");
            exit(2);
        }
        event.energy = *((double *)(&buf[0] + sizeof(unsigned int)));
        event.theta = *((double *)(&buf[0] + sizeof(unsigned int) + sizeof(double)));
        std::memcpy((void *)&event.times, (void *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(double)), NRECORDS*sizeof(float));
        std::memcpy((void *)&event.ePerp, (void *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(double) + NRECORDS*sizeof(float)), NRECORDS*sizeof(float));
        events.push_back(event);
    }
    binfile.close();
    
    printf("Read %lu Events!\n", events.size());
    
    double* randU01s = new double[events.size()*2*NRECORDS];
    double* randDeathTimes = new double[events.size()];
    
    for(unsigned long i = 0; i < events.size()*2*NRECORDS; i++) {
        randU01s[i] = nextU01();
    }
    for(unsigned long i = 0; i < events.size(); i++) {
        randDeathTimes[i] = -877.7*log(nextU01());
    }
    

    gsl_spline *spline;
    gsl_interp_accel *acc;
    createSplineQuantOxide(0.0, 4.4, &spline, &acc);
//    tStartDip = 500 + holdT + 20
//    tEndDip = 500 + holdT + 20 + 20
//    tEndInt = 500 + holdT + 20 + 100
    int holdt = 1000;
    measurement sumDip = sumCtsSpline(350+holdt+20, 350+holdt+20+20, 7, 1.4, 0.15, events, randU01s, randDeathTimes, spline, acc);
    measurement sumFull = sumCtsSpline(350+holdt+20, 350+holdt+20+100, 7, 1.4, 0.15, events, randU01s, randDeathTimes, spline, acc);
    printf("%f\n", sumDip.val/sumFull.val);
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
    

    delete[] randU01s;
    delete[] buf;

    return 0;
}

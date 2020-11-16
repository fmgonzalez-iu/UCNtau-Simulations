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

//std::vector<double> refHist = {24.0,276.0,2147.0,6357.0,6724.0,6445.0,6515.0,6303.0,5883.0,5937.0,5750.0,5862.0,5741.0,5499.0,5298.0,5118.0,5336.0,5208.0,4989.0,4914.0,4840.0,4773.0,4761.0,4603.0,4670.0,4585.0,4490.0,4504.0,4344.0,4234.0,4190.0,4215.0,4159.0,4099.0,4054.0,4042.0,3958.0,3875.0,3758.0,3770.0,6777.0,19591.0,26111.0,24937.0,23352.0,21985.0,21214.0,20137.0,19300.0,18771.0,17531.0,16790.0,16027.0,15480.0,15079.0,14571.0,13997.0,13510.0,13170.0,12590.0,18787.0,32668.0,30531.0,29039.0,27367.0,25942.0,24249.0,22440.0,21718.0,20811.0,19545.0,18862.0,17545.0,16773.0,15964.0,15590.0,14640.0,14114.0,13309.0,12663.0,19992.0,29095.0,27203.0,25432.0,24246.0,22596.0,21316.0,20074.0,18979.0,17967.0,16720.0,16213.0,15231.0,14874.0,13783.0,13126.0,12620.0,12152.0,11748.0,11099.0,19469.0,31086.0,29762.0,28316.0,26476.0,24893.0,23408.0,22080.0,21091.0,19686.0,18667.0,17434.0,16717.0,15647.0,15110.0,14155.0,13375.0,12639.0,12072.0,11831.0,18263.0,23110.0,22100.0,21373.0,20302.0,19436.0,18349.0,17475.0,16732.0,15687.0,14945.0,14201.0,13216.0,12934.0,12358.0,11613.0,11052.0,10341.0,9803.0,9612.0,15876.0,20286.0,19403.0,17943.0,16981.0,15999.0,14919.0,13821.0,12959.0,11977.0,11107.0,10590.0,9923.0,9119.0,8469.0,8032.0,7199.0,6634.0,6475.0,5908.0,8792.0,10150.0,8870.0,7618.0,6546.0,5978.0,5135.0,4541.0,3926.0,3486.0,3176.0,2779.0,2399.0,2159.0,1892.0,1670.0,1493.0,1387.0,1236.0,1038.0,959.0,787.0,739.0,663.0,589.0};
std::vector<double> refHist = {276.0,2147.0,6357.0,6724.0,6445.0,6515.0,6303.0,5883.0,5937.0,5750.0,5862.0,5741.0,5499.0,5298.0,5118.0,5336.0,5208.0,4989.0,4914.0,4840.0,4773.0,4761.0,4603.0,4670.0,4585.0,4490.0,4504.0,4344.0,4234.0,4190.0,4215.0,4159.0,4099.0,4054.0,4042.0,3958.0,3875.0,3758.0,3770.0,6777.0,19591.0,26111.0,24937.0,23352.0,21985.0,21214.0,20137.0,19300.0,18771.0,17531.0,16790.0,16027.0,15480.0,15079.0,14571.0,13997.0,13510.0,13170.0,12590.0,18787.0,32668.0,30531.0,29039.0,27367.0,25942.0,24249.0,22440.0,21718.0,20811.0,19545.0,18862.0,17545.0,16773.0,15964.0,15590.0,14640.0,14114.0,13309.0,12663.0,19992.0,29095.0,27203.0,25432.0,24246.0,22596.0,21316.0,20074.0,18979.0,17967.0,16720.0,16213.0,15231.0,14874.0,13783.0,13126.0,12620.0,12152.0,11748.0,11099.0,19469.0,31086.0,29762.0,28316.0,26476.0,24893.0,23408.0,22080.0,21091.0,19686.0,18667.0,17434.0,16717.0,15647.0,15110.0,14155.0,13375.0,12639.0,12072.0,11831.0,18263.0,23110.0,22100.0,21373.0,20302.0,19436.0,18349.0,17475.0,16732.0,15687.0,14945.0,14201.0,13216.0,12934.0,12358.0,11613.0,11052.0,10341.0,9803.0,9612.0,15876.0,20286.0,19403.0,17943.0,16981.0,15999.0,14919.0,13821.0,12959.0,11977.0,11107.0,10590.0,9923.0,9119.0,8469.0,8032.0,7199.0,6634.0,6475.0,5908.0,8792.0,10150.0,8870.0,7618.0,6546.0,5978.0,5135.0,4541.0,3926.0,3486.0,3176.0,2779.0,2399.0,2159.0,1892.0,1670.0,1493.0,1387.0,1236.0,1038.0,959.0,787.0,739.0,663.0,589.0};

typedef struct evt {
    double energy;
    float times[NRECORDS];
    float ePerp[NRECORDS];
} evt;

typedef struct weightedBin {
    double wgt;
    double wgtSqr;
    double num;
} weightedBin;

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
    const double voxide = (2*M_PI*(HBAR*HBAR)/MASS_N)*ABORON*NBORONB2O3 + (2*M_PI*(HBAR*HBAR)/MASS_N)*AOXYGEN*NOXYGENB2O3;
    const double woxide = (HBAR/2)*NBORONB2O3*SIGMABORON + (HBAR/2)*NOXYGENB2O3*SIGMAOXYGEN;
    const double vboron = (2*M_PI*(HBAR*HBAR)/MASS_N)*ABORON*NBORON;
    const double wboron = (HBAR/2)*NBORON*SIGMABORON;
    const double vzns = (2*M_PI*(HBAR*HBAR)/MASS_N)*AZINC*NZINC + (2*M_PI*(HBAR*HBAR)/MASS_N)*ASULFUR*NSULFUR;
    const double wzns = (HBAR/2)*NZINC*SIGMAZINC + (HBAR/2)*NSULFUR*SIGMASULFUR;
    
    std::vector<std::complex<double>> pots = {std::complex<double>(0, 0),
//                                              std::complex<double>(vcarbon, -wcarbon),
                                              std::complex<double>(voxide, -woxide),
                                              std::complex<double>(vboron, -wboron),
                                              std::complex<double>(vzns, -wzns)};
    std::vector<std::complex<double>> mbar = {std::complex<double>(1,0), std::complex<double>(0,0), std::complex<double>(0,0), std::complex<double>(1,0)};
    std::vector<double> zs = {0.0, thickOxide*1e-9, thickOxide*1e-9 + thickBoron*1e-9, 10000e-9};
    
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

std::vector<weightedBin> createHistQuantMultilayerEdE(double thickOxide, double thickBoron, double threshold, std::vector<evt>& events, double* randU01s, double* randDeathTimes) {
//void createHist(double absProb, double absCut, double threshold, double saturation, std::vector<evt>& events, double* randU01s) {
    std::vector<weightedBin> hist;
    weightedBin zero = {0.0, 0.0};
    hist.resize(184, zero);
    int numCount = 0;
    int numMiss = 0;
    for(unsigned long i = 0; i < events.size(); i++) {
        double weight = 0.0;
        if(events[i].energy*JTONEV < threshold) {
            continue;
        }
//        weight = events[i].energy*JTONEV/34.5;
        weight = (events[i].energy*JTONEV-threshold)/(34.5-threshold);

        for(int j = 0; j < NRECORDS; j++) {
            if(events[i].times[j] < 41) {
                continue;
            }
            if(events[i].times[j] - 41 > randDeathTimes[i]) {
                break;
            }
            if(events[i].times[j] >= 225) {
                break;
            }
            if(absorbMultilayer(events[i].ePerp[j], randU01s[i*2*NRECORDS + j], thickOxide, thickBoron)) {
//                if(events[i].times[j] - 40 > 160) {
//                    printf("%e\n", events[i].energy * JTONEV);
//                }
                if(int(events[i].times[j])-41 > 183) {
                    printf("Boo!\n");
                }
                hist[int(events[i].times[j])-41].wgt += weight;
                hist[int(events[i].times[j])-41].wgtSqr += weight*weight;
                hist[int(events[i].times[j])-41].num += 1;
                numCount += 1;
                break;
            }
            if(j == NRECORDS - 1) {
                numMiss += 1;
            }
        }
    }
    printf("%d %d\n", numCount, numMiss);
    return hist;
}

std::vector<weightedBin> createHistQuantMultilayerEdESpline(double thickOxide, double thickBoron, double threshold, std::vector<evt>& events, double* randU01s, double* randDeathTimes) {
    gsl_spline *spline;
    gsl_interp_accel *acc;
    createSplineQuantOxide(thickOxide, thickBoron, &spline, &acc);
    std::vector<weightedBin> hist;
    weightedBin zero = {0.0, 0.0};
    hist.resize(184, zero);
    int numCount = 0;
    for(unsigned long i = 0; i < events.size(); i++) {
        double weight = 0.0;
        if(events[i].energy*JTONEV < threshold) {
            continue;
        }
//        weight = events[i].energy*JTONEV/34.5;
        weight = (events[i].energy*JTONEV-threshold)/(34.5-threshold);

        for(int j = 0; j < NRECORDS; j++) {
            if(events[i].times[j] < 41) {
                continue;
            }
            if(events[i].times[j] - 41 > randDeathTimes[i]) {
                break;
            }
            if(events[i].times[j] >= 225) {
                break;
            }
//            if(absorbMultilayer(events[i].ePerp[j], randU01s[i*100 + j], thickOxide, thickBoron)) {
            if(absorbSpline(events[i].ePerp[j], randU01s[i*2*NRECORDS + j], spline, acc)) {
                if(int(events[i].times[j])-41 > 183) {
                    printf("Boo!\n");
                }
                hist[int(events[i].times[j])-41].wgt += weight;
                hist[int(events[i].times[j])-41].wgtSqr += weight*weight;
                hist[int(events[i].times[j])-41].num += 1;
                numCount += 1;
                break;
            }
        }
    }
    printf("%d\n", numCount);
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
    return hist;
}

TH1D* createHistQuantMultilayerEdESplineRoot(double thickOxide, double thickBoron, double threshold, std::vector<evt>& events, double* randU01s, double* randDeathTimes) {
    gsl_spline *spline;
    gsl_interp_accel *acc;
    createSplineQuantOxide(thickOxide, thickBoron, &spline, &acc);
    char hName[256];
    sprintf(hName, "MC%f%f%f", thickOxide, thickBoron, threshold);
    TH1D* hist = new TH1D(hName, "Monte Carlo", 184, 0, 184);
    for(unsigned long i = 0; i < events.size(); i++) {
        double weight = 0.0;
        if(events[i].energy*JTONEV < threshold) {
            continue;
        }
//        weight = events[i].energy*JTONEV/34.5;
        weight = (events[i].energy*JTONEV-threshold)/(34.5-threshold);

        for(int j = 0; j < NRECORDS; j++) {
            if(events[i].times[j] < 41) {
                continue;
            }
            if(events[i].times[j] - 41 > randDeathTimes[i]) {
                break;
            }
            if(events[i].times[j] >= 225) {
                break;
            }
//            if(absorbMultilayer(events[i].ePerp[j], randU01s[i*100 + j], thickOxide, thickBoron)) {
            if(absorbSpline(events[i].ePerp[j], randU01s[i*2*NRECORDS + j], spline, acc)) {
                if(int(events[i].times[j])-41 > 183) {
                    printf("Boo!\n");
                }
                hist->Fill(events[i].times[j]-41, weight);
                break;
            }
        }
    }
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
    return hist;
}

void fitTCQuantMultilayerEdESplineRoot(double thickOxide, double thickBoron, double threshold, std::vector<evt>& events, double* randU01s, double* randDeathTimes, double* tcs, double* tcsErr) {
    double fitOffsets[8] = {4, 42, 62, 82, 102, 122, 142, 162};
    double fitEnds[8] = {38, 58, 78, 98, 118, 138, 158, 178};
    gsl_spline *spline;
    gsl_interp_accel *acc;
    createSplineQuantOxide(thickOxide, thickBoron, &spline, &acc);
    char hName[256];
    sprintf(hName, "MC%f%f%f", thickOxide, thickBoron, threshold);
    TH1D* hist = new TH1D(hName, "Monte Carlo", 184, 0, 184);
    for(unsigned long i = 0; i < events.size(); i++) {
        double weight = 0.0;
        if(events[i].energy*JTONEV < threshold) {
            continue;
        }
//        weight = events[i].energy*JTONEV/34.5;
        weight = (events[i].energy*JTONEV-threshold)/(34.5-threshold);

        for(int j = 0; j < NRECORDS; j++) {
            if(events[i].times[j] < 41) {
                continue;
            }
            if(events[i].times[j] - 41 > randDeathTimes[i]) {
                break;
            }
            if(events[i].times[j] >= 225) {
                break;
            }
//            if(absorbMultilayer(events[i].ePerp[j], randU01s[i*100 + j], thickOxide, thickBoron)) {
            if(absorbSpline(events[i].ePerp[j], randU01s[i*2*NRECORDS + j], spline, acc)) {
                if(int(events[i].times[j])-41 > 183) {
                    printf("Boo!\n");
                }
                hist->Fill(events[i].times[j]-41, weight);
                break;
            }
        }
    }
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
    for(int i = 0; i < 8; i++) {
        TF1* fit = new TF1("mcFit", "[0]*exp(-(x-[2])/[1])", fitOffsets[i], fitEnds[i]);
        fit->SetParameter(2, fitOffsets[i]);
        fit->SetParLimits(2, fitOffsets[i], fitOffsets[i]);
        fit->SetParameter(0, 1000);
        fit->SetParameter(1, 20);
        hist->Fit(fit, "LQ", "", fitOffsets[i], fitEnds[i]);
        tcs[i] = fit->GetParameter(1);
        tcsErr[i] = fit->GetParError(1);
        delete fit;
    }
}

double calcChisq(std::vector<double>& hn, std::vector<weightedBin>& hm) {
    if(hn.size() != hm.size()) {
        return -1.0;
    }
    double cn = std::accumulate(hm.begin(), hm.end(), 0, [](double sum, weightedBin bin)->double{return sum + bin.wgt;});
    double cm = std::accumulate(hn.begin(), hn.end(), 0);
//    printf("%f %f\n", cn, cm);
    double chisqSum = 0.0;
    for(int i = 0; i < hn.size(); i++) {
        chisqSum += (1/(cn*cm))*pow(cn*hn[i] - cm*hm[i].wgt,2)/(hn[i]+hm[i].wgt);
    }
    return chisqSum;
}

double calcChisqNate(std::vector<double>& hn, std::vector<weightedBin>& hm) {
    if(hn.size() != hm.size()) {
        return -1.0;
    }
    double sumN = std::accumulate(hn.begin(), hn.end(), 0);
    double sumM = std::accumulate(hm.begin(), hm.end(), 0, [](double sum, weightedBin bin)->double{return sum + bin.wgt;});
    
    double chisqSum = 0.0;
    for(int i = 0; i < hn.size(); i++) {
        chisqSum += pow(hn[i]/sumN - hm[i].wgt/sumM,2)/(hn[i]/(sumN*sumN) + (hm[i].wgt*hm[i].wgt)/(sumM*sumM*hm[i].num));
    }
    return chisqSum;
}

double calcChisqRoot(TH1D* refHist, TH1D* mcHist) {
//    double chi2;
//    int ndf;
//    int igood;
//    refHist->Chi2TestX(mcHist, chi2, ndf, igood, "UW CHI2/NDF");
//    printf("igood: %d\n", igood);
//    if(igood) {
//        return -1;
//    }
//    return chi2/(double)ndf;
    return refHist->Chi2Test(mcHist, "UW CHI2/NDF");
}

int main(int argc, char** argv) {
    ROOT::EnableThreadSafety();
    initxorshift();
    //buff_len = 4 + 3*8 + 4
    const size_t buff_len = 4 + 1*8 + 2*NRECORDS*4 + 4;
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
        event.energy = *((double *)(&buf[0] + sizeof(unsigned int)));
        std::memcpy((void *)&event.times, (void *)(&buf[0] + sizeof(unsigned int) + sizeof(double)), NRECORDS*sizeof(float));
        std::memcpy((void *)&event.ePerp, (void *)(&buf[0] + sizeof(unsigned int) + sizeof(double) + NRECORDS*sizeof(float)), NRECORDS*sizeof(float));
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
    
    TH1D* refHistRoot = new TH1D("ref", "Reference", 184, 0, 184);
    for(int i = 0; i < 184; i++) {
        refHistRoot->SetBinContent(i+1, refHist[i]);
    }
    
    refHistRoot->Sumw2();
    
    int nBins = 20;
//    #pragma omp parallel for collapse(3)
  /*  for(int i = 0; i < nBins+1; i++) {
        for(int j = 0; j < nBins+1; j++) {
            for(int k = 0; k < nBins+1; k++) {
                double thresh = 3.0 + 6.0*i/(double)nBins;
                //double thresh = 0.0;
                double thickOxide = 0 + 30*j/(double)nBins;
                double thickBoron = 0 + 30*k/(double)nBins;
                if(thickOxide + thickBoron < 3) {
                    continue;
                }
//                std::vector<weightedBin> hist1 = createHistQuantMultilayerEdESpline(thickOxide, thickBoron, thresh, events, randU01s, randDeathTimes);
//                double chisqWgt = calcChisqWgt(refHist, hist1);
//                double chisqUnWgt = calcChisq(refHist, hist1);
//                double chisqNate = calcChisqNate(refHist, hist1);
//                printf("%f %f %f %f %f %f\n", thickOxide, thickBoron, thresh, chisqWgt/hist1.size(), chisqUnWgt/hist1.size(), chisqNate/hist1.size());
                TH1D* mcHistRoot = createHistQuantMultilayerEdESplineRoot(thickOxide, thickBoron, thresh, events, randU01s, randDeathTimes);
                double chisq = calcChisqRoot(refHistRoot, mcHistRoot);
                printf("%f %f %f %f\n", thickOxide, thickBoron, thresh, chisq);
                delete mcHistRoot;
                fflush(stdout);
            }
        }
    }*/
    
    /*double fitOffsets[8] = {4, 42, 62, 82, 102, 122, 142, 162};
    double fitEnds[8] = {38, 58, 78, 98, 118, 138, 158, 178};
    double refTCs[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double refTCsErr[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double mcTCs[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double mcTCsErr[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    for(int i = 0; i < 8; i++) {
        TF1* fit = new TF1("refFit", "[0]*exp(-(x-[2])/[1])", fitOffsets[i], fitEnds[i]);
        fit->SetParameter(2, fitOffsets[i]);
        fit->SetParLimits(2, fitOffsets[i], fitOffsets[i]);
        fit->SetParameter(0, 1000);
        fit->SetParameter(1, 20);
        refHistRoot->Fit(fit, "LQ", "", fitOffsets[i], fitEnds[i]);
        refTCs[i] = fit->GetParameter(1);
        refTCsErr[i] = fit->GetParError(1);
        delete fit;
    }
    int nBins = 20;
//    #pragma omp parallel for collapse(3)
    for(int i = 0; i < nBins+1; i++) {
        for(int j = 0; j < nBins+1; j++) {
            for(int k = 1; k < nBins+1; k++) {
                double thresh = 3.0 + 6.0*i/(double)nBins;
                double thickOxide = 0 + 30*j/(double)nBins;
                double thickBoron = 0 + 30*k/(double)nBins;
                if(thickOxide + thickBoron < 3) {
                    continue;
                }
                fitTCQuantMultilayerEdESplineRoot(thickOxide, thickBoron, thresh, events, randU01s, randDeathTimes, mcTCs, mcTCsErr);
                double chisq = 0.0;
                for(int l = 0; l < 8; l++) {
                    chisq += (refTCs[l]-mcTCs[l])*(refTCs[l]-mcTCs[l])/(refTCsErr[l]*refTCsErr[l] + mcTCsErr[l]*mcTCsErr[l]);
                    //printf("%f %f %f %f %f\n", refTCs[i], refTCsErr[i], mcTCs[i], mcTCsErr[i], refTCs[i]-mcTCs[i]);
                }
                printf("%f %f %f %f\n", thickOxide, thickBoron, thresh, chisq/7);
                fflush(stdout);
            }
        }
    }*/
    
//    double thresh = 5;
//    double thick = 4;
//    std::vector<weightedBin> hist1 = createHistQuantMultilayerEdE(thick, thresh, events, randU01s, randDeathTimes);
//    double chisqWgt = calcChisqWgt(refHist, hist1);
//    double chisqUnWgt = calcChisq(refHist, hist1);
//    double chisqNate = calcChisqNate(refHist, hist1);
////                printf("%f %f %f %f\n", GRAV*MASS_N*(0.01+0.05*i/10.0)*JTONEV, GRAV*MASS_N*(0.1+0.02*j/10.0)*JTONEV, 0.5+0.5*k/10.0, chisq/hist1.size());
//    printf("%f %f %f %f %f\n", thick, thresh, chisqWgt/hist1.size(), chisqUnWgt/hist1.size(), chisqNate/hist1.size());
//    fflush(stdout);
//    thresh = 6;
//    thick = 5;
//    hist1 = createHistQuantMultilayerEdE(thick, thresh, events, randU01s, randDeathTimes);
//    chisqWgt = calcChisqWgt(refHist, hist1);
//    chisqUnWgt = calcChisq(refHist, hist1);
//    chisqNate = calcChisqNate(refHist, hist1);
////                printf("%f %f %f %f\n", GRAV*MASS_N*(0.01+0.05*i/10.0)*JTONEV, GRAV*MASS_N*(0.1+0.02*j/10.0)*JTONEV, 0.5+0.5*k/10.0, chisq/hist1.size());
//    printf("%f %f %f %f %f\n", thick, thresh, chisqWgt/hist1.size(), chisqUnWgt/hist1.size(), chisqNate/hist1.size());
//    fflush(stdout);
    
    

//    std::vector<weightedBin> hist1 = createHist(0.155, 0.0, 15.377918, 17.428307, events, randU01s, randDeathTimes);
//    std::vector<weightedBin> hist1 = createHist(0.26, 2.1, 5.125973, 30.755837, events, randU01s, randDeathTimes);
//    std::vector<weightedBin> hist1 = createHist(0.225, 1.1, 7.176362, 22.9, events, randU01s, randDeathTimes); //Nate's chisq
//    std::vector<weightedBin> hist1 = createHist(0.27, 2.3, 5.125973, 34.5, events, randU01s, randDeathTimes); //Nate's chisq EdE
//      std::vector<weightedBin> hist1 = createHist(0.2551, 1.892, 5.701, 29.19, events, randU01s, randDeathTimes); //Prof II
//      std::vector<weightedBin> hist1 = createHist(0.285, 2.6, 4.100778, 34.5, events, randU01s, randDeathTimes); //Prof II
//      std::vector<weightedBin> hist1 = createHistQuant(0.1, 4.100778, 34.5, events, randU01s, randDeathTimes); //Test of quantum
//      std::vector<weightedBin> hist1 = createHistQuant(0.4, 4.100778, 34.5, events, randU01s, randDeathTimes); //Test of quantum
//      std::vector<weightedBin> hist1 = createHistQuant(0.444444, 3.8, 34.5, events, randU01s, randDeathTimes); //Test of quantum
//      std::vector<weightedBin> hist1 = createHistQuant(0.363636, 5.8, 34.5, events, randU01s, randDeathTimes); //Test of quantum frac
      //std::vector<weightedBin> hist1 = createHistQuant(0.35, 4.4, 34.5, events, randU01s, randDeathTimes); //quant frac zoom
//      std::vector<weightedBin> hist1 = createHistQuant(0.31, 6.0, 24.0, events, randU01s, randDeathTimes); //quant frac zoom
//      std::vector<weightedBin> hist1 = createHistQuantEdE(0.35, 0.0, events, randU01s, randDeathTimes); //quant frac zoom
//      std::vector<weightedBin> hist1 = createHistQuantEdE(0.32, 0.0, events, randU01s, randDeathTimes); //quant EdE opt, no cutoff
    //std::vector<weightedBin> hist1 = createHistQuantMultilayerEdE(5, 5, 5, events, randU01s, randDeathTimes);
//    std::vector<weightedBin> hist1 = createHistQuantMultilayerEdE(3, 4.5, 6, events, randU01s, randDeathTimes);
//    std::vector<weightedBin> hist1 = createHistQuantMultilayerEdE(6.0, 4.5, 8.7, events, randU01s, randDeathTimes);
//    std::vector<weightedBin> hist1 = createHistQuantMultilayerEdE(9.0, 6.0, 9.0, events, randU01s, randDeathTimes);
    std::vector<weightedBin> hist1 = createHistQuantMultilayerEdE(3.0, 6.0, 0.0, events, randU01s, randDeathTimes);
    //double chisq2 = calcChisq(refHist, hist1);
    //printf("%f\n\n", chisq2/hist1.size());
    for(auto it = hist1.begin(); it < hist1.end(); it++) {
      printf("%f,", it->wgt);
    }
    printf("\n");

    delete[] randU01s;
    delete[] buf;
    
    delete refHistRoot;

    return 0;
}

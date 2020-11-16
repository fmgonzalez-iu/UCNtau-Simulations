#include "../inc/chisq_spectrum.hpp"

/*----------------------------------------------------------------------
 * This file contains various methods of calculating the chi^2 of a
 * given unload.
 * 
 * All these (except the ROOT calculation) require an input reference
 * and an output calculation that's a "weightedBin" vector, formed by 
 * scaling each event by some scaling factor.
 * 
 * The Reference Histogram (refHist) is the rate of our unload on a 
 * second-to-second basis.
 *--------------------------------------------------------------------*/
 
double calcChisqGagunashvili(std::vector<double>& hn, std::vector<weightedBin>& hm) {
	// Gagunashvili algorithm (arXiv:0905.4221)
	// Chi^2 calculation assuming probabilistically generated histograms
	
    if(hn.size() != hm.size()) {
		printf("Size Mismatch!\n");
        return -1.0;
    }
    
    double sumN = std::accumulate(hn.begin(), hn.end(), 0);
    //double sumM_1 = std::accumulate(hm.begin(), hm.end(), 0, [](double sum, weightedBin bin)->double{return sum + bin.wgt;});
    
    double sumM = 0.0;
    for (auto h = hm.begin(); h < hm.end(); h++) {
		sumM += (*h).wgt;
	}
	printf("%f,%f\n",sumN,sumM);
    double chisqSum = 0.0;
    for(int i = 0; i < hn.size(); i++) {
		if (hm[i].wgt > 0.0 && hn[i] > 0.0) { // Require weighting to be non-zero
	        double phat = (sumM*hm[i].wgt - sumN*hm[i].wgtSqr
	                + sqrt(pow(sumM*hm[i].wgt - sumN*hm[i].wgtSqr, 2) + 4*sumM*sumM*hm[i].wgtSqr*hn[i]))
	                                /
	                          (2*sumM*sumM);
			chisqSum += pow(hn[i] - sumN*phat, 2)/(sumN*phat);
			chisqSum += pow(hm[i].wgt - sumM*phat, 2)/(hm[i].wgtSqr);
		}
    }
    
    return chisqSum;
}


double calcChisq(std::vector<double>& hn, std::vector<weightedBin>& hm) {
	// Normal chi squared algorithm
	
    if(hn.size() != hm.size()) {
        return -1.0;
    }
    double cn = std::accumulate(hm.begin(), hm.end(), 0, [](double sum, weightedBin bin)->double{return sum + bin.wgt;});
    double cm = std::accumulate(hn.begin(), hn.end(), 0);
//    printf("%f %f\n", cn, cm);
    double chisqSum = 0.0;
    for(int i = 0; i < hn.size(); i++) {
		if (hm[i].wgt > 0.0 && hn[i] > 0.0) { // Require weighting to be non-zero
			chisqSum += (1/(cn*cm))*pow(cn*hn[i] - cm*hm[i].wgt,2)/(hn[i]+hm[i].wgt);
		}
    }
    return chisqSum;
}

double calcChisqNate(std::vector<double>& hn, std::vector<weightedBin>& hm) {
	// An alternate chi squared algorithm that Nate wrote at some point
	
	
    if(hn.size() != hm.size()) {
        return -1.0;
    }
    double sumN = std::accumulate(hn.begin(), hn.end(), 0);
    double sumM = std::accumulate(hm.begin(), hm.end(), 0, [](double sum, weightedBin bin)->double{return sum + bin.wgt;});
    
    double chisqSum = 0.0;
    for(int i = 0; i < hn.size(); i++) {
		if (hm[i].wgt > 0.0 && hn[i] > 0.0) { // Require weighting to be non-zero
			chisqSum += pow(hn[i]/sumN - hm[i].wgt/sumM,2)/(hn[i]/(sumN*sumN) + (hm[i].wgt*hm[i].wgt)/(sumM*sumM*hm[i].num));
		}
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

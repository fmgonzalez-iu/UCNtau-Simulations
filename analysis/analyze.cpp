#include "TApplication.h"

//#include "./inc/arrivalTimes.hpp"
//#include "./inc/blockReader.hpp"

#include "./inc/2det_plotting.hpp"
#include "./inc/chisq_spectrum.hpp"
#include "./inc/file_reader.hpp"
#include "./inc/generate_weighting.hpp"
#include "./inc/quant_refl.hpp"
#include "./inc/stitcher.hpp"


extern "C" {
	#include "../general/setup.h"
	#include "../general/typedefs.h"
    //#include "./inc/xorshift.h"
}

#define START (20 + 41) // ??
#define END (20 + 41 + 184) // ??

/*----------------------------------------------------------------------
 * This is a code that does a bunch of weird analysis stuff based on the 
 * trajectory simulations for UCNtau.
 * 
 * Look in the folder "../general/" for setup.h and typedefs.h, both of 
 * which should be shared between the simulations and analysis software.
 * 
 * There are reader and stitcher functions, as well as generalized surface
 * reflection functions in the src and inc folders. Some specific block
 * reflection things and the general chi2 algorithms are there as well.
 * 
 * 
 * (Sorry for anyone who has to read this mess)
 * F
 * 
 *--------------------------------------------------------------------*/

int main(int argc, char** argv) {

	ROOT::EnableThreadSafety(); // ROOT multithreading
	
	//------------------------------------------------------------------
	// Command line arguments
	if(argc < 2) { // Check that we've input a file to analyze
		printf("Error! Usage: ./analyze fname hist fname2\n");
		return 2;
	}
	//------------------------------------------------------------------
	// Load events
	// Change this depending on what we're trying to run
	std::vector<std::vector<double>> positions = readFilePosition(argv[1]);
    //std::vector<fixedResult> events_tmp = readFileFixedRes(argv[1]); 
//    std::vector<noabsCleanDagResultZ> events_tmp = readFileNoAbsCleanDagResZ(argv[1]); 
    
	std::vector<double> refHist; // Initial reference histogram
	if (argc < 3) {
		refHist= {276.0,2147.0,6357.0,6724.0,6445.0,6515.0,6303.0,5883.0,5937.0,5750.0,5862.0,5741.0,5499.0,5298.0,5118.0,5336.0,5208.0,4989.0,4914.0,4840.0,4773.0,4761.0,4603.0,4670.0,4585.0,4490.0,4504.0,4344.0,4234.0,4190.0,4215.0,4159.0,4099.0,4054.0,4042.0,3958.0,3875.0,3758.0,3770.0,6777.0,19591.0,26111.0,24937.0,23352.0,21985.0,21214.0,20137.0,19300.0,18771.0,17531.0,16790.0,16027.0,15480.0,15079.0,14571.0,13997.0,13510.0,13170.0,12590.0,18787.0,32668.0,30531.0,29039.0,27367.0,25942.0,24249.0,22440.0,21718.0,20811.0,19545.0,18862.0,17545.0,16773.0,15964.0,15590.0,14640.0,14114.0,13309.0,12663.0,19992.0,29095.0,27203.0,25432.0,24246.0,22596.0,21316.0,20074.0,18979.0,17967.0,16720.0,16213.0,15231.0,14874.0,13783.0,13126.0,12620.0,12152.0,11748.0,11099.0,19469.0,31086.0,29762.0,28316.0,26476.0,24893.0,23408.0,22080.0,21091.0,19686.0,18667.0,17434.0,16717.0,15647.0,15110.0,14155.0,13375.0,12639.0,12072.0,11831.0,18263.0,23110.0,22100.0,21373.0,20302.0,19436.0,18349.0,17475.0,16732.0,15687.0,14945.0,14201.0,13216.0,12934.0,12358.0,11613.0,11052.0,10341.0,9803.0,9612.0,15876.0,20286.0,19403.0,17943.0,16981.0,15999.0,14919.0,13821.0,12959.0,11977.0,11107.0,10590.0,9923.0,9119.0,8469.0,8032.0,7199.0,6634.0,6475.0,5908.0,8792.0,10150.0,8870.0,7618.0,6546.0,5978.0,5135.0,4541.0,3926.0,3486.0,3176.0,2779.0,2399.0,2159.0,1892.0,1670.0,1493.0,1387.0,1236.0,1038.0,959.0,787.0,739.0,663.0,589.0};
	} else {
		refHist =refHistFromRoot(argv[2]);
	}
	
	// Load our events (second detector)	
	std::vector<noabsCleanResult> eventscleaner;
	if (argc >= 3) {
		eventscleaner = readFileNoAbsCleanRes(argv[3]);
	}
		
	std::vector<double> refHist2; // initial reference histogram for 2nd detector
	if (argc < 5) {
		refHist2 = {276.0,2147.0,6357.0,6724.0,6445.0,6515.0,6303.0,5883.0,5937.0,5750.0,5862.0,5741.0,5499.0,5298.0,5118.0,5336.0,5208.0,4989.0,4914.0,4840.0,4773.0,4761.0,4603.0,4670.0,4585.0,4490.0,4504.0,4344.0,4234.0,4190.0,4215.0,4159.0,4099.0,4054.0,4042.0,3958.0,3875.0,3758.0,3770.0,6777.0,19591.0,26111.0,24937.0,23352.0,21985.0,21214.0,20137.0,19300.0,18771.0,17531.0,16790.0,16027.0,15480.0,15079.0,14571.0,13997.0,13510.0,13170.0,12590.0,18787.0,32668.0,30531.0,29039.0,27367.0,25942.0,24249.0,22440.0,21718.0,20811.0,19545.0,18862.0,17545.0,16773.0,15964.0,15590.0,14640.0,14114.0,13309.0,12663.0,19992.0,29095.0,27203.0,25432.0,24246.0,22596.0,21316.0,20074.0,18979.0,17967.0,16720.0,16213.0,15231.0,14874.0,13783.0,13126.0,12620.0,12152.0,11748.0,11099.0,19469.0,31086.0,29762.0,28316.0,26476.0,24893.0,23408.0,22080.0,21091.0,19686.0,18667.0,17434.0,16717.0,15647.0,15110.0,14155.0,13375.0,12639.0,12072.0,11831.0,18263.0,23110.0,22100.0,21373.0,20302.0,19436.0,18349.0,17475.0,16732.0,15687.0,14945.0,14201.0,13216.0,12934.0,12358.0,11613.0,11052.0,10341.0,9803.0,9612.0,15876.0,20286.0,19403.0,17943.0,16981.0,15999.0,14919.0,13821.0,12959.0,11977.0,11107.0,10590.0,9923.0,9119.0,8469.0,8032.0,7199.0,6634.0,6475.0,5908.0,8792.0,10150.0,8870.0,7618.0,6546.0,5978.0,5135.0,4541.0,3926.0,3486.0,3176.0,2779.0,2399.0,2159.0,1892.0,1670.0,1493.0,1387.0,1236.0,1038.0,959.0,787.0,739.0,663.0,589.0};
	}else{
		refHist2 = refHistFromRoot(argv[4]);
	}
	//------------------------------------------------------------------
    
    
	// Change this for different file types (TODO: setup.h config?)
	
	
	//------------------------------------------------------------------
	
	//------------------------------------------------------------------
	// If we want to convert between file types, put a stitching function here
	//std::vector<noabsResult> events = stitch_noabs_to_noabsClean(events_tmp); 
	//TH1D* refHistRoot = stitch_weighted_vec_to_root(refHist);	
	//------------------------------------------------------------------
		
	//------------------------------------------------------------------
	// Initialize random numbers -- commented out since we're not doing temp events right now
	/*initxorshift(); // RNG seed
	double* randU01s = new double[events_tmp.size()*2*NRECORDS+eventscleaner.size()];
	for(unsigned long i = 0; i < events_tmp.size()*2*NRECORDS+eventscleaner.size(); i++) {
		randU01s[i] = nextU01();
	}
	double* randDeathTimes = new double[events_tmp.size()];
	for(unsigned long i = 0; i < events_tmp.size(); i++) {
		randDeathTimes[i] = -TAU_N*log(nextU01());
	}*/
    //------------------------------------------------------------------
    
    //------------------------------------------------------------------
	
	// ROOT reference histogram
	TH1D* refHistRoot = new TH1D("refHist", "Reference", (int)refHist.size(), 0, (int)refHist.size());
	for(int i = 0; i < refHist.size(); i++) {
		refHistRoot->SetBinContent(i+1, refHist[i]);
	}
	
	// Actual histogram we're making
	//TH1D* histRoot = stitch_fixed_result_to_root(events_tmp);
	TH2D* histRoot = create_root_2D(positions);
	histRoot->SaveAs("rootHist.root");
	
	// Root plotting for debug -- uncomment this to do debug plots
	//TApplication theApp("App",0,0); // This is so that we can plot ROOT stuff
			
	//TCanvas* c1 = new TCanvas("c1","c1",750,750);
	//c1->cd();
	//refHistRoot->Sumw2();
	//refHistRoot->Scale(1/refHistRoot->Integral(0,200));
	//refHistRoot->SetLineColor(kBlue);
	//refHistRoot->Draw();
	//refHistRoot->SaveAs("refHist.root");
	
	//TH1D* hist1Root = createHistQuantSplineRoot(0.0, 5.0, 7.0, 1.0, 0.25, events,refHist, randU01s, randDeathTimes);
	//std::vector<reducedEvent> test2d = reduceDag100(0.0,7.5,0.0,1.406250,0.25,0.00,events_tmp,randU01s);
	//std::vector<reducedEvent> test2d = reduceDag100(0.0,8.0,0.0,1.4,0.45,0.00,events_tmp,randU01s);
	//TH1D* hist1Root =  createHistReduced(test2d,refHist);
	//TCanvas* c3 = plotResiduals(hist1Root, refHistRoot);
	
	//std::vector<weightedBin> test2d = createHistQuantMult_2det(10.0,5.6, 0.0, 0.0, 0.0, 0.2, 0.0,0.013,events_tmp,eventscleaner,refHist.size(),randU01s,randDeathTimes);//,randU01s, randDeathTimes);
	//TH1D* hist1Root = stitch_weighted_vec_to_root(test2d);
	//hist1Root->Scale(1/hist1Root->Integral(0,200));
	//hist1Root->SetLineColor(kRed);
	//hist1Root->Draw("SAME");
	
	//std::vector<reducedEvent> test2d2 = reduceDag100(0.0,10.6,7.0,0.2,0.28,0.013,events_tmp,randU01s);
	//TH1D* hist2Root =  createHistReduced(test2d2,refHist);
	////std::vector<weightedBin> test2d2 = createHistQuantMult_2det(15.0,5.6, 0.0, 0.0, 0.0, 0.2, -0.28,0.013,events_tmp,eventscleaner,refHist.size(),randU01s,randDeathTimes);//,randU01s, randDeathTimes);
	////TH1D* hist2Root = stitch_weighted_vec_to_root(test2d2);
	//hist2Root->Scale(1/hist2Root->Integral(0,200));
	//hist2Root->SetLineColor(kGreen);
	//hist2Root->Draw("SAME");
	//hist2Root->SaveAs("guessHist.root");
	
	//std::vector<reducedEvent> test2d3 = reduceDag100(0.0,15.6,7.0,0.2,0.28,0.013,events_tmp,randU01s);
	//TH1D* hist3Root =  createHistReduced(test2d3,refHist);
	////std::vector<weightedBin> test2d3 = createHistQuantMult_2det(0.0,10.0, 0.0, 0.0, 7.0, 0.2, 0.28,0.013,events_tmp,eventscleaner,refHist.size(),randU01s,randDeathTimes);//,randU01s, randDeathTimes);
	////TH1D* hist3Root = stitch_weighted_vec_to_root(test2d3);
	//hist3Root->Scale(1/hist3Root->Integral(0,200));
	//hist3Root->SetLineColor(kBlack);
	//hist3Root->Draw("SAME");
	//c1->Update();
	
	//eventscleaner.erase(eventscleaner.begin(),eventscleaner.end()); // save memory by deleting cleaner (for now)
	
	//std::vector<reducedEvent> results = reduceDag100(0.0,8.0,0.0,1.4,0.45,0.00,events_tmp,randU01s);
	//TCanvas* c2 = plotSanityCheck(results);
	
	
	//theApp.Run();
	return 2;
	
	//hist1Root->SaveAs("Test.root");
	//double chisqR = calcChisqRoot(refHistRoot, hist1Root);
	//printf("%f\n",chisqR);
	
	// Calculate the chi^2 for a bunch of numbers across a parameter sweep
	/*int nBins = 0; // Might be nice to put this into a separate function?
	#pragma omp parallel for collapse(5) // parallelize a square loop to divide into 5 sections
	for(int i = 0; i < nBins+1; i++) {
		for(int j = 0; j < nBins+1; j++) {
			for(int k = 0; k < nBins+1; k++) {
				for(int l = 0; l < nBins+1; l++) {
				for(int m = 0; m < nBins+1; m++) {
					//for(int m = 0; m < nBins+1; m++) {
					// Define fitting parameters, relating to parts of the loop
					double thickOxide = 0.0;
					double thresh     = 0.0  + 14.0* i / (double)nBins; 
					double thickBoron = 5.0  + 10.0 * k / (double)nBins;
					double thickClean = 5.0  + 10.0 * k / (double)nBins; // For now doesn't do anything
					double power      = 0.75 + 1.75  * l / (double)nBins; // This is now zero-indexed
					double cosPower   = 0.0  + 0.5 * j / (double)nBins;
					double zeta       = 0.0  + 0.10 * m / (double)nBins;
//                    std::vector<weightedBin> hist1 = createHistQuantNoOxEPowdEThetaSpline_C(thickBoron, thresh, power, cosPower, events, randU01s, randDeathTimes);
					// Create a weighted histogram
					std::vector<weightedBin> hist1 = createHistQuantMult_2det(thickOxide, thickBoron, thickOxide, thickClean, thresh, power, cosPower,zeta,events_tmp,eventscleaner,refHist.size(),randU01s, randDeathTimes);
					//TH1D* hist1Root = stitch_weighted_vec_to_root(test2d);
					//hist1Root->SetLineColor(kRed);
					//hist1Root->Draw("SAME");
					//std::vector<weightedBin> hist1 = createHistQuantSpline(thickOxide, thickBoron, thresh, power, cosPower, events,refHist, randU01s, randDeathTimes);
					//double chisq = calcChisqGagunashvili(refHist, hist1)/(hist1.size()-1); // Calculate Chi^2 / NDF
					double chisq = calcChisq(refHist, hist1)/(hist1.size()-1); // Calculate Chi^2 / NDF
					//theApp.Run();
					printf("%f %f %f %f %f %f\n", thickBoron, thresh, power, cosPower,zeta, chisq); // Print
					fflush(stdout);
				
				}
				}
			}
		}
	}
	//8.75, 1.1875, 0.125,
	//7.5, 1.406250, 0.25
	//std::vector<weightedBin> hist1 = createHistQuantMultilayer(0.0,3.0, 6.0, 1.0,0.0, events, refHist,randU01s, randDeathTimes);
	// these values claim chi2 of 3.122
	std::vector<weightedBin> hist1 = createHistQuantMult_2det(0.0, 8.0, 0.0, 10.0, 0.0, 1.4, 0.45,0.00,events_tmp,eventscleaner,refHist.size(),randU01s, randDeathTimes);
	//for(auto it = hist1.begin(); it < hist1.end(); it++) {
	//	printf("%f,", it->wgt);
	//}
	printf("done\n");
	//theApp.Run();
	delete[] randU01s; // Free Memory by deleting vectors
	delete[] randDeathTimes;
	//delete[] buf;
	delete refHistRoot; // and also the root hist

	return 0; */
}
	
//-----------------------------------------------------------------------------------------------------------------------------
// Here there be dragons
//-----------------------------------------------------------------------------------------------------------------------------
    
//-----------------------------------------------------------------

    
   
////    #pragma omp parallel for collapse(3)
//for(int i = 0; i < nBins+1; i++) {
	//for(int j = 0; j < nBins+1; j++) {
		//for(int k = 0; k < nBins+1; k++) {
			//double thresh = 3.0 + 6.0*i/(double)nBins;
			////double thresh = 0.0;
			//double thickOxide = 0 + 30*j/(double)nBins;
			//double thickBoron = 0 + 30*k/(double)nBins;
			//if(thickOxide + thickBoron < 3) {
				//continue;
			//}
////                std::vector<weightedBin> hist1 = createHistQuantMultilayerEdESpline(thickOxide, thickBoron, thresh, events, randU01s, randDeathTimes);
////                double chisqWgt = calcChisqWgt(refHist, hist1);
////                double chisqUnWgt = calcChisq(refHist, hist1);
////                double chisqNate = calcChisqNate(refHist, hist1);
////                printf("%f %f %f %f %f %f\n", thickOxide, thickBoron, thresh, chisqWgt/hist1.size(), chisqUnWgt/hist1.size(), chisqNate/hist1.size());
			//TH1D* mcHistRoot = createHistQuantMultilayerEdESplineRoot(thickOxide, thickBoron, thresh, events, randU01s, randDeathTimes);
			//double chisq = calcChisqRoot(refHistRoot, mcHistRoot);
			//printf("%f %f %f %f\n", thickOxide, thickBoron, thresh, chisq);
			//delete mcHistRoot;
			//fflush(stdout);
		//}
	//}
//}

//double fitOffsets[8] = {4, 42, 62, 82, 102, 122, 142, 162};
//double fitEnds[8] = {38, 58, 78, 98, 118, 138, 158, 178};
//double refTCs[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
//double refTCsErr[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
//double mcTCs[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
//double mcTCsErr[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
//for(int i = 0; i < 8; i++) {
	//TF1* fit = new TF1("refFit", "[0]*exp(-(x-[2])/[1])", fitOffsets[i], fitEnds[i]);
	//fit->SetParameter(2, fitOffsets[i]);
	//fit->SetParLimits(2, fitOffsets[i], fitOffsets[i]);
	//fit->SetParameter(0, 1000);
	//fit->SetParameter(1, 20);
	//refHistRoot->Fit(fit, "LQ", "", fitOffsets[i], fitEnds[i]);
	//refTCs[i] = fit->GetParameter(1);
	//refTCsErr[i] = fit->GetParError(1);
	//delete fit;
//}
//int nBins = 20;
////    #pragma omp parallel for collapse(3)
//for(int i = 0; i < nBins+1; i++) {
	//for(int j = 0; j < nBins+1; j++) {
		//for(int k = 1; k < nBins+1; k++) {
			//double thresh = 3.0 + 6.0*i/(double)nBins;
			//double thickOxide = 0 + 30*j/(double)nBins;
			//double thickBoron = 0 + 30*k/(double)nBins;
			//if(thickOxide + thickBoron < 3) {
				//continue;
			//}
			//fitTCQuantMultilayerEdESplineRoot(thickOxide, thickBoron, thresh, events, randU01s, randDeathTimes, mcTCs, mcTCsErr);
			//double chisq = 0.0;
			//for(int l = 0; l < 8; l++) {
				//chisq += (refTCs[l]-mcTCs[l])*(refTCs[l]-mcTCs[l])/(refTCsErr[l]*refTCsErr[l] + mcTCsErr[l]*mcTCsErr[l]);
				////printf("%f %f %f %f %f\n", refTCs[i], refTCsErr[i], mcTCs[i], mcTCsErr[i], refTCs[i]-mcTCs[i]);
			//}
			//printf("%f %f %f %f\n", thickOxide, thickBoron, thresh, chisq/7);
			//fflush(stdout);
		//}
	//}
//}
//int nBins = 20;
    

    
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
    
        /*int nBins = 24;
    #pragma omp parallel for collapse(4)
    for(int i = 0; i < nBins+1; i++) {
        for(int j = 0; j < nBins+1; j++) {
            for(int k = 0; k < nBins+1; k++) {
                for(int l = 0; l < nBins+1; l++) {
                    double thresh = 5.5 + 1.5*i/(double)nBins;
                    double thickBoron = 4.4 + 0.4*k/(double)nBins;
                    double power = 1.0 + 0.4*l/(double)nBins;
                    double cosPower = 0.1 + 0.2*j/(double)nBins;
//                    std::vector<weightedBin> hist1 = createHistQuantNoOxEPowdEThetaSpline(thickBoron, thresh, power, cosPower, events, randU01s, randDeathTimes);
                    std::vector<weightedBin> hist1 = createHistQuantNoOxEPowdEThetaSpline_C(thickBoron, thresh, power, cosPower, events, randU01s, randDeathTimes);
                    double chisq = calcChisqGagunashvili(refHist, hist1)/(hist1.size()-1);
                    printf("%.15f %.15f %.15f %.15f %.15f\n", cosPower, thickBoron, thresh, power, chisq);
                    fflush(stdout);
                }
            }
        }
    }*/
    
//    int nBins = 100;
////    std::vector<double> param;
//    std::vector<double> minParam = {4.66666667, 5.875, 1.216666667, 0.25};
//    std::vector<double> rangeParam = {0.5, 5.875, 0.5, 0.25};
//    std::vector<std::vector<int>> pairs = {{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};
//    for(auto it = pairs.begin(); it < pairs.end(); it++) {
//        #pragma omp parallel for collapse(2)
//        for(int x = 0; x <= nBins; x++) {
//            for(int y = 0; y <= nBins; y++) {
//                std::vector<double> param;
//                param = minParam;
//                param[(*it)[0]] = param[(*it)[0]] - rangeParam[(*it)[0]] + 2*rangeParam[(*it)[0]]*(x/(double)nBins);
//                param[(*it)[1]] = param[(*it)[1]] - rangeParam[(*it)[1]] + 2*rangeParam[(*it)[1]]*(y/(double)nBins);
//                std::vector<weightedBin> hist1 = createHistQuantNoOxEPowdEThetaSpline_C(param[0], param[1], param[2], param[3], events, randU01s, randDeathTimes);
//                double chisq = calcChisqGagunashvili(refHist, hist1)/(hist1.size()-1);
//                printf("%.15f %.15f %.15f %.15f %.15f\n", param[3], param[0], param[1], param[2], chisq);
//                fflush(stdout);
//            }
//        }
//    }

//    std::vector<weightedBin> hist1 = createHistQuantMultilayerEPowdESpline(0.0, 4.6, 11.6, 1.1, events, randU01s, randDeathTimes);
//    std::vector<weightedBin> hist1 = createHistQuantNoOxEPowdEThetaSpline(4.6, 10.0, 1.1, 0.1, events, randU01s, randDeathTimes);
//    std::vector<weightedBin> hist1 = createHistQuantNoOxEPowdEThetaSpline(4.60606, 7.212121, 1.227273, 0.181818, events, randU01s, randDeathTimes);
//    std::vector<weightedBin> hist1 = createHistQuantNoOxEPowdEThetaSpline_C(4.60606, 7.212121, 1.227273-1, 0.181818, events, randU01s, randDeathTimes);
//    std::vector<weightedBin> hist1 = createHistQuantNoOxEPowdEThetaSpline_C(4.6, 5.875, 1.216666667, 0.25, events, randU01s, randDeathTimes);
////    double chisq2 = calcChisqGagunashvili(refHist, hist1);
////    printf("%f\n\n", chisq2/(hist1.size()-1));
////    for(auto it = hist1.begin(); it < hist1.end(); it++) {
//    for(int i = 0; i < hist1.size(); i++) {
//      printf("%f,", hist1[i].wgt);
//    }
//    printf("\n");
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



//#define START (300+20+41)
//#define END (300+20+41+184)

//#define START 0
//#define END 1000
//#define START 41
//#define END 225

//#define START (150 + 200 + 20 + 20)
//#define END (150 + 200 + 20 + 20 + 100)


//------------------------------------------------------------------
// This block of code does the readFile function
//------------------------------------------------------------------
//buff_len = 4 + 3*8 + 4
//const size_t buff_len = 4 + 1*8 + 2*NRECORDS*4 + 4;
//const size_t buff_len = 4 + 2*8 + 2*NRECORDS*4 + 4; // buff_len needs to be dynamic (?)
//char* buf = new char[buff_len];

//std::ifstream binfile(argv[1], std::ios::in | std::ios::binary);
//if(!binfile.is_open()) {
	//printf("Error! Could not open file %s\n", argv[1]);
	//return 1;
//}
//std::vector<evt> events;
//evt event;
//while(!binfile.eof()) {
	//binfile.read(buf, buff_len);
	//if(binfile.eof()) { //Breaks on last read of file (i.e. when 0 bytes are read and EOF bit is set)
		//break;
	//}
	//if(*((unsigned int *)&buf[0]) != 2*8 + 2*NRECORDS*4) {
		//fprintf(stderr, "Error! Aliased read on\n");
		//exit(2);
	//}
	//event.energy = *((double *)(&buf[0] + sizeof(unsigned int)));
	////event.theta = *((double *)(&buf[0] + sizeof(unsigned int) + sizeof(double)));
	//std::memcpy((void *)&event.times, (void *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(double)), NRECORDS*sizeof(float));
	//std::memcpy((void *)&event.ePerp, (void *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(double) + NRECORDS*sizeof(float)), NRECORDS*sizeof(float));
	////assert(*(int *)(&buf[0]) == buff_len-8);
	//events.push_back(event);
//}
//binfile.close();

//printf("Read %lu Events!\n", events.size());
//------------------------------------------------------------------


//double chisq2 = calcChisq(refHist, hist1);
//printf("%f\n\n", chisq2/hist1.size());

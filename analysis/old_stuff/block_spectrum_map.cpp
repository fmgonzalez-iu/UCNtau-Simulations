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

#define NRECORDS 50 
#define BLOCKSTATE 4

#define MASS_N 1.674927471e-27
#define GRAV 9.80665e0
#define JTONEV 6.2415091e27
#define HBAR 1.054571800e-34

/*  For the block, not sure if I need dagger spectrum fittings, keeping in for now. */
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

std::vector<double> refHist = {276.0,2147.0,6357.0,6724.0,6445.0,6515.0,6303.0,5883.0,5937.0,5750.0,5862.0,5741.0,5499.0,5298.0,5118.0,5336.0,5208.0,4989.0,4914.0,4840.0,4773.0,4761.0,4603.0,4670.0,4585.0,4490.0,4504.0,4344.0,4234.0,4190.0,4215.0,4159.0,4099.0,4054.0,4042.0,3958.0,3875.0,3758.0,3770.0,6777.0,19591.0,26111.0,24937.0,23352.0,21985.0,21214.0,20137.0,19300.0,18771.0,17531.0,16790.0,16027.0,15480.0,15079.0,14571.0,13997.0,13510.0,13170.0,12590.0,18787.0,32668.0,30531.0,29039.0,27367.0,25942.0,24249.0,22440.0,21718.0,20811.0,19545.0,18862.0,17545.0,16773.0,15964.0,15590.0,14640.0,14114.0,13309.0,12663.0,19992.0,29095.0,27203.0,25432.0,24246.0,22596.0,21316.0,20074.0,18979.0,17967.0,16720.0,16213.0,15231.0,14874.0,13783.0,13126.0,12620.0,12152.0,11748.0,11099.0,19469.0,31086.0,29762.0,28316.0,26476.0,24893.0,23408.0,22080.0,21091.0,19686.0,18667.0,17434.0,16717.0,15647.0,15110.0,14155.0,13375.0,12639.0,12072.0,11831.0,18263.0,23110.0,22100.0,21373.0,20302.0,19436.0,18349.0,17475.0,16732.0,15687.0,14945.0,14201.0,13216.0,12934.0,12358.0,11613.0,11052.0,10341.0,9803.0,9612.0,15876.0,20286.0,19403.0,17943.0,16981.0,15999.0,14919.0,13821.0,12959.0,11977.0,11107.0,10590.0,9923.0,9119.0,8469.0,8032.0,7199.0,6634.0,6475.0,5908.0,8792.0,10150.0,8870.0,7618.0,6546.0,5978.0,5135.0,4541.0,3926.0,3486.0,3176.0,2779.0,2399.0,2159.0,1892.0,1670.0,1493.0,1387.0,1236.0,1038.0,959.0,787.0,739.0,663.0,589.0};

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

/* Random number generator --
 * we use this for the block since we can't really initialize all the values */
double calcDecayTime(){
	double tauN, deathTime;
	tauN = 877.7;
	
	deathTime = (double) rand() / (RAND_MAX);
	deathTime = -tauN*log(deathTime);
	
	return deathTime;
}
/* Removed quantum calculations for abs prob of given B thickness. 
 * If we want to reweight the Boron thickness we'll need to put it back in. (and do the 50 hit thing on the dagger) */ 

double fieldstrength(double position[3]) {//-mu*mod(B)
        
	double B_h = 0.005; //holding field strength
	double B = 1.4; //remnant magnet field strength
	double d = 0.0254; //thickness of layer of PM array
	double L = 0.02; //characteristic spacing of magnets
	double mu = -9.662364e-27/1.674927351e-27; //mu in units where m=1
	double N = 3.0; //how far out to go in field ripple expansion
	double x = position[0];
	double y = position[1];
	double z = position[2];

	double R,r;
	double A = sqrt(8.0)*B/M_PI; //parameter related to B -- shows up in expansion

	if (x > 0.0)
	{
			R = 1.0;
			r = 0.5;
	}
	else
	{
			R = 0.5;
			r = 1.0;
	}

	double rho = sqrt(y*y+z*z);
	double r_zeta = sqrt((rho-R)*(rho-R)+x*x);
	double B_halbach = 0.0;

	if (z < -1.0 && r_zeta < r)
	{
			double zeta = r-r_zeta;
			double eta = r*atan(x/(rho-R));
			double Bsum = 0.0;
			double k_m,k_n,m,n;

			for (m = 1.0;m<=N;m+=1.0)
			{
					k_m = 2*M_PI*(4.0*m-3.0)/L;
					for (n = 1.0;n<=N;n+=1.0)
					{
							k_n = 2*M_PI*(4.0*n-3.0)/L;

							Bsum += pow(-1.0,m)*pow(-1.0,n)/(4.0*m-3.0)/(4.0*n-3.0)*(1-exp(-k_m*d))*(1-exp(-k_n*d))*exp(-(k_n+k_m)*zeta)*cos((k_n-k_m)*eta);
					}
			}

			B_halbach = sqrt(B_h*B_h*(r+R)*(r+R)/(y*y+z*z)+A*A*Bsum);

	}

	return B_halbach;
}

/* ---------------------- Here's the actual MAIN program!! -----------------------------------*/
int main(int argc, char** argv) {

	/* Initialize ROOT and our random number generator */
    ROOT::EnableThreadSafety();
    initxorshift();
    TApplication theApp("App",0,0);
    
    /*----------------------------------------------------------------*/
    /* !!!DECLARATION OF VARIABLES!!! */
    /*----------------------------------------------------------------*/
    /* Buffer declarations */
    const size_t buff_len = 4 + 8*8 + 4*4 + 4; // Change if MC data changes
	const size_t blok_len = (4 + 6*8 + 4*4 + 4);
    char* buf;
    char* buf_blok;
    buf = new char[buff_len];
	buf_blok = new char[blok_len];
	
    /* IO variables */    
    int test; // test debug parameter
    int length_dag, length_blok;
    double hold_t;
    
    /* Load events structures*/
    std::vector<dag_evt> events_dag;
    dag_evt event_dag, dag_buffer;
	std::vector<block_evt> events_blok;
    block_evt event_blok, block_buffer;
    
    /* Maximization and Minimization variables */
    double maxCoords[3] = {-DBL_MAX,-DBL_MAX,-DBL_MAX};
    double minCoords[3] = { DBL_MAX, DBL_MAX, DBL_MAX};
    
    /* ROOT plotting/filling variables*/ 
    int numHistos, nBins;
    double threshold, spacing, maxEnergy, minEnergy, maxTime, maxField;
    char* tBName = new char[10];
	char* tBTitle = new char[48];
	char* fBName = new char[10];
	char* fBTitle = new char[48];
	char* fPadName = new char[4];
	char* fPadTitle = new char[48];
	
    /* Analysis variables */
    bool upscatter = false;
    int brokenCounter = 0, doubleCounts = 0, blockCounter = 0;
    double angleTheta, energyCm, deathTime, yieldError, neutronYield= 0.0;
    double histoWeight, eneMin, eThMin, theMin, theCut, absMin, dagCut;
    double field_block, zeta, boronProb;
    double groupYield[3], groupError[3];
    double binBreak[2], binBuff, gBins[6];
    
    /* Chi2 variables */
    int numSteps, chi2Dim, minDex;
    char *rootName = new char[128];
    char *c2Filename = new char[128];
	FILE *chi2OutFile;
    double minChi2, chi2Min = DBL_MAX;
	
	/* Fitting variables */
	int fastIndex, slowIndex;
	double blockSlope[2], blockErr[2], blockPop[2], blockPopErr[2];
	double padMin, padMax;
	double sliceSlope[numHistos][2], sliceErr[numHistos][2], slicePop[numHistos][2], slicePopErr[numHistos][2];
	double drainFit[2], drainErr[2];
	
    /*----------------------------------------------------------------*/
    /* !!!BOOLEAN LOGIC to quickly change program analysis!!! */
    /*----------------------------------------------------------------*/
    bool debugMode = false; // Turns on debug mode
    bool quickLoad = false; // Only loads the first 10 dagger and block events
    bool thetaOn   = false; // Modify thetaOn to change plotting E and theta
    bool chiFit    = false; // Run Chi2 algorithm -- Automatically turns on?
    bool chiMap    = false; // Do we want to generate N-dimensional Chi2 Maps?
    bool decayOn   =  true; // Neutron decay (with lifetime 877.7)
    bool pseOn     = false; // Changes holding times with pse deep cleaning 
    bool dagEdge   =  true; // Turns on a tapered edge of the dagger
    
    /* Plotting and ROOT Analysis*/ 
	bool blockMapTot    = false;
	bool blockMapSlice  = false;
	bool drainingTimeOn = false;
	bool dataMCComp     =  true;
	bool saveMCFit      = false;
	bool bFieldSlice    = false;
	bool sanityCheck    =  true; // Shows the spectrum
	bool blockMap       = false;
	
	/*----------------------------------------------------------------*/
    /* !!!INITIAL PARAMETERS to quickly change analysis!!! */
    /*----------------------------------------------------------------*/
    /* ROOT plotting initialization */	
    int nBins1=500;
    nBins = 500;
	numHistos = 5; // Change this number if we want to change number of slices!
	threshold = 0.16; // Change this number if we want to change our slice threshold energy!
	maxEnergy=0.345; // Energies from spectrum generation
	minEnergy=0.018; 
	maxTime = 500.0;
	maxField = 0.2; // Few enough neutrons can penetrate past 0.2 that this should be set.
	spacing = (maxEnergy - threshold)/((numHistos-1));
	
	/* Yield finding initialization */
	binBreak[0] = 40.0; // Change these numbers to change the position of bin dividers (only 3dip for now)
	binBreak[1] = 60.0;                                                     
	binBuff = 0.0; // Change this number to modify the buffer space on either side of the bin
	gBins[0] = 0.0;
	gBins[5] = maxTime;
	for (int i = 0; i < 2; i++) { // this loop will calculate bin positions. Trust me!
		for (int j = 0; j < 2; j++) {
			gBins[(j+1)+2*i] = binBreak[i] + (double)(2*j-1)*binBuff;
		}
	}
	
	/* Chi2 Minimization guesses -- variables declared here since we need chi2Dim declared previously */
	chi2Dim = 6;
	numSteps = 5;
	double chi2Bounds[2*chi2Dim];
	double c2Mapper[chi2Dim];
	double chi2Vals[(int)pow(2,chi2Dim)+1];
	
	/* Chi2 fitting dimensions table: 
	 * 1 = Energy (power)
	 * 2 = Energy (threshold)
	 * 3 = Angle  (power) 
	 * 4 = Upscatter (prob) */
	//double c2Guess[chi2Dim] = {2.15, 0.088, 1.10, 1.3, 0.0,0.00}; //"Old," not as good numbers
	//double c2Guess[chi2Dim] = {2.0, 0.0, 1.0, 2.0, 0.0, 0.0}; // No weighting or truncation
	double c2Guess[chi2Dim] = {2.30448, 0.076071, 1.264079, 2.0,1.0, 0.0160763}; //"Nathan's suggested" numbers
	//double c2Step[chi2Dim] = {0.2, 0.02, 0.3, 0.1, 0.0,0.016};
	double c2Step[chi2Dim] = {0.0, 0.00, 0.0, 0.1, 0.0,0.0};
		
	/*----------------------------------------------------------------*/
    /* !!!LOAD DATA -- actual function begins here!!! */
    /*----------------------------------------------------------------*/
    /* Input/Output data */
    if((argc != 3) && (argc != 4) && (argc != 5)) {
        printf("Error! Usage: ./block_spectrum_map fname_dag fname_blok (hold_t) (fitting file) \n");
        return 0;
    }
	if (argc >= 4) {
		hold_t = std::stod(argv[3]);
		/*if (pseOn) {
			hold_t = hold_t+200.0;
		}*/
	} else {
		hold_t = 0.0;
		decayOn = false;
	}
	
	//nBins = 700+hold_t;
	printf("%d\n",nBins);
	if (argc >= 5) {
		chiFit = true; // Turn on chi^2 fitting if we have a fitting file.
	}
	
	/* Start by loading fname_dag */ 
    std::ifstream binfile(argv[1], std::ios::in | std::ios::binary);
    if(!binfile.is_open()) {
        printf("Error! Could not open file %s\n", argv[1]);
        return 1;
    }
    
    test=0;
    /* Read Dagger data in from buffer */
    while(!binfile.eof()) {
        binfile.read(buf, buff_len);
        if(binfile.eof()) {
            break;
        }   
		event_dag.time = *((double *)(&buf[0] + sizeof(unsigned int)));
        event_dag.energy = *((double *)(&buf[0] + sizeof(unsigned int)+sizeof(double)));
        event_dag.perp_energy = *((double *)(&buf[0] + sizeof(unsigned int)+2*sizeof(double)));
        for(int jj = 0; jj < 3; jj++) {
			event_dag.pos[jj] = *((double *)(&buf[0] + sizeof(unsigned int)+(3 + jj)*sizeof(double)));
		}
		event_dag.z_offset = *((double *)(&buf[0] + sizeof(unsigned int) + 6*sizeof(double)));
		for(int jj = 0; jj < 4; jj++) {
			event_dag.nHits[jj] = *((int *)(&buf[0] + sizeof(unsigned int) + 7*sizeof(double) + jj*sizeof(int)));
		}
		event_dag.theta = *((double *)(&buf[0] + sizeof(unsigned int) + 7*sizeof(double) + 4*sizeof(int)));
		if (debugMode) {
			if (test < 10) { 
				printf( " %e %e %e \n ", event_dag.time, event_dag.energy, event_dag.perp_energy);
				printf( " %e %e %e %e \n ", event_dag.pos[0], event_dag.pos[1], event_dag.pos[2], event_dag.z_offset);
				printf( " %u %u %u %u \n\n ", event_dag.nHits[0], event_dag.nHits[1], event_dag.nHits[2], event_dag.nHits[3]);
				test = test+1;
			} else if (quickLoad) {
				break;
			}
        }
        events_dag.push_back(event_dag);
    }
    binfile.close();
    length_dag = events_dag.size();
    printf("Read %u (dagger) Events!\n", length_dag);
    
    /* Now load fname_block events */
    std::ifstream binfile2(argv[2], std::ios::in | std::ios::binary);
    if(!binfile2.is_open()) {
        printf("Error! Could not open file %s\n", argv[2]);
        return 1;
    }
    
    test=0;
    /* Read Block data in from buffer */
    while(!binfile2.eof()) {
        binfile2.read(buf_blok, blok_len);
        if(binfile2.eof()) {
            break;
        }
        event_blok.time = *((double *)(&buf_blok[0] +sizeof(unsigned int)));
        event_blok.energy = *((double *)(&buf_blok[0] + +sizeof(unsigned int) + sizeof(double)));
        for(int jj = 0; jj < 4; jj++) {
			event_blok.nHits[jj] = *((int *)(&buf_blok[0] + sizeof(unsigned int) + 2*sizeof(double) + jj*sizeof(int)));
		}
		for(int jj = 0; jj < 3; jj++) {
			event_blok.position[jj] = *((double *)(&buf_blok[0] + sizeof(unsigned int) + (2 + jj)*sizeof(double) + 4*sizeof(int)));
		}
        event_blok.theta = *((double *)(&buf_blok[0] + sizeof(unsigned int) + 5*sizeof(double) + 4*sizeof(int)));               
        if (debugMode) {
	        if (test < 10) {
				printf( " %e %e \n ", event_blok.time, event_blok.energy);
				printf( " %u %u %u %u \n ", event_blok.nHits[0], event_blok.nHits[1], event_blok.nHits[2], event_blok.nHits[3]);
				printf( " %e %e %e \n ", event_blok.position[0], event_blok.position[1], event_blok.position[2]);
				printf( " %e \n ", event_blok.theta);
				test = test+1;
			} else if (quickLoad) {
				break;
			}
		}	
        events_blok.push_back(event_blok);
    }
    binfile2.close();
    length_blok = events_blok.size();
    printf("Read %u (block) Events!\n", length_blok);
    
    /* Have to declare these random number strings here (since they're variable length) */
    double* randDeathTimes = new double[length_dag+length_blok];
    double* randAbsProb = new double[2*length_blok];
    double* randDagProb = new double[length_dag];
    /* Create random numbers for weighting absorption spectrum */ 
    for(unsigned long i = 0; i < length_dag+length_blok; i++) {
		if (decayOn) {
			randDeathTimes[i] = -877.7*log(nextU01());
		} else {
			randDeathTimes[i] = DBL_MAX;
		}
    }
    for(unsigned long i = 0; i < 2*length_blok; i++) {
		randAbsProb[i] = nextU01();
		if (randAbsProb[i] < 0.0 || randAbsProb[i] > 1.0) {  // check and make sure random number generator doesn't crash
			i = i-1;
		}
	}
	for(unsigned long i = 0; i  < length_dag; i++) {
		randDagProb[i] = nextU01();
		//randDagProb[i] = 1.0;
		/*if (i > length_dag) {
			printf("%d,%f \n",i,randDagProb[i]);
		}*/
	}
    /*----------------------------------------------------------------*/
	/* !!!INITIALIZE all ROOT histograms!!! */
	/*----------------------------------------------------------------*/
	TH1D* timeDagHist = new TH1D("timeDag", "hit time (dagger)", nBins1,-150.0,maxTime); 
	TH1D* eneDagHist = new TH1D("eneDag", "energy (dagger)",300, 0.0, maxEnergy+minEnergy);
	TH1D* thetaDagRoot = new TH1D("thetaDag", "theta (dagger)",300, 0, 1.7); 
	TH2D* blockHistRoot = new TH2D("blockHist", "hit time vs. energy", 30, 0, maxTime, 30, minEnergy, maxEnergy+minEnergy);
	TH1D* timeBlockHist = new TH1D("timeBlock", "hit time (block)", nBins,0,maxTime);
	TH1D* eneBlockHist = new TH1D("eneBlock", "energy (block)",300, minEnergy, maxEnergy+minEnergy);
	TH1D* tBlockHist[numHistos];
	TH1D* fieldBHist[numHistos];
	
	/* Initialize slice histograms */
	/* Start with zeroth */
	sprintf(tBName,"tBlock0");
	sprintf(fBName,"fBlock0");
	sprintf(tBTitle,"hit time E > %f m", threshold);
	sprintf(fBTitle,"hit field E > %f m", threshold);
	tBlockHist[0] =  new TH1D(tBName, tBTitle, 300,0,maxTime+hold_t);
	fieldBHist[0] =  new TH1D(fBName, fBTitle, 300,0,maxField);
	/* Loop across the rest of the generated histograms */
	for (int i=1;i<numHistos;i++){
		sprintf(tBName,"tBlock%d",i);
		sprintf(fBName,"fBlock%d",i);
		sprintf(tBTitle,"hit time %f < E < %f m", (i-1)*spacing + threshold, i*spacing + threshold);
		sprintf(fBTitle,"hit field %f < E < %f m", (i-1)*spacing + threshold, i*spacing + threshold);
		tBlockHist[i] =  new TH1D(tBName, tBTitle, nBins,0,maxTime+hold_t);
		fieldBHist[i] =  new TH1D(fBName, fBTitle, 300,0,maxField);
	}
	/* Create hitpos and B field scatter plots*/
	TGraph2D* hitPos = new TGraph2D(length_blok);
		hitPos->SetName("hitPosition");
		hitPos->SetTitle("Hit Position");
	TGraph2D* hitB = new TGraph2D(length_blok);
		hitB->SetName("hitBField");
		hitB->SetTitle("Hit B Field");
	TMultiGraph* posMaps = new TMultiGraph();
	/* Define two-exponential decay function */
	//TF1* expo2 = new TF1("expo2","[0]*exp(-x/[1]) + [2]*exp(-x/[3])");
	TF1* expo2 = new TF1("expo2","[0]*exp(-x/[1])");
	//TF1* expo2 = new TF1("expo2","[0]*exp(-x/[1]) + [2]*exp(-x/877.7)");
		expo2->SetParameters(0.5, 800);//, 100, 100);
		expo2->SetParLimits(0,0.0,1.0);
		expo2->SetParLimits(1,0.0,2000.0);
		expo2->SetParName(0,"A");
		expo2->SetParName(1,"tau");
		//expo2->SetParLimits(2,0.0,10000);
		//expo2->SetParLimits(3,0.0,1000.0);	
	/* Chi2 ROOT initializations */
	TFile* rootFile;
	TH1D* dataHist;
	printf("Base ROOT histograms loaded\n");
	
	/*----------------------------------------------------------------*/
	/* !!!Run Chi^2 algorithm!!! */
	/*----------------------------------------------------------------*/
	if (chiFit) {
		printf("Chi^2 Fitting ON! \n");
		/* Load file from ROOT, and force it to fit inside our ROI */
		sprintf(rootName,"%s",argv[4]);
		printf("Loading %s\n",rootName);
		rootFile = TFile::Open(rootName);
		dataHist = (TH1D*)rootFile->Get("summedDipS");
		dataHist->GetXaxis()->SetRangeUser(42.0,220.0);
		
		/* IDEA: Step around in an N-D hypercube, and minimize differences in chi2 to find global minima */		
		for (int step=0;step<numSteps;step++) {
			/* Fill the chi2 bounds before doing anything else */
			for (int dim=0; dim<chi2Dim; dim++){
				chi2Bounds[2*dim]   = c2Guess[dim] - c2Step[dim];
				chi2Bounds[2*dim+1] = c2Guess[dim] + c2Step[dim];
			}
			
			/* Start by calculating middle of hypercube */
			eneMin = c2Guess[0];
			eThMin = c2Guess[1];
			theMin = c2Guess[2];
			theCut = c2Guess[3];
			absMin = c2Guess[4];
			dagCut = c2Guess[5];
			
			/* Load our relevant dagger events on a TH1D */
			/* Fill only one histogram at a time*/
			TH1D* timeDagFitting = new TH1D("timeDagFit", "hit time (dagger)", nBins,0,maxTime);
			blockCounter = 0;
			for (int i = 0; i<length_dag; i++) { 
				dag_buffer = events_dag[i];
				upscatter = false;
				while (dag_buffer.nHits[3] > 0) {
					if (randAbsProb[blockCounter] < absMin) {
						upscatter = true;
						blockCounter = blockCounter + dag_buffer.nHits[3];
						dag_buffer.nHits[3] = 0;
					} else {
						upscatter = false;
						blockCounter = blockCounter + 1;
						dag_buffer.nHits[3] = dag_buffer.nHits[3] - 1;
					}
				}
				/* check if we get absorbed by the bad part of the dagger */ 
			    if(dag_buffer.pos[0] > 0) {
			        zeta = 0.5 - sqrt(dag_buffer.pos[0]*dag_buffer.pos[0] + pow(fabs(dag_buffer.pos[2] - dag_buffer.z_offset) - 1.0, 2));
			    }
			    else {
					zeta = 1.0 - sqrt(dag_buffer.pos[0]*dag_buffer.pos[0] + pow(fabs(dag_buffer.pos[2] - dag_buffer.z_offset) - 0.5, 2));
			    }
			    boronProb = zeta > dagCut ? 1 : zeta/dagCut;
			    /* calculate weighting */
				if (!upscatter) {
					angleTheta = dag_buffer.theta;
					energyCm = dag_buffer.energy/(GRAV*MASS_N);
					if (thetaOn) {
						histoWeight = (pow(cos(angleTheta),theMin)*sin(angleTheta))/(cos(angleTheta)*sin(angleTheta)) * pow((energyCm),eneMin)/energyCm;
					} else {
						histoWeight = (pow(cos(angleTheta),theMin)*sin(angleTheta)) * pow((energyCm),eneMin)/energyCm;
					}
					if ((std::isfinite(histoWeight)) && (energyCm > eThMin) && (angleTheta < theCut) && (randDagProb[i] < boronProb)) {
						timeDagFitting->Fill(dag_buffer.time, histoWeight);
					}
				}
			}
			timeDagFitting->GetXaxis()->SetRangeUser(42.0,220.0);
			chi2Vals[(int)pow(2,chi2Dim)+1] = dataHist->Chi2Test(timeDagFitting,"UW CHI2/NDF");
			delete timeDagFitting;
			printf("%f %f %f %f %f \n",eneMin,eThMin,theMin,absMin,chi2Vals[(int)pow(2,chi2Dim)+1]);
			chi2Min = chi2Vals[(int)pow(2,chi2Dim)+1];
			minDex = pow(2,chi2Dim)+1;		
			
			/* Look at each endpoint. Potentially skip some? */
			for (int b=0;b<(pow(2,chi2Dim));b++) { 
							
				/* Calculate our fit parameters here. Doing binary math */
				eneMin = chi2Bounds[0 + ((int)floor(b/1) % 2)]; //0 or 1
				eThMin = chi2Bounds[2 + ((int)floor(b/2) % 2)]; //2 or 3
				theMin = chi2Bounds[4 + ((int)floor(b/4) % 2)]; //4 or 5
				theCut = chi2Bounds[6 + ((int)floor(b/8) % 2)]; //6 or 7
				absMin = chi2Bounds[8 + ((int)floor(b/16) % 2)]; //8 or 9
				dagCut = chi2Bounds[10 + ((int)floor(b/32) % 2)]; //10 or 11
				
				/* Skip if our scan range is 0 */
				for (int d=0; d < chi2Dim; d++ ) {
					if((chi2Bounds[2*d + ((int)floor(b/pow(2,d)) % 2)] == c2Guess[d]) && (((int)floor(b/pow(2,d)) % 2)==1)) {
						chi2Vals[b] = chi2Vals[(int)pow(2,chi2Dim)+1];
					} else {
						chi2Vals[b] = 0;
					}
				}
				if (chi2Vals[b] == 0) {
					/* Load our relevant dagger events on a TH1D */
					/* Fill only one histogram at a time*/
					TH1D* timeDagFitting = new TH1D("timeDagFit", "hit time (dagger)", nBins,0,maxTime);
					blockCounter = 0;
					for (int i = 0; i<length_dag; i++) { 
						dag_buffer = events_dag[i];
						/* Check if we upscatter off the block */
						upscatter = false;
						while (dag_buffer.nHits[3] > 0) {
							if (randAbsProb[blockCounter] < absMin) {
								upscatter = true;
								blockCounter = blockCounter + dag_buffer.nHits[3];
								dag_buffer.nHits[3] = 0;
							} else {
								upscatter = false;
								blockCounter = blockCounter + 1;
								dag_buffer.nHits[3] = dag_buffer.nHits[3] - 1;
							}
						}
						
						/* check if we get absorbed by the bad part of the dagger */ 
					    if(dag_buffer.pos[0] > 0) {
					        zeta = 0.5 - sqrt(dag_buffer.pos[0]*dag_buffer.pos[0] + pow(fabs(dag_buffer.pos[2] - dag_buffer.z_offset) - 1.0, 2));
					    }
					    else {
							zeta = 1.0 - sqrt(dag_buffer.pos[0]*dag_buffer.pos[0] + pow(fabs(dag_buffer.pos[2] - dag_buffer.z_offset) - 0.5, 2));
					    }
					    boronProb = zeta > dagCut ? 1 : zeta/dagCut;
					    /* Calculate weighting */
						if (!upscatter) {
							angleTheta = dag_buffer.theta;
							energyCm = dag_buffer.energy/(GRAV*MASS_N);
							if (thetaOn) {
								histoWeight = (pow(cos(angleTheta),theMin)*sin(angleTheta))/(cos(angleTheta)*sin(angleTheta)) * pow((energyCm),eneMin)/energyCm;
							} else {
								histoWeight = (pow(cos(angleTheta),theMin)*sin(angleTheta)) * pow((energyCm),eneMin)/energyCm;
							}
							if ((std::isfinite(histoWeight)) && (energyCm > eThMin) && (angleTheta < theCut) && (randDagProb[i] < boronProb)) {
								timeDagFitting->Fill(dag_buffer.time, histoWeight);
							}
						}
					}
					timeDagFitting->GetXaxis()->SetRangeUser(42.0,220.0);
					chi2Vals[b] = dataHist->Chi2Test(timeDagFitting,"UW CHI2/NDF");
					delete timeDagFitting;
					printf("%f %f %f %f %f \n",eneMin,eThMin,theMin,absMin,chi2Vals[b]);
				}
				/* Check where minimum is */
				if (chi2Vals[b] < chi2Min) {
					chi2Min = chi2Vals[b];
					minDex = b;
				}
			}
			/* Divide step size by 2 */
			for (int dim=0; dim<chi2Dim; dim++) {
				c2Step[dim] = c2Step[dim]/2;
				/* Move center (if necessary) */
				if (minDex != (pow(2,chi2Dim)+1)) {
					c2Guess[dim] = (c2Guess[dim] -(1 - 2*((int)floor(minDex/pow(2,dim)) % 2)) * c2Step[dim]);
				}
			}
		}
		printf("Min Chi^2: %f %f %f %f %f %f %f\n",eneMin,eThMin,theMin,theCut,absMin,dagCut,chi2Min);
		/*------------------------------------------------------------*/
		/* ----------------MAP CHI^2 TANGENT SURFACES---------------- */
		/*------------------------------------------------------------*/
		if (chiMap) {
			printf("Generating tangent surface plots! \n");			
			/* Increase scan size (increase our smallest step by dimensions) */
			for (int dim=0; dim<chi2Dim; dim++) {
				c2Step[dim] = pow(2,numSteps-1)*c2Step[dim];
				c2Mapper[dim] = c2Guess[dim];
			}
			
			/* Scan across our dimensions to produce our chi2 plots */			
			for (int dimx=0; dimx<chi2Dim; dimx++) {
				for (int dimy=0; dimy<dimx; dimy++) { // Scanning for only a given quadrant
					if (c2Step[dimx]!=0.0 && c2Step[dimy]!=0.0) { //Ignore if there's no step!
						
						/* Load our file for scanning */
						sprintf(c2Filename,"chi2Fit%d%d.csv",dimx,dimy);
						printf("Opening File %s\n",c2Filename);
						chi2OutFile = fopen(c2Filename,"w");
						
						/* Loop across our steps */
						for (int stepx=-numSteps; stepx<=numSteps; stepx++) {
						for (int stepy=-numSteps; stepy<=numSteps; stepy++) {
							
							/* Figure out what our positions are */
							for (int dim=0;dim<chi2Dim;dim++) {
								c2Mapper[dim] = c2Guess[dim];
							}
							c2Mapper[dimx] = c2Guess[dimx] + (double)stepx*c2Step[dimx];
							c2Mapper[dimy] = c2Guess[dimy] + (double)stepy*c2Step[dimy];
							
							eneMin = c2Mapper[0];
							eThMin = c2Mapper[1];
							theMin = c2Mapper[2];
							theCut = c2Mapper[3];
							absMin = c2Mapper[4];
							dagCut = c2Mapper[5];		
							
							/* Fill only one histogram at a time */
							TH1D* timeDagFitting = new TH1D("timeDagFit", "hit time (dagger)", nBins,0,maxTime);
							blockCounter = 0;
							for (int i = 0; i<length_dag; i++) { 
								dag_buffer = events_dag[i];
								upscatter = false;
								while (dag_buffer.nHits[3] > 0) {
									if (randAbsProb[blockCounter] < absMin) {
										upscatter = true;
										blockCounter = blockCounter + dag_buffer.nHits[3];
										dag_buffer.nHits[3] = 0;
									} else {
										upscatter = false;
										blockCounter = blockCounter + 1;
										dag_buffer.nHits[3] = dag_buffer.nHits[3] - 1;
									}
								}
								/* check if we get absorbed by the bad part of the dagger */ 
							    if(dag_buffer.pos[0] > 0) {
							        zeta = 0.5 - sqrt(dag_buffer.pos[0]*dag_buffer.pos[0] + pow(fabs(dag_buffer.pos[2] - dag_buffer.z_offset) - 1.0, 2));
							    }
							    else {
									zeta = 1.0 - sqrt(dag_buffer.pos[0]*dag_buffer.pos[0] + pow(fabs(dag_buffer.pos[2] - dag_buffer.z_offset) - 0.5, 2));
							    }
							    boronProb = zeta > dagCut ? 1 : zeta/dagCut;
							    /* Calculate weighting */
								if (!upscatter) {
									angleTheta = dag_buffer.theta;
									energyCm = dag_buffer.energy/(GRAV*MASS_N);
									if (thetaOn) {
										histoWeight = (pow(cos(angleTheta),theMin)*sin(angleTheta))/(cos(angleTheta)*sin(angleTheta)) * pow((energyCm),eneMin)/energyCm;
									} else {
										histoWeight = (pow(cos(angleTheta),theMin)*sin(angleTheta)) * pow((energyCm),eneMin)/energyCm;
									}
									if ((std::isfinite(histoWeight)) && (energyCm > eThMin) && (angleTheta < theCut) && (randDagProb[i] < boronProb)) {
										timeDagFitting->Fill(dag_buffer.time, histoWeight);
									}
								}
							}
							timeDagFitting->GetXaxis()->SetRangeUser(42.0,220.0);
							chi2Vals[0] = dataHist->Chi2Test(timeDagFitting,"UW CHI2/NDF");
							delete timeDagFitting;
							printf("%f %f %f %f %f \n",eneMin,eThMin,theMin,absMin,chi2Vals[0]);
							fprintf(chi2OutFile, "%lg %lg %lg \n",c2Mapper[dimx],c2Mapper[dimy],chi2Vals[0]);
						}
						}
						fclose(chi2OutFile);	
					}
				}
			}
		}
	}
	
	/* Load our energies from the minimum chi2guesses */
	eneMin = c2Guess[0];
	eThMin = c2Guess[1];
	theMin = c2Guess[2];
	theCut = c2Guess[3];
	absMin = c2Guess[4];
	dagCut = c2Guess[5];
	
	/*----------------------------------------------------------------*/
	/* !!!FILL all ROOT dagger histograms!!! */
	/*----------------------------------------------------------------*/
	printf("Filling Dagger histograms! \n");
	blockCounter = 0;
	for (int i=0; i<length_dag; i++) {
		dag_buffer=events_dag[i];
		
		if (pseOn) {
			dag_buffer.time = dag_buffer.time - (180.0 + hold_t);
		}
		
		/* Calculate out the weighting for theta */
		angleTheta = dag_buffer.theta;
		if (thetaOn) { 
			histoWeight = (pow(cos(angleTheta),theMin)*sin(angleTheta))/(cos(angleTheta)*sin(angleTheta)) * pow((dag_buffer.energy),eneMin)/dag_buffer.energy;
		} else {
			histoWeight = (pow(cos(angleTheta),theMin)*sin(angleTheta)) * pow((dag_buffer.energy),eneMin)/dag_buffer.energy;
		}
		/* check if we get absorbed by the bad part of the dagger */ 
	    if(dag_buffer.pos[0] > 0) {
	        zeta = 0.5 - sqrt(dag_buffer.pos[0]*dag_buffer.pos[0] + pow(fabs(dag_buffer.pos[2] - dag_buffer.z_offset) - 1.0, 2));
	    }
	    else {
			zeta = 1.0 - sqrt(dag_buffer.pos[0]*dag_buffer.pos[0] + pow(fabs(dag_buffer.pos[2] - dag_buffer.z_offset) - 0.5, 2));
	    }
	    boronProb = zeta > dagCut ? 1 : zeta/dagCut;
		/* Check if we're (a) loading a living neutron and (b) above threshold */
		if (((double)dag_buffer.time + (50.0 + hold_t)) < randDeathTimes[i] && (dag_buffer.energy/(GRAV*MASS_N) > eThMin) && (angleTheta < theCut) && dag_buffer.time > 0.0 && (randDagProb[i] < boronProb)){
			upscatter = false;
			/* Check if we've been killed by the block */
			while (dag_buffer.nHits[3] > 0) {
				if (randAbsProb[blockCounter] < absMin) {
					upscatter = true;
					blockCounter = blockCounter + dag_buffer.nHits[3];
					dag_buffer.nHits[3] = 0;
				} else {
					upscatter = false;
					blockCounter = blockCounter + 1;
					dag_buffer.nHits[3] = dag_buffer.nHits[3] - 1;
				}
			}
			
		    /* Calculate weighting */
			/* Fill Necessary histograms. */
			if (!upscatter) { 
				timeDagHist->Fill((double)dag_buffer.time, histoWeight);
				if (sanityCheck) {
					thetaDagRoot->Fill((double(dag_buffer.theta)),histoWeight);
				}
				
			}
			if (sanityCheck) {
				eneDagHist->Fill((dag_buffer.energy)/(GRAV*MASS_N),histoWeight);
			}
		}	
	}
	/* Calculate Yields based on histogram integrals */
	neutronYield = timeDagHist->IntegralAndError(0,nBins,yieldError,"");
	
	for (int i=0; i<sizeof(groupYield)/sizeof(*groupYield); i++) {
		groupYield[i] = timeDagHist->IntegralAndError(timeDagHist->GetXaxis()->FindBin(gBins[2*i]),timeDagHist->GetXaxis()->FindBin(gBins[2*i+1]),groupError[i],"");
	}
	
	printf("Filling Block histograms! \n");
	blockCounter = 0;
	/* Fill Necessary Block Data */
	for (int i=0; i<length_blok; i++) {
		block_buffer=events_blok[i];
		deathTime = randDeathTimes[i+length_dag];
		angleTheta = block_buffer.theta;
		upscatter = false;
		/* Check Block scattering */
		while (block_buffer.nHits[3] > 1) {
			if (randAbsProb[blockCounter] < absMin) {
				if (randAbsProb[blockCounter] < 0.0 || randAbsProb[blockCounter] > 1.0) { // fix buffer overflows or something?
					upscatter = true;
					block_buffer.nHits[3] = 1;
					blockCounter = 0;
				//	brokenCounter = brokenCounter + 1;
				} else {
					upscatter = true;
					blockCounter = blockCounter + block_buffer.nHits[3];
					block_buffer.nHits[3] = 1;
				}
				
			} else {
				upscatter = false;
				blockCounter = blockCounter + 1;
				block_buffer.nHits[3] = block_buffer.nHits[3] - 1;
			}
		}
		/* Load histograms that don't scatter and don't die */
		if ((!upscatter) && ((double)block_buffer.time < deathTime) && (block_buffer.position[2] < -1.0) && (block_buffer.energy/(GRAV*MASS_N) > eThMin) && (angleTheta < theCut)) {
			
			if (thetaOn) { 
				histoWeight = (pow(cos(angleTheta),theMin)*sin(angleTheta))/(cos(angleTheta)*sin(angleTheta)) * pow((block_buffer.energy),eneMin)/block_buffer.energy;
			} else {
				histoWeight = (pow(cos(angleTheta),theMin)*sin(angleTheta))*pow((block_buffer.energy),eneMin)/block_buffer.energy;
			}
			
			if (blockMapTot) {
				blockHistRoot->Fill(block_buffer.time, (block_buffer.energy)/(GRAV*MASS_N),histoWeight);
				timeBlockHist->Fill((double)block_buffer.time,histoWeight);
				eneBlockHist->Fill((block_buffer.energy)/(GRAV*MASS_N),histoWeight);
			}
			field_block = fieldstrength(block_buffer.position);
			if (blockMapSlice || bFieldSlice) {
				if ((block_buffer.energy)/(GRAV*MASS_N) <= threshold) {
						tBlockHist[0]->Fill((double)block_buffer.time,histoWeight);				
						fieldBHist[0]->Fill(field_block,histoWeight);
				} else {
					for (int j=1; j < numHistos; j++) {
						if (((block_buffer.energy)/(GRAV*MASS_N) > (((double)j-1.0)*spacing + threshold)) && ((block_buffer.energy)/(GRAV*MASS_N) <= (((double)j)*spacing + threshold))) {
							tBlockHist[j]->Fill((double)block_buffer.time,histoWeight);
							fieldBHist[j]->Fill(field_block,histoWeight);
						}
					}
				}			
			}
			if (blockMap) {
				hitPos->SetPoint(i,block_buffer.position[0],block_buffer.position[1],block_buffer.position[2]);
				if (block_buffer.position[2] > maxCoords[2]) {
					for (int q = 0; q<3; q++) {
						maxCoords[q] = block_buffer.position[q];
					}
				}
				if (block_buffer.position[2] < minCoords[2]) {
					for (int q = 0; q<3; q++) {
						minCoords[q] = block_buffer.position[q];
					}
				}
				hitB->SetPoint(i,block_buffer.position[0],block_buffer.position[1],field_block);
			}
		} else if (!upscatter && (double)block_buffer.time < deathTime && (block_buffer.energy/(GRAV*MASS_N) > eThMin) && (angleTheta < theCut)) { 
			/* Check if we have a block hit outside of the block */
			brokenCounter = brokenCounter + 1;
		} else if (upscatter  && (double)block_buffer.time < deathTime && (block_buffer.energy/(GRAV*MASS_N) > eThMin) && (angleTheta < theCut)) {
			doubleCounts = doubleCounts + 1;
		}
	}
	
	printf("Broken Trajectories: %d \nDouble Hits: %d\n", brokenCounter, doubleCounts);
	printf("\n---------------------------------------------------------------------\n");
	printf(" !!! OUTPUT DATA !!! ");
	printf("\n---------------------------------------------------------------------\n");
			
	/*----------------------------------------------------------------*/
	/* !!!PLOTTING AND ANALYSIS!!! */
	/*----------------------------------------------------------------*/
	/* CANVASSES: 
	 * 1: block mapping (block vs. hit time, block hits, etc.
	 * 2: block slice data
	 * 3: comparison between loaded data and MC simulation (on the dagger)
	 * 4: B field slice data
	 * 5: Chi^2 fitting region plot
	 * 6: energy/angle fitting sanity check
	 * 7: position on the block */
	 
	if (blockMapTot) {

		/* Generation of histograms on canvas */ 
		TCanvas* c1 = new TCanvas("c1","c1",750,500);
		TPad* pad1c1 = new TPad("p1","Energy vs Time",0.03,0.03,0.645,0.97,21);
		TPad* pad2c1 = new TPad("p2","Arrival Time",0.675,0.03,0.97,0.485,21);
		TPad* pad3c1 = new TPad("p3","Neutron Energy",0.675,0.515,0.97,0.97,21);
		pad1c1->Draw();
		pad2c1->Draw();
		pad3c1->Draw();	
		pad1c1->cd();
		blockHistRoot->Draw("CONT4Z");
		pad2c1->cd();
		timeBlockHist->Scale(1/timeBlockHist->Integral());
		timeBlockHist->Draw();
		timeBlockHist->Fit("expo2","WL M","",200.0,240.0+hold_t);
		/*timeBlockHist->Fit("expo2","QL","",200.0,240.0+hold_t);
		 
		blockPop[0]    = timeBlockHist->GetFunction("expo2")->GetParameter(0);
		blockPopErr[0] = timeBlockHist->GetFunction("expo2")->GetParError(0);
		blockSlope[0]  = timeBlockHist->GetFunction("expo2")->GetParameter(1);
		blockErr[0]    = timeBlockHist->GetFunction("expo2")->GetParError(1);
		blockPop[1]    = timeBlockHist->GetFunction("expo2")->GetParameter(2);
		blockPopErr[1] = timeBlockHist->GetFunction("expo2")->GetParError(2);
		blockSlope[1]  = timeBlockHist->GetFunction("expo2")->GetParameter(3);
		blockErr[1]    = timeBlockHist->GetFunction("expo2")->GetParError(3);*/
		pad3c1->cd();
		eneBlockHist->Scale(1/eneBlockHist->Integral());
		eneBlockHist->Draw();
		c1->Update();
		
		/* Block Map Scan Output data */
		if (blockSlope[0] < blockSlope[1]) {
			fastIndex = 0;
			slowIndex = 1;
		} else {
			fastIndex = 1;
			slowIndex = 0;
		}
		printf("For all neutrons:\n");
		blockPopErr[fastIndex] = blockPopErr[fastIndex]/blockPop[fastIndex]; // Convert to % error
		blockPop[fastIndex] = blockPop[fastIndex]/blockSlope[fastIndex]*(exp(-(90.0+hold_t)/blockSlope[fastIndex])-exp(-50.0*blockSlope[fastIndex]));
		blockPopErr[fastIndex] = blockPopErr[fastIndex]*blockPop[fastIndex];
		printf("   Fast Decay Constant: %e (%e)\n", blockSlope[fastIndex], blockErr[fastIndex]);
		printf("   Fast Population:     %e (%e)\n", blockPop[fastIndex], blockPopErr[fastIndex]);
		blockPopErr[slowIndex] = blockPopErr[slowIndex]/blockPop[slowIndex]; // Convert to % error
		blockPop[slowIndex] = blockPop[slowIndex]/blockSlope[slowIndex]*(exp(-(90.0+hold_t)/blockSlope[slowIndex])-exp(-50.0/blockSlope[slowIndex]));
		blockPopErr[slowIndex] = blockPopErr[slowIndex]*blockPop[slowIndex];
		printf("   Slow Decay Constant: %e (%e)\n", blockSlope[slowIndex], blockErr[slowIndex]);
		printf("   Slow Population:     %e (%e)\n", blockPop[slowIndex], blockPopErr[slowIndex]);
		printf("\n---------------------------------------------------------------------\n");
	}
	
	char* tPadName = new char[4];
	char* tPadTitle = new char[48];
	
	/* Slice data on the block */
	if (blockMapSlice) {
		TCanvas* c2 = new TCanvas("c2","c2",750,1000);
		c2->cd();
		TPad* padc2[numHistos];
					
		/* Start with zeroth pad */
		sprintf(tPadName,"p%d",0);
		sprintf(tPadTitle,"Arrival Time E < %f m", threshold);
		padMin = 0.01;
		padMax = (1.0/numHistos) - 0.01;
		padc2[0] = new TPad(tPadName,tPadTitle,0.03,padMin,0.97,padMax,21);
		padc2[0]->Draw();
		/* Create the remaining pads */
		for (int i=1;i<numHistos;i++) { 
			sprintf(tPadName,"p%d",i);
			sprintf(tPadTitle,"Arrival Time %f < E < %f m", (i-1)*spacing + threshold, i*spacing + threshold);
			padMax = ((i+1.0)/numHistos) - 0.01;
			padMin = (i+0.0)/numHistos + 0.01;
		
			padc2[i] = new TPad(tPadName,tPadTitle,0.03,padMin,0.97,padMax,21);
			padc2[i]->Draw();
		} 
		/* Fit our remaining histograms */
		for (int i=0; i<numHistos;i++) {
			padc2[i]->SetLogy();
			padc2[i]->cd();
			tBlockHist[i]->Draw();
			tBlockHist[i]->Fit("expo2","QL","",50.0,90.0+hold_t);
			slicePop[i][0]    = tBlockHist[i]->GetFunction("expo2")->GetParameter(0);
			slicePopErr[i][0] = tBlockHist[i]->GetFunction("expo2")->GetParError(0);
			sliceSlope[i][0]  = tBlockHist[i]->GetFunction("expo2")->GetParameter(1);
			sliceErr[i][0]    = tBlockHist[i]->GetFunction("expo2")->GetParError(1);
			slicePop[i][1]    = tBlockHist[i]->GetFunction("expo2")->GetParameter(2);
			slicePopErr[i][1] = tBlockHist[i]->GetFunction("expo2")->GetParError(2);
			sliceSlope[i][1]  = tBlockHist[i]->GetFunction("expo2")->GetParameter(3);
			sliceErr[i][1]    = tBlockHist[i]->GetFunction("expo2")->GetParError(3);
			if (sliceSlope[i][0] < sliceSlope[i][1]) {
				fastIndex = 0;
				slowIndex = 1;
			} else {
				fastIndex = 1;
				slowIndex = 0;
			}
			printf("For energy slice between %f and %f cm:\n",(i-1)*spacing + threshold, i*spacing + threshold);
					
			slicePopErr[i][fastIndex] = slicePopErr[i][fastIndex]/slicePop[i][fastIndex]; // Convert to % error
			slicePop[i][fastIndex] = slicePop[i][fastIndex]/sliceSlope[i][fastIndex]*(exp(sliceSlope[i][fastIndex]*(40.0+hold_t))-exp(sliceSlope[i][fastIndex]*50.0));
			slicePopErr[i][fastIndex] = slicePopErr[i][fastIndex]*slicePop[i][fastIndex];
			printf("   Fast Decay Constant: %e (%e)\n", sliceSlope[i][fastIndex], sliceErr[i][fastIndex]);
			printf("   Fast Population:     %e (%e)\n", slicePop[i][fastIndex], slicePopErr[i][fastIndex]);
						
			slicePopErr[i][slowIndex] = slicePopErr[i][slowIndex]/slicePop[i][slowIndex]; // Convert to % error
			slicePop[i][slowIndex] = slicePop[i][slowIndex]/sliceSlope[i][slowIndex]*(exp(sliceSlope[i][slowIndex]*(90.0+hold_t))-exp(sliceSlope[i][slowIndex]*50.0));
			slicePopErr[i][slowIndex] = slicePopErr[i][slowIndex]*slicePop[i][slowIndex];
			printf("   Slow Decay Constant: %e (%e)\n", sliceSlope[i][slowIndex], sliceErr[i][slowIndex]);
			printf("   Slow Population:     %e (%e)\n\n", slicePop[i][slowIndex], slicePopErr[i][slowIndex]);

			gStyle->SetOptFit(1);
		}
		c2->Update();
	}
	
	/* Data vs. MC simulations */
	if (dataMCComp) {
		
		TCanvas* c3 = new TCanvas("c3","c3",750,500);
		TPad* pad1c3 = new TPad("p1","Arrival Time (Dagger)",0.03,0.03,0.97,0.97,21);
		pad1c3->Draw();
		pad1c3->cd();
		timeDagHist->Scale(1/timeDagHist->Integral());
		timeDagHist->Draw();
		if (chiFit) {
			if (dataHist->GetSize() == timeDagHist->GetSize()){
				dataHist->Scale(1/dataHist->Integral());
				dataHist->SetLineColor(kRed);
				dataHist->Draw("SAME");
			}
		}	
		if (saveMCFit) {
			timeDagHist->SaveAs("./timeDagHist.root");
			printf("\n !!!SAVING DATA AS timeDagHist.root!!! \n");
		}

		if (drainingTimeOn) {timeDagHist->Fit("expo","Q","",45.0,59.0);
			drainFit[0] = timeDagHist->GetFunction("expo")->GetParameter(1);
			drainErr[0] = timeDagHist->GetFunction("expo")->GetParError(1);
			timeDagHist->Fit("expo","Q","",65.0,90.0);
			drainFit[1] = timeDagHist->GetFunction("expo")->GetParameter(1);
			drainErr[1] = timeDagHist->GetFunction("expo")->GetParError(1);
			printf("\n Draining Time: \n%e (Group 2)\n%e (Group 3)\n\n",drainFit[0],drainFit[1]);
		}		
		c3->Update();
	}
		
	if (bFieldSlice) {
		TCanvas* c4 = new TCanvas("c4","c4",750,1000);
		TPad* padc4[numHistos];
		
		sprintf(fPadName,"p%d",0);
		sprintf(fPadTitle,"B Field E < %f m", threshold);
		padMin = 0.01;
		padMax = 0.99;
		//padMax = (1.0/numHistos) - 0.01;
		padc4[0] = new TPad(fPadName,fPadTitle,0.03,padMin,0.97,padMax,21);
		padc4[0]->Draw();
		/*for (int i=1;i<numHistos;i++) { 
			sprintf(fPadName,"p%d",i);
			sprintf(fPadTitle,"B Field %f < E < %f m", (i-1)*spacing + threshold, i*spacing + threshold);
			padMax = ((i+1.0)/numHistos) - 0.01;
			padMin = (i+0.0)/numHistos + 0.01;
			padc4[i] = new TPad(fPadName,fPadTitle,0.03,padMin,0.97,padMax,21);
			padc4[i]->Draw();
		} */
		double normB = fieldBHist[3]->Integral();
		for (int i=0; i<numHistos;i++) {
			padc4[0]->cd();
			fieldBHist[i]->Scale(1/fieldBHist[i]->Integral());
			//fieldBHist[i]->Scale(1/normB);
			sprintf(fPadTitle,"B Field %f < E < %f m", (i-1)*spacing + threshold, i*spacing + threshold);
			fieldBHist[i]->SetTitle(fPadTitle);
			fieldBHist[i]->Draw("Same");		
		}
		c4->Update();
	}
	
	if (chiMap) {
		TCanvas* c5 = new TCanvas("c5","c5",750,750);
		
		/* Find how many chi2 pads we require */ 
		int nPads=0;
		for (int i=0;i<chi2Dim;i++) { 
			for (int j=0;j<i;j++) {
				if (c2Step[i]!=0 && c2Step[j]!=0) {
					nPads=nPads+1;
				}
			}
		}
		/* Create our pads */
		TPad* padc5[nPads];
		for (int i=0;i<nPads;i++) { 
			
				sprintf(tPadName,"p%d",i);
				sprintf(tPadTitle,"chi2 fit, variation %d", i);
				padMax = ((i+1.0)/nPads) - 0.01;
				padMin = (i+0.0)/nPads + 0.01;
				padc5[i] = new TPad(tPadName,tPadTitle,0.03,padMin,0.97,padMax,21);
				padc5[i]->Draw();
		} 
		/* Load our data */
		TGraph2D* chi2Map[nPads];
		int mapC=0;
		for (int i=0;i<chi2Dim;i++) { 
			for (int j=0;j<i;j++) {
				if (c2Step[i]!=0.0 && c2Step[j]!=0.0) {
					sprintf(c2Filename,"chi2Fit%d%d.csv",i,j);
					chi2Map[mapC] = new TGraph2D(c2Filename);
						chi2Map[mapC]->SetName(c2Filename);
						chi2Map[mapC]->SetTitle("Chi^2 Fit");
					padc5[mapC]->cd();
					chi2Map[mapC]->Draw("CONT4Z");
					mapC=mapC+1;
				}
			}
		}
		c5->Update();
	}
	
	if (sanityCheck) {
		TCanvas* c6 = new TCanvas("c6","c6",750,750);
		TPad* pad1c6 = new TPad("p1","Energy (all neutrons)",0.03,0.03,0.97,0.485,21);	
		TPad* pad2c6 = new TPad("p2","Theta (dagger)",0.03,0.515,0.97,0.97,21);
		pad1c6->Draw();
		pad2c6->Draw();
		pad1c6->cd();	
		eneDagHist->Scale(1/eneDagHist->Integral());
		eneDagHist->SaveAs("./eneDagHist.root");
		eneDagHist->Draw();
		pad2c6->cd();
		thetaDagRoot->Scale(1/thetaDagRoot->Integral());
		thetaDagRoot->SaveAs("./thetaDagHist.root");
		thetaDagRoot->Draw();
		c6->Update();
	}
	
	if (blockMap) {
		TCanvas* c7 = new TCanvas("c7","c7",500,500);
		TPad* pad1c7 = new TPad("p1","Block Hit Position",0.03,0.03,0.97,0.97,21);
		pad1c7->Draw();
		pad1c7->cd();
		hitPos->GetHistogram()->SetMinimum(-1.46);
		hitPos->GetHistogram()->SetMaximum(-1.42);
		hitPos->GetYaxis()->SetRangeUser(0.10,0.15);
		hitPos->GetXaxis()->SetRangeUser(0.20,0.25);
	 	hitPos->Draw("PCOL");
	 	c7->Update();
	 	printf("MAX: (%f, %f, %f);\nMIN: (%f, %f, %f).\n",maxCoords[0],maxCoords[1],maxCoords[2],minCoords[0],minCoords[1],minCoords[2]);
	}
	
	/*
	TCanvas* c8 = new TCanvas("c8","c8",750,750);
	double totalHits, percentage;
	TH1D* percentDag = new TH1D("percentDag", "% of neutrons hitting block",300, minEnergy, maxEnergy + minEnergy);
	
	/*for (int ii = 0; ii < length_dag; ii++)	{
		totalHits = eneDagHist->GetBinContent(ii)+eneDagOnBlock->GetBinContent(ii);
		if(totalHits != 0.0){
			percentage = eneDagHist->GetBinContent(ii)/totalHits;
		} else{ 
			percentage = 1.0;
		}
		percentDag->SetBinContent(ii, 1.0-percentage);
	}
	percentDag->Draw();*/
	//eneDagHist->Add(eneBlockHist);
			
	printf("The neutron yield is: %e (+- %e) !\n", neutronYield, yieldError);
	printf("This includes a yield of: \n");
	printf("%e (+- %e) (%f +- %f percent) (group 1) \n", groupYield[0],groupError[0], groupYield[0]/neutronYield, groupError[0]/neutronYield);
	printf("%e (+- %e) (%f +- %f percent) (group 2) \n", groupYield[1],groupError[1], groupYield[1]/neutronYield, groupError[1]/neutronYield);
	printf("%e (+- %e) (%f +- %f percent) (group 3) \n", groupYield[2],groupError[2], groupYield[2]/neutronYield, groupError[2]/neutronYield);
	
	
	theApp.Run();
    
    delete[] randDeathTimes;
    delete[] randAbsProb;
    delete[] buf;
	
    return 0;
}


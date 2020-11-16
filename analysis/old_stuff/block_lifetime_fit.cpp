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
/*#include "TH1D.h"
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
*/
extern "C" {
    #include "xorshift.h"
}

#define MASS_N 1.674927471e-27
#define GRAV 9.80665e0
#define JTONEV 6.2415091e27
#define HBAR 1.054571800e-34

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

/* ---------------------- Here's the actual MAIN program!! -----------------------------------*/
int main(int argc, char** argv) {
  //  ROOT::EnableThreadSafety();
    initxorshift();
    //TApplication theApp("App",0,0);
    
    /*----------------------------------------------------------------*/
    /* !!!DECLARATION OF VARIABLES!!! */
    /*----------------------------------------------------------------*/
    /* Buffer declarations */
    const size_t buff_len = 4 + 8*8 + 4*4 + 4; // Change if MC data changes
	const size_t blok_len = (4 + 6*8 + 4*4 + 4);
	
	/* Load events structures*/
    std::vector<dag_evt> events_dag;
    dag_evt event_dag, dag_buffer;
	std::vector<block_evt> events_blok;
    block_evt event_blok, block_buffer;
		
    /* IO variables */    
    int test; // test debug parameter
    int length_dag, length_blok;
    
    if(argc <= 3) {
        printf("Error! Usage: ./block_lifetime_fit folder t1 t2 t... \n");
        return 0;
    }
    int numArgs = argc - 2; // number of arguments -- find now 
    double hold_t[numArgs];
	
    /* Analysis and weighting variables */
    bool upscatter = false;
    int blockCounter, length_od, length_ob = 0;
    double angleTheta, zeta, boronProb;
    double neutronYield[numArgs] = {0.0};
    double neutronError[numArgs] = {0.0};
    double histoWeight, eneMin, eThMin, theMin, theCut, absMin, dagCut;
    double groupYield[3][numArgs] = {0.0};
    double groupError[3][numArgs] = {0.0};
    double binBreak[2];
    
    /* Yield finding initialization */
	binBreak[0] = 40.0; // Change these numbers to change the position of bin dividers (only 3dip for now)
	binBreak[1] = 60.0;
		
	/* Chi2 fitting dimensions table: 
	 * 1 = Energy (power)
	 * 2 = Energy (threshold)
	 * 3 = Angle  (power) 
	 * 4 = Angle  (cutoff)
	 * 5 = Upscatter (prob)
	 * 6 = B10 Bad Region (m) */
	int chi2Dim = 6;
	//double c2Guess[chi2Dim] = {2.15, 0.088, 1.10, 2.0, 1.0, 0.0}; //"Old," not as good numbers
	//double c2Guess[chi2Dim] = {2.0, 0.0, 1.0, 2.0, 1.0, 0.0}; // No weighting or truncation
	double c2Guess[chi2Dim] = {2.30448, 0.076071, 1.264079, 2.0, 1.0e-3, 0.0160763}; //"Nathan's suggested" numbers
	//double c2Guess[chi2Dim] = {3.504480, 0.048571, 2.564079, 1.0, 1.0, 0.052076}; //I found these numbers!	 -- They're super unphysical!
    char dagFile[48];
    char blokFile[48];
        
    /*----------------------------------------------------------------*/
    /* !!!BOOLEAN LOGIC to quickly change program analysis!!! */
    /*----------------------------------------------------------------*/
    bool debugMode = false; // Turns on debug mode
    bool quickLoad = false; // Only loads the first 10 dagger and block events
    bool thetaOn   = false; // Modify thetaOn to change plotting E and theta
    	
	//printf("%d\n",numArgs);
	eneMin = c2Guess[0];
	eThMin = c2Guess[1];
	theMin = c2Guess[2];
	theCut = c2Guess[3];
	absMin = c2Guess[4];
	dagCut = c2Guess[5];
		
	for (int ii=0; ii<numArgs; ii++) {
		
		char* buf;
		char* buf_blok;
		buf = new char[buff_len];
		buf_blok = new char[blok_len];

		sprintf(dagFile,"%sdag.out.%ss", argv[1],argv[ii+2],'c');

		hold_t[ii]= std::stod(argv[ii+2]);
		std::ifstream binfile(dagFile, std::ios::in | std::ios::binary);
		
		if(!binfile.is_open()) {
			printf("Error! Could not open file %s\n", dagFile);
			return 1;
		}

		int test=1;  		
		while(!binfile.eof()) {
			binfile.read(buf, buff_len);
			if(binfile.eof()) {
				break;
			}
       
			event_dag.time = *((double *)(&buf[0] + sizeof(unsigned int)));
			event_dag.energy = *((double *)(&buf[0] + sizeof(unsigned int)+sizeof(double)));
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
	
	    length_dag = events_dag.size() - length_od;
	    printf("Read %u (dagger) Events in file %s!\n", length_dag, dagFile);
	   
		sprintf(blokFile,"%sblok.out.%ss", argv[1],argv[ii+2],'c');	    
	    std::ifstream binfile2(blokFile, std::ios::in | std::ios::binary);
		if(!binfile2.is_open()) {
			printf("Error! Could not open file %s\n", blokFile);
			return 1;
		}	
	    test = 1;
	       	    
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
	   	    
	    length_blok = events_blok.size() - length_ob;
	    printf("Read %u (block) Events in file %s!\n", length_blok, blokFile);
	    
	    /* Create random numbers for weighting absorption spectrum (dagger) */ 
	    double* randDeathTimes = new double[length_dag];
		double* randAbsProb = new double[2*length_blok];
	    double* randDagProb = new double[length_dag];
	    for(unsigned long i = 0; i < length_dag; i++) {
	        //randDeathTimes[i] = -877.7*log(nextU01());
	        randDeathTimes[i] = 99999999999;
	    }
	    for(unsigned long i = 0; i < 2 * length_blok; i++) {
			randAbsProb[i] = nextU01();
			if (randAbsProb[i] < 0.0 || randAbsProb[i] > 1.0) {  // check and make sure random number generator doesn't crash
				i = i-1;
			}
		}
		for(unsigned long i = 0; i  < length_dag; i++) {
			randDagProb[i] = nextU01();
		}
		
		/* Calculate Neutron weighting */				
		printf("Calculating Weighting!\n");
		blockCounter = 0;
		for (int jj=0; jj<length_dag; jj++) {
			dag_buffer=events_dag[jj];
			angleTheta = dag_buffer.theta;
			/* check if we get absorbed by the bad part of the dagger */ 
		    if(dag_buffer.pos[0] > 0) {
		        zeta = 0.5 - sqrt(dag_buffer.pos[0]*dag_buffer.pos[0] + pow(fabs(dag_buffer.pos[2] - dag_buffer.z_offset) - 1.0, 2));
		    }
		    else {
				zeta = 1.0 - sqrt(dag_buffer.pos[0]*dag_buffer.pos[0] + pow(fabs(dag_buffer.pos[2] - dag_buffer.z_offset) - 0.5, 2));
		    }
		    boronProb = zeta > dagCut ? 1 : zeta/dagCut;
		    
		    /* load our dagger things */
			if ((((double)dag_buffer.time + (50.0 + hold_t[ii])) < randDeathTimes[jj]) && (dag_buffer.energy/(GRAV*MASS_N) > eThMin) && (angleTheta < theCut) && dag_buffer.time > 0.0 && (randDagProb[jj] < boronProb)){
				upscatter = false;
				while (dag_buffer.nHits[3] > 0) {
					if (blockCounter > 2*length_blok) {
						printf("Error with blockCounter! Resetting random number string!\n");
						for(unsigned long i = 0; i < 2 * length_blok; i++) {
							randAbsProb[i] = nextU01();
						}
						blockCounter = 0;
					}
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
				if (!upscatter) {
					
					if (thetaOn) { 
						histoWeight = (pow(cos(angleTheta),theMin)*sin(angleTheta))/(cos(angleTheta)*sin(angleTheta)) * pow((dag_buffer.energy),eneMin)/dag_buffer.energy;
					} else {
						histoWeight = (pow(cos(angleTheta),theMin)*sin(angleTheta)) * pow((dag_buffer.energy),eneMin)/dag_buffer.energy;
					}
					if (std::isfinite(histoWeight)) {
						neutronYield[ii] = neutronYield[ii] + histoWeight;
						neutronError[ii] = neutronError[ii] + histoWeight*histoWeight;
						if (dag_buffer.time < binBreak[0]) { 
							groupYield[0][ii] = groupYield[0][ii] + histoWeight;
							groupError[0][ii] = groupError[0][ii] + histoWeight*histoWeight;
						} else if (dag_buffer.time >= binBreak[0] && dag_buffer.time < binBreak[1]) {
							groupYield[1][ii] = groupYield[1][ii] + histoWeight;
							groupError[1][ii] = groupError[1][ii] + histoWeight*histoWeight;
						} else if (dag_buffer.time >=binBreak[1]) {
							groupYield[2][ii] = groupYield[2][ii] + histoWeight;
							groupError[2][ii] = groupError[2][ii] + histoWeight*histoWeight;
						}
						
					}
		
				}
			}
		}
				
		/* normalize the errors */
		neutronError[ii] = sqrt(neutronError[ii]);
		for (int jj=0; jj<3; jj++) { 
			groupError[jj][ii] = sqrt(groupError[jj][ii]);
		}
		printf("Done Weighting!\n");
		
		/* Clear memory cache */
		events_dag.erase(events_dag.begin(),events_dag.end());
		events_blok.erase(events_blok.begin(),events_blok.end());
		delete[] randDeathTimes;
		delete[] randAbsProb;
		delete[] randDagProb;
		delete[] buf;
		delete[] buf_blok;

	}
		
	printf("Normalizing group weights!\n");
	double normalization = neutronYield[0];
	for(int ii =0 ;ii<numArgs;ii++){
		neutronYield[ii] = neutronYield[ii]/normalization;
		neutronError[ii] = neutronError[ii]/normalization;
		for(int jj=0; jj<3; jj++) {
			groupYield[jj][ii] = groupYield[jj][ii]/normalization;
			groupError[jj][ii] = groupError[jj][ii]/normalization;
		}
		
		
		printf("\n%e, %e, %e, %e\n", neutronYield[ii], groupYield[0][ii],groupYield[1][ii],groupYield[2][ii]);
		printf("%e, %e, %e, %e\n", neutronError[ii], groupError[0][ii],groupError[1][ii],groupError[2][ii]);
	}
			
	if (numArgs > 2) {
		printf("Saving output file! \n");
		
		FILE* outfile;
		outfile = fopen("block_yields.csv","w");
		for (int ii = 0; ii<numArgs; ii++){
			fprintf(outfile,"%f, ",hold_t[ii]);
			fprintf(outfile, "%f, %f, %f, %f, ",     neutronYield[ii], neutronError[ii], groupYield[0][ii], groupError[0][ii]);
			fprintf(outfile, "%f, %f, %f, %f, \n", groupYield[1][ii], groupError[1][ii], groupYield[2][ii] ,groupError[2][ii]);
		}
		fclose(outfile);
	} else {
		     
	    double totalLifetime;
	    double blockLifetime;
	    double totalError;
	    
	    double groupLife[3];
	    double dT = hold_t[0] - hold_t[1];
	    totalLifetime = (-dT)/log(1/(exp(dT/877.7)*neutronYield[1]/neutronYield[0]));
		blockLifetime = totalLifetime - 877.7;
		//blockLifetime = 1/(1/totalLifetime - 1/877.7);
		totalError = -dT / (neutronError[1]*neutronYield[1]/neutronYield[0]) / pow(log(1/(exp(-dT/877.7)*(neutronError[1]*neutronYield[1]/neutronYield[0]))),2);
		//totalError = -dT*sqrt((neutronError[1]/neutronYield[1])*(neutronError[1]/neutronYield[1]) + (neutronError[0]/neutronYield[0])*(neutronError[0]/neutronYield[0]))/(log(neutronYield[0]/neutronYield[1])*log(neutronYield[0]/neutronYield[1]));
		
		for (int ii=0;ii<3;ii++){
			groupLife[ii] = -dT/log(1/(exp(dT/877.7)*groupYield[ii][1]/groupYield[ii][0]));
		}
		
		printf("The total lifetime is: %f +- %f !\n", totalLifetime, totalError);
		printf("This corresponds to a block shift of %f !\n", blockLifetime);
		printf("Individual group lifetimes:\n   %f (group 1) \n   %f (group 2) \n   %f (group 3)\n", groupLife[0], groupLife[1], groupLife[2]);
	}
	
    return 0;
}


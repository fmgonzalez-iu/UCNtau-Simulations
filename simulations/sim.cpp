#include "../general/setup.h"
#include "../general/typedefs.h"

#include <vector>
#include <stdio.h>
#include <cmath>
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <getopt.h>
#include <cstring>

#include "./inc/data_io.hpp"
#include "./inc/track_gen.hpp"
#include "./inc/trackUCN.hpp"
#include "./inc/fields_nate.h"
#include "./inc/lyap.hpp"
#include "./inc/geometry.hpp"

extern "C" {
    #include "./inc/xorshift.h"
}

/*typedef struct fixedResult {
    double energy;
    double theta;
    double t;
    double ePerp;
    double x;
    double y;
    double z;
    double zOff;
    int nHit;
    int nHitHouseLow;
    int nHitHouseHigh;
    double eStart;
    double deathTime;
} fixedResult;*/



int main(int argc, char** argv) {
	// Main function of simulation

	//------------------------------------------------------------------
	int ierr = MPI_Init(&argc, &argv); // Multithreading initialization
	int nproc;
	int rank;

	ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	ierr = MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	//------------------------------------------------------------------

	int c; // Load options from command line
	char fName[1024];
	fName[0] = '\0';
	double dt = 0.0;
	double nTraj = 0;
	double holdT = HOLDTIME;

	while (1) {
		static struct option long_options[] = {
			{"file", required_argument, 0, 'f'},
			{"dt", required_argument, 0, 'd'},
			{"ntraj", required_argument, 0, 'n'},
			{"holdt", required_argument, 0, 'h'},// added this
			{0, 0, 0, 0},
		};
		// getopt_long stores the option index here. 
		int option_index = 0;

		c = getopt_long(argc, argv, "f:d:n:h:", long_options, &option_index);

		// Detect the end of the options.
		if(c == -1) {
			break;
		}

		switch(c) {
			case 'f':
				strncpy(fName, optarg, 1024-1);
				fName[1024-1] = '\0';
				break;
			case 'd':
				dt = atof(optarg);
				break;
			case 'n':
				nTraj = atoi(optarg);
				break;
			case 'h':
				holdT = atoi(optarg);
				break;
			case '?':
				// getopt_long already printed an error message.
				break;
			default:
				exit(1);
		}
	}
	
	if(fName[0] == '\0' || dt == 0 || nTraj == 0 || holdT == 0) {
		fprintf(stderr, "Error! Usage: ./arrival_time_fixed_eff --file=fName --dt=timestep --ntraj=N --holdt=time\n");
		exit(1);
	}
	//------------------------------------------------------------------
	// Load files
	
//	trace tr = readTrace("./xvals.bin", "./yvals.bin", "./zvals.bin");
	trace tr = readTrace(XFNAME, YFNAME, ZFNAME); // Load heating file
	if(tr.x == NULL || tr.y == NULL || tr.z == NULL) {
		fprintf(stderr, "Bad trace\n");
		return 2;
	}
	printf("num trace bins: %d\n", tr.num);

	char fNameRank[1024];
	snprintf(fNameRank, 1024, "%s%d", fName, rank); // Load binary outfile
	std::ofstream binfile(fNameRank, std::ios::out | std::ios::binary);

	initxorshift(INITBLOCK); // Load seed

//	printf("%d - %d\n", rank*(nTraj/nproc), (rank+1)*(nTraj/nproc));

	//------------------------------------------------------------------
	// Load trajectories
	std::vector<std::vector<double>> traj;
	traj.reserve(nTraj/nproc);
	for(int i = 0; i < nTraj; i++) {
		if(i >= rank*(nTraj/nproc) && i < (rank+1)*(nTraj/nproc)) {
			traj.push_back(TRACKGENERATOR(tr));
		} else {
			TRACKGENERATOR(tr);
		}
	}
	//printf("Running %f trajectories\n",nTraj);
	for(int i = 0; i < rank; i++) {
		jump();
	}

	//------------------------------------------------------------------
	// Actual computation
	
	if(TRACKANDPRINT) { // If running a tracking movie maker
		for(size_t ii = 0; ii < traj.size(); ii++) {
			trackAndPrint(traj[ii], dt, tr, ii, binfile);
		}
		binfile.close();
		ierr = MPI_Finalize();
		return 0;
	}
#ifdef MULTIDETSAVE // MULTIDETSAVE lets us run the classes of trajectory sims with 2 output files
	char fNameRank2[1024];
	snprintf(fNameRank2, 1024, "%s_2_%d", fName, rank); // Load binary outfile (2nd one)
	std::ofstream binfile2(fNameRank2, std::ios::out | std::ios::binary);

	//printf("%s\n", fNameRank);
	//printf("%s\n", fNameRank2);
	for(auto it = traj.begin(); it < traj.end(); it++) { // Run our trajectories
		// double lyap = calcLyap(*it, dt, tr, 0.0);
		// printf("%f\n", lyap);
		auto res = TRACKER(*it, dt, tr, holdT, binfile2);
		// writeFixedRes(binfile, res);
		WRITER(binfile, res);
	}
	//binfile2.close();
#else // This lets us run the trajectory with only one output file
	//printf("%s\n", fNameRank);
	for(auto it = traj.begin(); it < traj.end(); it++) { // Run our trajectories
		// double lyap = calcLyap(*it, dt, tr, 0.0);
		// printf("%f\n", lyap);
		auto res = TRACKER(*it, dt, tr, holdT);
		// writeFixedRes(binfile, res);
		WRITER(binfile, res);
	}
#endif
	
	binfile.close();
	ierr = MPI_Finalize();
    
	return 0;
}

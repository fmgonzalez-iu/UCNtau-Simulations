#pragma once

typedef struct measurement { // Results structure
    double val;
    double err;
} measurement;

typedef struct bounceResult {
	double energy;
	double theta;
	float times[NRECORDS];
	int refl[NRECORDS];
} bounceResult;

typedef struct noabsResult { // For finding dagger effective thickness
    double energy;
    double theta;
    float times[NRECORDS];
    float ePerps[NRECORDS];
    float zetas[NRECORDS];
} noabsResult;


//----------------------------------------------------------------------
typedef struct noabsCleanDagResult {
	double energy;
	double theta;
	float  times[NRECORDS];
	float  ePerps[NRECORDS];
	int    nClean[NRECORDS];
} noabsCleanDagResult;

typedef struct noabsCleanDagResultZ {
	double energy;
	double theta;
	float  times[NRECORDS];
	float  ePerps[NRECORDS];
	float  zetas[NRECORDS];
	int    nClean[NRECORDS];
} noabsCleanDagResultZ;

typedef struct noabsCleanResult {
	float energy;
	float theta;
	float t;
	float ePerp;
	float x;
	float y;
	float z;
}noabsCleanResult;

//----------------------------------------------------------------------

typedef struct reducedEvent {
	double energy;
	double theta;
	double time;
	double x;
	double y;
	double z;
	double weight;
		
} reducedEvent;


typedef struct fixedResultBlock { // Legacy structure for block analysis
	double time;
	double energy;
	//int nHits[4];
	int nHit;
    int nHitHouseLow;
    int nHitHouseHigh;
    int nHitBlock;
	//double position[3];
	double x;
	double y;
	double z;
	double theta;
} fixedResultBlock;

typedef struct fixedResultDag { // Legacy structure for block analysis
	double time;
    double energy;
    double perp_energy;
    //double pos[3];
    double x;
    double y;
    double z;
    double z_offset;
    double theta;
    //int nHits[4];
    int nHits;
    int nHitHouseLow;
    int nHitHouseHigh;
    int nHitBlock;
} fixedResultDag;

typedef struct fixedResult { // Assuming dagger thickness is fixed
    double energy;
    double theta;
    double t;
    double settlingT;
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
} fixedResult;

typedef struct cleanResult { // For looking at UCN hitting the cleaner
    float energy;
    float theta;
    float t;
    float x;
    float y;
    float z;
    int code;
} cleanResult;

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


typedef struct trace {
    double* x;
    double* y;
    double* z;
    int num;
} trace;
typedef struct lyapResult {
    float eStart;
    float eEnd;
    float theta;
    float lce;
} lyapResult;



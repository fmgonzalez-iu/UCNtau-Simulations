#include "../inc/file_reader.hpp"

/*----------------------------------------------------------------------
 * This code converts various binary outputs of the simulation files
 * into c++ vectors, which can then be parsed with calculations done.
 *
 * The various structures referenced here should be referenced in the 
 * "../general/typedefs.h" so that they are consistent between data 
 * structure types.
 *--------------------------------------------------------------------*/

std::vector<evt> readFile(const char* fName) {
	// Read file of NRECORDS bounces (for calibration of dagger efficiency)
	
	std::vector<evt> events;
	
	// Check that the file is readable
	std::ifstream binfile(fName, std::ios::in | std::ios::binary);
	if(!binfile.is_open()) { // Open file with given name
		fprintf(stderr, "Error! Could not open file %s\n", fName);
		return events;
	}

	evt event; // Parse file in chunks of the structure's size
	const size_t buff_len = 4 + 1*sizeof(double) + 2*NRECORDS*sizeof(float) + 4;
	char* buf = new char[buff_len];
	
	while(!binfile.eof()) {
		binfile.read(buf, buff_len);
		if(binfile.eof()) { //Breaks on last read of file (i.e. when 0 bytes are read and EOF bit is set)
			break;
		}
		// Load binary into structural form
		event.energy = *((double *)(&buf[0] + sizeof(unsigned int)));
		std::memcpy((void *)&event.times, (void *)(&buf[0] + sizeof(unsigned int) + sizeof(double)), NRECORDS*sizeof(float)); // copy vector
		std::memcpy((void *)&event.ePerp, (void *)(&buf[0] + sizeof(unsigned int) + sizeof(double) + NRECORDS*sizeof(float)), NRECORDS*sizeof(float));
		// Load single structure into event
		events.push_back(event);
	}
	binfile.close();
	printf("Read %lu Events!\n", events.size());

	delete[] buf;
	return events;
}

std::vector<fixedResultDag> readFileDagRes(const char* fName) { 
	// Old fixed result DAG from 2017 simulations
	
	std::vector<fixedResultDag> events;
	
	// Check that the file is readable
    std::ifstream binfile(fName, std::ios::in | std::ios::binary);
	if(!binfile.is_open()) { // Open file with given name
		fprintf(stderr, "Error! Could not open file %s\n", fName);
		return events;
	}
	
	fixedResultDag event; // Parse file in chunks of the structure's size
	const size_t buff_len =  4 + 8 * sizeof(double) + 4 * sizeof(int) + 4;
	char* buf = new char[buff_len];
		
	while(!binfile.eof()) {
		binfile.read(buf, buff_len);
		if(binfile.eof()) { // Breaks on last read of file (i.e. when 0 bytes are read and EOF bit is set)
			break;
		}

		// Load binary into structural form
		event.time        = *((double *)(&buf[0] + sizeof(unsigned int)));
		event.energy      = *((double *)(&buf[0] + sizeof(unsigned int) + 1*sizeof(double)));
		event.perp_energy = *((double *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(double)));
		event.x           = *((double *)(&buf[0] + sizeof(unsigned int) + 3*sizeof(double)));
		event.y           = *((double *)(&buf[0] + sizeof(unsigned int) + 4*sizeof(double)));
		event.z           = *((double *)(&buf[0] + sizeof(unsigned int) + 5*sizeof(double)));
		event.z_offset    = *((double *)(&buf[0] + sizeof(unsigned int) + 6*sizeof(double)));
		event.nHits          = *((int *)(&buf[0] + sizeof(unsigned int) + 7*sizeof(double)));
		event.nHitHouseLow   = *((int *)(&buf[0] + sizeof(unsigned int) + 7*sizeof(double) + 1*sizeof(int)));
		event.nHitHouseHigh  = *((int *)(&buf[0] + sizeof(unsigned int) + 7*sizeof(double) + 2*sizeof(int)));
		event.nHitBlock      = *((int *)(&buf[0] + sizeof(unsigned int) + 7*sizeof(double) + 3*sizeof(int)));
		event.theta       = *((double *)(&buf[0] + sizeof(unsigned int) + 7*sizeof(double) + 4*sizeof(int)));
		// Load single structure into event
		events.push_back(event);
	}
	binfile.close();
	printf("Read %lu Events!\n", events.size());

	delete[] buf; 
	return events;
}

std::vector<fixedResultBlock> readFileBlockRes(const char* fName) {
	// Block hit results
	
	std::vector<fixedResultBlock> events;
	
	// Check that the file is readable
	std::ifstream binfile(fName, std::ios::in | std::ios::binary);
	if(!binfile.is_open()) { // Open file with given name
		fprintf(stderr, "Error! Could not open file %s\n", fName);
		return events;
	}

	fixedResultBlock event; // Parse file in chunks of the structure's size
	const size_t buff_len =  4 + 6 * sizeof(double) + 2 * sizeof(int) + 4;
	char* buf = new char[buff_len];
	
	while(!binfile.eof()) {
		binfile.read(buf, buff_len);
		if(binfile.eof()) { // Breaks on last read of file (i.e. when 0 bytes are read and EOF bit is set)
			break;
		}

		// Load binary into structural form
		event.time = *((double *)(&buf[0] +sizeof(unsigned int)));
		event.energy = *((double *)(&buf[0] + +sizeof(unsigned int) + sizeof(double)));
		event.nHit   = *((int *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(double)));
		event.nHitHouseLow  = *((int *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(double) + 1*sizeof(int)));
		event.nHitHouseHigh = *((int *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(double) + 2*sizeof(int)));
		event.nHitBlock     = *((int *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(double) + 3*sizeof(int)));
		/*for(int jj = 0; jj < 4; jj++) {
			event.nHits[jj] = *((int *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(double) + jj*sizeof(int)));
		}*/
		event.x = *((double *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(double) + 4*sizeof(int)));
		event.y = *((double *)(&buf[0] + sizeof(unsigned int) + 3*sizeof(double) + 4*sizeof(int)));
		event.z = *((double *)(&buf[0] + sizeof(unsigned int) + 4*sizeof(double) + 4*sizeof(int)));
		/*for(int jj = 0; jj < 3; jj++) {
			event.position[jj] = *((double *)(&buf[0] + sizeof(unsigned int) + (2 + jj)*sizeof(double) + 4*sizeof(int)));
		}*/	        
		event.theta = *((double *)(&buf[0] + sizeof(unsigned int) + 5*sizeof(double) + 4*sizeof(int)));
		// Load single structure into event
		events.push_back(event);
	}
	binfile.close();
	printf("Read %lu Events!\n", events.size());

	delete[] buf; 
	return events;
}

std::vector<noabsResult> readFileNoAbsRes(const char* fName) {
	// Read noabsResult output files
    
	std::vector<noabsResult> events;
	
	// Check that the file is readable
	std::ifstream binfile(fName, std::ios::in | std::ios::binary);
	if(!binfile.is_open()) { // Open file with given name
		fprintf(stderr, "Error! Could not open file %s\n", fName);
		return events;
	}

	noabsResult event; // Parse file in chunks of the structure's size
	const size_t buff_len = 4 + 2 * sizeof(double) + 3 * NRECORDS * sizeof(float) + 4; // buff len has 4 bits on front/back
	char* buf = new char[buff_len]; 

	while(!binfile.eof()) { 
		binfile.read(buf, buff_len);
		if(binfile.eof()) { // Breaks on last read of file (i.e. when 0 bytes are read and EOF bit is set)
			break;
		}

		// Load binary into structural form
		event.energy = *((double *)(&buf[0] + sizeof(unsigned int)));
		event.theta = *((double *)(&buf[0] + sizeof(unsigned int) + sizeof(double)));
		std::memcpy((void *)&event.times, (void *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(double)), NRECORDS*sizeof(float));
		std::memcpy((void *)&event.ePerps, (void *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(double) + NRECORDS*sizeof(float)), NRECORDS*sizeof(float));
		std::memcpy((void *)&event.zetas, (void *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(double) + 2*NRECORDS*sizeof(float)), NRECORDS*sizeof(float));
		// Load single structure into event
		events.push_back(event);
	}
	binfile.close();
	printf("Read %lu Events!\n", events.size());

    delete[] buf;
    return events;
}

std::vector<noabsCleanDagResultZ> readFileNoAbsCleanDagResZ(const char* fName) {
// Read noabsResult output files
    
	std::vector<noabsCleanDagResultZ> events;
	
	// Check that the file is readable
	std::ifstream binfile(fName, std::ios::in | std::ios::binary);
	if(!binfile.is_open()) { // Open file with given name
		fprintf(stderr, "Error! Could not open file %s\n", fName);
		return events;
	}

	noabsCleanDagResultZ event; // Parse file in chunks of the structure's size
	const size_t buff_len = 4 + 2*sizeof(double) + 3*NRECORDS*sizeof(float) +NRECORDS*sizeof(int)+ 4; // buff len has 4 bits on front/back
	char* buf = new char[buff_len]; 

	while(!binfile.eof()) { 
		binfile.read(buf, buff_len);
		if(binfile.eof()) { // Breaks on last read of file (i.e. when 0 bytes are read and EOF bit is set)
			break;
		}

		// Load binary into structural form
		event.energy = *((double *)(&buf[0] + sizeof(unsigned int)));
		event.theta = *((double *)(&buf[0] + sizeof(unsigned int) + sizeof(double)));
		std::memcpy((void *)&event.times, (void *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(double)), NRECORDS*sizeof(float));
		std::memcpy((void *)&event.ePerps, (void *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(double) + NRECORDS*sizeof(float)), NRECORDS*sizeof(float));
		std::memcpy((void *)&event.zetas, (void *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(double) + 2*NRECORDS*sizeof(float)), NRECORDS*sizeof(float));
		std::memcpy((void *)&event.nClean, (void *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(double) + 3*NRECORDS*sizeof(float)), NRECORDS*sizeof(int));
		// Load single structure into event
		events.push_back(event);
	}
	binfile.close();
	printf("Read %lu Events!\n", events.size());

    delete[] buf;
    return events;
    
}

std::vector<noabsCleanDagResult> readFileNoAbsCleanDagRes(const char* fName) {
// Read noabsResult output files
    
	std::vector<noabsCleanDagResult> events;
	
	// Check that the file is readable
	std::ifstream binfile(fName, std::ios::in | std::ios::binary);
	if(!binfile.is_open()) { // Open file with given name
		fprintf(stderr, "Error! Could not open file %s\n", fName);
		return events;
	}

	noabsCleanDagResult event; // Parse file in chunks of the structure's size
	const size_t buff_len = 4 + 2*sizeof(double) + 2*NRECORDS*sizeof(float) +NRECORDS*sizeof(int)+ 4; // buff len has 4 bits on front/back
	char* buf = new char[buff_len]; 

	while(!binfile.eof()) { 
		binfile.read(buf, buff_len);
		if(binfile.eof()) { // Breaks on last read of file (i.e. when 0 bytes are read and EOF bit is set)
			break;
		}

		// Load binary into structural form
		event.energy = *((double *)(&buf[0] + sizeof(unsigned int)));
		event.theta = *((double *)(&buf[0] + sizeof(unsigned int) + sizeof(double)));
		std::memcpy((void *)&event.times, (void *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(double)), NRECORDS*sizeof(float));
		std::memcpy((void *)&event.ePerps, (void *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(double) + NRECORDS*sizeof(float)), NRECORDS*sizeof(float));
		std::memcpy((void *)&event.nClean, (void *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(double) + 2*NRECORDS*sizeof(float)), NRECORDS*sizeof(int));
		// Load single structure into event
		events.push_back(event);
	}
	binfile.close();
	printf("Read %lu Events!\n", events.size());

    delete[] buf;
    return events;
    
}

std::vector<fixedResult> readFileFixedRes(const char* fName) {
    // Read fixedResult output files
    
	std::vector<fixedResult> events;

	// Check that the file is readable
	std::ifstream binfile(fName, std::ios::in | std::ios::binary);
	if(!binfile.is_open()) { // Open file with given name
		fprintf(stderr, "Error! Could not open file %s\n", fName);
		return events;
	}

	fixedResult event; // Parse file in chunks of the structure's size	
	const size_t buff_len = sizeof(unsigned int) + 9*sizeof(double) + 3*sizeof(int) + 2*sizeof(double) + sizeof(unsigned int); // buff len has 4 bits on front/back 
	char* buf = new char[buff_len];
	
	while(!binfile.eof()) {
		binfile.read(buf, buff_len);
		if(binfile.eof()) { //Breaks on last read of file (i.e. when 0 bytes are read and EOF bit is set)
			break;
		}

		// Load binary into structural form
		event.energy        = *((double *)(&buf[0] + sizeof(unsigned int)));
		event.theta         = *((double *)(&buf[0] + sizeof(unsigned int) + 1*sizeof(double)));
		event.t             = *((double *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(double)));
		event.settlingT     = *((double *)(&buf[0] + sizeof(unsigned int) + 3*sizeof(double)));
		event.ePerp         = *((double *)(&buf[0] + sizeof(unsigned int) + 4*sizeof(double)));
		event.x             = *((double *)(&buf[0] + sizeof(unsigned int) + 5*sizeof(double)));
		event.y             = *((double *)(&buf[0] + sizeof(unsigned int) + 6*sizeof(double)));
		event.z             = *((double *)(&buf[0] + sizeof(unsigned int) + 7*sizeof(double)));
		event.zOff          = *((double *)(&buf[0] + sizeof(unsigned int) + 8*sizeof(double)));
		event.nHit          = *((int *)(&buf[0] + sizeof(unsigned int) + 9*sizeof(double)));
		event.nHitHouseLow  = *((int *)(&buf[0] + sizeof(unsigned int) + 9*sizeof(double) + 1*sizeof(int)));
		event.nHitHouseHigh = *((int *)(&buf[0] + sizeof(unsigned int) + 9*sizeof(double) + 2*sizeof(int)));
		event.eStart        = *((double *)(&buf[0] + sizeof(unsigned int) + 9*sizeof(double) + 3*sizeof(int)));
		event.deathTime     = *((double *)(&buf[0] + sizeof(unsigned int) + 9*sizeof(double) + 3*sizeof(int) + sizeof(double)));    
		// Load single structure into event
		events.push_back(event);

	}
	binfile.close();
	printf("Read %lu Events!\n", events.size());

	delete[] buf;
	return events;
}

std::vector<cleanResult> readFileCleanRes(const char* fName) {
	// Read fixedResult output files

	std::vector<cleanResult> events;

	// Check that the file is readable
	std::ifstream binfile(fName, std::ios::in | std::ios::binary);
	if(!binfile.is_open()) { // Open file with given name
		fprintf(stderr, "Error! Could not open file %s\n", fName);
		return events;
	}

	cleanResult event; // Parse file in chunks of the structure's size	
	const size_t buff_len = 4 + 6 * sizeof(double) + 1 * sizeof(int) + 4; // buff len has 4 bits on front/back
	char* buf = new char[buff_len];
	
	while(!binfile.eof()) {
		binfile.read(buf, buff_len);
		if(binfile.eof()) { //Breaks on last read of file (i.e. when 0 bytes are read and EOF bit is set)
			break;
		}
		
		// Load binary into structural form
		event.energy   = *((float *)(&buf[0] + sizeof(unsigned int)));
		event.theta    = *((float *)(&buf[0] + sizeof(unsigned int) + sizeof(float)));
		event.t        = *((float *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(float)));
		event.x        = *((float *)(&buf[0] + sizeof(unsigned int) + 3*sizeof(float)));
		event.y        = *((float *)(&buf[0] + sizeof(unsigned int) + 4*sizeof(float)));
		event.z        = *((float *)(&buf[0] + sizeof(unsigned int) + 5*sizeof(float)));
		event.code     = *((int *)(&buf[0] + sizeof(unsigned int) + 6*sizeof(double)));
		// Load single structure into event
		events.push_back(event);
		
	}
	binfile.close();
    printf("Read %lu Events!\n", events.size());

	delete[] buf;
	return events;
}

std::vector<noabsCleanResult> readFileNoAbsCleanRes(const char* fName) {
	// Read fixedResult output files

	std::vector<noabsCleanResult> events;

	// Check that the file is readable
	std::ifstream binfile(fName, std::ios::in | std::ios::binary);
	if(!binfile.is_open()) { // Open file with given name
		fprintf(stderr, "Error! Could not open file %s\n", fName);
		return events;
	}

	noabsCleanResult event; // Parse file in chunks of the structure's size	
	const size_t buff_len = 4 + 7*sizeof(float) + 4; // buff len has 4 bits on front/back
	char* buf = new char[buff_len];
	
	while(!binfile.eof()) {
		binfile.read(buf, buff_len);
		if(binfile.eof()) { //Breaks on last read of file (i.e. when 0 bytes are read and EOF bit is set)
			break;
		}
		
		// Load binary into structural form
		event.energy   = *((float *)(&buf[0] + sizeof(unsigned int)));
		event.theta    = *((float *)(&buf[0] + sizeof(unsigned int) + sizeof(float)));
		event.t        = *((float *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(float)));
		event.ePerp    = *((float *)(&buf[0] + sizeof(unsigned int) + 3*sizeof(float)));
		event.x        = *((float *)(&buf[0] + sizeof(unsigned int) + 4*sizeof(float)));
		event.y        = *((float *)(&buf[0] + sizeof(unsigned int) + 5*sizeof(float)));
		event.z        = *((float *)(&buf[0] + sizeof(unsigned int) + 6*sizeof(float)));
		
		// Load single structure into event
		events.push_back(event);
		//if (events.size() % 1000000 == 0){printf("Read %lu events\n", events.size());}
	}
	binfile.close();
    printf("Read %lu Events!\n", events.size());

	delete[] buf;
	return events;
}

std::vector<double> refHistFromRoot(const char* fName) {
	
	std::vector<double> outHist;
	
	TFile* rootFile = TFile::Open(fName);
	if (!rootFile->GetListOfKeys()->Contains("summedDipS")) {
		return outHist;
	}
	
	TH1D* rootHist = (TH1D*)rootFile->Get("summedDipS"); // Now we load our histogram
	
	// Want to scan across the range of DIPTIMES global.
	double dipEnds[NDIPS] = DIPTIMES;
	double minT = 0.0;
	if (NDIPS > 2) {
		minT = dipEnds[1]; // Cut off first (zeroth) dip
	}
	double maxT = minT+200.;//dipEnds[NDIPS-1]; // hardcoded 200s refHist
	if (maxT <= 0) {
		printf("ERROR: Unable to load refHist!");
		return outHist;
	}
	
	// Now we sum the histogram into 1s bins
	int time = 0; // time counter ( to divide into s bins)
	double binVal = 0.0; // Value of bin
	for (int i=0; i < (int)rootHist->GetNbinsX(); i++) {
		double lowedge = rootHist->GetXaxis()->GetBinLowEdge(i);
		if (lowedge < minT || lowedge > maxT) {
			continue;
		}
		binVal += fabs(rootHist->GetBinContent(i)); // sum up binValue
		if ((int)lowedge > time) {
			outHist.push_back(binVal); // Create a c++ vector
			time = (int)lowedge;
			binVal = 0.0;
		}
	}
	
	delete rootHist;
	rootFile->Close();
	//for (int i=0; i<outHist.size(); i++) {
	//	printf("%f,", outHist.at(i));
	//}
	return outHist;
}
	
std::vector<std::vector<double>> readFilePosition(const char* fName) {
	
	// Read Position (traj, time, x,y,z) output files. 

	std::vector<std::vector<double>> events;

	// Check that the file is readable
	std::ifstream binfile(fName, std::ios::in | std::ios::binary);
	if(!binfile.is_open()) { // Open file with given name
		fprintf(stderr, "Error! Could not open file %s\n", fName);
		return events;
	}

	std::vector<double> event; // Parse file in chunks of the structure's size	
	for (size_t ii = 0; ii < 4; ii++) {
		event.push_back(0); // Want the vector to contain 4 events
	}
	event.resize(4);
	const size_t buff_len = sizeof(unsigned int) + sizeof(long) + 4*sizeof(double) + sizeof(unsigned int); // buff len has 4 bits on front/back
	char* buf = new char[buff_len];
	long nTraj = 0;
	while(!binfile.eof()) {
		binfile.read(buf, buff_len);
		if(binfile.eof()) { //Breaks on last read of file (i.e. when 0 bytes are read and EOF bit is set)
			break;
		}
		
		nTraj += 1;
		// Load event data
		event[0] = *((double *)(&buf[0] + sizeof(unsigned int) + sizeof(long)));
		event[1] = *((double *)(&buf[0] + sizeof(unsigned int) + sizeof(long) + 1*sizeof(double)));
		event[2] = *((double *)(&buf[0] + sizeof(unsigned int) + sizeof(long) + 2*sizeof(double)));
		event[2] = *((double *)(&buf[0] + sizeof(unsigned int) + sizeof(long) + 3*sizeof(double)));
		
		// Load single structure into event
		events.push_back(event);
	}
	
	binfile.close();
    printf("Read %lu Events!\n", events.size());

	delete[] buf;
	return events;
	
	
	
	
}

	
	//dataHist = (TH1D*)rootFile->Get("summedDipS");
	//dataHist->GetXaxis()->SetRangeUser(60.0,140.0);
	
	


/*
std::vector<evt> readFile(const char* fName) {
    std::vector<evt> events;
    
    const size_t buff_len = 4 + 8*8 + 3*4 + 4;
    char* buf = new char[buff_len];
    std::ifstream binfile(fName, std::ios::in | std::ios::binary);
    if(!binfile.is_open()) {
        fprintf(stderr, "Error! Could not open file %s\n", fName);
        return events;
    }
    
    evt event;
    while(!binfile.eof()) {
        binfile.read(buf, buff_len);
        if(binfile.eof()) { //Breaks on last read of file (i.e. when 0 bytes are read and EOF bit is set)
            break;
        }
        if(*((unsigned int *)&buf[0]) != buff_len - 8) {
			fprintf(stderr, "Error! Wrong binary format!\n");
			exit(2);
		}
		
		
        event.energy = *((double *)(&buf[0] + sizeof(unsigned int)));
        event.theta = *((double *)(&buf[0] + sizeof(unsigned int) + sizeof(double)));
        event.time = *((double *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(double)));
        event.eperp = *((double *)(&buf[0] + sizeof(unsigned int) + 3*sizeof(double)));
        event.x = *((double *)(&buf[0] + sizeof(unsigned int) + 4*sizeof(double)));
        event.y = *((double *)(&buf[0] + sizeof(unsigned int) + 5*sizeof(double)));
        event.z = *((double *)(&buf[0] + sizeof(unsigned int) + 6*sizeof(double)));
        event.zoff = *((double *)(&buf[0] + sizeof(unsigned int) + 7*sizeof(double)));
        event.nhit = *((int *)(&buf[0] + sizeof(unsigned int) + 8*sizeof(double)));
        event.nhitBot = *((int *)(&buf[0] + sizeof(unsigned int) + 7*sizeof(double) + sizeof(int)));
        event.nhitTop = *((int *)(&buf[0] + sizeof(unsigned int) + 7*sizeof(double) + 2*sizeof(int)));        
        events.push_back(event);
    }
    binfile.close();
    
    printf("Read %lu Events!\n", events.size());

    delete[] buf;
    return events;
}
*/


/*int main(int argc, char** argv) {
    if(argc != 2) {
        printf("Error! Usage: ./minimal fname\n");
        return 1;
    }

    std::vector<evt> data = readFile(argv[1]);

    return 0;
}*/


/*typedef struct evt {
    double energy;
    double theta;
    float times[NRECORDS];
    float ePerp[NRECORDS];
} evt;

typedef struct noabsResult { // For finding dagger effective thickness
    double energy;
    double theta;
    float times[NRECORDS];
    float ePerps[NRECORDS];
    float zetas[NRECORDS];
} noabsResult;

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

typedef struct evt_in {
    double energy;
    float times[NRECORDS];
    float ePerp[NRECORDS];
} evt_in;

typedef struct evt_out {
    double energy;
    double theta;
    float times[NRECORDS];
    float ePerp[NRECORDS];
} evt_out;

typedef struct evt_th {
    double energy;
    double theta;
} evt_th;
*/

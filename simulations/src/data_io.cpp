#include "../inc/data_io.hpp"

//======================================================================
void writePos(std::ofstream &binfile, long traj, double t, std::vector<double> state) {
	// Write out the x-y-z position for a given trajectory/time.
	const size_t buff_len = sizeof(unsigned int) + sizeof(long) + 4*sizeof(double) + sizeof(unsigned int);
	char buf[buff_len];
	if(!binfile.is_open()) {
		fprintf(stderr, "Error! file closed\n");
		return;
	}
	*((unsigned int *)(&buf[0])) = buff_len - 2*sizeof(unsigned int);
	*((long *)(&buf[0] + sizeof(unsigned int))) = traj;
	*((double *)(&buf[0] + sizeof(unsigned int) + sizeof(long))) = t;
	*((double *)(&buf[0] + sizeof(unsigned int) + sizeof(long) + 1*sizeof(double))) = state[0];
	*((double *)(&buf[0] + sizeof(unsigned int) + sizeof(long) + 2*sizeof(double))) = state[1];
	*((double *)(&buf[0] + sizeof(unsigned int) + sizeof(long) + 3*sizeof(double))) = state[2];
	*((unsigned int *)(&buf[0] + sizeof(unsigned int) + sizeof(long) + 4*sizeof(double))) = buff_len - 2*sizeof(unsigned int);
	binfile.write(buf, buff_len);
	
}

void writeFixedRes(std::ofstream &binfile, fixedResult res) {
	// Save a result structure into a binary file
	
    const size_t buff_len = sizeof(unsigned int) + 9*sizeof(double) + 3*sizeof(int) + 2*sizeof(double) + sizeof(unsigned int);
    char buf[buff_len];
    if(!binfile.is_open()) {
        fprintf(stderr, "Error! file closed\n");
        return;
    }
    *((unsigned int *)(&buf[0])) = buff_len - 2*sizeof(unsigned int);
    *((double *)(&buf[0] + sizeof(unsigned int))) = res.energy;
    *((double *)(&buf[0] + sizeof(unsigned int) + 1*sizeof(double))) = res.theta;
    *((double *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(double))) = res.t;
    *((double *)(&buf[0] + sizeof(unsigned int) + 3*sizeof(double))) = res.settlingT;
    *((double *)(&buf[0] + sizeof(unsigned int) + 4*sizeof(double))) = res.ePerp;
    *((double *)(&buf[0] + sizeof(unsigned int) + 5*sizeof(double))) = res.x;
    *((double *)(&buf[0] + sizeof(unsigned int) + 6*sizeof(double))) = res.y;
    *((double *)(&buf[0] + sizeof(unsigned int) + 7*sizeof(double))) = res.z;
    *((double *)(&buf[0] + sizeof(unsigned int) + 8*sizeof(double))) = res.zOff;
    *((int *)(&buf[0] + sizeof(unsigned int) + 9*sizeof(double))) = res.nHit;
    *((int *)(&buf[0] + sizeof(unsigned int) + 9*sizeof(double) + 1*sizeof(int))) = res.nHitHouseLow;
    *((int *)(&buf[0] + sizeof(unsigned int) + 9*sizeof(double) + 2*sizeof(int))) = res.nHitHouseHigh;
    *((double *)(&buf[0] + sizeof(unsigned int) + 9*sizeof(double) + 3*sizeof(int))) = res.eStart;
    *((double *)(&buf[0] + sizeof(unsigned int) + 9*sizeof(double) + 3*sizeof(int) + 1*sizeof(double))) = res.deathTime;
    *((unsigned int *)(&buf[0] + sizeof(unsigned int) + 9*sizeof(double) + 3*sizeof(int) + 2*sizeof(double))) = buff_len - 2*sizeof(unsigned int);
    binfile.write(buf, buff_len);
}

void writeBounceRes(std::ofstream &binfile, bounceResult res) {
	// Save a result structure into a binary file
	
	const size_t buff_len = sizeof(unsigned int) + 2*sizeof(double) + NRECORDS*sizeof(float) + NRECORDS*sizeof(int) + sizeof(unsigned int);
    char buf[buff_len];
    if(!binfile.is_open()) {
        fprintf(stderr, "Error! file closed\n");
        return;
    }
    *((unsigned int *)(&buf[0])) = buff_len - 2*sizeof(unsigned int);
    *((double *)(&buf[0] + sizeof(unsigned int))) = res.energy;
    *((double *)(&buf[0] + sizeof(unsigned int) + 1*sizeof(double))) = res.theta;
    std::memcpy((void *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(double)), (void *)&(res.times[0]), NRECORDS*sizeof(float));
    std::memcpy((void *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(double) + NRECORDS*sizeof(float)), (void *)&(res.refl[0]), NRECORDS*sizeof(int));
    *((unsigned int *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(double) + NRECORDS*sizeof(float) + NRECORDS*sizeof(int))) = buff_len - 2*sizeof(unsigned int);
    binfile.write(buf, buff_len);
}
void writeNoAbsCleanDagRes(std::ofstream &binfile, noabsCleanDagResult res) {
	// Save a result structure into binary
	
	if(res.energy < 0.0) {
        return;
    }
    const size_t buff_len = sizeof(unsigned int) + 2*sizeof(double) + 2*NRECORDS*sizeof(float) +1*NRECORDS*sizeof(int)+ sizeof(unsigned int);
    char buf[buff_len];
    if(!binfile.is_open()) {
        fprintf(stderr, "Error! file closed\n");
        return;
    }
    *((unsigned int *)(&buf[0])) = buff_len - 2*sizeof(unsigned int);
    *((double *)(&buf[0] + sizeof(unsigned int))) = res.energy;
    *((double *)(&buf[0] + sizeof(unsigned int) + 1*sizeof(double))) = res.theta;
    std::memcpy((void *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(double)), (void *)&(res.times[0]), NRECORDS*sizeof(float));
    std::memcpy((void *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(double) + NRECORDS*sizeof(float)), (void *)&(res.ePerps[0]), NRECORDS*sizeof(float));
    std::memcpy((void *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(double) + 2*NRECORDS*sizeof(float)), (void *)&(res.nClean[0]), NRECORDS*sizeof(int));
    *((unsigned int *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(double) + 2*NRECORDS*sizeof(float)+1*NRECORDS*sizeof(int))) = buff_len - 2*sizeof(unsigned int);
    binfile.write(buf, buff_len);
	
}
void writeNoAbsCleanDagResZ(std::ofstream &binfile, noabsCleanDagResultZ res) {
	// Save a result structure into binary
	
	if(res.energy < 0.0) {
        return;
    }
    const size_t buff_len = sizeof(unsigned int) + 2*sizeof(double) + 3*NRECORDS*sizeof(float) +1*NRECORDS*sizeof(int)+ sizeof(unsigned int);
    char buf[buff_len];
    if(!binfile.is_open()) {
        fprintf(stderr, "Error! file closed\n");
        return;
    }
    *((unsigned int *)(&buf[0])) = buff_len - 2*sizeof(unsigned int);
    *((double *)(&buf[0] + sizeof(unsigned int))) = res.energy;
    *((double *)(&buf[0] + sizeof(unsigned int) + 1*sizeof(double))) = res.theta;
    std::memcpy((void *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(double)), (void *)&(res.times[0]), NRECORDS*sizeof(float));
    std::memcpy((void *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(double) + NRECORDS*sizeof(float)), (void *)&(res.ePerps[0]), NRECORDS*sizeof(float));
    std::memcpy((void *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(double) + 2*NRECORDS*sizeof(float)), (void *)&(res.zetas[0]), NRECORDS*sizeof(float));
    std::memcpy((void *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(double) + 3*NRECORDS*sizeof(float)), (void *)&(res.nClean[0]), NRECORDS*sizeof(int));
    *((unsigned int *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(double) + 3*NRECORDS*sizeof(float)+1*NRECORDS*sizeof(int))) = buff_len - 2*sizeof(unsigned int);
    binfile.write(buf, buff_len);
	
}


void writenoabsCleanRes(std::ofstream &binfile, noabsCleanResult res) {
	// Save a result structure into a binary file
	
    const size_t buff_len = sizeof(unsigned int) + 7*sizeof(float) + sizeof(unsigned int);
    char buf[buff_len];
    if(!binfile.is_open()) {
        fprintf(stderr, "Error! file closed\n");
        return;
    }
    *((unsigned int *)(&buf[0])) = buff_len - 2*sizeof(unsigned int);
    *((float *)(&buf[0] + sizeof(unsigned int))) = res.energy;
    *((float *)(&buf[0] + sizeof(unsigned int) + 1*sizeof(float))) = res.theta;
    *((float *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(float))) = res.t;
    *((float *)(&buf[0] + sizeof(unsigned int) + 3*sizeof(float))) = res.ePerp;
    *((float *)(&buf[0] + sizeof(unsigned int) + 4*sizeof(float))) = res.x;
    *((float *)(&buf[0] + sizeof(unsigned int) + 5*sizeof(float))) = res.y;
    *((float *)(&buf[0] + sizeof(unsigned int) + 6*sizeof(float))) = res.z;
    *((unsigned int *)(&buf[0] + sizeof(unsigned int) + 7*sizeof(float))) = buff_len - 2*sizeof(unsigned int);
    binfile.write(buf, buff_len);
}


void writeNoabsRes(std::ofstream &binfile, noabsResult res) {
	// Save a result structure into a binary file
	
    if(res.energy < 0.0) {
        return;
    }
    const size_t buff_len = sizeof(unsigned int) + 2*sizeof(double) + 3*NRECORDS*sizeof(float) + sizeof(unsigned int);
    char buf[buff_len];
    if(!binfile.is_open()) {
        fprintf(stderr, "Error! file closed\n");
        return;
    }
    *((unsigned int *)(&buf[0])) = buff_len - 2*sizeof(unsigned int);
    *((double *)(&buf[0] + sizeof(unsigned int))) = res.energy;
    *((double *)(&buf[0] + sizeof(unsigned int) + 1*sizeof(double))) = res.theta;
    std::memcpy((void *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(double)), (void *)&(res.times[0]), NRECORDS*sizeof(float));
    std::memcpy((void *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(double) + NRECORDS*sizeof(float)), (void *)&(res.ePerps[0]), NRECORDS*sizeof(float));
    std::memcpy((void *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(double) + 2*NRECORDS*sizeof(float)), (void *)&(res.zetas[0]), NRECORDS*sizeof(float));
    *((unsigned int *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(double) + 3*NRECORDS*sizeof(float))) = buff_len - 2*sizeof(unsigned int);
    binfile.write(buf, buff_len);
}

void writeCleanRes(std::ofstream &binfile, cleanResult res) {
	// Save a result structure into a binary file
	
    const size_t buff_len = sizeof(unsigned int) + 6*sizeof(float) + 1*sizeof(int) + sizeof(unsigned int);
    char buf[buff_len];
    if(!binfile.is_open()) {
        fprintf(stderr, "Error! file closed\n");
        return;
    }
    *((unsigned int *)(&buf[0])) = buff_len - 2*sizeof(unsigned int);
    *((float *)(&buf[0] + sizeof(unsigned int))) = res.energy;
    *((float *)(&buf[0] + sizeof(unsigned int) + 1*sizeof(float))) = res.theta;
    *((float *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(float))) = res.t;
    *((float *)(&buf[0] + sizeof(unsigned int) + 3*sizeof(float))) = res.x;
    *((float *)(&buf[0] + sizeof(unsigned int) + 4*sizeof(float))) = res.y;
    *((float *)(&buf[0] + sizeof(unsigned int) + 5*sizeof(float))) = res.z;
    *((int *)(&buf[0] + sizeof(unsigned int) + 6*sizeof(float))) = res.code;
    *((unsigned int *)(&buf[0] + sizeof(unsigned int) + 6*sizeof(float) + 1*sizeof(int))) = buff_len - 2*sizeof(unsigned int);
    binfile.write(buf, buff_len);
}

void writeLyapRes(std::ofstream &binfile, lyapResult res) {
	// Save a result structure into a binary file
	
    const size_t buff_len = sizeof(unsigned int) + 4*sizeof(float) + sizeof(unsigned int);
    char buf[buff_len];
    if(!binfile.is_open()) {
        fprintf(stderr, "Error! file closed\n");
        return;
    }
    *((unsigned int *)(&buf[0])) = buff_len - 2*sizeof(unsigned int);
    *((float *)(&buf[0] + sizeof(unsigned int))) = res.eStart;
    *((float *)(&buf[0] + sizeof(unsigned int) + 1*sizeof(float))) = res.eEnd;
    *((float *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(float))) = res.theta;
    *((float *)(&buf[0] + sizeof(unsigned int) + 3*sizeof(float))) = res.lce;
    *((unsigned int *)(&buf[0] + sizeof(unsigned int) + 4*sizeof(float))) = buff_len - 2*sizeof(unsigned int);
    binfile.write(buf, buff_len);
}

trace readTrace(const char *xfile, const char *yfile, const char *zfile) {
	// Read a trace file for heating
	
    trace t;
    t.x = NULL;
    t.y = NULL;
    t.z = NULL;
    t.num = -1;
    
    const size_t buff_len = 1*8;
    char* buf = new char[buff_len];
    
    std::ifstream binfileX(xfile, std::ios::in | std::ios::binary);
    std::ifstream binfileY(yfile, std::ios::in | std::ios::binary);
    std::ifstream binfileZ(zfile, std::ios::in | std::ios::binary);
    if(!binfileX.is_open() || !binfileY.is_open() || !binfileZ.is_open()) {
        fprintf(stderr, "Error! Could not open files!\n");
        // In this case set x,y,z to zeros
        t.x = new double[2];
        t.y = new double[2];
        t.z = new double[2];
        for (int i=0;i<2;i++) {
			t.x[i] = 0;
			t.y[i] = 0;
			t.z[i] = 0;
        }
        return t;
    }
    
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
    
    while(!binfileX.eof()) {
        binfileX.read(buf, buff_len);
        if(binfileX.eof()) { //Breaks on last read of file (i.e. when 0 bytes are read and EOF bit is set)
            break;
        }
        x.push_back(*(double *)&buf[0]);
    }
    binfileX.close();
    
    while(!binfileY.eof()) {
        binfileY.read(buf, buff_len);
        if(binfileY.eof()) { //Breaks on last read of file (i.e. when 0 bytes are read and EOF bit is set)
            break;
        }
        y.push_back(*(double *)&buf[0]);
    }
    binfileY.close();
    
    while(!binfileZ.eof()) {
        binfileZ.read(buf, buff_len);
        if(binfileZ.eof()) { //Breaks on last read of file (i.e. when 0 bytes are read and EOF bit is set)
            break;
        }
        z.push_back(*(double *)&buf[0]);
    }
    binfileZ.close();
    
    if(z.size() != x.size() || z.size() != y.size()) {
        fprintf(stderr, "Error! Sample length mismatch!\n");
        return t;
    }
    
    t.x = new double[x.size()];
    t.y = new double[y.size()];
    t.z = new double[z.size()];
    
    for(int i = 0; i < x.size(); i++) {
        t.x[i] = x[i];
        t.y[i] = y[i];
        t.z[i] = z[i];
    }
    
    t.num = x.size();
    
    return t;
}

//======================================================================

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>

#define NRECORDS 50

typedef struct evt {
    double energy;
    double theta;
    float times[NRECORDS];
    float ePerp[NRECORDS];
} evt;

std::vector<evt> readFile(const char* fName) {
    std::vector<evt> events;
    
    const size_t buff_len = 4 + 1*8 + 1*8 + 2*NRECORDS*4 + 4;
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
        event.energy = *((double *)(&buf[0] + sizeof(unsigned int)));
        event.theta = *((double *)(&buf[0] + sizeof(unsigned int) + sizeof(double)));
        std::memcpy((void *)&event.times, (void *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(double)), NRECORDS*sizeof(float));
        std::memcpy((void *)&event.ePerp, (void *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(double) + NRECORDS*sizeof(float)), NRECORDS*sizeof(float));
        events.push_back(event);
    }
    binfile.close();
    
    printf("Read %lu Events!\n", events.size());

    delete[] buf;
    return events;
}

int main(int argc, char** argv) {
    if(argc != 2) {
        printf("Error! Usage: ./minimal fname\n");
        return 1;
    }

    std::vector<evt> data = readFile(argv[1]);

    return 0;
}

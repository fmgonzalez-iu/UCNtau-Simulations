#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <numeric>
#include <cstring>
#include <cmath>
#include <getopt.h>

extern "C" {
    #include "xorshift.h"
}

std::vector<double> refHist = {276.0,2147.0,6357.0,6724.0,6445.0,6515.0,6303.0,5883.0,5937.0,5750.0,5862.0,5741.0,5499.0,5298.0,5118.0,5336.0,5208.0,4989.0,4914.0,4840.0,4773.0,4761.0,4603.0,4670.0,4585.0,4490.0,4504.0,4344.0,4234.0,4190.0,4215.0,4159.0,4099.0,4054.0,4042.0,3958.0,3875.0,3758.0,3770.0,6777.0,19591.0,26111.0,24937.0,23352.0,21985.0,21214.0,20137.0,19300.0,18771.0,17531.0,16790.0,16027.0,15480.0,15079.0,14571.0,13997.0,13510.0,13170.0,12590.0,18787.0,32668.0,30531.0,29039.0,27367.0,25942.0,24249.0,22440.0,21718.0,20811.0,19545.0,18862.0,17545.0,16773.0,15964.0,15590.0,14640.0,14114.0,13309.0,12663.0,19992.0,29095.0,27203.0,25432.0,24246.0,22596.0,21316.0,20074.0,18979.0,17967.0,16720.0,16213.0,15231.0,14874.0,13783.0,13126.0,12620.0,12152.0,11748.0,11099.0,19469.0,31086.0,29762.0,28316.0,26476.0,24893.0,23408.0,22080.0,21091.0,19686.0,18667.0,17434.0,16717.0,15647.0,15110.0,14155.0,13375.0,12639.0,12072.0,11831.0,18263.0,23110.0,22100.0,21373.0,20302.0,19436.0,18349.0,17475.0,16732.0,15687.0,14945.0,14201.0,13216.0,12934.0,12358.0,11613.0,11052.0,10341.0,9803.0,9612.0,15876.0,20286.0,19403.0,17943.0,16981.0,15999.0,14919.0,13821.0,12959.0,11977.0,11107.0,10590.0,9923.0,9119.0,8469.0,8032.0,7199.0,6634.0,6475.0,5908.0,8792.0,10150.0,8870.0,7618.0,6546.0,5978.0,5135.0,4541.0,3926.0,3486.0,3176.0,2779.0,2399.0,2159.0,1892.0,1670.0,1493.0,1387.0,1236.0,1038.0,959.0,787.0,739.0,663.0,589.0};

typedef struct evt {
    double energy;
    double theta;
    double time;
    double eperp;
    double x;
    double y;
    double z;
    double zoff;
    int nhit;
    int nhitBot;
    int nhitTop;
} evt;

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

int main(int argc, char** argv) {
    int c;
    
    char fName[256];
    fName[0] = '\0';
    double begin = -1;
    double end = -1;
    int nbins = -1;

    while (1) {
        static struct option long_options[] = {
            {"file", required_argument, 0, 'f'},
            {"begin", required_argument, 0, 'b'},
            {"end", required_argument, 0, 'e'},
            {"nbins", required_argument, 0, 'n'},
            {0, 0, 0, 0},
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long(argc, argv, "f:b:e:n:", long_options, &option_index);

        /* Detect the end of the options. */
        if(c == -1) {
            break;
        }

        switch(c) {
            case 'f':
                strncpy(fName, optarg, 256-1);
                fName[256-1] = '\0';
                break;
            case 'b':
                begin = atof(optarg);
                break;
            case 'e':
                end = atof(optarg);
                break;
            case 'n':
                nbins = atoi(optarg);
                break;
            case '?':
                /* getopt_long already printed an error message. */
                break;
            default:
                exit(1);
        }
    }
    
    if(fName[0] == '\0' || begin == -1 || end == -1 || nbins == -1) {
        fprintf(stderr, "Error! Usage: ./arrival_time_fixed_eff --file=fName --begin=start_t --end=end_t --nbins=N\n");
        exit(1);
    }
    
//    printf("%f %f %d\n", begin, end, nbins);

    std::vector<evt> data = readFile(fName);
    std::vector<double> mcHist;
    
    mcHist.resize(nbins, 0.0);
    
    for(auto it = data.begin(); it < data.end(); it++) {
        double time = it->time;
//        printf("%f\n", time);
        if(time < begin || time >= end) {
            continue;
        }
        int bin = floor(nbins*(time - begin)/(end-begin));
//        printf("%d\n", bin);
        mcHist[bin] += 1;
    }
    
    for(auto it = mcHist.begin(); it < mcHist.end(); it++) {
        printf("%f,", *it);
    }
    printf("\n");

    return 0;
}

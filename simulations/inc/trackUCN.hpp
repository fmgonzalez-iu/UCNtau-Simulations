#ifndef TRACKUCN_H
#define TRACKUCN_H
#pragma once
#include <vector>

extern "C" {
    #include "../../general/setup.h"
    #include "../../general/typedefs.h"    
}

fixedResult fixedEffDaggerHitTime(std::vector<double> state, double dt, trace tr, double holdT);
fixedResult fixedEffDaggerHitTime_reinsert(std::vector<double> state, double dt, trace tr, double holdT);
fixedResult fixedEffDaggerHitTime_PSE(std::vector<double> state, double dt, trace tr, double holdT);
cleanResult cleanTime(std::vector<double> state, double dt, trace tr, double holdT);
noabsResult noabsCleanTime(std::vector<double> state, double dt, trace tr, double holdT);
noabsResult daggerHitTimes(std::vector<double> state, double dt, trace tr, double holdT);
noabsCleanDagResultZ noAbsCleanTrack_multidata(std::vector<double> state, double dt, trace tr, double holdT, std::ofstream &binfile);

bounceResult timeToBounce(std::vector<double> state, double dt, trace tr, double holdT);
void trackAndPrint(std::vector<double> state, double dt, trace tr, long traj, std::ofstream &binfile);

#endif /* TRACKUCN_H */

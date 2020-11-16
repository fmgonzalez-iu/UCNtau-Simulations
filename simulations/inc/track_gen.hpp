#ifndef TRACK_GEN_H
#define TRACK_GEN_H
#pragma once
#include <vector>
extern "C" {
    #include "../../general/setup.h"
    #include "../../general/typedefs.h"
}

std::vector<double> randomPointTrapOptimum(trace tr);
std::vector<double> randomPointTrapEdE(trace tr);
std::vector<double> randomPointTrapOptimumCleanable(trace tr);
std::vector<double> randomPointTrapEdECleanable(trace tr);
std::vector<double> randomPointTrapOptimumOnlyCleanable(trace tr);
std::vector<double> randomPointTrap(trace tr);

#endif /* TRACK_GEN_H */

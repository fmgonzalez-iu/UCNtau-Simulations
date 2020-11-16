#ifndef SYMPLECTIC_H
#define SYMPLECTIC_H
#pragma once
#include <vector>

extern "C" {
    #include "../../general/setup.h"
    #include "../../general/typedefs.h"
}

void symplecticStep(std::vector<double> &state, double deltaT, double &energy, double t, trace tr);

#endif /* SYMPLECTIC_H */

#ifndef LYAP_H
#define LYAP_H

#include <vector>
#pragma once

extern "C" {
    #include "../../general/setup.h"
    #include "../../general/typedefs.h"
    #include "./fields_nate.h"
}

lyapResult calcLyap(std::vector<double> ref, double dt, trace tr, double tStart = 0);

#endif /* LYAP_H */

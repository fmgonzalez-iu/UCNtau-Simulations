#ifndef STITCHER_H
#define STITCHER_H

#pragma once

#include <vector>
#include <stdio.h>
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"

extern "C" {
	#include "../../general/setup.h"
	#include "../../general/typedefs.h"
}

/*----------------------------------------------------------------------
 * This file converts between various types of data so we don't need to 
 * rewrite a lot of stuff.
 * 
 * -------------------------------------------------------------------*/
std::vector<noabsResult> stitch_noabs_to_noabsClean(std::vector<noabsCleanDagResult> evts_in);
TH1D* stitch_weighted_vec_to_root(std::vector<weightedBin> evts_in);
TH1D* stitch_fixed_result_to_root(std::vector<fixedResult> evts_in);
TH2D* create_root_2D(std::vector<std::vector<double>> state);
TH3D* create_root_3D(std::vector<std::vector<double>> state);


#endif /* STITCHER_H */

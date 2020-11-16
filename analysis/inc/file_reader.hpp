#pragma once
#ifndef FILEREADER_H
#define FILEREADER_H

#include <cstring>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include "TH1D.h"
#include "TF1.h"
#include "TFile.h"

extern "C" {
	#include "../../general/setup.h"
	#include "../../general/typedefs.h"
}
/*----------------------------------------------------------------------
 * Include file for file reading objects
 *--------------------------------------------------------------------*/ 

std::vector<evt> readFile(const char* fName);
std::vector<noabsResult> readFileNoAbsRes(const char* fName);
std::vector<noabsCleanDagResult> readFileNoAbsCleanDagRes(const char* fName);
std::vector<noabsCleanDagResultZ> readFileNoAbsCleanDagResZ(const char* fName);
std::vector<noabsCleanResult> readFileNoAbsCleanRes(const char* fName);
std::vector<fixedResult> readFileFixedRes(const char* fName);
std::vector<cleanResult> readFileCleanRes(const char* fName);
std::vector<fixedResultBlock> readFileBlockRes(const char* fName);
std::vector<fixedResultDag> readFileDagRes(const char* fName);
std::vector<double> refHistFromRoot(const char* fName);
std::vector<std::vector<double>> readFilePosition(const char* fName);


#endif /* FILEREADER_H */

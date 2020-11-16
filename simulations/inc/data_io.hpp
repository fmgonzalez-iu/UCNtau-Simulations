#ifndef DATAIO_H
#define DATAIO_H

#pragma once
#include <vector>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cstring>
#include "../../general/setup.h"
#include "../../general/typedefs.h"

void writePos(std::ofstream &binfile, long traj, double t, std::vector<double> state);
void writeFixedRes(std::ofstream &binfile, fixedResult res);
void writeBounceRes(std::ofstream &binfile, bounceResult res);
void writenoabsCleanRes(std::ofstream &binfile, noabsCleanResult res);
void writeNoabsRes(std::ofstream &binfile, noabsResult res);
void writeCleanRes(std::ofstream &binfile, cleanResult res);
void writeLyapRes(std::ofstream &binfile, lyapResult res);
void writeNoAbsCleanDagRes(std::ofstream &binfile, noabsCleanDagResult res);
void writeNoAbsCleanDagResZ(std::ofstream &binfile, noabsCleanDagResultZ res);
trace readTrace(const char *xfile, const char *yfile, const char *zfile);

#endif /* DATAIO_H */

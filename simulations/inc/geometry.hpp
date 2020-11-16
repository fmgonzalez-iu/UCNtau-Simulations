#ifndef GEOMETRY_H
#define GEOMETRY_H
#pragma once

#include <vector>

std::vector<double> cross(std::vector<double> a, std::vector<double> b);
void normalize(std::vector<double> &a);
double zOffDipCalc(double t);
double calcCleanHeight(double t, double h1, double h2);
void reflect(std::vector<double> &state, std::vector<double> norm, std::vector<double> tang);
bool checkDagHit(double x, double y, double z, double zOff);
int checkClean(std::vector<double> state, std::vector<double> prevState, double cleanHeight);
bool checkBlockHit(std::vector<double> &state, std::vector<double> prevState, double &ePerp);
double calcDagZeta(double x, double y, double z, double zOff);
bool checkMagReflection(std::vector<double> s1, std::vector<double> s2);
bool hitTD(std::vector<double> s1, std::vector<double> s2);
std::vector<double> convertToTrap(std::vector<double> state);
bool checkHouseHitLow(double x, double y, double z, double zOff);
bool checkHouseHitHigh(double x, double y, double z, double zOff);
std::vector<double> initializeLyapState(std::vector<double> ref);
void resetStates(std::vector<double> ref, std::vector<double> &pair);
double distance(std::vector<double> ref, std::vector<double> pair);
//void shiftTD(double *x, double *y, double *z, double t);

#endif /* GEOMETRY_H */

#ifndef UCNTFIELDS_F_H
#define UCNTFIELDS_F_H

#pragma once
#include "../../general/setup.h"
#include "../../general/typedefs.h"

#include <math.h>
#include <stdio.h>

void shift(double *x, double *y, double *z, double t, trace* tr);
void force(double *x_in, double *y_in, double *z_in, double *fx, double *fy, double *fz, double *totalU, double* t, trace* tr);
void fieldstrength(double *x_in, double *y_in, double *z_in, double *totalB, double* t, trace* tr);
void potential(double *x_in, double *y_in, double *z_in, double *totalU, double* t, trace* tr);

#endif /* UCNTFIELDS_F_H */

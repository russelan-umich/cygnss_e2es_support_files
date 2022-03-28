#ifndef GEOM_SPECULAR
#define GEOM_SPECULAR

// This header file was created to allow geomSpecular.c
// to be compiled outside of the end-to-end simulator
// project.

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "sp_math.h"


#define DBL_MAX 9.99e99;
#define pi 3.14159265359

#define fatalError(...) { char str[1000]; sprintf(str, __VA_ARGS__); printf( "%s (%s in %s, line %d)\n", str, __FUNCTION__, __FILE__, __LINE__); exit(EXIT_FAILURE); }


void solveSpecularPt(double rx_pos_ecef[3], double tx_pos_ecef[3], double sx_pos_ecef[3],
                     double rx_vel_ecef[3], double tx_vel_ecef[3], double sx_vel_ecef[3], 
                     double *sx_valid, double r_init_m );

int solveSpecularPtPosition(double rx_pos_ecef[3], double tx_pos_ecef[3], 
                     double sx_pos_ecef[3], double r_init_m, int maxIterations);

int solveSpecularPt_sphericalEarth(double rx_pos_ecef[3], double tx_pos_ecef[3], 
                                    double sx_pos_ecef_unit[3], double r);

int determineLOSSatelliteVis( double rx_pos_ecef[3], double tx_pos_ecef[3] );

double get_earthRadius_m( double x_ecef[3] );

double minPositiveRoot( double polynomCoefs[5] );

#endif // GEOM_SPECULAR

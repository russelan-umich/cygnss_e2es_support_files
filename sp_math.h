#ifndef SP_MATH
#define SP_MATH

// This header file was created to allow sp_math.c
//  to be compiled outside the end-to-end simulator
//	project.
// The math.c file was renamed to sp_math.c so its
//  header file was not named math.h. 

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define R2D 57.295779513082320876798154814105  // radians to degrees

double cot(double z);

double sec(double z);

double csc(double z);

double db2mag( double db );

double mag2db( double mag );

double linear_interp( double a, double b, int direction, double time_01 );

void cubic_interpolation( double f0, double f1, double df0, double df1, 
                         double t, double *ft, double *dft, double timeInterval_s);

void cubic_interpolation_3vector( double f0[3], double f1[3], double df0[3], double df1[3], 
                         double t, double ft[3], double dft[3], double timeInterval_s);

void vector_orthoNorm( double a[3], double b[3] );

void vector_scale( double a[3], double b[3], double s );

double  vector_dot_product (double *a, double *b);

void vector_cross_product(double a[3], double b[3],double c[3]);

double vector_norm(double v[3]);

void vector_unit(double *x, double *x_unit);

void vector_add(double a[3], double b[3],double c[3]);

void vector_subtract(double a[3], double b[3],double c[3]);

void vector_constrainToPlane( double a[3], double b[3], double c[3] );

void matrix_multiply(unsigned rowsA,unsigned columnsA,unsigned columnsB,double *A,double *B,double *C);

void matrix_transpose(unsigned rows,unsigned columns,double *A,double *B);

void matrixVectorMult3x3( double M[9], double x[3], double y[3] );

void matrix_form3x3( double row1[3], double row2[3], double row3[3], double *M );

double matrix_det_3x3(double m[3][3]);

void matrix_scaleAdjoint_3x3(double a[3][3],double s,double m[3][3]);

void matrix_invert_3x3(double A[9], double invA[9]);

double uniformRandf( void );

double randn( void );

void randn_2(double *n1, double *n2, double sigma);
    
void cart2sph( double x[3], double y[3] );

void wgslla2xyz(double x_llh[3], double x_ecef[3]);

void wgsxyz2enu(double p_ecef[3] ,double x_llh[3], double pe_enu[3]);

void wgsxyz2lla( double *x_ecef, double *x_lla );

#endif // SP_MATH


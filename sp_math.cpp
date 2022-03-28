//***************************************************************************
// This file is part of the CYGNSS E2ES.
// 
// CYGNSS E2ES is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// CYGNSS E2ES is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with CYGNSS E2ES.  If not, see <http://www.gnu.org/licenses/>.
//
//---------------------------------------------------------------------------
//
// This file contains various vector, matrix, and math routines used
// throughout the E2ES, including random number generation and coordinate
// system transforms
//
//****************************************************************************/


// This is the only change required to compile
//  this file outside the end-to-end simulator
//  project. The original name of this file
//  was math.c, but was renamed so its new
//  header file was not named math.h.
//#include "gnssr.h"
#include "sp_math.h"


double cot(double z){ return 1.0 / tan(z); }
double sec(double z){ return 1.0 / cos(z); }
double csc(double z){ return 1.0 / sin(z); }

double db2mag( double db ){  return pow(10,(db/20.0)); }
double mag2db( double mag ){  return 20.0 * log10(mag); }

//-----------------------------------------------------------------------------------

double linear_interp( double a, double b, int direction, double time_01 ){
    if( direction == 0 )
        return (1-time_01) * a + time_01 * b;
    else
        return (1-time_01) * b + time_01 * a;
}

void cubic_interpolation( double f0, double f1, double df0, double df1, 
                         double t, double *ft, double *dft, double timeInterval_s){

    // t must be between 0 and 1
    if( (t < 0) || (t > 1) ){
        printf("Time error in cubic interp!\n");
        exit(0);
    }
    
    // special cases when no interpolation needed
    if(t == 0){ *ft  = f0; *dft = df0; return; }
    if(t == 1){ *ft  = f1; *dft = df1; return; }
    
    df0 = df0 * timeInterval_s;
    df1 = df1 * timeInterval_s;
    
    // solve for coef of thrid deg polynom
    const double a = 2*f0 - 2*f1 + df0 + df1;
    const double b = -3*f0 + 3*f1 - 2*df0 - df1;
    const double c = df0;
    const double d = f0;
    
    // values of polynom and derivative at t
    *ft  = a*pow(t,3) + b*pow(t,2) + c*t + d;
    *dft = 3*a*pow(t,2) + 2*b*t + c;
    
    *dft = *dft * ( 1.0 / timeInterval_s );
}


void cubic_interpolation_3vector( double f0[3], double f1[3], double df0[3], double df1[3], 
                         double t, double ft[3], double dft[3], double timeInterval_s){

    // cubic interpolate a 3d vector
    cubic_interpolation(f0[0], f1[0], df0[0], df1[0], t, &(ft[0]), &(dft[0]), timeInterval_s);
    cubic_interpolation(f0[1], f1[1], df0[1], df1[1], t, &(ft[1]), &(dft[1]), timeInterval_s);
    cubic_interpolation(f0[2], f1[2], df0[2], df1[2], t, &(ft[2]), &(dft[2]), timeInterval_s);
    
}


//-----------------------------------------------------------------------------------

void vector_orthoNorm( double a[3], double b[3] ){
    // take vector b, remove component parallel to a
    // and then normalize length to 1
    double s = vector_dot_product(a,b);
    b[0] = b[0] - s*a[0];
    b[1] = b[1] - s*a[1];
    b[2] = b[2] - s*a[2];
    vector_unit(b,b); 
}

void vector_scale( double a[3], double b[3], double s ){
    // b = a * s
    b[0] = a[0]*s;
    b[1] = a[1]*s;
    b[2] = a[2]*s;
}

double  vector_dot_product (double *a, double *b) {
    return (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
}

void vector_cross_product(double a[3], double b[3],double c[3]){
    // c = a x b
    c[0] = a[1]*b[2]-a[2]*b[1];
    c[1] = a[2]*b[0]-a[0]*b[2];
    c[2] = a[0]*b[1]-a[1]*b[0];
}

double vector_norm(double v[3]) {
    return(sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]));
}

void vector_unit(double *x, double *x_unit) {
    // takes x and produces a unit vector in the
    // same direction (x and x_unit can be the same
    // when calling this function)
    double X = vector_norm(x);
    x_unit[0] = x[0] / X;
    x_unit[1] = x[1] / X;
    x_unit[2] = x[2] / X;
}

void vector_add(double a[3], double b[3],double c[3]) {
    c[0] = a[0] + b[0];
    c[1] = a[1] + b[1];
    c[2] = a[2] + b[2];
}

void vector_subtract(double a[3], double b[3],double c[3]) {
    c[0] = a[0] - b[0];
    c[1] = a[1] - b[1];
    c[2] = a[2] - b[2];
}

void vector_constrainToPlane( double a[3], double b[3], double c[3] ){
    double tempa[3], tempb[3];
    vector_unit(a, tempa);
    vector_scale(b, tempb, 1);
    vector_orthoNorm(tempa, tempb);
    double c1 = vector_dot_product(tempa, c);
    double c2 = vector_dot_product(tempb, c);
    
    c[0] = tempa[0] * c1 + tempb[0] * c2;
    c[1] = tempa[1] * c1 + tempb[1] * c2;
    c[2] = tempa[2] * c1 + tempb[2] * c2;
}

//-----------------------------------------------------------------------------------

void matrix_multiply(unsigned rowsA,unsigned columnsA,unsigned columnsB,double *A,double *B,double *C) {
    unsigned i, j, k;

    for(i = 0; i < rowsA; i++) {
        for(j = 0; j < columnsB; j++) {
            C[i*columnsB +j] = 0.0;
            for(k = 0; k < columnsA; k++) {
                C[i*columnsB +j] +=  A[i*columnsA+k]*B[k*columnsB+j];
            }
        }
    }
}

void matrix_transpose(unsigned rows,unsigned columns,double *A,double *B) {
    unsigned i, j;
    for(i = 0; i < columns; i++) {
        for(j = 0; j < rows; j++) {
            B[i*rows+j] = A[j*columns+i];
        }
    }
}

void matrixVectorMult3x3( double M[9], double x[3], double y[3] ){
    // y = M*x
    y[0] = M[0] * x[0] + M[1] * x[1] + M[2] * x[2];
    y[1] = M[3] * x[0] + M[4] * x[1] + M[5] * x[2];
    y[2] = M[6] * x[0] + M[7] * x[1] + M[8] * x[2];
}

void matrix_form3x3( double row1[3], double row2[3], double row3[3], double *M ){
    // form a matrix M by combining three row vectors
    for(int i=0;i<3;i++){
        M[i]   = row1[i];
        M[i+3] = row2[i];
        M[i+6] = row3[i];
    }
}

double matrix_det_3x3(double m[3][3]){
    return   m[0][0] * (m[1][1]*m[2][2] - m[1][2] * m[2][1]) 
           - m[0][1] * (m[1][0]*m[2][2] - m[1][2] * m[2][0]) 
           + m[0][2] * (m[1][0]*m[2][1] - m[1][1] * m[2][0]); 
}

void matrix_scaleAdjoint_3x3(double a[3][3],double s,double m[3][3]){
    a[0][0] = (s) * (m[1][1] * m[2][2] - m[1][2] * m[2][1]);
    a[1][0] = (s) * (m[1][2] * m[2][0] - m[1][0] * m[2][2]);
    a[2][0] = (s) * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    a[0][1] = (s) * (m[0][2] * m[2][1] - m[0][1] * m[2][2]);
    a[1][1] = (s) * (m[0][0] * m[2][2] - m[0][2] * m[2][0]);
    a[2][1] = (s) * (m[0][1] * m[2][0] - m[0][0] * m[2][1]);
    a[0][2] = (s) * (m[0][1] * m[1][2] - m[0][2] * m[1][1]);
    a[1][2] = (s) * (m[0][2] * m[1][0] - m[0][0] * m[1][2]);
    a[2][2] = (s) * (m[0][0] * m[1][1] - m[0][1] * m[1][0]);
}


void matrix_invert_3x3(double A[9], double invA[9]){
    
    double invM[3][3], M[3][3];
    for(int row=0; row<3; row++){
        for(int col=0; col<3; col++){
            M[row][col] = A[col + row*3];
        }
    }
    
    matrix_scaleAdjoint_3x3(invM, (1.0 / matrix_det_3x3(M)), M);
    
    for(int row=0; row<3; row++){
        for(int col=0; col<3; col++){
            invA[col + row*3] = invM[row][col];
        }
    }
    
}



/****************************************************************************/
// uniform random number (0,1]. never = 0 so that we can safely take log of it

double uniformRandf( void ) { return ((double)rand() + 1.0)/((double)RAND_MAX + 1.0); }

/****************************************************************************/
//  Generate a Gausian random variables with zero mean and std dev = 1

double randn( void ){
    double n1, n2;
    randn_2(&n1, &n2, 1);
    return(n1);
}

/****************************************************************************/
//  Generate two Gausian random variables with zero mean and std dev sigma
//  ( standard random number generates create two realizations rather than
//    one )

void randn_2(double *n1, double *n2, double sigma){
    
#if (RAND_GEN_TYPE == 1) // standard Box-Muller Method
    double U1, U2, R, T1, T2;
    U1 = uniformRandf();
    U2 = uniformRandf();
    R  = sigma*sqrt(-2.0*log(U1))*(1/sqrt(2));
    T1 = R*cos(2*pi*U2);
    T2 = R*sin(2*pi*U2);
    *n1 = T1;
    *n2 = T2;
    return();
#endif
    
#if (RAND_GEN_TYPE == 2) // polar form (supposedly faster ...)
    double S,U1, U2, R, T1, T2;
    do {
        U1=2 * uniformRandf() - 1; /* U1=[-1,1] */
        U2=2 * uniformRandf() - 1; /* U2=[-1,1] */
        S=U1 * U1 + U2 * U2;
    } while(S >= 1);
    R  = sigma*sqrt(-2.0*log(S) / S)*(1/sqrt(2));
    T1 = R*U1;
    T2 = R*U2;
    *n1 = T1;
    *n2 = T2;
    return();
#endif
    
}

/****************************************************************************/
// Cartesian to spherical coordinates

void cart2sph( double x[3], double y[3] ){
    // Matlab / mathematicians definition cartesian to spherical coords
    y[0] = sqrt (pow(x[0],2) + pow(x[1],2) + pow(x[2],2)); //R
    y[1] = atan2(x[1], x[0]); // Theta  (azimuth)
    y[2] = atan2(x[2], sqrt(pow(x[0],2) + pow(x[1],2))); //Phi (elevation)
}


/****************************************************************************/
// WGS-84 Lat-Lon-Alt to ECEF X-Y-Z
//**************************************************************
// Replaced by new code from Andrew 2015-03-13
//**************************************************************
//void wgslla2xyz(double x_llh[3], double x_ecef[3]){
//    // LLH to ECEF
    
//    const double f = 1/298.257223563;        //   WGS-84 Flattening.
//    const double e = sqrt(f*(2 - f));        //   Eccentricity.
//    const double R_0 = 6378137;              //   WGS-84 equatorial radius (m).
    
//    //   Compute East-West Radius of curvature at current position
//    const double R_E = R_0 / sqrt(1 - pow(e*sin(x_llh[0]),2));
    
//    //   Compute ECEF coordinates
//    x_ecef[0] = (R_E + x_llh[2])*cos(x_llh[0])*cos(x_llh[1]);
//    x_ecef[1] = (R_E + x_llh[2])*cos(x_llh[0])*sin(x_llh[1]);
//    x_ecef[2] = ((1 - pow(e,2))*R_E + x_llh[2])*sin(x_llh[0]);
//}



/****************************************************************************/
// WGS-84 ECEF X-Y-Z to East-North-Up coordinates at a given location on Earth
void wgsxyz2enu(double p_ecef[3] ,double x_llh[3], double pe_enu[3]){
    // ECEF to ENU
    
    double delta_xyz[3],x_ecef[3];
    
    // Calculate the relative position vector
    wgslla2xyz(x_llh, x_ecef);
    vector_subtract(p_ecef, x_ecef, delta_xyz);
    
    // Calculate ENU coordinates
    pe_enu[2] = cos(x_llh[0])*cos(x_llh[1])*delta_xyz[0] +
    cos(x_llh[0])*sin(x_llh[1])*delta_xyz[1] +
    sin(x_llh[0])*delta_xyz[2];
    pe_enu[0] = -sin(x_llh[1])*delta_xyz[0] + cos(x_llh[1])*delta_xyz[1];
    pe_enu[1] = -sin(x_llh[0])*cos(x_llh[1])*delta_xyz[0] -
    sin(x_llh[0])*sin(x_llh[1])*delta_xyz[1] +
    cos(x_llh[0])*delta_xyz[2];
}


/****************************************************************************/
// WGS-84 ECEF X-Y-Z to Lat-Lon-Alt
//**************************************************************
// Replaced by new code from Andrew 2015-03-13
//**************************************************************
//void wgsxyz2lla( double *x_ecef, double *x_lla ) {
//    // ECEF to LLH
    
//    double f,e,omega_ie,R_0,R_P,mu_E;
//    double x,y,z,lon;
//    double p,E,F,G,c,s,P,Q,k_1,k_2,k_3,k_4,r_0,k_5,U,V;
//    double alt,z_0,e_p,lat;
    
//    // paramters describing the WGS-84 reference ellipsoid
//    f = 1/298.257223563;        //   WGS-84 Flattening.
//    e = sqrt(f*(2 - f));        //   Eccentricity.
//    omega_ie = 7.292115e5;      //   WGS-84 Earth rate (rad/s).
//    R_0 = 6378137;              //   WGS-84 equatorial radius (m).
//    R_P = R_0*(1 - f);          //   Polar radius (m).
//    mu_E = 3.986004418e14;      //   WGS-84 Earth's gravitational
    
//    //  Place holder for output and temporary variables
//    x = x_ecef[0];
//    y = x_ecef[1];
//    z = x_ecef[2];
    
//    //  Calculate longitude
//    lon = atan2(y,x)*R2D;
//    if (lon < 0.0) lon = lon + 360.0;
    
//    p = sqrt(pow(x,2) + pow(y,2));
//    E = sqrt(pow(R_0,2) - pow(R_P,2));
//    F = 54*pow((R_P*z),2);
//    G = pow(p,2) + (1 - pow(e,2))*pow(z,2) - pow((e*E),2);
//    c = pow(e,4)*F*pow(p,2)/pow(G,3);
//    s = pow((1 + c + sqrt(pow(c,2) + 2*c)),(1/3));
//    P = (F/(3*pow(G,2)))/(pow((s + (1/s) + 1),2));
//    Q = sqrt(1 + 2*pow(e,4)*P);
//    k_1 = -P*pow(e,2)*p/(1 + Q);
//    k_2 = 0.5*pow(R_0,2)*(1 + 1/Q);
//    k_3 = -P*(1 - pow(e,2))*pow(z,2)/(Q*(1 + Q));
//    k_4 = -0.5*P*pow(p,2);
//    r_0 = k_1 + sqrt(k_2 + k_3 + k_4);
//    k_5 = (p - pow(e,2)*r_0);
//    U = sqrt(pow(k_5,2) + pow(z,2));
//    V = sqrt(pow(k_5,2) + (1 - pow(e,2))*pow(z,2));
//    alt = U*(1 - (pow(R_P,2)/(R_0*V)));
    
//    //   Compute additional values required for computing
//    //   latitude
//    z_0 = (pow(R_P,2)*z)/(R_0*V);
//    e_p = (R_0/R_P)*e;
//    lat = atan((z + z_0*pow(e_p,2))/p) * R2D;
    
//    x_lla[0] = lat;
//    x_lla[1] = lon;
//    x_lla[2] = alt;
//}


//**************************************************************
// Updated by Andrew 2015-03-13
//**************************************************************
void wgslla2xyz(double x_llh[3], double x_ecef[3]){
    // LLH to ECEF

    const double f = 1/298.257223563;        //   WGS-84 Flattening.
    const double e = sqrt(f*(2 - f));        //   Eccentricity.
    const double R_0 = 6378137;              //   WGS-84 equatorial radius (m).

    //   Compute East-West Radius of curvature at current position
    const double R_E = R_0 / sqrt(1 - pow(e*sin(x_llh[0]),2));

    //   Compute ECEF coordinates
    x_ecef[0] = (R_E + x_llh[2])*cos(x_llh[0])*cos(x_llh[1]);
    x_ecef[1] = (R_E + x_llh[2])*cos(x_llh[0])*sin(x_llh[1]);
    x_ecef[2] = ((1 - pow(e,2))*R_E + x_llh[2])*sin(x_llh[0]);
}

//**************************************************************
// Updated by Andrew 2015-03-13
//**************************************************************
void wgsxyz2lla( double *x_ecef, double *x_lla ) {
   // ECEF to LLH

   double f,e,omega_ie,R_0,R_P,mu_E;
   double x,y,z,lon;
   double p,E,F,G,c,s,P,Q,k_1,k_2,k_3,k_4,r_0,k_5,U,V;
   double alt,z_0,e_p,lat;

   // paramters describing the WGS-84 reference ellipsoid
   f = 1/298.257223563;        //   WGS-84 Flattening.
   e = sqrt(f*(2 - f));        //   Eccentricity.
   omega_ie = 7.292115e5;      //   WGS-84 Earth rate (rad/s).
   R_0 = 6378137;              //   WGS-84 equatorial radius (m).
   R_P = R_0*(1 - f);          //   Polar radius (m).
   mu_E = 3.986004418e14;      //   WGS-84 Earth's gravitational

   //  Place holder for output and temporary variables
   x = x_ecef[0];
   y = x_ecef[1];
   z = x_ecef[2];

   //  Calculate longitude
   lon = atan2(y,x)*R2D;

   p = sqrt(pow(x,2) + pow(y,2));
   E = sqrt(pow(R_0,2) - pow(R_P,2));
   F = 54*pow((R_P*z),2);
   G = pow(p,2) + (1 - pow(e,2))*pow(z,2) - pow((e*E),2);
   c = pow(e,4)*F*pow(p,2)/pow(G,3);
   s = pow((1 + c + sqrt(pow(c,2) + 2*c)),(1/3));
   P = (F/(3*pow(G,2)))/(pow((s + (1/s) + 1),2));
   Q = sqrt(1 + 2*pow(e,4)*P);
   k_1 = -P*pow(e,2)*p/(1 + Q);
   k_2 = 0.5*pow(R_0,2)*(1 + 1/Q);
   k_3 = -P*(1 - pow(e,2))*pow(z,2)/(Q*(1 + Q));
   k_4 = -0.5*P*pow(p,2);
   r_0 = k_1 + sqrt(k_2 + k_3 + k_4);
   k_5 = (p - pow(e,2)*r_0);
   U = sqrt(pow(k_5,2) + pow(z,2));
   V = sqrt(pow(k_5,2) + (1 - pow(e,2))*pow(z,2));
   alt = U*(1 - (pow(R_P,2)/(R_0*V)));

   //   Compute additional values required for computing
   //   latitude
   z_0 = (pow(R_P,2)*z)/(R_0*V);
   e_p = (R_0/R_P)*e;
   lat = atan((z + z_0*pow(e_p,2))/p) * R2D;

   x_lla[0] = lat;
   x_lla[1] = lon;
   x_lla[2] = alt;
}



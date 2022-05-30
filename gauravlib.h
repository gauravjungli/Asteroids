#ifndef GAURAV_LIB
#define GAURAV_LIB

#include <iostream>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <string>
#include <iomanip>
#include <vector>
#include "parameters.h"

using namespace std;

//To write files
void write(const double,const double t );
void write( vector<CV>&,const double t );


//To read files
void read ( vector<double>& );

//Boundary conditions
void bc(vector<CV>& );


//To catch errors
void error(string , string );


//Class for storing a 2D gravity field
class grav{
    public:
        double X, Y, Z;
        grav(): X(0), Y(0), Z(0)
        {}
};

//For the gravity of the sphere
 void grav_sph(vector<grav>&);


//To store conserved variables
class CV{
    public:
        double p, q, r,h,u,v;
        CV(): p(0), q(0), r(0),h(0),u(0),v(0){}
        CV(double h, double u, double v,double x){}
};

//To compute source terms
double S1(int , const CV& , grav& );
double S2(int , const CV& , grav& );
double S3(int , const CV& , grav& );


//To compute characteristic speeds
double ax(int j );


//To incorporate topography
double xi(double , double ); 
double dbdx(double , double );
double dbdy(double , double );


//Calculate minmod limiter
double minmod(double , double , double ); 


//To initialize the simulation
void uniform_IC (vector<CV> & , const CV& , double& );
void grid(vector<double>&);


//To be used in the solver terms
double Hx(int , int , vector<CV>& ); 


//To calculate flux
double flux(CV ,  int  );


//To calculate eigen value
double maxeigenx(double , int );


//To be used for the limiters
double derivative( int , int , vector<CV>& );

//For the boundary conditions
double boundary_condt();


//To check cfl condition
void cfl();


//To reconstruct variables
void reconstruction(CV&  );


//To be used in the time marching
void predictor_corrector(vector<CV>& ,  vector<CV>& , double , vector<double> , double );

//To update omega
void omega(double , double ,double , double , double , double , double );
void Ang_mom(const vector<CV>&, const double&,const vector<double>&,double,double);

#endif
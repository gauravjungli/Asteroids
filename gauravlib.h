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
void write(const vector<CV> &, double  );
void write( const double &, double  );


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
        void modify(double p, double q, double r);
};

//To compute source terms
CV Source(int , const CV& , grav& );


//To incorporate topography
double xi(double , double ); 
double dbdx(double , double );
double dbdy(double , double );


//Calculate minmod limiter
double minmod(double , double , double ); 


//To initialize the simulation
void uniform_IC (vector<CV> & , vector<double> & , double& );
void grid(vector<double>&);


//To be used in the solver terms
CV Hx(CV , CV , double ); 


//To calculate flux
CV flux(CV , double  );


//To calculate eigen value and chareacteristics speed
double eigen(CV, double );
double ax(CV , CV , double );

//To be used for the limiters
CV derivative( CV , CV , CV );

//To check cfl condition
void cfl();

//To be used in the time marching
void march (vector<CV>& , double & , vector<double> & , vector<grav>& ,double);
void predictor(vector<CV>& ,  vector<CV>& , double , vector<double>& , double, vector<grav>& );
void corrector(vector<CV>& ,  vector<CV>& ,vector<CV>&, double , vector<double>& , double, vector<grav>& );

//To update omega
void Ang_mom(const vector<CV>&, const vector<CV>&, double& ,const vector<double>&,double,double);

//mass shedding and pressure
void shed(vector<CV>&,vector<double>&);
double phi(CV&,double&);

//values at edges
void edge(vector<CV>& , vector<CV>& , vector<CV>& );
#endif
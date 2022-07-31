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
void write(const vector<CV> &,const double  );
void write( const double , const double  );

//To read files
void read ( vector<double>&, string file );

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
//spherical gravity
void grav_sph(vector<grav> &);

//To store conserved variables
class CV
{
    public:
        double p, q, r,h,u,v,b,x,psi,lambda;
        grav g;
        CV(double h, double u, double v, double b, grav g, double x, double psi ){}
        void modify(double p, double q, double r);
};

class FS
{
    public:
        double p, q, r;
        FS(): p(0), q(0), r(0){}      
};



//To compute source terms
FS Source(const CV& w);



//To incorporate topography
double xi(double , double ); 
double dbdx(double , double );
double dbdy(double , double );
void Base(vector<double> & base);


//Calculate minmod limiter
double minmod(double a, double b, double c); 


//To initialize the simulation

void uniform_IC (vector<CV> & w, vector<double> & x, vector<grav>& g, double & Om );

void grid(vector<double> & x);


//To be used in the solver terms

FS Hx( CV wl, CV wr);


//To calculate flux
FS flux( CV w );

//To calculate eigen value and chareacteristics speed
double eigen(CV w );
double ax(CV wl, CV wr );

//To be used for the limiters
double derivative( double w1, double w2, double w3);

//To check cfl condition
void cfl(vector<CV>& wl,vector<CV>& wr, double & dt);


//To be used in the time marching
void march (vector<CV>& w, double & Om, double finalt);
void predictor(vector<CV>& w,  vector<CV>& wl, vector<CV>& wr, double Om, double dt);
void corrector(vector<CV>& w,  vector<CV>& wl, vector<CV>& wr, vector<CV>& wtemphat, double Om, double dt);

//To update omega
void Ang_mom(const vector<CV>& w,const vector<CV>& wtemphat, double& Om , double dt,double momincb);

//mass shedding and pressure
void shed(vector<CV>& w);
double Psi();

//values at edges
void edge(vector<CV>& w, vector<CV>& wl, vector<CV>& wr);

void time_step(vector <CV>& wl, vector <CV>& wr, double & dt, double & t, int & timesteps);


#endif
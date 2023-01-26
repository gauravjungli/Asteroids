#ifndef GAURAV_LIB
#define GAURAV_LIB

#include <iostream>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <string>
#include <iomanip>
#include <vector>
#include <sstream>
#include <filesystem>
#include <map>
//#include "parameters.h"
#include <sys/stat.h>
#include <sys/types.h>
using namespace std;

//----------------------------------------------------------------------------
/*#define res (int) round(par["res"])
#define PI 3.14159265
#define dump int(par["dump"])
#define offset par["offset"]
#define xmax   par["xmax"]
#define xmin par["xmin"]
#define weight par["weight"] // used in 2nd  order RK
#define uni_h par["uni_h"]
#define finalt par["finalt"]
#define delta par["delta"] // angle of friction in degrees
#define theta par["theta"] // used in min-mod limiter
#define slides par["slides"]
#define epsilon par["epsilon"] // shallowness parameter
#define omega par["omega"]
#define dx par["dx"]*/
extern std::map <std::string, double> par;

extern  int res;
extern double PI;
extern int dump;
extern double offset;
extern double xmax;
extern double xmin;
extern double weight;
extern double uni_h;
extern double finalt;
extern double delta;
extern double theta;
extern double slides;
extern double epsilon;
extern double omega;
extern double dx;




//------------------------------------------------------------------------------

//Class for storing a 2D gravity field

class Grav{
    public:
        double X1, X2, X3;
        Grav(): X1(0), X2(0), X3(0)
        {}
        
};

class AMB{
    public:
        std::vector<double> ang_mom;
        std::vector<double> inertia;
        AMB(double mom, double inr)
        {
            ang_mom=vector<double>(5,mom);
            inertia=vector<double>(5,inr);
        }
        
};

//To store conserved variables
class CV
{
    public:
        double p, q, r,h,u,v,b,grad_b,x,psi,lambda;
        Grav g;
        CV(double h, double u, double v, double b, double grad_b, Grav g, double x, double om );
        void Modify(double p, double q, double r, double om);
        
        double H(double om);
        double U(double om);
        double V(double om);
};

class FS
{
    public:
        double p, q, r;
        FS(): p(0), q(0), r(0){} 
        FS(CV w);
        FS operator+ (FS w); 
        FS operator- (FS w);  
        FS operator* (double w);   
        FS operator/ (double w);    
};


//---------------------------------------------------------------------------------

/////   I/O

//To write files 
void Write(const vector<CV> & w, string file  );
void Write( const double om, string file );
void Write( const double om, const double t, string file );

//To read files
void Read ( vector<double>&, string file );
bool Parameters();

//To catch errors
void Error(string , string );

//to delete files
void deleteDirectoryContents(const std::string& dir_path);

//----------------------------------------------------------------------------------------
///// bc.cpp
//Boundary conditions
void BC(vector<CV>& w, double om );

//---------------------------------------------------------------------------------

//////// gavity.cpp
//spherical gravity
void Grav_sph(vector<Grav> & g);

//--------------------------------------------------------------------------------

///////// TVD.cpp

//To be used for the limiters
double Derivative( double w1, double w2, double w3);

//Calculate minmod limiter
double Minmod(double a, double b, double c); 
FS Minmod(FS w, FS v); 
//--------------------------------------------------------------------------
//////////  IC.cpp

//To initialize the simulation
void Uniform_IC (vector<CV> & w, vector<double> & x, vector<Grav>& g, double om );
void Grid(vector<double> & x);
//To incorporate topography
void Base(vector<CV>& w,vector<double> & b);

//-------------------------------------------------------------------------

////////////// solver.cpp

//To be used in the solver terms

FS Hx( CV wl, CV wr, double om);


//------------------------------------------------------------------------------------------

////////// cv_flux-source.cpp

//To calculate flux
FS Flux( CV w, double om );

//To compute source terms
FS Source( CV w, double om, double alpha);
FS Eigen(CV w );
//---------------------------------------------------------------------------------

///////// characteristics.cpp

//To calculate eigen value and chareacteristics speed

double Ax(CV wl, CV wr, string s);
//values at edges
void Edge(vector<CV>& w, vector<CV>& wl, vector<CV>& wr, double om);
void Reconstruct(CV& wl, CV w1, CV w2, CV w3, int sign );


//--------------------------------------------------------------------------------------------------------

/////// march.cpp

//To be used in the time marching
void March (vector<CV>& w, double & om, double final_t);
void Predictor(vector<CV>& w,  vector<CV>& wl, vector<CV>& wr, double& om, double dt, AMB amb);
void Corrector(vector<CV>& w,  vector<CV>& wl, vector<CV>& wr, vector<CV>& w_init, double& om, double om_init, double dt, AMB amb);
void Time_step(vector <CV>& wl, vector <CV>& wr, double & dt, double & t, int & timesteps);
void CFL(vector<CV>& wl,vector<CV>& wr, double & dt);

//------------------------------------------------------------------------------------------------------------

///////// Angular_mom.cpp
//To update omega
void Ang_mom(AMB& amb, double integral1, double integral2);
double Alpha( double& om , double dt, AMB amb );
void Integrals(double& integral1, double& integral2, vector<CV>& w, vector<CV>& wl, vector<CV>& wr);
double Int_B(CV w);
double Int_A(CV w);
double Diff(vector<double> w, double dt);


//--------------------------------------------------------------------------------------------------------------

///////// pressure_shed.cpp

//mass shedding and pressure
void Shed(vector<CV>& w, double om);
double Psi(CV w, double om);





#endif
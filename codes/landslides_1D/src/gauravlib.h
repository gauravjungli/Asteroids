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
#include <regex>
#include <sys/stat.h>
#include <sys/types.h>

using namespace std;
namespace fs = std::filesystem;

extern std::map <std::string, string> par;

extern const  int res;
extern const double PI;
extern const int dump;
extern const double offset;
extern const double xmax;
extern const double xmin;
extern const double weight;
extern const double finalt;
extern const double Delta;
extern const double theta;
extern const double slides;
extern const double epsilon;
extern const double omega;
extern const double dx;
extern const double past_time;
extern const double dia;
extern const double min_h;
extern const double min_u;
extern double delta;
extern const string par_add;
extern const string fric_type;
extern const string Output_folder;
extern const string verbose_dir;
extern const string verbose;
extern const double seismic_time;
extern double mass_shed;
extern double k_d;
//------------------------------------------------------------------------------

//Class for storing a 2D gravity field

class Grav{
    public:
        double X1, X2, X3;
        Grav(): X1(-1), X2(0), X3(0)
        {}
        Grav operator+ (Grav w);
        Grav operator/ (double w);
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
        double w, p, q, r,h,u,u_c,v,V,b,db,ddb,x,psi,metric,theta,phi_norm,phi_tan,phi,J;
        Grav g;
        CV(double h, double u,double u_c, double v, double b, double db, double ddb, Grav g, double x);
      //  CV(const CV& temp);
      //  CV& operator=(const CV& temp);
        void Modify(double p, double q, double r);
     
};

 template <typename T>
int sign(T val) {
    return (T(0) < val) - (val < T(0));
}

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
void Write_data(const vector<CV> & w, string file  );
void Write_base ( const vector<CV>& w, string file);
void Write( const double om, string file );
void Write( const double om, const double t, string file );
void Write_par ( std::map <std::string, string> par);

//To read files
void Read ( vector<double>&, string file );
bool Parameters();
void Read_grav( vector<Grav>& g, const string& file);
void Read_data ( vector<double>& x,vector<double>& b,vector<double>& db,vector<double>& ddb,vector<double>& h,string file);
//To catch errors
void Error(string , string );

//to delete files
void deleteDirectoryContents(const std::string& dir_path);
std::string to_string(double value, int precision);

//----------------------------------------------------------------------------------------
///// bc.cpp
//Boundary conditions
void BC(vector<CV>& w1 ,vector<CV>& w2);

//---------------------------------------------------------------------------------

//////// gavity.cpp
void Init_grav(vector<Grav> & g, string file);

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
void Initial_Condition (vector<CV> & w, vector<CV> & wl, vector<CV> & wr,vector<Grav>& g);
void Grid(vector<double> & x);

//-------------------------------------------------------------------------

////////////// solver.cpp

//To be used in the solver terms

FS Hx( CV wl, CV wr);


//------------------------------------------------------------------------------------------

////////// cv_flux-source.cpp

//To calculate flux
FS Flux( CV w );

//To compute source terms
FS Source( CV w, CV w1, CV w2, CV w3, CV w4);
FS Eigen(CV w );

FS Friction (CV w);

FS Body_force (CV w, CV w1, CV w2, CV w3, CV w4);


double Ang_mom_reg (CV w);

double Jinertia1_reg (CV w);

double Jinertia1_ast (CV w);

double Jinertia2_reg (CV w);

double Jinertia2_ast (CV w);

//---------------------------------------------------------------------------------

///////// characteristics.cpp

//To calculate eigen value and chareacteristics speed

double Ax(CV wl, CV wr, string s);
//values at edges
void Edge(vector<CV>& w, vector<CV>& wl, vector<CV>& wr);
void Reconstruct(CV& wl, CV& wr, CV w1, CV w2, CV w3 );
void Balancing (CV& w, CV& wl, CV& wr);
//--------------------------------------------------------------------------------------------------------

/////// march.cpp

//To be used in the time marching
void March (vector<CV>& w, vector<CV>& wl, vector<CV>& wr, double& Ang_Shed);
void Predictor(vector<CV>& w,  vector<CV>& wl, vector<CV>& wr, double dt);
void Corrector(vector<CV>& w,  vector<CV>& wl, vector<CV>& wr, vector<CV>& w_init, double dt);
void Time_step(vector <CV>& wl, vector <CV>& wr, double & dt, double & t, int & timesteps);
void CFL(vector<CV>& wl,vector<CV>& wr, double & dt);

//------------------------------------------------------------------------------------------------------------

///////// Omega.cpp
//To update omega
double Inertia(vector<CV>& w, int no);
double Reg_Inertia(vector<CV>& w);
//--------------------------------------------------------------------------------------------------------------

///////// pressure_shed.cpp

//mass shedding and pressure
void Shed(vector<CV>& w, double& Ang_Shed);
double Psi(CV w);
double Psi_basal(CV w);
double J_Psi(CV w);


////////
//double sin(double x);


#endif
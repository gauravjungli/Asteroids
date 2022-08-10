
#ifndef PARAMETERS
#define PARAMETERS

const int res=200;
const double PI=3.14159265; 
const int dump=50; //dumping interval, no. of time steps
const double offset= 0.05;
const double xmax=PI;
const double xmin=0;
const double weight = 0.5; // used in 2nd  order RK
const double uni_h=1.0;
const double finalt = 1;

// parameters
const double delta = 20 * PI / 180; // angle of friction in degrees
const double theta = 1.0; // used in min-mod limiter



const double epsilon = 0.01; // shallowness parameter

///donot change beyond this

const double dx = (xmax-xmin-2*offset)/res;

const double omega=1.5;
#endif


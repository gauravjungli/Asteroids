const int res=200;
const double PI=3.14159265; 
const int dump=50; //dumping interval, no. of time steps
const int xmax=2*PI;
const int xmin=0;
const double weight = 0.5; // used in 2nd  order RK
const double uni_h=1.0;
const double finalt = 10;

// parameters
const double delta = 20 * PI / 180; // angle of friction in degrees
const double theta = 1.0; // used in min-mod limiter


const double epsilon = 0.005; // shallowness parameter

///donot change beyond this

const double dx = (xmax-xmin) / res;

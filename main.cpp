
#include "gauravlib.h"

int main()
{	
	double Om=0;
    vector<grav> g(res);
	grav_sph(g);
	vector<double> x(res);
	grid(x);
	vector<CV> w;
	uniform_IC(w,x,g,Om);
	march(w,Om,finalt);
	return 0;
	
}
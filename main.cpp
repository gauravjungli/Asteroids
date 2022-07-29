
#include "gauravlib.h"

int main()
{	
	double Om=0;
    vector<grav> g(res);
	grav_sph(g);
	vector<double> x;
	grid(x);
	vector<CV> w;

	uniform_IC(w,x,Om);
	march(w,Om,x,g,finalt);
	
	return 0;
}

#include "gauravlib.h"

int main()
{	
	double dt = dx / 4;
	double t = 0;
	double Om=0;
    vector<grav> g(res);
	grav_sph(g);
	vector<CV> w;
	vector<CV> wl;
	vector<CV> wr;
	CV ic(uni_h,0,0);
	vector<double> x;
	grid(x);
	uniform_IC(w,ic,Om);
	write(w,t);
	write(Om,t);
	vector<CV> wtemp(w);
	bc(wtemp);
	return 0;
}
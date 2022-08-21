
#include "gauravlib.h"

int main()
{	

	double om=omega;
    vector<Grav> g(res);
	Grav_sph(g);
	vector<double> x(res);
	Grid(x);
	vector<CV> w;
	Uniform_IC(w,x,g,om);
	March(w,om,finalt);
	return 0;
	
}
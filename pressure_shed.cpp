#include "gauravlib.h"
void shed(vector<CV>& w, double Om)
{
for ( int j = 2; j <= res; j++)
{

	if (Psi(w[j],Om) > 0) // -1e-6
	{
		//cout << "mass shed";
		w[j].modify(  pow(10, -8),
			w[j].p * w[j].u,
		 w[j].p * w[j].v,Om);	
	}
}
}
double Psi(CV w, double Om)
{
return	Om*Om*sin(w.x)*sin(w.x)+2*Om*w.v*sin(w.x)+w.g.X+w.u*w.u+w.v*w.v;//change for the sphere
}
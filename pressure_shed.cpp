#include "gauravlib.h"
void Shed(vector<CV>& w)
{
for ( int j = 2; j < res-2; j++)
{

 	if (Psi(w[j]) < 0) 
	{
		w[j].Modify(  pow(10, -8),
			pow(10, -8) * w[j].u,
		 pow(10, -8) * w[j].v);	
	} 
}
}
double Psi(CV w)
{
return	-(omega*omega*sin(w.x)*sin(w.x)+2*omega*w.v*sin(w.x)+w.g.X1+w.u*w.u+w.v*w.v);//change for the sphere
}
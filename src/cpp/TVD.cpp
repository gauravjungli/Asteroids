#include "gauravlib.h"

double Minmod(double a, double b, double c) // calculate minmod
{

	if (a < 0 && b < 0 && c < 0)
		return std::max({a, b, c});
	else if (a > 0 && b > 0 && c > 0)
		return std::min({a,b,c});
	else
		return 0;
}

FS Minmod(FS w, FS v) // calculate minmod
{
	FS temp;
	temp.p=Minmod(w.p,v.p,v.p);
	temp.q=Minmod(w.q,v.q,v.q);
	temp.r=Minmod(w.r,v.r,v.r);
	return temp;
}


double Derivative( double w1, double w2, double w3)
{
	double w;	
	w= Minmod(theta * (w2 - w1) / dx, (w3 - w1) / 2 / dx, theta * (w3 - w2) / dx);
	return w;		 
}

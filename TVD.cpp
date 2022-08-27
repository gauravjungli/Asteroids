#include "gauravlib.h"

double Minmod(double a, double b, double c) // calculate minmod
{

	if (a < 0 && b < 0 && c < 0)
		return std::max(a, std::max(b, c));
	else if (a > 0 && b > 0 && c > 0)
		return std::min(a, std::min(b, c));
	else
		return 0;
}


double Derivative( double w1, double w2, double w3)
{
	double w;	
	w= Minmod(theta * (w2 - w1) / dx, (w3 - w1) / 2 / dx, theta * (w3 - w2) / dx);
	return w;		 
}

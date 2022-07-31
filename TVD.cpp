#include "gauravlib.h"

double minmod(double a, double b, double c) // calculate minmod
{
	double temps;
	//temps = std::min(abs(a), std::min(abs(b), abs(c)));

	//return temps;

	if (a < 0 && b < 0 && c < 0)
		return std::max(a, std::max(b, c));
	else if (a > 0 && b > 0 && c > 0)
		return std::min(a, std::min(b, c));
	else
		return 0;
}


double derivative( double w1, double w2, double w3)
{
	double w;	
	w= minmod(theta * (w2 - w1) / dx, (w3 - w1) / 2 / dx, theta * (w3 - w2) / dx);
	return w;		 
}

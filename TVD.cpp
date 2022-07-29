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

CV derivative( CV w1, CV w2, CV w3)
{
	CV w;	
	w.modify( minmod(theta * (w2.p - w1.p) / dx, (w3.p - w1.p) / 2 / dx, theta * (w3.p - w2.p) / dx),
			 minmod(theta * (w2.q - w1.q) / dx, (w3.q - w1.q) / 2 / dx, theta * (w3.q - w2.q) / dx),
			 minmod(theta * (w2.r - w1.r) / dx, (w3.r - w1.r) / 2 / dx, theta * (w3.r - w2.r) / dx));

	return w;		 
}

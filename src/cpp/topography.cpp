#include "gauravlib.h"

double xi(double xx, double yy) // topography parameter
{
	return 0 * PI / 180.0;
}

double dbdx(double xx, double yy)
{
	return -tan(xi(xx, yy)); // b = -tan(xi(x,y))*xx + 1
}
double dbdy(double xx, double yy)
{
	return 0;
}
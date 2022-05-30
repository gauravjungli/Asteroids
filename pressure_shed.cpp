#include "gauravlib.h"

void shed(vector<CV>& w,vector<double> x)
{
for ( int j = 2; j <= res; j++)
			{
				w[j].h = w[j].p;
				w[j].u = w[j].q / w[j].p;
				w[j].v = w[j].r / w[j].q;
				if (phi > 0) // -1e-6
				{
					//cout << "mass shed";
					w[j].p = pow(10, -8);
					w[j].q = w[j].p * w[j].u;
					w[j].r = w[j].p * w[j].v;		
				}
			}
}
double phi()
{
return	0;//change for the sphere
}
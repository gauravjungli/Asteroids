#include "gauravlib.h"

void grid(vector<double> & x)
{
    for (int i=0;i<res;i++)
        {
        double z=(xmax-xmin) * i / res + (xmax-xmin) / 2 / res;
        x.push_back(z);
        }
}

void uniform_IC (vector<CV> & w, const CV& v, double& Om )
{   
    Om = 1.5;
	for (int j = 0; j < res; j++)
	    w.push_back(v);
}
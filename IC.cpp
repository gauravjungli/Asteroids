#include "gauravlib.h"

void grid(vector<double> & x)
{
    for (int i=0;i<res;i++)
        {
        double z=(xmax-xmin) * i / res + (xmax-xmin) / 2 / res;
        x.push_back(z);
        }
}

void uniform_IC (vector<CV> & w, vector<double> & x, double & Om )
{   
    Om = omega;
    
	for (int j = 0; j < res; j++)
    {   CV temp(uni_h,0,0,x[j]);
        w.push_back(temp);
    }
	    
}
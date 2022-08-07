#include "gauravlib.h"

void grid(vector<double> & x)
{
    for (int i=0;i<res;i++)
        {
        double z=(xmax-xmin) * i / res + (xmax-xmin) / 2 / res;
        x[i]=z;
        }
}


void uniform_IC (vector<CV> & w, vector<double> & x, vector<grav>& g, double & Om )
{   
    Om = omega;
    vector<double>  b(res);
    
    Base(b);
	for (int j = 0; j < res; j++)
    {
        CV temp(uni_h,0,0,b[j],g[j],x[j],Om);
        w.push_back(temp);
    }
	    
}

void Base(vector<double> & b)
{
    b=vector<double>(res,0);
}
#include "gauravlib.h"

void Grid(vector<double> & x)
{
    for (int i=0;i<res;i++)
        {
        x[i]=offset+dx * (i+0.5);
        }
}


void Uniform_IC (vector<CV> & w, vector<double> & x, vector<Grav>& g, double om )
{   
    vector<double>  b(res);
    
    Base(w,b);
    if (w.empty())
    {
	for (int j = 0; j < res; j++)
    {
        CV temp(uni_h,0,0,b[j],g[j],x[j],om);
        w.push_back(temp);
    }
    }
    else
    {
        for (int j = 0; j < res; j++)
            w[j]=CV(uni_h,0,0,b[j],g[j],x[j],om);
    }
	    
}

void Base(vector<CV>& w, vector<double> & b)
{   
    if (w.empty())
        b=vector<double>(res,5);
    else
    {
        for (int j = 0; j < res; j++) 
            b[j]=w[j].b+w[j].h-uni_h;
    }
        
}
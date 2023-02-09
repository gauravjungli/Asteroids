#include "gauravlib.h"

std::map <std::string, double> par;
bool set_parameter=Parameters();
const int res= (int) round(par["res"]);
const double PI= 3.14159265;
const int dump= int(par["dump"]);
const double offset= par["offset"];
const double xmax=   par["xmax"];
const double xmin= par["xmin"];
const double weight= par["weight"]; 
const double uni_h= par["uni_h"];
const double finalt= par["finalt"];
const double delta= par["delta"]; 
const double theta= par["theta"]; 
const double slides= par["slides"];
const double epsilon= par["epsilon"]; 
const double omega= par["omega"];
const double dx= (xmax-xmin-2*offset)/res;

void Grid(vector<double> & x)
{
    for (int i=0;i<res;i++)
        {
            x[i]=offset+dx * (i+0.5);
        }
}


void Uniform_IC (vector<CV> & w, vector<double> & x, vector<Grav>& g)
{   
    vector<double>  b;
     vector <double> h(res), u(res),v(res);

    for (int j=0;j<res;j++)
    {
        h[j]=uni_h-pow((0.75*j/res),2); u[j]=0;v[j]=0;
        
    }
    Base(w,b,h);
    if (w.empty())
    {
	    for (int j = 0; j < res; j++)
        {
            CV temp(h[j],u[j],v[j],b[j],g[j],x[j]);
            w.push_back(temp);
        } 
    }
    else
    {
        for (int j = 1; j < res-1; j++)
       { 
            w[j]=CV(h[j],u[j],v[j],b[j],g[j],x[j]);
       }
    } 
	    
}

void Base(vector<CV>& w, vector<double> & b, vector<double>& h)
{   
    if (w.empty())
    {
        b=vector<double>(res,0);
        for (int j=0;j<res;j++)
            b[j]=uni_h-h[j];
    }
    else
    {
        for (int j = 1; j < res-1; j++) 
         {   
            b[j]=w[j].w-h[j];
         }
    }
        
}
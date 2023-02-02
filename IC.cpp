#include "gauravlib.h"

std::map <std::string, double> par;
bool set_parameter=Parameters();
int res= (int) round(par["res"]);
double PI= 3.14159265;
int dump= int(par["dump"]);
double offset= par["offset"];
double xmax=   par["xmax"];
double xmin= par["xmin"];
double weight= par["weight"]; 
double uni_h= par["uni_h"];
double finalt= par["finalt"];
double delta= par["delta"]; 
double theta= par["theta"]; 
double slides= par["slides"];
double epsilon= par["epsilon"]; 
double omega= par["omega"];
double dx= (xmax-xmin-2*offset)/res;

void Grid(vector<double> & x)
{
    for (int i=0;i<res;i++)
        {
            x[i]=offset+dx * (i+0.5);
        }
}


void Uniform_IC (vector<CV> & w, vector<double> & x, vector<Grav>& g, double om )
{   
    vector<double>  b;
     vector <double> h(res), u(res);

    for (int j=0;j<res;j++)
    {
        if (j<res/2)
        {h[j]=uni_h; u[j]=0;}
        else
        {h[j]=uni_h; u[j]=0;}
    }
    Base(w,b,h);
    if (w.empty())
    {
	    for (int j = 0; j < res; j++)
        {
            CV temp(h[j],u[j],0,b[j],g[j],x[j],om);
            w.push_back(temp);
        } 
    }
    else
    {
        for (int j = 1; j < res-1; j++)
       { 
            w[j]=CV(h[j],u[j],0,b[j],g[j],x[j],om);
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
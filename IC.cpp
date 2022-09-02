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
       { 
            w[j]=CV(uni_h,0,0,b[j],g[j],x[j],om);
       }
    }
	    
}

void Base(vector<CV>& w, vector<double> & b)
{   
    if (w.empty())
        b=vector<double>(res,0);
    else
    {
        for (int j = 0; j < res; j++) 
            b[j]=w[j].b+w[j].h-uni_h;
    }
        
}
#include "gauravlib.h"

std::map <std::string, double> par;
bool set_parameter=Parameters();
const int res= (int) round(par["res"]);
const double PI= M_PI;
const int dump= int(par["dump"]);
const double offset= par["offset"];
const double xmax=   PI;
const double xmin=   0;
const double weight= par["weight"]; 
const double uni_h= par["uni_h"];
const double finalt= par["finalt"];
const double delta= par["delta"]; 
const double theta= par["theta"]; 
const double slides= par["slides"];
const double epsilon= par["epsilon"]; 
const double omega= par["omega"];
const double omega_initial=par["omega_initial"];
const double dx= (xmax-xmin-2*offset)/res;
const int layers= par["layers"];
const double radius=par["radius"];
const double past_time=par["time"];
const double dia=par["dia"];
const double min_h=pow(dx,4);

void Grid(vector<double> & x)
{
    for (int i=0;i<res;i++)
        {
            x[i]=offset+dx * (i+0.5);
        }
}


void Uniform_IC (vector<CV> & w, vector<double> & x, vector<Grav>& g, string file)
{   
    vector<double>  b(res,0);
    vector <double> h(res,uni_h), u(res,0),v(res,0);
    
//Uncomment this one only if you want special initial conditions
  /*  for (int j=0;j<res;j++)
    {
        
        b[j]=sin(x[j]);
        h[j]=uni_h-b[j];
        
    } */
    

    if (slides>1)
        Base(b,file);
    
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

void Base ( vector<double>& b,string file)
{
ifstream myfile(file);

string line;
int i=2;
while(getline(myfile,line))
{
	istringstream iss(line);
    string word1,temp,word2;
	getline(iss, temp, ',');  
    getline(iss, word1, ',');
    getline(iss, word2, ',');       

	b[i]=(1+epsilon*(stod(word1)+stod(word2))-radius)/epsilon;
	i++;
}
myfile.close();

}
#include "gauravlib.h"

std::map <std::string, string> par;
bool set_parameter=Parameters();
const int res= (int) round(stod(par["res"]));
const double PI= M_PI;
const int dump= int(stod(par["dump"]));
const double offset= stod(par["offset"]);
const double xmax=   PI;
const double xmin=   0;
const double weight= stod(par["weight"]); 
const double uni_h= stod(par["uni_h"]);
const double finalt= stod(par["finalt"]);
const double delta= stod(par["delta"]); 
const double theta= stod(par["theta"]); 
const double slides= stod(par["slides"]);
const double epsilon= stod(par["epsilon"]); 
const double omega= stod(par["omega"]);
const double omega_initial=stod(par["omega_in"]);
const double dx= (xmax-xmin-2*offset)/res;
const double past_time=stod(par["time"]);
const double dia=stod(par["dia"]);
const double min_h=pow(dx,4);
const double gamma=stod(par["gamma"]);

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
    
    Base(b,x,file);

    for (int j=0;j<res;j++)
    {
        
        // if (b[j]<-1)
        // {
        //     h[j]=min_h;
        //     b[j]=b[j]+epsilon/gamma*uni_h;
        // }

        // if (b[j]>1)
        // {
        //     h[j]=gamma/epsilon*(b[j]-1);
        //     b[j]=1;
        // }
        
    } 
    
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

void Base ( vector<double>& b,vector<double>& x, string file)
{
ifstream myfile(file+"/base.txt");

string line;
int i=0;
while(getline(myfile,line))
{
	istringstream iss(line);
    string word1,word2;
    getline(iss, word1, ',');
    getline(iss, word2, ',');       

    x[i]=stod(word1);
	b[i]=stod(word2)-epsilon/gamma*uni_h;
	i++;
}
myfile.close();

}
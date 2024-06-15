#include "gauravlib.h"


void Grid(vector<double> & x)
{
    for (int i=0;i<res;i++)
        {
            x[i]=offset+dx * (i+0.5);
        }
}


void Uniform_IC (vector<CV> & w, vector<double> & x, vector<Grav>& g)
{   
    vector<double>  b(res,0);
    vector <double> h(res,1), u(res,0),v(res,0);
    
//Uncomment this one only if you want special initial conditions

    Base(b,h,x);

    // for (int j=0;j<res;j++)
    // {
    //     double max_h=1;
    //     if (b[j]<-max_h)
    //     {
    //         h[j]=min_h;
    //         b[j]=b[j]+epsilon/Gamma*uni_h;
    //     }

    //     if (b[j]>max_h)
    //     { 
    //         h[j]=h[j]+Gamma/epsilon*(b[j]-max_h);
    //         b[j]=max_h;
    //     }
        
    // } 
    
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

/* void Base(vector<CV>& w, vector<double> & b, vector<double>& h)
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
        
} */

void Base ( vector<double>& b,vector<double>& h,vector<double>& x)
{
    std::string  file1=folder +"/base.txt";
    ifstream myfile(file1);

    std::string line;
    int i=0;
    while(getline(myfile,line))
    {
        istringstream iss(line);
        string word1,word2,word3;
        getline(iss, word1, ',');
        getline(iss, word2, ','); 
        getline(iss, word3, ',');        

        x[i]=stod(word1);
        b[i]=stod(word2);
        h[i]=stod(word3);
        i++;
    }
    myfile.close();

}
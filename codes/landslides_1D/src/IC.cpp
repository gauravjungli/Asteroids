#include "gauravlib.h"


void Grid(vector<double> & x)
{
    for (int i=0;i<res;i++)
        {
            x[i]=offset+dx * (i+0.5);
        }
}


void Initial_Condition (vector<CV> & w, vector<CV> & wl, vector<CV> & wr,vector<Grav>& g)
{   
    vector<double>  b(res,1),x(res,0),db(res,0),ddb(res,0);
    vector <double> h(res,1), u(res,min_u),v(res,min_u);

    Read_data(x,b,db,ddb,h,"base.txt");

    
    if (w.empty())
    {
	    for (int j = 0; j < res; j++)
        {   
            double dbase = 0;
            double ddbase = 0;
            double base =1;

           // h[j] = 1; // + exp(-pow((x[j]-x[500]),2)/(2*0.5*0.5)); 

            CV temp(h[j],u[j],v[j],base,dbase ,ddbase, g[j],x[j]);

            if (j>0 and j<res-1)
            {
             dbase = db[j];
             ddbase = ddb[j];
             base = b[j];
             temp = CV(h[j],u[j],v[j],base,dbase,ddbase,g[j],x[j]);
            }
            
            w.push_back(temp);

               if (j>0)
            {
            base = (b[j-1]+b[j])/2;
            dbase = (db[j-1]+db[j])/2;
            ddbase = (ddb[j-1]+ddb[j])/2;          
            temp = CV(h[j],u[j],v[j],base,dbase,ddbase,(g[j-1]+g[j])/2,(x[j-1]+x[j])/2);
            }

            wl.push_back(temp);


             if ( j<res-1)
            {
            base = (b[j+1]+b[j])/2;
            dbase = (db[j+1]+db[j])/2;
            ddbase = (ddb[j+1]+ddb[j])/2;
            temp = CV(h[j],u[j],v[j],base, dbase, ddbase, (g[j+1]+g[j])/2,(x[j+1]+x[j])/2);
            }  
            
            wr.push_back(temp); 
        }  
    }
    else
    {
       cout<<"There already appears to be elements in the w vector"<<endl;
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


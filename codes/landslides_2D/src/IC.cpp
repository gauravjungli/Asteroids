#include "gauravlib.h"


void Grid(vector<double> & x, vector <double> & y)
{
    for (int i=0;i<rows;i++)
    for (int j=0;j<cols;j++)
    {
        {   int i1 = index(i,j);
            x[i1]= offset+dx * (i+0.5);
            y[i1]= dy * (j+0.5);
        }
    }
}


void Initial_Condition (vector<CV> & w, vector<CV> & wl, vector<CV> & wr,vector<CV> & wb, vector<CV> & wt, vector<Grav>& g)
{   
    vector<double>  b(rows*cols,1),x(rows*cols,0),y(rows*cols,0),Ra(rows*cols,0),dR(rows*cols,0),ddR(rows*cols,0) ;
    vector <double> h(rows*cols,1), u(rows*cols,0),v(rows*cols,0);
    

    Read_data(x,y,b,h,Ra,dR,ddR,"base.txt");
    
     if (w.empty())
    {
	    for(int i=0;i<rows;i++)
    {
	    for (int j = 0; j < cols; j++)
        {
            int i0 = index(i,j);
            int i_1 = index(i-1,j);
		    int i1 = index(i+1,j);
            int j0 = index(i,j);
		    int j_1 = index(i,j-1);
		    int j1 = index(i,j+1);
            
            double dRadius = 0;
            double ddRadius = 0;
            double base = 0;
            double Radius =1;
      

            CV temp(h[i0],u[i0],v[i0],base, Radius, dRadius ,ddRadius, g[i0],x[i0],y[i0]);

            if (i>=0 and i<=rows-1)
            {
             dRadius = dR[i0];
             ddRadius = ddR[i0];
             base = b[i0];
             Radius = Ra[i0];
             temp = CV(h[i0],u[i0],v[i0],base,Radius,dRadius,ddRadius,g[i0],x[i0],y[i0]);
            }
            
            w.push_back(temp);

            if (i>0)
            {
              
            base = (b[i_1]+b[i0])/2;
            dRadius = (dR[i_1]+dR[i0])/2;
            ddRadius = (ddR[i_1]+ddR[i0])/2; 
            Radius = (Ra[i_1]+Ra[i0])/2;         
            temp = CV(h[i0],u[i0],v[i0],base,Radius,dRadius,ddRadius,(g[i_1]+g[i0])/2,(x[i_1]+x[i0])/2,y[i0]);
            }

            wl.push_back(temp);


             if (i<rows-1)
            {
            
            base = (b[i1]+b[i0])/2;
            dRadius = (dR[i1]+dR[i0])/2;
            ddRadius = (ddR[i1]+ddR[i0])/2;
            Radius = (Ra[i1]+Ra[i0])/2; 
            temp = CV(h[i0],u[i0],v[i0],base,Radius,dRadius,ddRadius, (g[i1]+g[i0])/2,(x[i1]+x[i0])/2,y[i0]);
            }  
            
            wr.push_back(temp); 

            base = (b[j1]+b[j0])/2;
            dRadius = (dR[j1]+dR[j0])/2;
            ddRadius = (ddR[j1]+ddR[j0])/2;
            Radius = (Ra[j1]+Ra[j0])/2; 
            temp = CV(h[j0],u[j0],v[j0],base,Radius, dRadius, ddRadius, (g[j1]+g[j0])/2,x[j0],(y[j1]+y[j0])/2);

            wt.push_back(temp);

            base = (b[j_1]+b[j0])/2;
            dRadius = (dR[j_1]+dR[j0])/2;
            ddRadius = (ddR[j_1]+ddR[j0])/2;
            Radius = (Ra[j_1]+Ra[j0])/2; 
            temp = CV(h[j0],u[j0],v[j0],base,Radius, dRadius, ddRadius, (g[j_1]+g[j0])/2,x[j0],(y[j_1]+y[j0])/2);

            wb.push_back(temp);

            
        } 
    } 
    }
    else
    {
       cout<<"There already appears to be elements in the w vector"<<endl;
    } 
	    
}



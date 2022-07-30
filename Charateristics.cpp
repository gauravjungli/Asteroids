#include "gauravlib.h"


double eigen(CV w,double x)

{
	double e1, e2, e3; //eigenvalues
	double root;
		root = 0;//change for the sphere
		if (root < 0)
		{//cout<<"non-hperbolic";
			root = 0;
		}
		else
		{
			root = 0;//chnage for the sphere
		}
		e1 = abs(w.q / w.p);
		e2 = abs(e1 - root);
		e3 = abs(e1 + root);

		return (std::max(e1, std::max(e2, e3)));
	
}

double ax(CV wl, CV wr, double x)
{    	

	return (max(eigen(wl, x), eigen(wr, x))); 
}

void edge(vector<CV>& w, vector<CV>& wl, vector<CV>& wr)
{   
    for(int i=1; i<res-1;i++)
    {
        CV deriv= derivative(w[i-1],w[i],w[i+1]);
        wr[i].modify(  w[i].p - dx * deriv.p / 2,
                    w[i].q - dx * deriv.q / 2,
                    w[i].r - dx * deriv.r/ 2);
        deriv= derivative(w[i-2],w[i-1],w[i]);            
        wl[i].modify(  w[i-1].p + dx * deriv.p / 2,
                    w[i-1].q + dx * deriv.q / 2,
                    w[i-1].r + dx * deriv.r / 2);
    }
    

}
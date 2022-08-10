#include "gauravlib.h"


double eigen(CV w)

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

double ax(CV wl, CV wr)
{    	

	return (max(eigen(wl), eigen(wr))); 
}

void edge(vector<CV>& w, vector<CV>& wl, vector<CV>& wr, double Om)
{   
    for(int i=1; i<res-1;i++)
    {	
        reconstruct(wr[i],w[i-1],w[i],w[i+1],-1);
		wr[i].modify(  wr[i].p,wr[i].q,wr[i].r,Om);

        reconstruct(wl[i],w[i-2],w[i-1],w[i],1);            
        wl[i].modify( wl[i].p,wl[i].q,wl[i].r,Om);
    }
}

  void reconstruct(CV& wl, CV w1, CV w2, CV w3, int sign )
  {
	wl.p=w2.p+sign*dx*derivative(w1.p,w2.p,w3.p)/2;
	wl.q=w2.q+sign*dx*derivative(w1.q,w2.q,w3.q)/2;
	wl.r=w2.r+sign*dx*derivative(w1.r,w2.r,w3.r)/2;
	wl.b=w2.b+sign*dx*derivative(w1.b,w2.b,w3.b)/2;
	wl.x=w2.x+sign*dx*derivative(w1.x,w2.x,w3.x)/2;
	wl.g.X1=w2.g.X1+sign*dx*derivative(w1.g.X1,w2.g.X1,w3.g.X1)/2;
	wl.g.X2=w2.g.X2+sign*dx*derivative(w1.g.X2,w2.g.X2,w3.g.X2)/2;
	wl.g.X3=w2.g.X3+sign*dx*derivative(w1.g.X3,w2.g.X3,w3.g.X3)/2;

  } 
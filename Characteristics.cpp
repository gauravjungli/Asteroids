#include "gauravlib.h"




double Ax(CV wl, CV wr)
{    	

	return (max(Eigen(wl), Eigen(wr))); 
}

void Edge(vector<CV>& w, vector<CV>& wl, vector<CV>& wr, double om)
{   
    for(int i=1; i<res-1;i++)
    {	
        Reconstruct(wr[i],w[i-1],w[i],w[i+1],1);
		wr[i].Modify(  wr[i].p,wr[i].q,wr[i].r,om);

        Reconstruct(wl[i],w[i-1],w[i],w[i+1],-1);            
        wl[i].Modify( wl[i].p,wl[i].q,wl[i].r,om);
    }
}

  void Reconstruct(CV& w, CV w1, CV w2, CV w3, int sign )
  {
	w.p=w2.p+sign*dx*Derivative(w1.p,w2.p,w3.p)/2;
	w.q=w2.q+sign*dx*Derivative(w1.q,w2.q,w3.q)/2;
	w.r=w2.r+sign*dx*Derivative(w1.r,w2.r,w3.r)/2;
	w.b=w2.b+sign*dx*Derivative(w1.b,w2.b,w3.b)/2;
	w.x=w2.x+sign*dx*Derivative(w1.x,w2.x,w3.x)/2;
	w.g.X1=w2.g.X1+sign*dx*Derivative(w1.g.X1,w2.g.X1,w3.g.X1)/2;
	w.g.X2=w2.g.X2+sign*dx*Derivative(w1.g.X2,w2.g.X2,w3.g.X2)/2;
	w.g.X3=w2.g.X3+sign*dx*Derivative(w1.g.X3,w2.g.X3,w3.g.X3)/2;

  } 
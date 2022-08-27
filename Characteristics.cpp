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
		wr[i]=CV(  wr[i].h,wr[i].u,wr[i].v,wr[i].b,wr[i].g,wr[i].x,om);

        Reconstruct(wl[i],w[i-1],w[i],w[i+1],-1);            
        wl[i]=CV(  wl[i].h,wl[i].u,wl[i].v,wl[i].b,wl[i].g,wl[i].x,om);
    }
}

  void Reconstruct(CV& w, CV w1, CV w2, CV w3, int sign )
  {
	w.h=w2.h+sign*dx*Derivative(w1.h,w2.h,w3.h)/2;
	w.u=w2.u+sign*dx*Derivative(w1.u,w2.u,w3.u)/2;
	w.v=w2.v+sign*dx*Derivative(w1.v,w2.v,w3.v)/2;
	w.b=w2.b+sign*dx*Derivative(w1.b,w2.b,w3.b)/2;
	w.x=w2.x+sign*dx*Derivative(w1.x,w2.x,w3.x)/2;
	w.g.X1=w2.g.X1+sign*dx*Derivative(w1.g.X1,w2.g.X1,w3.g.X1)/2;
	w.g.X2=w2.g.X2+sign*dx*Derivative(w1.g.X2,w2.g.X2,w3.g.X2)/2;
	w.g.X3=w2.g.X3+sign*dx*Derivative(w1.g.X3,w2.g.X3,w3.g.X3)/2;

  } 
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

void edge(vector<CV>& w, vector<CV>& wl, vector<CV>& wr)
{   
    for(int i=1; i<res-1;i++)
    {	
        LE(wr[i],w[i-1],w[i],w[i+1]);
		wr[i].modify(  wr[i].p,wr[i].q,wr[i].r);

        RE(wl[i],w[i-2],w[i-1],w[i]);            
        wl[i].modify( wl[i].p,wl[i].q,wl[i].r);
    }
}
void LE(CV& wr,CV w1, CV w2, CV w3 )
  {
	wr.p=w2.p-dx*derivative(w1.p,w2.p,w3.p)/2;
	wr.q=w2.q-dx*derivative(w1.q,w2.q,w3.q)/2;
	wr.r=w2.r-dx*derivative(w1.r,w2.r,w3.r)/2;
	wr.b=w2.b-dx*derivative(w1.b,w2.b,w3.b)/2;
	wr.x=w2.x-dx*derivative(w1.x,w2.x,w3.x)/2;
	wr.g.X=w2.g.X-dx*derivative(w1.g.X,w2.g.X,w3.g.X)/2;
	wr.g.Y=w2.g.Y-dx*derivative(w1.g.Y,w2.g.Y,w3.g.Y)/2;
	wr.g.Z=w2.g.Z-dx*derivative(w1.g.Z,w2.g.Z,w3.g.Z)/2;

  } 

  void RE(CV& wl, CV w1, CV w2, CV w3 )
  {
	wl.p=w2.p+dx*derivative(w1.p,w2.p,w3.p)/2;
	wl.q=w2.q+dx*derivative(w1.q,w2.q,w3.q)/2;
	wl.r=w2.r+dx*derivative(w1.r,w2.r,w3.r)/2;
	wl.b=w2.b+dx*derivative(w1.b,w2.b,w3.b)/2;
	wl.x=w2.x+dx*derivative(w1.x,w2.x,w3.x)/2;
	wl.g.X=w2.g.X+dx*derivative(w1.g.X,w2.g.X,w3.g.X)/2;
	wl.g.Y=w2.g.Y+dx*derivative(w1.g.Y,w2.g.Y,w3.g.Y)/2;
	wl.g.Z=w2.g.Z+dx*derivative(w1.g.Z,w2.g.Z,w3.g.Z)/2;

  } 
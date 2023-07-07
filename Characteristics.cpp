#include "gauravlib.h"

double max( FS w)
{
	return std::max({w.p,w.q,w.r});
}
double min( FS w)
{
	return std::min({w.p,w.q,w.r});
}


double Ax(CV wl, CV wr,string s)
{    	
if (s=="max")
	return max(max(Eigen(wl)),max( Eigen(wr))); 
else 
	return min(min(Eigen(wl)), min(Eigen(wr)));
}



void Edge(vector<CV>& w, vector<CV>& wl, vector<CV>& wr)
{   
    for(int i=2; i<res-1;i++)
    {	
		Reconstruct(wl[i],w[i-1],w[i],w[i+1],-1);            
        

        Reconstruct(wr[i-1],w[i-2],w[i-1],w[i],1);
		
		wl[i].b=wr[i-1].b=(wl[i].b+wr[i-1].b)/2;

	//	Balancing(w,wl,wr,i);
	

		wl[i]=CV(  wl[i].h,wl[i].u,wl[i].v,wl[i].b,wl[i].g,wl[i].x);
	//	wl[i].Modify(wl[i].p,wl[i].q,wl[i].r);  

		wr[i-1]=CV(  wr[i-1].h,wr[i-1].u,wr[i-1].v,wr[i-1].b,wr[i-1].g,wr[i-1].x);
	//	wr[i].Modify(wr[i].p,wr[i].q,wr[i].r);	
	//if (w[i-1].w<wr[i-1].b || w[i-1].w<wl[i-1].b)
	//		cout<<"Partially filled cells "<< i-1<<  endl;

	 }
}

 void Reconstruct(CV& w, CV w1, CV w2, CV w3, int sign )
  {
	w.h=w2.h+sign*dx*Derivative(w1.h,w2.h,w3.h)/2;
	w.w=w2.w+sign*dx*Derivative(w1.w,w2.w,w3.w)/2;
	w.b=w2.b+sign*dx*Derivative(w1.b,w2.b,w3.b)/2;
	w.u=w2.u+sign*dx*Derivative(w1.u,w2.u,w3.u)/2;
	w.v=w2.v+sign*dx*Derivative(w1.v,w2.v,w3.v)/2;
	w.x=w2.x+sign*dx/2.0;//can add theta here
	w.g.X1=w2.g.X1+sign*dx*Derivative(w1.g.X1,w2.g.X1,w3.g.X1)/2;
	w.g.X2=w2.g.X2+sign*dx*Derivative(w1.g.X2,w2.g.X2,w3.g.X2)/2;
	w.g.X3=w2.g.X3+sign*dx*Derivative(w1.g.X3,w2.g.X3,w3.g.X3)/2;

  }



  /*void Reconstruct(CV& w, CV w1, CV w2, CV w3, int sign )
  {
	w.p=w2.p+sign*dx*Derivative(w1.p,w2.p,w3.p)/2;
	w.q=w2.q+sign*dx*Derivative(w1.q,w2.q,w3.q)/2;
	w.r=w2.r+sign*dx*Derivative(w1.r,w2.r,w3.r)/2;

  }*/

  void Balancing (vector<CV>& w, vector<CV>& wl, vector<CV>& wr, int i)
  {	
		
			
		if (wr[i-1].w<wr[i-1].b)
		{
			wr[i-1].w=wr[i-1].b;
			wl[i-1].w=2*w[i-1].w-wr[i-1].b;
		}
		if (wl[i].w<wl[i].b)
		{
			wr[i].w=2*w[i].w-wl[i].b;
			wl[i].w=wl[i].b;
		}

		wl[i].h=wl[i].w-wl[i].b;
		wr[i-1].h=wr[i-1].w-wr[i-1].b; 
	
  }


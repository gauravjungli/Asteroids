#include "gauravlib.h"

double max( FS w)
{
	return std::max({w.p,w.q,w.r});
}
double min( FS w)
{
	return std::min({w.p,w.q,w.r});
}


double Ax(CV wl, CV wr, string s)
{    	
if (s=="max")
	return max({max(Eigenx(wl)),max( Eigenx(wr))}); 
else 
	return min({min(Eigenx(wl)), min(Eigenx(wr))});
}

double Ay(CV wb, CV wt, string s)
{    	
if (s=="max")
	return max({max(Eigeny(wb)),max( Eigeny(wt))}); 
else 
	return min({min(Eigeny(wb)), min(Eigeny(wt))});
}

void Edge(vector<CV>& w, vector<CV>& wl, vector<CV>& wr,  vector<CV>& wb, vector<CV>& wt)
{   
    for(int i=2; i<rows-1;i++)
    {	
		for (int j=0; j<cols;j++)
		{
			int i0 = index(i,j);
			int i_1 = index(i-1,j);
			int i_2 = index(i-2,j);
			int i1 = index(i+1,j);

			int j0 = index(i,j);
			int j_1 = index(i,j-1);
			int j_2 = index(i,j-2);
			int j1 = index(i,j+1);


			Reconstruct(wl[i0],w[i_1],w[i0],w[i1],-1);            
			Reconstruct(wb[j0],w[j_1],w[j0],w[j1],-1);  

        	Reconstruct(wr[i_1],w[i_2],w[i_1],w[i0],1);
			Reconstruct(wt[j_1],w[j_2],w[j_1],w[j0],1);
		
			wl[i0].b=wr[i_1].b=(wl[i0].b+wr[i_1].b)/2;

	//	Balancing(w,wl,wr,i);
	

			wl[i0]=CV(  wl[i0].h,wl[i0].u,wl[i0].v,wl[i0].b,wl[i0].g,wl[i0].x,wl[i0].y);
			wb[j0]=CV(  wb[j0].h,wb[j0].u,wb[j0].v,wb[j0].b,wb[j0].g,wb[j0].x,wb[j0].y);
	//	wl[i0].Modify(wl[i0].p,wl[i0].q,wl[i0].r);  

			wr[i_1]=CV(  wr[i_1].h,wr[i_1].u,wr[i_1].v,wr[i_1].b,wr[i_1].g,wr[i_1].x,wr[i_1].y);
			wt[j_1]=CV(  wt[j_1].h,wt[j_1].u,wt[j_1].v,wt[j_1].b,wt[j_1].g,wt[j_1].x,wt[j_1].y);
	//	wr[i0].Modify(wr[i0].p,wr[i0].q,wr[i0].r);	
	//if (w[i_1].w<wr[i_1].b || w[i_1].w<wl[i_1].b)
	//		std::cout<<"Partially filled cells "<< i_1<<  endl;
		}
	 }
}

 void Reconstruct(CV& w, CV w1, CV w2, CV w3, int sign )
  {
	
	w.h=w2.h+sign*dx*Derivative(w1.h,w2.h,w3.h)/2;//currently only using first order scheme. So commented most of it. Uncomment for higher order schemes
	w.w=w2.w+sign*dx*Derivative(w1.w,w2.w,w3.w)/2;
	w.b=w2.b+sign*dx*Derivative(w1.b,w2.b,w3.b)/2;
	w.u=w2.u+sign*dx*Derivative(w1.u,w2.u,w3.u)/2;
	w.v=w2.v+sign*dx*Derivative(w1.v,w2.v,w3.v)/2;
	w.x=w2.x+sign*dx/2.0;//can add theta here 
	w.y=w2.y+sign*dy/2.0;
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
  {	for (int j=0; j<cols;j++)
	{
	int i0 = index(i,j);
	int i_1 = index(i-1,j);
	int i_2 = index(i-2,j);
	int i1 = index(i+1,j);
			
		if (wr[i_1].w<wr[i_1].b)
		{
			wr[i_1].w=wr[i_1].b;
			wl[i_1].w=2*w[i_1].w-wr[i_1].b;
		}
		if (wl[i0].w<wl[i0].b)
		{
			wr[i0].w=2*w[i0].w-wl[i0].b;
			wl[i0].w=wl[i0].b;
		}

		wl[i0].h=wl[i0].w-wl[i0].b;
		wr[i_1].h=wr[i_1].w-wr[i_1].b; 
	}
  }


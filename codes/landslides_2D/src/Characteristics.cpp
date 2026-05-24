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
    for(int i=2; i<rows-2;i++)
    {	
		for (int j=0; j<cols;j++)
		{
			int i0 = index(i,j);
			int i_1 = index(i-1,j);
			int i1 = index(i+1,j);

			int j0 = index(i,j);
			int j_1 = index(i,j-1);
			int j1 = index(i,j+1);

		 if (reconst=="Cartesian conserved variables")
		{
			Reconstruct_U(wl[i0],w[i_1],w[i0],w[i1],-1,dx);            
			Reconstruct_U(wb[j0],w[j_1],w[j0],w[j1],-1,dy);  
        	Reconstruct_U(wr[i0],w[i_1],w[i0],w[i1],1,dx);
			Reconstruct_U(wt[j0],w[j_1],w[j0],w[j1],1,dy);

			wl[i0].Modify_U(wl[i0].P,wl[i0].Q,wl[i0].R);  		
			wr[i0].Modify_U(wr[i0].P,wr[i0].Q,wr[i0].R);
			wb[j0].Modify_U(wb[j0].P,wb[j0].Q,wb[j0].R);  		
			wt[j0].Modify_U(wt[j0].P,wt[j0].Q,wt[j0].R);
		}

		else
		{
			cout<<"No valid reconstruction scheme selected. No reconstruction done."<<endl;
			return ;
		}



		}
	 }
	 	BC(wb,wb);
	 	BC(wt,wt);
	 	BC(wl,wr);
	 	BC(wr,wl);

}

//  void Reconstruct(CV& w, CV w1, CV w2, CV w3, int sign )
//   {
	
// 	w.h=w2.h+sign*dx*Derivative(w1.h,w2.h,w3.h)/2;//currently only using first order scheme. So commented most of it. Uncomment for higher order schemes
// 	w.w=w2.w+sign*dx*Derivative(w1.w,w2.w,w3.w)/2;
// 	w.b=w2.b+sign*dx*Derivative(w1.b,w2.b,w3.b)/2;
// 	w.u=w2.u+sign*dx*Derivative(w1.u,w2.u,w3.u)/2;
// 	w.v=w2.v+sign*dx*Derivative(w1.v,w2.v,w3.v)/2;
// 	w.x=w2.x+sign*dx/2.0;//can add theta here 
// 	w.y=w2.y+sign*dy/2.0;
// 	w.g.X1=w2.g.X1+sign*dx*Derivative(w1.g.X1,w2.g.X1,w3.g.X1)/2;
// 	w.g.X2=w2.g.X2+sign*dx*Derivative(w1.g.X2,w2.g.X2,w3.g.X2)/2;
// 	w.g.X3=w2.g.X3+sign*dx*Derivative(w1.g.X3,w2.g.X3,w3.g.X3)/2;

//   }

   void Reconstruct_U(CV& w, CV w1, CV w2, CV w3, int sign, double dsize )
  {
	if (w2.h<epsilon*epsilon)
		sign = 0;
	w.P=w2.P+sign*dsize*Derivative(w1.P,w2.P,w3.P,dsize)/2;
	w.Q=w2.Q+sign*dsize*Derivative(w1.Q,w2.Q,w3.Q,dsize)/2;
	w.R=w2.R+sign*dsize*Derivative(w1.R,w2.R,w3.R,dsize)/2;

  } 

    void Reconstruct(CV& w, CV w1, CV w2, CV w3, int sign, double dsize )
  {
	w.p=w2.p+sign*dx*Derivative(w1.p,w2.p,w3.p,dsize)/2;
	w.q=w2.q+sign*dx*Derivative(w1.q,w2.q,w3.q,dsize)/2;
	w.r=w2.r+sign*dx*Derivative(w1.r,w2.r,w3.r,dsize)/2;

  } 

  void Balancing (vector<CV>& w, vector<CV>& wl, vector<CV>& wr,vector<CV>& wb, vector<CV>& wt,int i)
  {	
	for (int j=0; j<cols;j++)
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


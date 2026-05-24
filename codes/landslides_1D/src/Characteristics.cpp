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
    for(int i=2; i<res-2;i++)
    {	
		if (reconst == "Primitive variables")
		{
		Reconstruct( wl[i], wr[i], w[i-1],  w[i], w[i+1] );   
		wl[i]=CV( wl[i].h, wl[i].u, wl[i].v, wl[i].b, wl[i].db, wl[i].ddb, wl[i].g, wl[i].x );
		wr[i]=CV( wr[i].h, wr[i].u, wr[i].v, wr[i].b, wr[i].db, wr[i].ddb, wr[i].g, wr[i].x );
		}

		else if (reconst=="Cartesian conserved variables")
		{
			Reconstruct_U( wr[i], w[i-1],  w[i], w[i+1], 1 );      
			Reconstruct_U( wl[i], w[i-1],  w[i], w[i+1], -1 );
			wl[i].Modify_U(wl[i].P,wl[i].Q,wl[i].R);  		
			wr[i].Modify_U(wr[i].P,wr[i].Q,wr[i].R);
		}

		else if (reconst == "Conserved variables")
		{
			Reconstruct( wr[i], w[i-1],  w[i], w[i+1], 1 );      
			Reconstruct( wl[i], w[i-1],  w[i], w[i+1], -1 );
			wl[i].Modify(wl[i].p,wl[i].q,wl[i].r);  
			wr[i].Modify(wr[i].p,wr[i].q,wr[i].r);
		}

				else
		{
			cout<<"No valid reconstruction scheme selected. No reconstruction done."<<endl;
			return ;
		}

	//if (w[i-1].w<wr[i-1].b || w[i-1].w<wl[i-1].b)
	//		std::cout<<"Partially filled cells "<< i-1<<  endl;
		
	 }
	 BC(wl,wr);
	 BC(wr,wl);
}

 void Reconstruct(CV& wl, CV& wr, CV w1, CV w2, CV w3 )
  { int sign;

	wl.h=w2.h + sign*dx*Derivative(w1.h,w2.h,w3.h)/2;
	wl.u=w2.u + sign*dx*Derivative(w1.u,w2.u,w3.u)/2;
	wl.v=w2.v + sign*dx*Derivative(w1.v,w2.v,w3.v)/2;


	wr.h=w2.h + sign*dx*Derivative(w1.h,w2.h,w3.h)/2;
	wr.u=w2.u + sign*dx*Derivative(w1.u,w2.u,w3.u)/2;
	wr.v=w2.v + sign*dx*Derivative(w1.v,w2.v,w3.v)/2;


  }



   void Reconstruct_U(CV& w, CV w1, CV w2, CV w3, int sign )
  {
	if (w2.h<epsilon*epsilon)
		sign = 0;
	w.P=w2.P+sign*dx*Derivative(w1.P,w2.P,w3.P)/2;
	w.Q=w2.Q+sign*dx*Derivative(w1.Q,w2.Q,w3.Q)/2;
	w.R=w2.R+sign*dx*Derivative(w1.R,w2.R,w3.R)/2;

  } 

   void Reconstruct(CV& w, CV w1, CV w2, CV w3, int sign )
  {
	w.p=w2.p+sign*dx*Derivative(w1.p,w2.p,w3.p)/2;
	w.q=w2.q+sign*dx*Derivative(w1.q,w2.q,w3.q)/2;
	w.r=w2.r+sign*dx*Derivative(w1.r,w2.r,w3.r)/2;

  } 

  void Balancing (CV& w, CV& wl, CV& wr)
  {	
		
			
		if (wr.w<wr.b)
		{	cout <<"Dry region appeared" << endl;
			wr.w=wr.b;
			wl.w=2*w.w-wr.b;
		}
		if (wl.w<wl.b)
		{	cout <<"Dry region appeared" << endl;
			wr.w=2*w.w-wl.b;
			wl.w=wl.b;
		}

		wl.h= wl.w - wl.b;
		wr.h= wr.w - wr.b; 
	
  }


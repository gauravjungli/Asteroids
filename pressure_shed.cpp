#include "gauravlib.h"
void Shed(vector<CV>& w)
{
for ( int j = 2; j < res-2; j++)
{

 	if (Psi(w[j]) < 0) 
	{
		w[j]=CV(  min_h,0,0,w[j].b,w[j].g,w[j].x  );
		 cout<<"Mass Shedding Happens"<<endl;
	} 
}
}
double Psi(CV w)
{
return	-(omega*omega*sin(w.x)*sin(w.x)+2*omega*w.v*sin(w.x)+w.g.X1+w.u*w.u+w.v*w.v);
}
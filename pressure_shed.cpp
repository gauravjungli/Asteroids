#include "gauravlib.h"
void Shed(vector<CV>& w, ofstream& myfile, double& Ang_Shed)
{

for ( int j = 2; j < res-2; j++)
{

 	if (Psi(w[j]) <= 1e-12) 
	{	
		Ang_Shed=Ang_Shed+(PI/2*(w[j].v*(pow(1+gamma*w[j].b+epsilon*w[j].h,4)-pow(1+gamma*w[j].b,4)))+2*PI/5*(
					omega*(pow(1+gamma*w[j].b+epsilon*w[j].h,5)-pow(1+gamma*w[j].b,4))))*dx;
		w[j]=CV(  min_h,0,0,w[j].b,w[j].g,w[j].x  );
		 
	} 
}
}
double Psi(CV w)
{
return	-(omega*omega*sin(w.x)*sin(w.x)+2*omega*w.v*sin(w.x)+w.g.X1+w.u*w.u+w.v*w.v);
}
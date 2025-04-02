#include "gauravlib.h"
void Shed(vector<CV>& w, double& Ang_Shed)
{
for (int i =2;i<rows-2;i++)
{

for ( int j = 0; j < cols; j++)
{
	int i0 = index(i,j);
	int i_1 = index(i-1,j);
	int i1 = index(i+1,j);
	//change
//	w[i0].psi=w[i0].psi/(1-tan(delta/180*PI)*Gamma*((w[i0].u)/(w[i0].u*w[i0].u+w[i0].v*w[i0].v)*max((w[i0].b-w[i_1].b),(w[i1].b-w[i0].b))/dx+
//				(w[i0].v)/(w[i0].u*w[i0].u+w[i0].v*w[i0].v)*max((w[i1].b-w[i0].b),(w[i0].b-w[i_1].b))/dy));
 	if (w[i0].psi <= pow(epsilon,1) ) 
	{	
	//	Ang_Shed=Ang_Shed+(PI/2*(w[j].v*(pow(1+Gamma*w[j].b+epsilon*w[j].h,4)-pow(1+Gamma*w[j].b,4)))+2*PI/5*(
	//				omega*(pow(1+Gamma*w[j].b+epsilon*w[j].h,5)-pow(1+Gamma*w[j].b,5))))*dx;
		
	//	Ang_Shed=Ang_Shed-2*PI/5*(omega*(pow(1+Gamma*w[j].b+epsilon*min_h,5)-pow(1+Gamma*w[j].b,5)))*dx;
		w[i0]=CV(  min_h,0,0,w[i0].b,w[i0].g,w[i0].x ,w[i0].y );
		cout<<"mass shedding is happening"<<"  "<<i<<"   "<<j<<endl;
	} 

}}

}
double Psi(CV w)
{
return	-(omega*omega*sin(w.x)*sin(w.x)+2*omega*w.v*sin(w.x)+w.g.X1+(w.u*w.u+w.v*w.v)/w.lambda);
}
#include "gauravlib.h"

double Inertia(vector<CV>& w, int no)
{
	double integral1=8*PI/15.0; 
	double integral2=8*PI/15.0;
	for (int i=1;i<res;i=i+2)
	{
	integral1+=2*PI/5*2*dx/6*((pow(1+gamma*w[i-1].b+epsilon*w[i-1].h,5)-1)*pow(sin(w[i-1].x),3)+
	4*(pow(1+gamma*w[i].b+epsilon*w[i].h,5)-1)*pow(sin(w[i].x),3)+
	(pow(1+gamma*w[i+1].b+epsilon*w[i+1].h,5)-1)*pow(sin(w[i+1].x),3)); 
	
	integral2+=PI/5*(pow(1+gamma*w[i].b+epsilon*w[i].h,5)-1)*(2-pow(sin(w[i].x),2))*sin(w[i].x)*dx; 
	}

	if (no==1)
		return integral1;
	else
		return integral2;
}

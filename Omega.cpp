#include "gauravlib.h"

double Inertia(vector<CV>& w, int no)
{
	double integral1=8*PI/15.0; 
	double integral2=8*PI/15.0;
	for (int i=0;i<res;i++)
	{
	integral1+=2*PI/5*(pow(1+gamma*w[i].b+epsilon*w[i].h,5)-1)*pow(sin(w[i].x),3)*dx; 
	
	integral2+=PI/5*(pow(1+gamma*w[i].b+epsilon*w[i].h,5)-1)*(2-pow(sin(w[i].x),2))*sin(w[i].x)*dx; 
	}

	if (no==1)
		return integral1;
	else
		return integral2;
}

#include "gauravlib.h"

double Alpha( double& om , double dt, AMB amb )
{  	
	return -(1.0/2*Diff(amb.ang_mom,dt)+2.0/5.0*om*Diff(amb.inertia,dt))/(2.0/5*amb.ang_mom[4]+8.0/15);
}

void Integrals(double& integral1, double& integral2, vector<CV>& w, vector<CV>& wl, vector<CV>& wr)
{	integral1=integral2=0;
	for (int i=2;i<res-2;i++)
	{
	integral1+=dx/6*(Int_B(wl[i])+4*Int_B(w[i])+Int_B(wr[i]));
	integral2+=dx/6*(Int_A(wl[i])+4*Int_A(w[i])+Int_A(wr[i]));
	}
}

void Ang_mom(AMB& amb, double integral1, double integral2)
{
		amb.ang_mom.erase(amb.ang_mom.begin());
		amb.ang_mom.push_back(integral1);
		amb.inertia.erase(amb.inertia.begin());
		amb.inertia.push_back(integral1);
}

double Int_B(CV w)
{
	return w.v*(pow((1+epsilon*(w.b+w.h)),4)-pow((1+epsilon*(w.b)),4))*sin(w.x);
}
double Int_A(CV w)
{
	return (pow((1+epsilon*(w.b+w.h)),5)-pow((1),5))*pow(sin(w.x),3);
}

double Diff(vector<double> w, double dt)
{
	return (w[0]-16/3.0*w[1]+12*w[2]-16*w[3]+25/3.0*w[4])/(4*dt);
}

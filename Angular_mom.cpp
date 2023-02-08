#include "gauravlib.h"

double Alpha( double& om , double dt, AMB amb, double inertia )
{  	
	return -15.0/4.0*(Diff(amb.ang_mom,dt)+omega*(1+epsilon*om)*Diff(amb.inertia,dt))*inertia;
}

void Integrals(double& integral1, double& integral2,  vector<CV>& w, vector<CV>& wl, vector<CV>& wr)
{	integral1=integral2=0;
	for (int i=2;i<res-2;i++)
	{
	integral1+=dx/6*(Int_B(wl[i])+4*Int_B(w[i])+Int_B(wr[i]));
	integral2+=dx/6*(Int_A(w[i-1],wl[i],w[i])+4*Int_A(wl[i],w[i],wr[i])+Int_A(w[i],wr[i],w[i+1]));
	}
	
}

void Ang_mom(AMB& amb, double integral1, double integral2, bool replace)
{
	if (replace)
	{
		amb.ang_mom.erase(amb.ang_mom.begin());
		amb.ang_mom.push_back(integral1);
		amb.inertia.erase(amb.inertia.begin());
		amb.inertia.push_back(integral2);
	}
	else 
	{
		amb.ang_mom[4]=integral1;
		amb.inertia[4]=integral2;
	}
}

double Int_B(CV w)
{
	return (w.v+epsilon*omega*w.h*sin(w.x))*w.h*(1+3*epsilon*w.lambda)*pow(sin(w.x),2);
}
double Int_A(CV wl, CV w, CV wr)
{
	return w.h*(1+2*epsilon*w.lambda)*pow(sin(w.x),1)*((sin(w.x-dx/2)*sin(w.x-dx/2)+sin(w.x+dx/2)*sin(w.x+dx/2))/2
			+2*epsilon*(wl.b*wl.b*sin(w.x-dx/2)*sin(w.x-dx/2)+wr.b*wr.b*sin(w.x+dx/2)*sin(w.x+dx/2))/2);
}

double Diff(vector<double> w, double dt)
{
		//cout<<(w[4]-w[3])<<endl;

	return (w[0]-16/3.0*w[1]+12*w[2]-16*w[3]+25/3.0*w[4])/(4*dt);
}

double AMB_corrector (vector<CV>& w, vector<CV>& wl, vector<CV>& wr)
{
	double multiplier=0;
	for (int i=2;i<res-2;i++)
	{
	multiplier+=dx/6*((wl[i].h+wl[i].b)*pow(sin(wl[i].x),3)+4*(w[i].h+w[i].b)*pow(sin(w[i].x),3)+(wr[i].h+wr[i].b)*pow(sin(wr[i].x),3));
	}
	multiplier=1-15.0/4.0*epsilon*multiplier;
	return multiplier;
}

#include "gauravlib.h"

void Ang_mom(const vector<CV>& w,const vector<CV>& wtemphat, double& Om ,const vector<double>& x, double dt,double momincb)
{  
	double integral1_o=0,integral2_o=0,integral1=0,integral2=0;
   double mass_total_o = 0;
	for (int j = 2; j <= res + 1; j++)
	{
		mass_total_o = mass_total_o + 2 * PI * w[j].p * dx;
		integral1 = integral1 ;//complete
		integral2 = integral2 ; //complete 
		integral1_o=integral1_o;
	}
    double dom = -(Om * (integral1 - integral1_o) + (integral2 - integral2_o)) / (momincb + integral1);
	double alpha = dom / dt;	// angular acceleration
    double	alpha_momin = -(Om * (integral1 - integral1_o)) / (momincb + integral1) / dt;
	double	alpha_angmom = -(integral2 - integral2_o) / (momincb + integral1) / dt;
}

#include "gauravlib.h"

void Ang_mom(const vector<CV>& w, const double& Om ,const vector<double>& x,double integral1,double integral2)
{
   double mass_total_o = 0;
			for (int j = 2; j <= res + 1; j++)
			{
				mass_total_o = mass_total_o + 2 * PI * w[j].p * dx;
				integral1 = integral1 ;//complete
				integral2 = integral2 ; //complete 
			}
}
void omega(double Om, double integral1,double integral2, double integral1_o, double integral2_o, double momincb, double dt )
{
    double dom = -(Om * (integral1 - integral1_o) + (integral2 - integral2_o)) / (momincb + integral1);
	double alpha = dom / dt;	// angular acceleration
    double	alpha_momin = -(Om * (integral1 - integral1_o)) / (momincb + integral1) / dt;
	double	alpha_angmom = -(integral2 - integral2_o) / (momincb + integral1) / dt;
}

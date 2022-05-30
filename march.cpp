#include "gauravlib.h"
void predictor_corrector(vector<CV>& w,  vector<CV>& wtemp, double Om, vector<double> x, double dt, vector<grav>& g)
{   double momincb=0;//chnage for the sphere
    double integral1_o=0;
    double integral2_o=0;
    Ang_mom( w, Om , x, integral1_o, integral2_o);

	// predictor step for interior
    vector<CV> wtemphat(w);
	for (int j = 0; j < res; j++)
	{
		w[j].p = wtemp[j].p - ((Hx(j + 1, 1, wtemp) - Hx(j, 1,wtemp)) / dx - S1(j,wtemp[j],g[j])) * dt;
		w[j].q = wtemp[j].q - ((Hx(j + 1, 2,wtemp) - Hx(j, 2,wtemp)) / dx - S2(j,wtemp[j],g[j])) * dt;
    	w[j].r = wtemp[j].r - ((Hx(j + 1, 3,wtemp) - Hx(j, 3,wtemp)) / dx - S3(j,wtemp[j], g[j])) * dt;
	}

	// corrector step for interior
    wtemp=w;
	for (int j=0; j < res; j++)
	{
	w[j].p = wtemphat[j].r * weight + (1 - weight) * (wtemp[j].p - ((Hx(j + 1, 1,wtemp) - Hx(j, 1,wtemp)) / dx - S1(j,wtemp[j],g[j])) * dt);
	w[j].q = wtemphat[j].r * weight + (1 - weight) * (wtemp[j].q - ((Hx(j + 1, 2,wtemp) - Hx(j, 2,wtemp)) / dx - S2(j,wtemp[j],g[j])) * dt);
	w[j].r = wtemphat[j].r * weight + (1 - weight) * (wtemp[j].r - ((Hx(j + 1, 3,wtemp) - Hx(j, 3,wtemp)) / dx - S3(j,wtemp[j],g[j])) * dt);
	}
    double integral1=0;
    double integral2=0;
    Ang_mom( w, Om , x, integral1, integral2);
    omega( Om, integral1, integral2, integral1_o, integral2_o, momincb,dt );
}
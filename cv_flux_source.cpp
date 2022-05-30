#include <gauravlib.h>

 CV::CV(double h, double u, double v,double x){
            this->h=h;this->h=u;this->v=v;
            this->p=h*sin(x); //change for specific case
			this->q=h*u*sin(x)+epsilon*sin(x)*(h*h*u)/2; // change for specific case
			this->r=h*v*sin(x)+epsilon*sin(x)*(h*h*v)/2; // change for specific case
        }


double flux( CV w, int j )
// no is 1,2,3 for f1,f2,f3 and fg 1,2 for f,g
{

	double f1 = w.;//change for specific case
	double f2 = w.q*w.q / w.p; //change for specific case
	double f3 = w.q * w.r / w.p;// change for the specific case

}



double S1(int j, const CV& wtemp, grav& g)
{
	return 0;//change for the specific case
}

double S2(int j, const CV& wtemp, const grav& g)//change for the specific case
{
	double ht = wtemp.h;
	double ut = wtemp.u;
	double vt = wtemp.v;
	int sgn = (ut > 0) ? 1 : ((ut < 0) ? -1 : 0);
	double Gt = 0;//change for the sphere
	return 0;
}

double S3(int j, const CV& wtemp, grav& g)//change for the specfic case
{
	double ht = wtemp.h;
	double ut = wtemp.u;
	double vt = wtemp.v;
	int sgn = (vt > 0) ? 1 : ((vt < 0) ? -1 : 0);
	return 0 ;//change for the sphere
}

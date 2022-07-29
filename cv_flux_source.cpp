#include <gauravlib.h>

CV::CV(): p(0), q(0), r(0),h(0),u(0),v(0){}

CV::CV(double h, double u, double v,double x)
{
    this->h=h;this->h=u;this->v=v;
    this->p=h*sin(x); //change for specific case
	this->q=h*u*sin(x)+epsilon*sin(x)*(h*h*u)/2; // change for specific case
	this->r=h*v*sin(x)+epsilon*sin(x)*(h*h*v)/2; // change for specific case
}

void CV::modify(double p, double q, double r)
{
	this->p=p; this->q=q; this->r=r;
	this->h=0;//change for specific case
	this->u=0;//change for specific case
	this->v=0;//change for the specific case
}

CV flux( CV w, double x )

{
	CV f;
	double f1 = w.p;//change for specific case
	double f2 = w.q*w.q / w.p; //change for specific case
	double f3 = w.q * w.r / w.p;// change for the specific case
	f.modify(f1,f2,f3);
	return f;
}



CV Source(int j, const CV& wtemp, grav& g)
{
	 CV source;
	double ht = wtemp.h;
	double ut = wtemp.u;
	double vt = wtemp.v;
	int sgn_u = (ut > 0) ? 1 : ((ut < 0) ? -1 : 0);
	double Gt = 0;//change for the sphere
	int sgn_v = (vt > 0) ? 1 : ((vt < 0) ? -1 : 0);
	source.modify(0,0,0) ;//change for the sphere
	return source;
}

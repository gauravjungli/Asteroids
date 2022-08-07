#include "gauravlib.h"

CV::CV(double h, double u, double v, double b, grav g, double x, double Om )
{	
    this->h=h;this->h=u;this->v=v;this->b=b;this->g=g;this->x=x;this->psi=Psi(*this,Om);
	lambda= 1+epsilon*(b+h/2.0);
    this->p=h*sin(x)*(1+2*epsilon*lambda); //change for specific case
	this->q=h*u*sin(x)*(1+3*epsilon*lambda); // change for specific case
	this->r=h*v*sin(x)*(1+3*epsilon*lambda); // change for specific case
}

void CV::modify(double p, double q, double r,double Om)
{	h=(-(1+2*epsilon*b)+pow(pow((1+2*epsilon*b),2)+2*epsilon*p/sin(x),0.5))/epsilon;
	lambda= 1+epsilon*(b+h/2.0);
	this->p=p; this->q=q; this->r=r;
	u=q/p*(1+2*epsilon*lambda)/(1+3*epsilon*lambda);//change for specific case
	v=r/p*(1+2*epsilon*lambda)/(1+3*epsilon*lambda);//change for specific case
	psi=Psi(*this, Om);
}

  

FS flux( CV w )

{
	FS f;
	 f.p = w.u*w.h*(1+epsilon*w.lambda)*sin(w.x);//change for specific case
	 f.q = (w.u*w.u+epsilon*w.psi*w.h)*w.h*(1+2*epsilon*w.lambda)*sin(w.x); //change for specific case
	 f.r = (w.u*w.v)*w.h*(1+2*epsilon*w.lambda)*sin(w.x);// change for the specific case
	return f;

}


FS Source(const CV& w, double Om)
{
	FS source;
	source.p=0;
	double bf=Om*Om*cos(w.x)*sin(w.x)*(1+4*epsilon*w.lambda)/(1+3*epsilon*w.lambda)+2*Om*cos(w.x)*w.v;
	double mu=0;
	double fr=-mu*w.psi*w.h*w.u/(pow(pow(w.u,2)+pow(w.v,2),0.5));
	source.q = (w.v*w.v+epsilon*w.psi*w.h)*w.h*(1+2*epsilon*w.lambda)+bf*w.h*sin(w.x)*(1+3*epsilon*w.lambda)+fr*sin(w.x)*(1+epsilon*w.b);
	fr=-mu*w.psi*w.h*w.v/(pow(pow(w.u,2)+pow(w.v,2),0.5));
	bf=-2*Om*cos(w.x)*w.u;
	source.r=(w.u*w.v)*w.h*(1+epsilon*w.lambda)+bf*w.h*sin(w.x)*(1+3*epsilon*w.lambda)+fr*sin(w.x)*(1+epsilon*w.b);;//change for the sphere
	return source;
}

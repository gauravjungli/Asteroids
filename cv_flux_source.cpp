#include "gauravlib.h"

CV::CV(double h, double u, double v, double b, Grav g, double x, double om )
{	
    this->h=h;this->u=u;this->v=v;this->b=b;this->g=g;this->x=x;this->psi=Psi(*this,om);
	lambda= (b+h/2.0);
    this->p=h*sin(x)*(1+2*epsilon*lambda); 
	this->q=h*u*sin(x)*(1+3*epsilon*lambda); 
	this->r=h*(v+epsilon*om*h)*sin(x)*(1+3*epsilon*lambda); 
}

void CV::Modify(double p, double q, double r,double om)
{	h=0.5*(-2*b*epsilon-1+pow((pow(2*epsilon*b+1,2)+4*epsilon*p/sin(x)),0.5))/epsilon;
	this->p=p; this->q=q; this->r=r;
	//h=H(om);
	lambda= (b+h/2.0);
	u=q/(h*(1+3*epsilon*lambda)*sin(x));
	//u=U(om);
	v=r/(h*sin(x)*(1+3*epsilon*lambda))-epsilon*om*h;
	//v=V(om);
	psi=Psi(*this, om);
}

double CV::H( double om)
{
	return -p*(2*b*epsilon*sin(x)+epsilon*p-sin(x))/pow(sin(x),2);
}
double CV::U( double om)
{
	return -q*(2*b*epsilon*sin(x)+epsilon*p-2*sin(x))/(2*sin(x)*p);
}

double CV::V( double om)
{
	return -r*(2*b*epsilon*sin(x)+epsilon*p-2*sin(x))/(2*sin(x)*p)-epsilon*om*p/sin(x);
}


FS Flux( CV w, double om )

{
	FS f;
	 f.p = w.u*w.h*(1+epsilon*w.lambda)*sin(w.x);
	 f.q = (w.u*w.u-epsilon*w.psi*w.h/2)*w.h*(1+2*epsilon*w.lambda)*sin(w.x); 
	 f.r = (w.u*(w.v+epsilon*om*w.h))*w.h*(1+2*epsilon*w.lambda)*sin(w.x);
	//	f.p=-w.q*(2*w.b*epsilon*sin(w.x)+epsilon*w.p-sin(w.x))/sin(w.x);
	//	f.q=-(4*w.b*epsilon*sin(w.x)*w.q*w.q+2*epsilon*w.p*w.q*w.q-2*sin(w.x)*w.q*w.q-epsilon*pow(w.p,3)*w.psi)/(2*sin(w.x)*w.p);
	//	f.r=-w.q*w.r*(2*w.b*epsilon*sin(w.x)+epsilon*w.p-sin(w.x))/(sin(w.x)*w.p);
	return f;

}


FS Source( CV w, double om, double alpha)
{
	FS source;
	source.p=0;
	double mu=tan(delta* PI / 180);
	double bf=(om*om*cos(w.x)*sin(w.x)+2*om*cos(w.x)*w.v+w.g.X2);
	double fr=0;
	if (pow(pow(w.u,2)+pow(w.v,2),0.5)>1e-10) fr=mu*(w.psi/2.0)*w.h*w.u/pow(pow(w.u,2)+pow(w.v,2),0.5);
	else fr=-bf;
	source.q = (w.v*w.v-epsilon*w.psi*w.h/2)*w.h*(1+2*epsilon*w.lambda)*cos(w.x)+bf*w.h*sin(w.x)*(1+3*epsilon*w.lambda)+fr*sin(w.x)*(1+3*epsilon*w.b);
	if (pow(pow(w.u,2)+pow(w.v,2),0.5)>1e-10) fr=mu*(w.psi*w.h/2)*w.v/pow(pow(w.u,2)+pow(w.v,2),0.5);
	else fr=-bf;
	bf=-2*om*cos(w.x)*w.u-alpha*sin(w.x)+w.g.X3;
	source.r=-(w.u*w.v)*w.h*(1+2*epsilon*w.lambda)*cos(w.x)+bf*w.h*sin(w.x)*(1+3*epsilon*w.lambda)+fr*sin(w.x)*(1+3*epsilon*w.b);
	return source;
}

double Eigen(CV w)

{
	double e1, e2, e3; //eigenvalues
	double root,base;
		root =sqrt(-epsilon*w.p*w.p*(8*sin(w.x)*w.b*epsilon*w.p*w.psi+4*epsilon*w.p*w.p*w.psi
				-4*sin(w.x)*w.p*w.psi-epsilon*w.q*w.q));
		base= 4*sin(w.x)*w.b*epsilon*w.q+3*epsilon*w.p*w.q-2*w.q*sin(w.x);
		e1 = abs(-w.q / (w.p*sin(w.x))*(2*w.b*epsilon*sin(w.x)+epsilon*w.p-sin(w.x)));
		e2 = abs(-1.0/(2*w.p*sin(w.x))*(base+root));
		e3 = abs(-1.0/(2*w.p*sin(w.x))*(base-root));

		return (std::max(e1, std::max(e2, e3)));
	
}

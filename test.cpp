#include "gauravlib.h"

CV::CV(double h, double u, double v, double b, Grav g, double x, double om )
{	
    this->h=h;this->u=u;this->v=v;this->b=b;
	this->g=g;this->x=x;this->psi=-g.X1;//(*this,om);
	lambda= 0;
    this->p=h; 
	this->q=h*u; 
	this->r=h*v; 
}

void CV::Modify(double p, double q, double r,double om)
{	h=p;
	this->p=p; this->q=q; this->r=r;
	//h=H(om);
	lambda= 0;
	u=q/h;
	//u=U(om);
	v=r/h;
	
}

double CV::H( double om)
{
	return p;
}
double CV::U( double om)
{
	return q/p;
}

double CV::V( double om)
{
	return r/p;
}


  


FS Flux( CV w, double om )

{
	FS f;
	 f.p = w.u*w.h;
	 f.q = (w.u*w.u-w.g.X1*w.h/2)*w.h; 
	 f.r = 0;
	return f;

}


FS Source( CV w, double om, double alpha)
{
	FS source;
	source.p=0;

	source.q = 0;
	source.r=0;
	return source;
}

FS Eigen(CV w)

{
	double root,base;
	FS e;
		root =sqrt(-w.g.X1*w.h);
		base= w.u;
		e.p = 0;
		e.q = (base+root);
		e.r = (base-root);

		return e;
	
}

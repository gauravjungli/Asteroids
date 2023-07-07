#include "gauravlib.h"

CV::CV(double h, double u, double v, double b, Grav g, double x )
{	
    this->h=h;this->u=u;this->v=v;this->b=b;
	this->g=g;this->x=x;this->psi=Psi(*this); this->w=h+b;
	lambda= (b+h/2.0);
    this->p=h; 
	this->q=h*u; 
	this->r=h*v; 
}


void CV::Modify(double p, double q, double r)
{	h=p;
	this->p=p; this->q=q; this->r=r; w=h+b;
	lambda= (b+h/2.0);
	u=q/h;
	v=r/h;
	psi=Psi(*this);
}


FS Flux( CV w )

{
	FS f;
	 f.p = w.u*w.h;
	 f.q = (w.u*w.u+w.psi*w.h/2)*w.h; 
	 f.r = w.u*w.v*w.h;
	//	f.p=-w.q*(2*w.b*epsilon*sin(w.x)+epsilon*w.p-sin(w.x))/sin(w.x);
	//	f.q=-(4*w.b*epsilon*sin(w.x)*w.q*w.q+2*epsilon*w.p*w.q*w.q-2*sin(w.x)*w.q*w.q-epsilon*pow(w.p,3)*w.psi)/(2*sin(w.x)*w.p);
	//	f.r=-w.q*w.r*(2*w.b*epsilon*sin(w.x)+epsilon*w.p-sin(w.x))/(sin(w.x)*w.p);
	return f;

}


FS Source( CV w, CV w1, CV w2, CV w3, CV w4,  double om, double alpha)
{
	FS source;
	FS bf=Body_force( w,  w1,  w2, w3,  w4, om, alpha);
	FS fr=Friction(w,bf);
	

	source.p=0;


	double grad_b= epsilon*((w3.b-w1.b)/dx)*((w1.h+w2.h+w3.h+w4.h)/4)*((w1.psi+w2.psi+w3.psi+w4.psi)/4)*(sin(w3.x)+sin(w1.x))/2;

	double pressure=(((pow((w1.v+w2.v)/2,2)+pow((w3.v+w4.v)/2,2))/2)*(1+2*epsilon*(w4.lambda+w3.lambda+w2.lambda+w1.lambda)/4)
					*((w1.h+w2.h+w3.h+w4.h)/4)+epsilon*((w1.psi+w2.psi+w3.psi+w4.psi)/4)*(pow((w1.h+w2.h)/2,2)+pow((w3.h+w4.h)/2,2))/4)
					*(sin(w3.x)-sin(w1.x))/dx;

	source.q = 0;//pressure+(bf.q-fr.q)*w.h*sin(w.x)-grad_b;

	source.r=0;//bf.r-fr.r;

	return source;
}

FS Friction (CV w, FS bf)

{
	FS fr;
	double mu=tan(delta* PI / 180);

	if (pow(pow(w.u,2)+pow(w.v,2),0.5)>1e-4) 
	{
		//cout<<"Inside the if of q"<<endl;
		fr.q=(mu*w.u/pow(pow(w.u,2)+pow(w.v,2),0.5))*w.psi*(1+3*epsilon*w.b);
	}
	else 
	{
		fr.q=(bf.q>0?1:-1)*min(abs(bf.q),mu*w.psi*(1+3*epsilon*w.b));
		//fr.q=4e+3*omega*w.u*mu*w.psi*(1+3*epsilon*w.b);
		//fr.q=mu*w.psi*(1+3*epsilon*w.b);
	}

	if (pow(pow(w.u,2)+pow(w.v,2),0.5)>1e-4) 
	{	
		//cout<<"Inside the if of r"<<endl;
		fr.r=mu*(w.psi)*w.v/pow(pow(w.u,2)+pow(w.v,2),0.5)*(1+3*epsilon*w.b)*w.h*pow(sin(w.x),2);
	}
	else 
	{	
		fr.r=(bf.r>0?1:-1)*min(abs(bf.r),mu*w.psi*(1+3*epsilon*w.b)*w.h*pow(sin(w.x),2));
		//fr.r=4e+3*omega*w.v*mu*w.psi*(1+3*epsilon*w.b)*w.h*pow(sin(w.x),2);
		//fr.r=mu*w.psi*(1+3*epsilon*w.b)*w.h*pow(sin(w.x),2);
	}

	return fr;
}

FS Body_force (CV w, CV w1, CV w2, CV w3, CV w4,  double om, double alpha)
{
	FS bf;
	bf.q=(omega*omega*sin(w.x)*cos(w.x))*(1+4*epsilon*w.lambda+2*epsilon*om)+(2*omega*cos(w.x)*w.v+w.g.X2)*(1+3*epsilon*w.lambda+epsilon*om);

	FS hl= Hx(w1,w2);
	FS hr= Hx(w3,w4);
	bf.r=-omega*(1+epsilon*om)*(hl.p+hr.p)/2*((sin(w3.x)*sin(w3.x)-sin(w1.x)*sin(w1.x))/dx
				+2*epsilon*(w3.b*sin(w3.x)*sin(w3.x)-w1.b*sin(w1.x)*sin(w1.x))/dx)
				-epsilon*alpha*pow(sin(w.x),3)*w.h; 
	return bf;
}

FS Eigen(CV w)

{
	double root,base;
	FS e;
		
		e.p = w.u;
		e.q = w.u+sqrt(w.h);
		e.r = w.u-sqrt(w.h);

		return e;
	
}

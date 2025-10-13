#include "gauravlib.h"

CV::CV(double h, double u, double v, double b, double db, double ddb, Grav g, double x )
{	
	if (h<=min_h) 
	{  if (h<0) std::cout<<"Much smaller values encountered "<<h<<"  "<<x<< endl;
	h=min_h;
	u=min_u; v=min_u;
	}
    this->h=h; this->u=u; this->v=v; this->b=b; this->db=db; this->ddb=ddb;
	lambda= (1+Gamma*b);
	this->g=g;this->x=x;this->psi=Psi(*this); this->w=h+Gamma/epsilon*b;
    this->p=h*sin(x)*pow(lambda,2); 
	this->q=h*u*pow(lambda,3); 
	this->r=h*v*sin(x)*sin(x)*pow(lambda,3); 
}


void CV::Modify(double p, double q, double r)
{	h=p/sin(x)/pow(lambda,2);
	if (h>100) std::cout<<"Much bigger values encountered in Modify "<<h<<"  "<<x<<"   "<<p<<"    "<<q<<"    "<<r<<"   "<< endl;
	if (h<min_h)
	{  if (h<0) std::cout<<"Much smaller values encountered in Modify "<<h<<"  "<<x<<"   "<<p<<"    "<<q<<"    "<<r<<"   "<< endl;
	h = min_h; u = min_u; v = min_u;
	p=h*sin(x)*pow(lambda,2);
	q=h*u*pow(lambda,3);
	r=h*v*sin(x)*sin(x)*pow(lambda,3);
	}
	this->p=p; this->q=q; this->r=r; w=h+Gamma/epsilon*b;
	
	double u_temp=q*sin(x)/(lambda*p);

	//if (u_temp*u < 0) 
	//{
	//	u= u_temp/abs(u_temp) * min_u;
	//	this->q=h*u*pow(lambda,3);
	//}
	//else
		u = u_temp;
	
	 double v_temp =r/(p*lambda*sin(x));
	 
	// if (v_temp*v < 0) 
	// {
	//	 v = v_temp/abs(v_temp) * min_u;
	//	 this->r=h*v*sin(x)*sin(x)*pow(lambda,3); 
	// }

 	//else
		 v = v_temp;


	psi=Psi(*this);
}


FS Flux( CV w )

{
	FS f;
	 f.p =  w.u*w.h*w.lambda*sin(w.x); 
	 f.q = (w.u*w.u*pow(w.lambda,2)+epsilon*w.psi*w.h/2)*w.h;  
	 f.r = w.u*w.v*w.h*pow(w.lambda,2)*sin(w.x)*sin(w.x); 
	//	f.p=-w.q*(2*w.b*epsilon*sin(w.x)+epsilon*w.p-sin(w.x))/sin(w.x);
	//	f.q=-(4*w.b*epsilon*sin(w.x)*w.q*w.q+2*epsilon*w.p*w.q*w.q-2*sin(w.x)*w.q*w.q-epsilon*pow(w.p,3)*w.psi)/(2*sin(w.x)*w.p);
	//	f.r=-w.q*w.r*(2*w.b*epsilon*sin(w.x)+epsilon*w.p-sin(w.x))/(sin(w.x)*w.p);
	return f;

}


FS Source( CV w, CV w1, CV w2, CV w3, CV w4) 
{
	FS source;
	FS bf=Body_force( w,  w1,  w2, w3,  w4);
	FS fr=Friction(w,bf);
	
	source.p=0;

	double grad_b = Gamma*w.db*(w.u*w.u+w.v*w.v) +(w.v*w.v+w.u*w.u)/w.lambda*cos(w.x)/sin(w.x);


	source.q = (bf.q-fr.q+grad_b)*pow(w.lambda,3)*w.h;
	source.r = (bf.r-fr.r)*pow(w.lambda,3)*w.h*pow(sin(w.x),2);

	return source;
}

FS Friction (CV w, FS bf)

{
	FS fr;
	double mu=tan(delta* PI / 180);

	if (pow(pow(w.u,2)+pow(w.v,2),0.5)>min_u*dx) 
	{
	
		fr.q=(mu*w.u/pow(pow(w.u,2)+pow(w.v,2),0.5))*w.psi;
	}
	else 
	{	
		fr.q=(bf.q>0?1:-1)*min(abs(bf.q),mu*w.psi);
		//fr.q=4e+3*omega*w.u*mu*w.psi*(1+3*epsilon*w.b);
		//fr.q=mu*w.psi*(1+3*epsilon*w.b);
	}

	if (pow(pow(w.u,2)+pow(w.v,2),0.5)>min_u*dx) 
	{	
		//std::cout<<"Inside the if of r"<<endl;
		fr.r=mu*(w.psi)*w.v/pow(pow(w.u,2)+pow(w.v,2),0.5);
	}
	else 
	{	
		fr.r=(bf.r>0?1:-1)*min(abs(bf.r),mu*w.psi);
		//fr.r=4e+3*omega*w.v*mu*w.psi*(1+3*epsilon*w.b)*w.h*pow(sin(w.x),2);
		//fr.r=mu*w.psi*(1+3*epsilon*w.b)*w.h*pow(sin(w.x),2);
	}

	return fr;
}

FS Body_force (CV w, CV w1, CV w2, CV w3, CV w4)
{
	FS bf;
	bf.q=(omega*omega*sin(w.x)*cos(w.x))*w.lambda + (2*omega*cos(w.x)*w.v+w.g.X2)+Gamma*(w.db)*omega*sin(w.x)*(omega*sin(w.x)+2*w.v);
	bf.r  = -2*omega*w.u*(cos(w.x)+Gamma*(w.db)*sin(w.x));
	return bf;
}

FS Eigen(CV w)
{
	double root,base;
	FS e;
		root =sqrt(epsilon*std::max(w.psi,0.0)*w.p/sin(w.x))/pow(w.lambda,3);
		base= w.q*sin(w.x)/(w.p*pow(w.lambda,2));
		e.p = base;
		e.q = (base+root);
		e.r = (base-root);

		return e;
}


double Psi(CV w)
{
return	-(omega*omega*w.lambda*sin(w.x)*sin(w.x)+2*omega*w.v*sin(w.x)+w.g.X1+(w.u*w.u+w.v*w.v)/w.lambda-Gamma*w.db*cos(w.x)*(omega*omega*sin(w.x)+2*w.v*omega+w.v*w.v/sin(w.x))
		-Gamma*w.ddb*w.u*w.u);
}
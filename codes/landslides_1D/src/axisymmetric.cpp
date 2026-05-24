#include "gauravlib.h"

CV::CV(double h, double u, double v, double b, double db, double ddb, Grav g, double x )
{	
 	if (h<=min_h) 
	{  if (h<0) std::cout<<"Much smaller values encountered "<<h<<"  "<<x<< endl;
	h=min_h; u=min_u; v=min_u;
	} 

    this->h=h; this->u=u; this->v=v; this->b=b;
	this->db=db; this->ddb=ddb;
	
    metric = sqrt(b*b +db*db);
    theta = (b*b + 2*db*db - b*ddb)/pow(metric,3);

	phi =  1/(b*sin(x));
    phi_norm = (b*sin(x)-db*cos(x))/(b*metric*sin(x));
		if (epsilon*theta*h<-1)
		cout<<"Very large value of theta "<< theta<< " at " << x<< endl;
    phi_tan = (b*cos(x)+ db*sin(x))/(b*metric*sin(x));
	J =  metric/phi; 
	V = v + omega/phi;
	this->g= g;
	this->x=x;
	psi=Psi(*this);
    this->p= J*h; 
	this->q= J*u*h; 
	this->r= J*V*h/phi; 
	this->P= h; 
	this->Q= u*h; 
	this->R= V*h; 
}


void CV::Modify(double p, double q, double r)
{	
    h=p/J;//(-J_basal + sqrt(2*J_basal*epsilon*p*theta + pow(J_basal,2)))/(J_basal*epsilon*theta);
	 if (h>100) std::cout<<"Much bigger values encountered in Modify "<<h<<"  "<<x<<"   "<<p<<"    "<<q<<"    "<<r<<"   "<< endl;
	if (h<min_h)
	{  
    std::cout<<"Much smaller values encountered in Modify "<<h<<"  "<<x<<"   "<<p<<"    "<<q<<"    "<<r<<"   "<< endl;
	h = min_h; u = sign(q)*min_u; v = sign(r*phi/p-omega/phi)*min_u; V = v+omega/phi;

	p=J*h;
	q=J*u*h; 
	r=J*V*h/phi; 
    } 


	this->p=p; this->q=q; this->r=r; 
	
	u = q/p;
	V = r*phi/p;
	v = V - omega/phi;
	psi=Psi(*this);
	P = h; Q =h*u; R = h*V;

}

void CV::Modify_U(double P, double Q, double R)
{	
    h=P;
	this->P=P; this->Q=Q; this->R=R; 
	
	u = Q/P;
	V = R/P;
	v = V-omega/phi;
	p = J*P; q= J*Q; r= J*R/phi;
	psi = Psi(*this);

}


FS Flux( CV w )

{
	FS f;
	 f.p =  w.u*w.h/w.phi; 
	 f.q = (w.u*w.u + epsilon*w.psi*w.h/2)*w.h/w.phi;
	 f.r = w.u*(w.V)*w.h/(w.phi*w.phi); 
	return f;
}


FS Source( CV w, CV w1, CV w2, CV w3, CV w4) 
{
	FS source;
	FS bf=Body_force( w, w1, w2, w3, w4);
	FS fr=Friction(w);
	source.p = 0;
	double pressure = (pow(w.V,2) + epsilon*w.psi*w.h/2)*(1/w3.phi+1/w4.phi-1/w1.phi-1/w2.phi)/(2*dx) +
						epsilon*w.h/2*(w3.g.X1+w4.g.X1-w1.g.X1-w2.g.X1)/(2*dx*w.phi); 
	source.q = (w.g.X2*w.J+pressure)*w.h ;// - epsilon*(w3.g.X1+w4.g.X1-w1.g.X1-w2.g.X1)/(2*dx*w.phi)*w.h*w.h/2;
	source.r = 0;
	return source;
}

FS Friction (CV w)

{
	FS fr;
	double mu=tan(delta* PI / 180);

	if (pow(w.u,2)+pow(w.v,2)>0)

	{
		
		fr.q=(mu*w.u/pow(pow(w.u,2)+pow(w.v,2),0.5))*w.psi*w.J*w.h;
	
		fr.r=mu*w.v/pow(pow(w.u,2)+pow(w.v,2),0.5)*w.psi*w.J*w.h/w.phi;
	}
	return fr;
}

FS Body_force (CV w, CV w1, CV w2, CV w3, CV w4)
{
	FS bf;
	bf.q= (omega/w.phi*(omega/w.phi+2*w.v)*(1/w4.phi+1/w3.phi-1/w1.phi-1/w2.phi)/(2*dx) +w.g.X2*w.J);
	FS hl = Hx(w1,w2);
	FS hr = Hx(w3,w4);
	bf.r  = 0;//((w4.u*w4.h/w4.phi)+(w3.u*w3.h/w3.phi)-(w2.u*w2.h/w2.phi)-(w1.u*w1.h/w1.phi))/(2*dx);
	return bf;
}

FS Eigen(CV w)
{
double root,base;
	FS e;
	 double b = epsilon*w.h*w.theta;
		root = sqrt((4*w.psi +epsilon*w.h*pow(w.theta,2))*epsilon*w.h); 
			if (root<0)
			{
				root =0;
				cout <<"The root is turning out to be negative. Making it zero" << endl;
			}
		e.p = (w.u)/w.J/w.phi;
		e.q = (2*w.u-epsilon*w.h*w.theta*w.u+root)/(2*w.J*w.phi);
		e.r = (2*w.u-epsilon*w.h*w.theta*w.u-root)/(2*w.J*w.phi);

		return e;
}


double Psi(CV w)
{
	
return	-( w.u*w.u*w.theta+ pow(w.v+omega/w.phi,2)*w.phi_norm+ 
		w.g.X1);
}



double Ang_mom_reg (CV w)
{
	return  2*PI*epsilon*w.v/w.phi*w.p*(1);//+(w.theta/2+w.phi_norm)*epsilon*w.h+pow(epsilon*w.h,2)/3*w.phi_norm*(2*w.theta+w.phi_norm));

}
double Jinertia1_reg (CV w)
{
	return 2*PI*w.p/pow(w.phi,2)*epsilon*(1);//+epsilon/2*(3*w.phi_norm +w.theta)*w.h+pow(epsilon*w.h,2)*w.phi_norm*(w.theta+w.phi_norm));
}
double Jinertia1_ast (CV w)
{
	return 2*PI/5*pow(w.b,2)/pow(w.phi,3);
}
double Jinertia2_reg (CV w)
{
	return PI*pow(w.b,3)*w.metric*epsilon*w.h*sin(w.x)*(1+pow(cos(w.x),2));
}
double Jinertia2_ast (CV w)
{
	return PI/5*(-pow(sin(w.x),3)+2*sin(w.x))*pow(w.b,5);
}
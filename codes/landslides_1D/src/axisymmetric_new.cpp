#include "gauravlib.h"

CV::CV(double h, double u, double v, double b, double db, double ddb, Grav g, double x )
{	
 	if (h<=min_h) 
	{  
        if (h<0) std::cout<<"Much smaller values encountered "<<h<<"  "<<x<< endl;

	    h=min_h; u=min_u; v=min_u;
	} 

	if (abs(u)>10)
			cout <<"Much higher velocity"<< endl;

    this->h=h; this->u=u; this->v=v; 

    this->b=b; this->db=db; this->ddb=ddb;
	
    metric = sqrt(b*b +db*db);

    theta_b = (b*b + 2*db*db - b*ddb)/pow(metric,3);

	phi_b =  1/(b*sin(x));

    R_phi_b = 1/phi_b;

    J_b = metric*R_phi_b;

    phi_norm_b = (b*sin(x)-db*cos(x))/(b*metric*sin(x));

    phi_tan_b = (b*cos(x)+ db*sin(x))/(b*metric*sin(x));

    J = J_b*(1+epsilon*(theta_b+phi_norm_b)*h/2);

    R_phi = R_phi_b*(1+epsilon*phi_norm_b*h/2); 

    phi = 1/R_phi; 

    theta = theta_b*(1-epsilon*theta_b*h/2); 

    phi_norm = phi_norm_b*(1-epsilon*phi_norm_b*h/2);

    phi_tan = phi_tan_b*(1-epsilon*phi_norm_b*h/2);

	if (J<0)
		cout<<"Very large value of theta "<< theta<< " at " << x<<"  bcause of h value becoming "<<h<< endl;
    
	V = v + omega/phi;
	this->g= g;
	this->x= x;
    this->p= J*h; 
	this->q= J*J*phi*u*h; 
	this->r= J*R_phi*V*h; 
	this->P= h; 
	this->Q= u*h; 
	this->R= V*h; 
	psi=Psi(*this);
}


void CV::Modify(double p, double q, double r)
{	
    h=-(1 - sqrt(2*epsilon*p*(theta_b + phi_norm_b)/J_b + 1))/(epsilon*(theta_b + phi_norm_b));
	 if (h>100) 
	 std::cout<<"Much bigger values encountered in Modify "<<h<<"  "<<x<<"   "<<p<<"    "<<q<<"    "<<r<<"   "<< endl;
	if (h<0)
	{  
    std::cout<<"Much smaller values encountered in Modify "<<h<<"  "<<x<<"   "<<p<<"    "<<q<<"    "<<r<<"   "<< endl;
	h = min_h; u = sign(q)*min_u; v = sign(r*phi/p-omega*R_phi)*min_u; V = v+omega/phi;

	p=J_b*h;
	q=J_b*J_b*phi_b*u*h; 
	r=J_b*R_phi_b*V*h; 
    } 

    J = J_b*(1+epsilon*(theta_b+phi_norm_b)*h/2);

    R_phi = R_phi_b*(1+epsilon*phi_norm_b*h/2); 

    phi = 1/R_phi; 

    theta = theta_b*(1-epsilon*theta_b*h/2); 

    phi_norm = phi_norm_b*(1-epsilon*phi_norm_b*h/2);

    phi_tan = phi_tan_b*(1-epsilon*phi_norm_b*h/2);

	this->p=p; this->q=q; this->r=r; 
	
	u = q/p/J/phi;

	if (abs(u)>10)
		cout <<"Much higher velocity"<< endl;

	V = r*phi/p;
	v = V - omega*R_phi;
	P = h; Q =h*u; R = h*V;
	psi = Psi(*this);

}

void CV::Modify_U(double P, double Q, double R)
{	
    h=P;

	 if (h>100 || Q/P>100 || R/P>100) std::cout<<"Much bigger values encountered in Modify "<<x<<"   "<<P<<"    "<<Q<<"    "<<R<<"   "<< endl;
	if (h<min_h)
	{  
   // std::cout<<"Much smaller values encountered in Modify "<<h<<"  "<<x<<"   "<<p<<"    "<<q<<"    "<<r<<"   "<< endl;
	h = min_h; u = sign(q)*min_u; 
	v = sign(r*phi/p-omega*R_phi)*min_u; V = v+omega/phi;
	P = h;
	Q = u*h; 
	R = V*h; 
    }

	this->P=P; this->Q=Q; this->R=R; 
	
    J = J_b*(1+epsilon*(theta_b+phi_norm_b)*h/2);

    R_phi = R_phi_b*(1+epsilon*phi_norm_b*h/2); 

    phi = 1/R_phi; 

    theta = theta_b*(1-epsilon*theta_b*h/2); 

    phi_norm = phi_norm_b*(1-epsilon*phi_norm_b*h/2);

    phi_tan = phi_tan_b*(1-epsilon*phi_norm_b*h/2);

	u = Q/P;

	if (abs(u)>10)
		cout <<"Much higher velocity"<< endl;

	V = R/P;
	v = V-omega/phi;
	p = J*P; q= J*J*phi*Q; r= J*R_phi*R;
	psi = Psi(*this);


}


FS Flux( CV w )

{
	FS f;
	 f.p =  w.u*w.h*w.R_phi; 
	 f.q = (w.u*w.u + epsilon*w.psi*w.h/2)*w.h*w.J ;
	 f.r = w.u*(w.V)*w.h*(w.R_phi*w.R_phi);
	return f;
}


FS Source( CV w, CV w1, CV w2, CV w3, CV w4) 
{
	FS source;
	FS fr=Friction(w);
	source.p = 0;
	double term1 = ((pow(w.V,2) + epsilon*w.psi*w.h/2)*w.phi_tan + w.g.X2 + 
                    epsilon/(w.J*w.phi)*(w3.g.X1+w4.g.X1-w1.g.X1-w2.g.X1)/(2*dx)*w.h/2)*w.J*w.J*w.phi*w.h;
    
    double term2 = (pow(w.u,2) + epsilon*w.psi*w.h/2)*w.R_phi*w.h*((w3.J_b*w3.phi_b+w4.J_b*w4.phi_b
                - w2.J_b*w2.phi_b-w1.J_b*w1.phi_b)/(2*dx) +
                epsilon* (w3.J_b*w3.phi_b*w3.theta_b + w4.J_b*w4.phi_b*w4.theta_b
                - w2.J_b*w2.phi_b*w2.theta_b-w1.J_b*w1.phi_b*w1.theta_b)/(2*dx)*w.h/2 );

 
	source.q = term1 + term2;
	source.r = 0 ;
	return source;
}

FS Friction (CV w)

{
	FS fr;

	
		fr.q = (mu*w.u/pow(pow(w.u,2)+pow(w.v,2)+1e-16,0.5))*w.psi*w.J_b*w.J_b*w.phi_b*w.h;
	
		fr.r = mu*w.v/pow(pow(w.u,2)+pow(w.v,2)+1e-16,0.5)*w.psi*w.J_b*w.R_phi_b*w.h;
	
	return fr;
}

FS Eigen(CV w)
{
double root,base;
	FS e;
	 
		root = epsilon*pow(w.p,3)*w.J*pow(w.phi,2)*w.psi; 

			if (root<0)
			{
				root =0;
			//	cout <<"The root is turning out to be negative. Making it zero" << endl;
			}
		e.p = w.q/pow(w.J*w.phi,2)/w.p;
		e.q = (w.q*(1-epsilon*w.theta_b*w.h) + sqrt(root))/(pow(w.J,2)*w.p*pow(w.phi,2));
		e.r = (w.q*(1-epsilon*w.theta_b*w.h) - sqrt(root))/(pow(w.J,2)*w.p*pow(w.phi,2));

		return e;
}


double Psi(CV w)
{
	
	double psi =	-( w.u*w.u*w.theta+ pow(w.V,2)*w.phi_norm + w.g.X1);

	return std::max(0.0,psi);
}


double Ang_mom_reg (CV w)
{
	return  2*PI*epsilon*w.v/w.phi*w.p*(1);
}

double Jinertia1_reg (CV w)
{
	return 2*PI*w.p/pow(w.phi,2)*epsilon;
}

double Jinertia1_ast (CV w)
{
	return 2*PI/5*pow(w.b,2)/pow(w.phi_b,3);
}

double Jinertia2_reg (CV w)
{
	return PI*pow(w.b,3)*w.metric*epsilon*w.h*sin(w.x)*(1+pow(cos(w.x),2));
}

double Jinertia2_ast (CV w)
{
	return PI/5*(-pow(sin(w.x),3)+2*sin(w.x))*pow(w.b,5);
}
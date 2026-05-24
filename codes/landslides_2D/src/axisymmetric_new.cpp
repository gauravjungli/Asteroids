#include "gauravlib.h"

CV::CV(double h, double u, double v, double b, double Ra, double dR, double ddR, Grav g, double x, double y )
{	
 	if (h<=min_h) 
	{  
        if (h<0) std::cout<<"Much smaller values encountered "<<h<<"  "<<x<< endl;

	    h=min_h; u=min_u; v=min_u;
	} 

    this->h = h; this->b = b;
    
    this->u = u; this->v = v; 
 
    this->Ra = Ra; this->dR = dR; this->ddR = ddR;
	
    N_0 =  sqrt(Ra*Ra +dR*dR);

    theta_0 =  (Ra*Ra + 2*dR*dR - Ra*ddR)/pow(N_0,3);

	phi_0 =  1/(Ra*sin(x));

    R_phi_0 =  1/phi_0;

    J_0 = N_0*R_phi_0;

    phi_norm_0 =  (Ra*sin(x)-dR*cos(x))/(J_0);

    phi_tan_0 =  (Ra*cos(x)+ dR*sin(x))/(J_0);

    J_b = J_0*(1 + epsilon*(theta_0 + phi_norm_0)*b);

    N_b = N_0*(1 + epsilon*theta_0*b);

    R_phi_b = R_phi_0*(1 + epsilon*phi_norm_0*b); 

    phi_b = phi_0*(1 - epsilon*phi_norm_0*b); 

    theta_b = theta_0*(1 - epsilon*theta_0*b); 

    phi_norm_b = phi_norm_0*(1 - epsilon*phi_norm_0*b);

    phi_tan_b = phi_tan_0*(1 - epsilon*phi_norm_0*b);

    J = J_b*(1 + epsilon*(theta_0 + phi_norm_0)*h/2);

    N = N_b*(1 + epsilon*theta_0*h/2);

    R_phi = R_phi_b*(1 + epsilon*phi_norm_0*h/2); 

    phi = phi_b*(1 - epsilon*phi_norm_0*h/2); 

    theta = theta_b*(1 - epsilon*theta_0*h/2); 

    phi_norm = phi_norm_b*(1 - epsilon*phi_norm_0*h/2);

    phi_tan = phi_tan_b*(1 - epsilon*phi_norm_0*h/2);

	if (epsilon*theta*h<-1)
		cout<<"Very large value of theta "<< theta<< " at " << x<< endl;
    
	V = v + omega*R_phi;
	this->g= g;
	this->x=x;
    this->y=y;
	psi=Psi(*this);

    this->p= J*h; 
	this->q= J*N*u*h; 
	this->r= J*R_phi*V*h; 
	this->P= h; 
	this->Q= u*h; 
	this->R= V*h; 
}


void CV::Modify(double p, double q, double r)
{	
    h= -(J_b - sqrt(2*J_b*epsilon*p*theta_b + 2*J_b*epsilon*p*phi_norm_b + J_b*J_b))/(J_b*epsilon*(theta_b + phi_norm_b));
	 if (h>100) std::cout<<"Much bigger values encountered in Modify "<<h<<"  "<<x<<"   "<<p<<"    "<<q<<"    "<<r<<"   "<< endl;
	if (h<min_h)
	{  
   // std::cout<<"Much smaller values encountered in Modify "<<h<<"  "<<x<<"   "<<p<<"    "<<q<<"    "<<r<<"   "<< endl;
	h = min_h; u = sign(q)*min_u; v = sign(r/p/R_phi-omega*R_phi)*min_u; V = v+omega*R_phi;
	p=J_b*h;
	q=J_b*N_b*u*h; 
	r=J_b*R_phi_b*V*h; 
    } 

    J = J_b*(1+epsilon*(theta_0+phi_norm_0)*h/2);

    R_phi = R_phi_b*(1+epsilon*phi_norm_0*h/2); 

    phi = phi_b*(1-epsilon*phi_norm_0*h/2); 

    theta = theta_b*(1-epsilon*theta_0*h/2); 

    phi_norm = phi_norm_b*(1-epsilon*phi_norm_0*h/2);

    phi_tan = phi_tan_b*(1-epsilon*phi_norm_0*h/2);
	
	N = N_b*(1 + epsilon*theta_0*h/2);

	this->p=p; this->q=q; this->r=r; 
	
	u = q/p/N;
	V = r/p/R_phi;
	v = V - omega*R_phi;
	psi=Psi(*this);
	P = h; Q =h*u; R = h*V;

}

void CV::Modify_U(double P, double Q, double R)
{	
    h=P;

	 if (h>100 || Q/P>100 || R/P>100) std::cout<<"Much bigger values encountered in Modify "<<x<<"   "<<P<<"    "<<Q<<"    "<<R<<"   "<< endl;
	if (h<min_h)
	{  
    if (h<0) std::cout<<"Much smaller values encountered in Modify "<<h<<"  "<<x<<"   "<<p<<"    "<<q<<"    "<<r<<"   "<< endl;
	h = min_h; u = sign(q)*min_u; 
	v = sign(r/p/R_phi-omega*R_phi)*min_u; V = v + omega*R_phi;
	P = h;
	Q = u*h; 
	R = V*h; 
    }

	this->P=P; this->Q=Q; this->R=R; 
	
    J = J_b*(1+epsilon*(theta_0+phi_norm_0)*h/2); 

    R_phi = R_phi_b*(1+epsilon*phi_norm_0*h/2); 

    phi = phi_b*(1-epsilon*phi_norm_0*h/2); 

    theta = theta_b*(1-epsilon*theta_0*h/2); 

    phi_norm = phi_norm_b*(1-epsilon*phi_norm_0*h/2);

    phi_tan = phi_tan_b*(1-epsilon*phi_norm_0*h/2);

	N = N_b*(1 + epsilon*theta_0*h/2);

	u = Q/P;
	V = R/P;
	v = V-omega*R_phi;
	p = J*P; q= J*N*Q; r= J*R_phi*R;
	psi = Psi(*this);

}


FS Flux_x( CV w )

{
	FS f;
	 f.p =    w.u*w.h*w.R_phi ; 
	 f.q =   (w.u*w.u + epsilon*w.psi*w.h/2)*w.h*w.J ;
	 f.r =   w.u*w.V*w.h*pow(w.R_phi, 2) ;
	return f;
}


FS Flux_y( CV w )

{
	FS g;
	 g.p =  w.v*w.h*w.N ; 
	 g.q =  w.u*w.v*w.h*pow(w.N,2) ;
	 g.r =  (w.v*w.V + epsilon*w.psi*w.h/2)*w.h*w.J ;
	return g;
}


FS Source( CV w, CV w1, CV w2, CV w3, CV w4, CV w5, CV w6, CV w7, CV w8) 
{
	FS source;
	FS bf=Body_force( w, w1, w2, w3, w4);
	FS fr=Friction(w);
	source.p = 0;
	double term1 = (pow(w.V,2)*w.phi_tan + w.g.X2 + 
                    epsilon/(w.N)*(w3.g.X1+w4.g.X1-w1.g.X1-w2.g.X1)/(2*dx)*w.h/2)*w.J*w.N*w.h;
    
    double term2 = pow(w.u,2)*w.R_phi*w.h*((w3.N_b+w4.N_b
                - w2.N_b-w1.N_b)/(2*dx) +
                epsilon* (w3.N_b*w3.theta_b + w4.N_b*w4.theta_b
                - w2.N_b*w2.theta_b-w1.N_b*w1.theta_b)/(2*dx)*w.h/2 );
    
    double term3 = epsilon*w.psi*w.h*((w3.J + w4.J - w2.J -w1.J)/(2*dx)*w.h/2 - 
                    w.J_b*(w3.b + w4.b - w2.b -w1.b)/(2*dx));
 
	source.q =   term1 + term2 + term3;
	source.r =   -epsilon*w.psi*w.h*w.J_b*(w7.b + w8.b - w5.b -w6.b)/(2*dy);
	return source;
}

FS Friction (CV w)

{
	FS fr;

	if (pow(w.u,2)+pow(w.v,2)>0)

	{
		
		fr.q=(mu*w.u/pow(pow(w.u,2)+pow(w.v,2),0.5))*w.psi*w.J_b*w.N_b*w.h;
	
		fr.r=mu*w.v/pow(pow(w.u,2)+pow(w.v,2),0.5)*w.psi*w.J_b*w.R_phi_b*w.h;
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

FS Eigenx(CV w)
{
double root,base;
	FS e;

	root = sqrt(-8*w.p*epsilon*(pow(w.phi_b,2)*(pow(w.r*w.phi_b,2)*w.phi_norm_b + pow(w.p,2)*w.g.X1)*pow(w.J_b,3) 
        - w.p*epsilon*pow(w.phi_b,2)*w.theta_b*(pow(w.r*w.phi_b,2)*w.phi_norm_b + pow(w.p,2)*w.g.X1)*pow(w.J_b,2) 
        + w.J_b*pow(w.q,2)*w.theta_b - (9*epsilon*w.p*pow(w.q*w.theta_b,2))/8)); 

	if (root<0)
			{
				root =0;
				cout <<"The root is turning out to be negative. Making it zero" << endl;
			}
	e.p =  w.q*(-epsilon*w.p*w.theta_b + w.J_b)/(pow(w.J_b,3)*w.p*pow(w.phi_b,2));
	e.q =  (-5*epsilon*w.p*w.q*w.theta_b + 2*w.J_b*w.q + root)/(2*pow(w.J_b,3)*w.p*pow(w.phi_b,2));
	e.r =  (-5*epsilon*w.p*w.q*w.theta_b + 2*w.J_b*w.q - root)/(2*pow(w.J_b,3)*w.p*pow(w.phi_b,2));

	return e;
}

FS Eigeny(CV w)
{
double root,base;
	FS e;

	root = sqrt(-epsilon*w.p*(pow(w.J_b*w.r,2)*pow(w.phi_b,4)*w.phi_norm_b + pow(w.J_b*w.p*w.phi_b,2)*w.g.X1
			 + w.q*w.q*w.theta_b)*(-epsilon*w.p*w.phi_norm_b + w.J_b));
	if (root<0)
			{
				root =0;
				cout <<"The root is turning out to be negative. Making it zero" << endl;
			}
	e.p =  (w.r*(-epsilon*w.p*w.phi_norm_b + w.J_b)*pow(w.phi_b,2)-w.J_b*omega*w.p)/(w.J_b*w.p);
	e.q = ((w.J_b*w.r*(-2*epsilon*w.p*w.phi_norm_b + w.J_b)*pow(w.phi_b,2)-w.J_b*w.J_b*omega*w.p) + root)/(w.J_b*w.J_b*w.p);
	e.r = ((w.J_b*w.r*(-2*epsilon*w.p*w.phi_norm_b + w.J_b)*pow(w.phi_b,2)-w.J_b*w.J_b*omega*w.p) - root)/(w.J_b*w.J_b*w.p);

	return e;
}


double Psi(CV w)
{
	
return	-( w.u*w.u*w.theta+ pow(w.V,2)*w.phi_norm+ w.g.X1);
}



double Ang_mom_reg (CV w)
{
	return  epsilon*w.v*w.R_phi*w.p*(1);//+(w.theta/2+w.phi_norm)*epsilon*w.h+pow(epsilon*w.h,2)/3*w.phi_norm*(2*w.theta+w.phi_norm));

}
double Jinertia1_reg (CV w)
{
	return w.p*pow(w.R_phi,2)*epsilon*(1);//+epsilon/2*(3*w.phi_norm +w.theta)*w.h+pow(epsilon*w.h,2)*w.phi_norm*(w.theta+w.phi_norm));
}
double Jinertia1_ast (CV w)
{
	return 2*PI/5*pow(w.R,2)*pow(w.R_phi,3) +w.J_b*w.b*pow(w.R_phi,2)*epsilon;
}
double Jinertia2_reg (CV w)
{
	return PI*pow(w.b,3)*w.N_b*epsilon*w.h*sin(w.x)*(1+pow(cos(w.x),2));
}
double Jinertia2_ast (CV w)
{
	return PI/5*(-pow(sin(w.x),3)+2*sin(w.x))*pow(w.b,5);
}
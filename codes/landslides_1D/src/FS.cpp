#include "gauravlib.h"

FS FS::operator+ (FS w){
    FS temp;
    temp.p=this->p+w.p;
    temp.q=this->q+w.q;
	temp.r=this->r+w.r;
    return temp;
  }

  FS FS::operator- (FS w){
    FS temp;
    temp.p=this->p-w.p;
    temp.q=this->q-w.q;
	temp.r=this->r-w.r;
    return temp;
  }

  FS FS::operator* (double w){
    FS temp;
    temp.p=w*this->p;
    temp.q=w*this->q;
	temp.r=w*this->r;
    return temp;
  }
  FS FS::operator/ (double w){
    FS temp;
    temp.p=this->p/w;
    temp.q=this->q/w;
	temp.r=this->r/w;
    return temp;
  }
  FS::FS(CV w)
  {
	this->p=w.p; this->q=w.q; this->r=w.r;
  }

  Grav Grav::operator+ (Grav w){
    Grav temp;
    temp.X1=this->X1+w.X1;
    temp.X2=this->X2+w.X2;
	  temp.X3=this->X3+w.X3;
    return temp;
  }

  Grav Grav::operator/ (double w){
    Grav temp;
    temp.X1=this->X1/w;
    temp.X2=this->X2/w;
	  temp.X3=this->X3/w;
    return temp;
  }

/*   CV::CV(const CV& temp)
  {
    this->w =temp.w; this->p =temp.p; this->q =temp.q;
    this->r =temp.r; this->u =temp.u; this->v =temp.v; 
    this->h =temp.h; this->b =temp.b; this->db =temp.db; 
    this->ddb =temp.ddb; this->x =temp.x; this->psi =temp.psi;  this->psi_basal =temp.psi_basal;
    this->metric =temp.metric; this->theta =temp.theta; this->phi_norm =temp.phi_norm; 
    this->phi_tan =temp.phi_tan; this->phi =temp.phi; this->J =temp.J; this->J_basal =temp.J_basal;
  } */

/* CV CV::operator= (const  temp&){
    Grav temp;
    temp.X1=this->X1/w;
    temp.X2=this->X2/w;
	  temp.X3=this->X3/w;
    return temp;
  } */

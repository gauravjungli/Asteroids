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

 /*  double sin(double x)
  {
    return 1;
  } */
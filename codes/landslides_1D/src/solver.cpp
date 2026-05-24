#include "gauravlib.h"

FS Hx( CV wl, CV wr)
{

	FS fr = Flux(wr); 
	FS fl = Flux(wl);
	double el=Ax(wl,wr,"min");
	double er=Ax(wl,wr,"max");

	
	FS Wr(wr),Wl(wl);
	FS w;
	 
	double a_plus = max(er,0.0);
	double a_minus = min(el,0.0);

	if ((a_plus-a_minus)>0)

	w = (fl*a_plus - fr*a_minus + (Wr-Wl)*a_plus*a_minus)/(a_plus-a_minus);

	else

	w=(fl+fr)/2-(Wr-Wl)*max(abs(el),abs(er))/2; 

	
	return w;
}




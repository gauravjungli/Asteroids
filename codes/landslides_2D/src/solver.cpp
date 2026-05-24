#include "gauravlib.h"

FS Hx( CV wl, CV wr)
{

	FS fr = Flux_x(wr); 
	FS fl = Flux_x(wl);

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

FS Hy( CV wb, CV wt)
{

	FS ft = Flux_y(wt); 
	FS fb = Flux_y(wb);
	double eb = Ay(wb,wt,"min");
	double et = Ay(wb,wt,"max");

	FS Wb(wb),Wt(wt);
	FS w;
	 
	double a_plus = max(et,0.0); 
	double a_minus = min(eb,0.0);

	if ((a_plus-a_minus)>0)

		w = (fb*a_plus - ft*a_minus + (Wt-Wb)*a_plus*a_minus)/(a_plus-a_minus);

	else

		w = (ft+fb)/2-(Wt-Wb)*max(abs(et),abs(eb))/2; 

	return w;
}




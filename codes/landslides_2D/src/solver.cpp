#include "gauravlib.h"

FS Hx( CV wl, CV wr)
{

	FS fr = Flux_x(wr); 
	FS fl = Flux_x(wl);
	double el=Ax(wl,wr,"min");
	double er=Ax(wl,wr,"max");
	FS delta_w(wl);
	
	FS Wr(wr),Wl(wl);
	FS W_star= (Wr*er-Wl*el-(fr-fl))/(er-el);
	delta_w=Minmod(Wr-W_star,W_star-Wl);

	FS w;
	  
	//	w= (fl*er-fr*el)/(er-el)+(Wr-Wl-delta_w)*(er*el)/(er-el);
		w=(fl+fr)/2-(Wr-Wl)*max(abs(el),abs(er))/2; 
	
	return w;
}

FS Hy( CV wb, CV wt)
{

	FS ft = Flux_y(wt); 
	FS fb = Flux_y(wb);
	double el=Ay(wb,wt,"min");
	double er=Ay(wb,wt,"max");
	FS delta_w(wb);
	
	FS Wb(wb),Wt(wt);
	//FS W_star= (Wr*er-Wl*el-(fr-fl))/(er-el);
	//delta_w=Minmod(Wr-W_star,W_star-Wl);

	FS w;
	  
	//	w= (fl*er-fr*el)/(er-el)+(Wr-Wl-delta_w)*(er*el)/(er-el);
		w=(fb+ft)/2-(Wt-Wb)*max(abs(el),abs(er))/2; 
	
	return w;
}




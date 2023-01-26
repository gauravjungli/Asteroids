#include "gauravlib.h"

FS Hx( CV wl, CV wr, double om)
{

	FS fr = Flux(wr,om); 
	FS fl = Flux(wl,om);
	double el=Ax(wl,wr,"min");
	double er=Ax(wl,wr,"max");
	FS delta_w(wl);
	
	FS Wr(wr),Wl(wl);
	FS W_star= (Wr*er-Wl*el-(fr-fl))/(er-el);
	delta_w=Minmod(Wr-W_star,W_star-Wl);

	FS w;
	  
		w= (fl*er-fr*el)/(er-el)+(Wr-Wl-delta_w)*(er*el)/(er-el);
	return w;
}



/*FS Hx( CV wl, CV wr, double om)
{

	FS fp = Flux(wr,om); 
	FS fm = Flux(wl,om);

	FS w;
	  
		w.p= ((fp.p + fm.p) / 2 - Ax(wl,wr) * (wr.p - wl.p) / 2);
		w.q= ((fp.q + fm.q) / 2 - Ax(wl,wr) * (wr.q - wl.q) / 2);
	 	w.r= ((fp.r + fm.r) / 2 - Ax(wl,wr) * (wr.r - wl.r) / 2);
	return w;
}*/


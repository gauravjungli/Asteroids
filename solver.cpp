#include "gauravlib.h"

FS Hx( CV wl, CV wr)
{

	FS fp = Flux(wr); 
	FS fm = Flux(wl);

	FS w;
	  
		w.p= ((fp.p + fm.p) / 2 - Ax(wl,wr) * (wr.p - wl.p) / 2);
		w.q= ((fp.q + fm.q) / 2 - Ax(wl,wr) * (wr.q - wl.q) / 2);
	 	w.r= ((fp.r + fm.r) / 2 - Ax(wl,wr) * (wr.r - wl.r) / 2);
	return w;
}


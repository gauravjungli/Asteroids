#include "gauravlib.h"

FS Hx( CV wl, CV wr)
{

	FS fp = flux(wr); 
	FS fm = flux(wl);

	FS w;
	  
		w.p= ((fp.p + fm.p) / 2 - ax(wl,wr) * (wr.p - wl.p) / 2);
		w.q= ((fp.q + fm.q) / 2 - ax(wl,wr) * (wr.q - wl.q) / 2);
	 	w.r= ((fp.r + fm.r) / 2 - ax(wl,wr) * (wr.r - wl.r) / 2);
	return w;
}


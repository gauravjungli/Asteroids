#include "gauravlib.h"

CV Hx( CV wl, CV wr,double x) // this xx, yy corresponds to boundaries of cells
// no is 1,2,3 for denoting p,q,r
{


	CV fp = flux(wr,x); // first 1 is for f and second 1 is for plus
	CV fm = flux(wl,x);

	CV w;
	  
		w.modify( ((fp.p + fm.p) / 2 - ax(wl,wr,x) * (wr.p - wl.p) / 2),
		 ((fp.q + fm.q) / 2 - ax(wl,wr,x) * (wr.q - wl.q) / 2),
	 ((fp.r + fm.r) / 2 - ax(wl,wr,x) * (wr.r - wl.r) / 2));
	return w;
}




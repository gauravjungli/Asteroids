#include "gauravlib.h"

void BC(vector<CV>& w, double om )
{
	//wall boundary conditions at pole
	w[0].Modify(w[3].p,-w[3].q,w[3].r,om);
	w[1].Modify(w[2].p,-w[2].q,w[2].r,om);
	w[res-2].Modify(w[res-3].p,-w[res-3].q,w[res-3].r,om);
	w[res-1].Modify(w[res-4].p,-w[res-4].q,w[res-4].r,om);
	
}
#include "gauravlib.h"

void BC(vector<CV>& w, double om )
{
	//wall boundary conditions at pole
	w[0]=CV(w[3].h,-w[3].u,w[3].v,w[0].b,w[0].grad_b,w[0].g,w[0].x,om);
	w[1]=CV(w[2].h,-w[2].u,w[2].v,w[1].b,w[1].grad_b,w[1].g,w[1].x,om);
	w[res-2]=CV(w[res-3].h,-w[res-3].u,w[res-3].v,w[res-2].b,w[res-2].grad_b,w[res-2].g,w[res-2].x,om);
	w[res-1]=CV(w[res-4].h,-w[res-4].u,w[res-4].v,w[res-1].b,w[res-1].grad_b,w[res-1].g,w[res-1].x,om);
	
}
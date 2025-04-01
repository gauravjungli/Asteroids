#include "gauravlib.h"


//Currently working best
 void BC(vector<CV>& w )
{
	//wall boundary conditions at pole
	for (int j = 0; j < cols; ++j) {
		int i0 = index(0,j);
		int i1 = index(1,j);
		int i2 = index(2,j);
		int i3 = index(3,j);
		
		int i_1 = index(rows-1,j);
		int i_2 = index(rows-2,j);
		int i_3 = index(rows-3,j);
		int i_4 = index(rows-4,j);

		w[i0] = CV(w[i3].h,-w[i3].u,w[i3].v,w[i3].b,w[i3].g,w[i0].x, w[i0].y);

		w[i1] = CV(w[i2].h,-w[i2].u,w[i2].v,w[i2].b,w[i2].g,w[i1].x, w[i1].y);

		w[i_2] = CV(w[i_3].h,-w[i_3].u,w[i_3].v,w[i_3].b,w[i_3].g,w[i_2].x, w[i_2].y);

		w[i_1] = CV(w[i_4].h,-w[i_4].u,w[i_4].v,w[i_4].b,w[i_4].g,w[i_1].x, w[i_1].y);
	}
}

 /*void BC(vector<CV>& w )
{
	//wall boundary conditions at pole
	w[0]=w[3];w[0].q=-w[3].q;w[0].u=-w[3].u;
	w[1]=w[2];w[1].q=-w[2].q;w[1].u=-w[2].u;
	//w[0].Modify(w[3].p,-w[3].q,w[3].r); 
	//w[1].Modify(w[2].p,-w[2].q,w[2].r); 
	w[res-2]=w[res-3];w[res-2].q=-w[res-3].q;w[res-2].u=-w[res-3].u;
	w[res-1]=w[res-4];w[res-1].q=-w[res-4].q;w[res-1].u=-w[res-4].u;
	//w[res-2].Modify(w[res-3].p,-w[res-3].q,w[res-3].r);
//w[res-1].Modify(w[res-4].p,-w[res-4].q,w[res-4].r);
	
} */

/*void BC(vector<CV>& w )
{
	//wall boundary conditions at pole
	
	w[1]=CV(w[2].h-theta*(w[3].h-w[2].h), -w[2].u,w[2].v, w[2].b-theta*(w[3].b-w[2].b), w[1].g, w[1].x);
	w[0]=CV(w[1].h-theta*(w[2].h-w[1].h), -w[3].u,w[3].v, w[1].b-theta*(w[2].b-w[1].b), w[0].g, w[0].x);
	w[res-2]=CV(w[res-3].h+theta*(w[res-3].h-w[res-4].h), -w[res-3].u,w[res-3].v, w[res-3].b+theta*(w[res-3].b-w[res-4].b), w[res-2].g, w[res-2].x);
	w[res-1]=CV(w[res-2].h+theta*(w[res-2].h-w[res-3].h), -w[res-4].u,w[res-4].v, w[res-2].b+theta*(w[res-2].b-w[res-3].b), w[res-1].g, w[res-1].x);
	
}*/
#include "gauravlib.h"


//Currently working best
 void BC(vector<CV>& w1 ,vector<CV>& w2)
{
	//wall boundary conditions at pole


	w1[0]=CV(w2[3].h,w2[3].u,-w2[3].u_c, w2[3].v, w2[3].b,w2[3].db,w2[3].ddb, w2[3].g, w2[0].x);

	w1[1]=CV(w2[2].h,w2[2].u,-w2[2].u_c, w2[2].v,w2[2].b,w2[2].db,w2[2].ddb, w2[2].g, w2[1].x);

	w1[res-2]=CV(w2[res-3].h,w2[res-3].u, -w2[res-3].u_c, w2[res-3].v,w2[res-3].b,w2[res-3].db,w2[res-3].ddb,w2[res-3].g,w2[res-3].x);

	w1[res-1]=CV(w2[res-4].h,w2[res-4].u, -w2[res-4].u_c,w2[res-4].v,w2[res-4].b,w2[res-4].db,w2[res-4].ddb,w2[res-4].g,w2[res-4].x);
	
}

 /*void BC(vector<CV>& w )
{
	//wall boundary conditions at pole
	w2[0]=w2[3];w2[0].q=-w2[3].q;w2[0].u=-w2[3].u;
	w2[1]=w2[2];w2[1].q=-w2[2].q;w2[1].u=-w2[2].u;
	//w2[0].Modify(w2[3].p,-w2[3].q,w2[3].r); 
	//w2[1].Modify(w2[2].p,-w2[2].q,w2[2].r); 
	w2[res-2]=w2[res-3];w2[res-2].q=-w2[res-3].q;w2[res-2].u=-w2[res-3].u;
	w2[res-1]=w2[res-4];w2[res-1].q=-w2[res-4].q;w2[res-1].u=-w2[res-4].u;
	//w2[res-2].Modify(w2[res-3].p,-w2[res-3].q,w2[res-3].r);
//w2[res-1].Modify(w2[res-4].p,-w2[res-4].q,w2[res-4].r);
	
} */

/*void BC(vector<CV>& w )
{
	//wall boundary conditions at pole
	
	w2[1]=CV(w2[2].h-theta*(w2[3].h-w2[2].h), -w2[2].u,w2[2].v, w2[2].b-theta*(w2[3].b-w2[2].b), w2[1].g, w2[1].x);
	w2[0]=CV(w2[1].h-theta*(w2[2].h-w2[1].h), -w2[3].u,w2[3].v, w2[1].b-theta*(w2[2].b-w2[1].b), w2[0].g, w2[0].x);
	w2[res-2]=CV(w2[res-3].h+theta*(w2[res-3].h-w2[res-4].h), -w2[res-3].u,w2[res-3].v, w2[res-3].b+theta*(w2[res-3].b-w2[res-4].b), w2[res-2].g, w2[res-2].x);
	w2[res-1]=CV(w2[res-2].h+theta*(w2[res-2].h-w2[res-3].h), -w2[res-4].u,w2[res-4].v, w2[res-2].b+theta*(w2[res-2].b-w2[res-3].b), w2[res-1].g, w2[res-1].x);
	
}*/
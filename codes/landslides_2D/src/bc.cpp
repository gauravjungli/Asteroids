#include "gauravlib.h"


//Currently working best
 void BC(vector<CV>& w1 ,vector<CV>& w2)
{
	//wall boundary conditions at pole
	for (int j = 0; j < cols; ++j) 
	{
		int i0 = index(0,j);
		int i1 = index(1,j);
		int i2 = index(2,j);
		int i3 = index(3,j);
		
		int i_1 = index(rows-1,j);
		int i_2 = index(rows-2,j);
		int i_3 = index(rows-3,j);
		int i_4 = index(rows-4,j);

		w1[i0] = CV(w2[i3].h, -w2[i3].u,  w2[i3].v, w2[i3].b, w2[i3].Ra, w2[i3].dR, w2[i3].ddR, w2[i3].g, w2[i3].x,w2[i3].y);

		w1[i1] = CV(w2[i2].h, -w2[i2].u,  w2[i2].v, w2[i2].b, w2[i2].Ra, w2[i2].dR, w2[i2].ddR, w2[i2].g, w2[i2].x,w2[i2].y);

		w1[i_2] = CV(w2[i_3].h, -w2[i_3].u, w2[i_3].v, w2[i_3].b, w2[i_3].Ra, w2[i_3].dR, w2[i_3].ddR, w2[i_3].g, w2[i_3].x,w2[i_3].y);

		w1[i_1] = CV(w2[i_4].h, -w2[i_4].u, w2[i_4].v, w2[i_4].b, w2[i_4].Ra, w2[i_4].dR, w2[i_4].ddR, w2[i_4].g, w2[i_4].x,w2[i_4].y);
	}
}

//Currently working best
 void BC(vector<CV>& w )
{
	//wall boundary conditions at pole
	for (int j = 0; j < cols; ++j) 
	{
		int i0 = index(0,j);
		int i1 = index(1,j);
		int i2 = index(2,j);
		int i3 = index(3,j);
		
		int i_1 = index(rows-1,j);
		int i_2 = index(rows-2,j);
		int i_3 = index(rows-3,j);
		int i_4 = index(rows-4,j);

	w[i0]= CV(w[i3].h, -w[i3].u, w[i3].v, w[i3].b, w[i3].Ra, w[i3].dR, w[i3].ddR, w[i3].g, w[i0].x,w[i0].y);

	w[i1]= CV(w[i2].h, -w[i2].u, w[i2].v, w[i2].b, w[i2].Ra, w[i2].dR, w[i2].ddR, w[i2].g, w[i1].x,w[i1].y);

	w[i_2]= CV(w[i_3].h, -w[i_3].u, w[i_3].v, w[i_3].b, w[i_3].Ra, w[i_3].dR, w[i_3].ddR, w[i_3].g, w[i_2].x,w[i_2].y);

	w[i_1]= CV(w[i_4].h, -w[i_4].u, w[i_4].v, w[i_4].b, w[i_4].Ra, w[i_4].dR, w[i_4].ddR, w[i_4].g, w[i_1].x,w[i_1].y);
	}
}


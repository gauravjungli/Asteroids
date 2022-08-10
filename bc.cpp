#include "gauravlib.h"

void bc(vector<CV>& wtemp )
{           double factor=1;
	//wall boundary conditions at pole
	wtemp[0].p = wtemp[3].p * factor;wtemp[0].q = -wtemp[3].q;wtemp[0].r = wtemp[3].r;
	wtemp[1].p = wtemp[2].p * factor;wtemp[1].q = -wtemp[2].q;wtemp[1].r = wtemp[2].r;

	//wall boundary conditions at equator
	wtemp[res-2].p = wtemp[res-3].p;wtemp[res-2].q = -wtemp[res-3].q;wtemp[res-2].r = wtemp[res-3].r;
	wtemp[res-1].p = wtemp[res-4].p;    wtemp[res-1].q = -wtemp[res-4].q;    wtemp[res-1].r = wtemp[res-4].r;

    //ghost cell update 
    wtemp[0].p = wtemp[3].p;wtemp[0].q = -wtemp[3].q;wtemp[0].r = wtemp[3].r;
    wtemp[0].h = wtemp[3].h;wtemp[0].u = -wtemp[3].u;wtemp[0].v = (1 + 2 * theta) * wtemp[2].v - 2 * theta * wtemp[3].v;
	wtemp[1].p = wtemp[2].p;wtemp[1].q = -wtemp[2].q;wtemp[1].r = wtemp[2].r;
    wtemp[1].h = wtemp[2].h;wtemp[1].u = -wtemp[2].u;wtemp[1].v = (1 + theta) * wtemp[2].v - theta * wtemp[3].v;

	wtemp[res -2].p = wtemp[res-3].p;wtemp[res-2].q = -wtemp[res-3].q;wtemp[res-2].r = wtemp[res-3].r;
    wtemp[res -2].h = wtemp[res-3].h;wtemp[res-2].u = -wtemp[res-3].u;wtemp[res-2].v = (1 + theta) * wtemp[res-3].v - theta * wtemp[res].v;
	wtemp[res -1].p = wtemp[res-4].p;wtemp[res-1].q = -wtemp[res-4].q;wtemp[res-1].r = wtemp[res-4].r;
    wtemp[res -1].h = wtemp[res-4].h;wtemp[res-1].u = -wtemp[res-4].u;wtemp[res-1].v = (1 + 2 * theta) * wtemp[res-3].v - 2 * theta * wtemp[res-4].v;
}
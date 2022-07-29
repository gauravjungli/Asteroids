#include "gauravlib.h"
void grav_sph (vector<grav>& v )
{
    for(int i=0;i<v.size();i++)
    {
        v[i].X=0;v[i].Y=0;v[i].Z=1;
    }
}

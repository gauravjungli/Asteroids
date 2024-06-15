#include "gauravlib.h"
void Init_grav (vector<Grav>& g, string file )
{   
    file=file+"/grav.txt";
    Read_grav(g,file);
}



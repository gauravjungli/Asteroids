#include "gauravlib.h"
void Init_grav (vector<Grav>& g, string file )
{   
    fs::path base_path = file;

	fs::path file_name= "grav.txt";
	fs::path full_path = base_path / file_name;
	string file1=	full_path.string();
    Read_grav(g,file1);
}



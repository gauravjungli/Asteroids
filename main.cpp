
#include "gauravlib.h"



int main()
{	
	
	double om=omega;
    vector<Grav> g(res);
	Grav_sph(g);
	vector<double> x(res);
	Grid(x);
	vector<CV> w;
	
	string file=string("output/files_")+to_string(delta)+string("_")+to_string(om);
	if(filesystem::exists(file))
		deleteDirectoryContents(file);
	filesystem::create_directory(file);
	for (int i=0;i<=slides;i++)
	{	
		Uniform_IC(w,x,g,om);
		BC(w,om);
		string file1=file+string("/field_")+to_string(i)+string(".csv");
		if(i==0){Write(w,file1);continue;}
		
		March(w,om,i*finalt);
		Write(w,file1);
		cout<<i<<endl;
	}
	
	return 0;
	
}
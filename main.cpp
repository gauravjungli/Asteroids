
#include "gauravlib.h"



int main()
{	
	vector<Grav> g(res);
	Grav_sph(g);
	
    
	vector<double> x(res);
	Grid(x);
	vector<CV> w;

	std::string file="output/files_"+to_string(delta)+"_"+to_string(omega);
	if(filesystem::exists(file))
		deleteDirectoryContents(file);
	filesystem::create_directory(file); 
 	for (int i=0;i<=slides;i++)
	{	
		Uniform_IC(w,x,g);
		string file1=file+string("/field_")+to_string(i)+string(".csv");
		if(i==0){Write(w,file1);i++;}
		
		double final_om=March(w,i*finalt);
		Write(w,file1);
		Write (final_om,"output/omega.txt");
		cout<<i<<endl;
	} 
	
	return 0;
	
}
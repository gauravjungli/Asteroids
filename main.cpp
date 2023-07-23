
#include "gauravlib.h"
#include <chrono>


int main()
{	
	chrono::steady_clock sc;
	auto start = sc.now();

	std::string file="output/files_"+to_string(delta)+"_"+to_string(omega_initial);
	cout<<file<<endl;

	vector<Grav> g(res);
	Init_grav(g,file);
	
	vector<double> x(res);
	//Grid(x);
	vector<CV> w;

	//uncomment only for the solo run

//	if(filesystem::exists(file))
//		deleteDirectoryContents(file);
//	filesystem::create_directory(file); 

	string file1=file+string("/field_")+to_string(int(slides))+string(".csv");		
	string file2=file+string("/log.txt");	
	Uniform_IC(w,x,g,file);
	
	double Ang_mom=stod(par["jinertia"])*omega;
	March(w,finalt);
	par["jinertia"]=to_string(Inertia(w,1),15);
	par["jinertia1"]=to_string(Inertia(w,2),15);
	par["omega"]=to_string(Ang_mom/stod(par["jinertia"]),15);
	Write(w,file1);
	Write(x,w,file);
	Write (par,"parameters");

	ofstream myfile(file2,std::ofstream::app);
	if (!myfile) Error("Can't open the file","log.txt");
	
	myfile<<"Impact number "<<slides<<endl;

	auto end = sc.now();
	auto time_span = static_cast<chrono::duration<double>>(end - start);   // measure time span between start & end
   	myfile<<"Operation took: "<<time_span.count()<<" seconds !!!";
	return 0;
}
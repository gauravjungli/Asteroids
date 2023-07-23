
#include "gauravlib.h"
#include <chrono>


int main()
{	
	chrono::steady_clock sc;
	auto start = sc.now();
	double Ang_Shed=0;
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
	ofstream myfile(file2,std::ofstream::app);

	myfile<<"Impact number "<<slides<<endl;
	myfile<<"Initial omega --> "<<omega<<endl<<" Initial Inertia from the parameters--> "<<par["jinertia"]<<endl ;

	Uniform_IC(w,x,g,file);

	par["jinertia"]=to_string(Inertia(w,1),15);
	par["jinertia1"]=to_string(Inertia(w,2),15);
	double Ang_Mom=stod(par["jinertia"])*omega;

	myfile<<" Initial Inertia from the calculation--> "<<par["jinertia"]<<endl ;


	March(w,myfile,Ang_Shed);
	
	par["jinertia"]=to_string(Inertia(w,1),15);
	par["jinertia1"]=to_string(Inertia(w,2),15);
	par["omega"]=to_string((Ang_Mom-Ang_Shed)/stod(par["jinertia"]),15);

	Write(w,file1);
	Write(x,w,file);
	Write (par,"parameters");
	myfile<<"Initial Angular Momentum --> "<<stod(par["jinertia"])*omega<<endl<<" Total Angular Momentum Shed --> "<<Ang_Shed<<endl ;
	myfile<<"Final omega --> "<<par["omega"]<<endl<<" Final Inertia--> "<<par["jinertia"]<<endl ;

	auto end = sc.now();
	auto time_span = static_cast<chrono::duration<double>>(end - start);   // measure time span between start & end
   	myfile<<"Operation took: "<<time_span.count()<<" seconds !!! "<<endl;
	myfile<<"----------------------------------------------------------------------------"<<endl;
	myfile<<"----------------------------------------------------------------------------"<<endl;
	return 0;
}
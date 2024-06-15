
#include "gauravlib.h"
#include <chrono>

std::map <std::string, string> par;
const string par_add = "parameters";
bool set_parameter=Parameters();


const int res= (int) round(stod(par["Resolution"]));
const double PI= M_PI;
const int dump= int(stod(par["dump"]));
const double offset= stod(par["offset"]);
const double xmax=   PI;
const double xmin=   0;
const double weight= stod(par["weight"]); 
const double finalt= stod(par["Landslide simulation period"]);
const double Delta= stod(par["Friction angle"]); 
const double theta= stod(par["theta"]); 
const double slides= stod(par["slides"]);
const double epsilon= stod(par["epsilon"]); 
const double omega= stod(par["omega"]);
const double dx= (xmax-xmin-2*offset)/res;
const double past_time=stod(par["time"]);
const double dia=stod(par["Diameter"]);
const double min_h=pow(dx,4);
const double Gamma=stod(par["Gamma"]);
double delta=Delta;
const string fric_type=par["Friction type"];
const string folder=par["folder"];
const string verbose_dir=par["verbose_dir"];
const string verbose=par["verbose"];

int main()
{

	chrono::steady_clock sc;
	auto start = sc.now();
	double Ang_Shed=0;
	std::string file=folder;
	std::cout<<file<<endl;
	
	vector<Grav> g(res);
	
	Init_grav(g,file);
	
	vector<double> x(res);
	//Grid(x);
	vector<CV> w;

	//uncomment only for the solo run

	//if(filesystem::exists(file))
	//	deleteDirectoryContents(file);
	//filesystem::create_directory(file); 
	

	fs::path base_path = file;
	fs::path file_path = "data";

	fs::path file_name= string("field_")+to_string(int(slides))+string(".csv");
	fs::path full_path = base_path / file_path/file_name;
	string file1=	full_path.string();

	file_name= string("log.txt");
	full_path = base_path / file_path/file_name;
	string file2=full_path.string();	
	ofstream myfile(file2,std::ofstream::app);

	myfile<<"Impact number "<<slides<<endl;
	myfile<<"Impact time "<<past_time<<endl;
	myfile<<"Initial omega --> "<<omega<<endl<<" Initial Inertia --> "<<par["jinertia"]<<endl ;


	Uniform_IC(w,x,g);

	par["jinertia"]=to_string(Inertia(w,1),15);
	par["jinertia1"]=to_string(Inertia(w,2),15);
	double Ang_Mom=stod(par["jinertia"])*omega;

	March(w,Ang_Shed);
	
	par["jinertia"]=to_string(Inertia(w,1),15);
	par["jinertia1"]=to_string(Inertia(w,2),15);
	par["omega"]=to_string((Ang_Mom-Ang_Shed)/stod(par["jinertia"]),15);

	Write(w,file1);
	Write(x,w,file);
	Write (par);
	myfile<<"Initial Angular Momentum --> "<<Ang_Mom<<endl<<" Total Angular Momentum Shed --> "<<Ang_Shed<<endl ;
	myfile<<"Final omega --> "<<par["omega"]<<endl<<" Final Inertia--> "<<par["jinertia"]<<endl ;

	auto end = sc.now();
	auto time_span = static_cast<chrono::duration<double>>(end - start);   // measure time span between start & end
   	myfile<<"Operation took: "<<time_span.count()<<" seconds !!! "<<endl;
	myfile<<"----------------------------------------------------------------------------"<<endl;
	myfile<<"----------------------------------------------------------------------------"<<endl;
	return 0;
}


#include "gauravlib.h"
#include <chrono>

std::map <std::string, string> par;
//for the debugging mode only

const string par_add = "parameters";
bool set_parameter = Parameters();

const int res = (int) round(stod(par["Resolution"]));
const double PI = M_PI;
const int dump = int(stod(par["dump"]));
const double offset = stod(par["offset"]);
const double xmax =   PI;
const double xmin =   0;
const double weight = stod(par["Correction weight"]); 
const double finalt = stod(par["Landslide simulation period"]);
const double Delta = stod(par["Friction angle"]); 
const double theta = stod(par["Minmod Limiter"]); 
const double slides = stod(par["slides"]);
const double epsilon = stod(par["epsilon"]); 
const double omega = stod(par["omega"]);
const double dx = stod(par["dx"]);
const double past_time = stod(par["time"]);
const double dia = stod(par["Current diameter"]);
const double min_h = 1E-15;
const double min_u = 1E-15;
const double seismic_time =  stod(par["Seismic_shaking_time"]); 
double delta = Delta;
const string fric_type = par["Friction type"];
const string Output_folder = par["Data folder"];
const string verbose_dir = par["verbose_dir"];
const string verbose = par["verbose"];
double mass_shed = stod(par["Mass shed"]);
double k_d = stod(par["k_d"]);
int main()
{

	chrono::steady_clock sc;
	auto start = sc.now();
	double Ang_Shed=0;
	std::string file=Output_folder;
	
	std::ofstream outfile("cpp_output.txt",std::ofstream::app);  // Create or open output file

	std::streambuf *coutbuf = std::cout.rdbuf(); 
    std::streambuf *cerrbuf = std::cerr.rdbuf(); 
    if (outfile.is_open()) {
        std::cout.rdbuf(outfile.rdbuf()); // Redirect cout	
        std::cerr.rdbuf(outfile.rdbuf()); // Redirect cerr
	}
	vector<Grav> g(res+1);
	fs::path base_path = file;
	Init_grav(g,base_path.parent_path());
	
	vector<CV> w,wl,wr;
	

	fs::path file_name= string("field_")+to_string(int(slides))+string(".csv");
	fs::path full_path = base_path / file_name;
	string file1=	full_path.string();

	file_name= string("log.txt");
	full_path = base_path / file_name;
	string file2=full_path.string();	
	ofstream myfile(file2,std::ofstream::app);

	file_name= string("dia.txt");
	full_path = base_path / file_name;
	file2=full_path.string();	
	ofstream dia_file(file2,std::ofstream::app);

	dia_file<< past_time <<"\t"<<std::setprecision(18)<<"\t"<< dia << "\t" << epsilon <<"\t"<<mass_shed<< endl;

	myfile<<"Impact number -->  "<<slides<<endl;
	myfile<<"Impact time -->  "<<past_time<<endl<<"Seismic shaking duration --> "<< seismic_time<<endl;
	myfile<<" Initial Inertia as per python module --> "<<par["jinertia"]<<endl ;



	Initial_Condition(w,wl,wr,g);

	par["jinertia"]=to_string(Inertia(w,1),15);
	par["jinertia1"]=to_string(Inertia(w,2),15);
	double Ang_Mom=stod(par["jinertia"])*omega;

	myfile<<"Initial omega --> "<<omega<<endl<<" Initial Inertia --> "<<par["jinertia"]<<endl ;

	March(w,wl,wr,Ang_Shed);
	
	par["jinertia"]=to_string(Inertia(w,1),15);
	par["jinertia1"]=to_string(Inertia(w,2),15);
	par["omega"]=to_string((Ang_Mom-Ang_Shed)/stod(par["jinertia"]),15);
	par["Mass shed"] = to_string(mass_shed,15);
	
	Write_data(w,file1);
	Write_base(w,base_path.parent_path()); 
	Write_par(par); 

	myfile<<"Initial Angular Momentum --> "<<Ang_Mom<<endl<<" Total Angular Momentum Shed --> "<<Ang_Shed<<endl<<" Total Mass Shed --> "<<mass_shed<<endl ;
	myfile<<"Final omega --> "<<par["omega"]<<endl<<" Final Inertia--> "<<par["jinertia"]<<endl ;

	auto end = sc.now();
	auto time_span = static_cast<chrono::duration<double>>(end - start);   // measure time span between start & end
   	cout<<"Operation took: "<<time_span.count()<<" seconds !!! "<<endl;
	myfile<<"----------------------------------------------------------------------------"<<endl;
	myfile<<"----------------------------------------------------------------------------"<<endl;
	myfile.close();
	std::cout.rdbuf(coutbuf); 
    std::cerr.rdbuf(cerrbuf);
	outfile.close();
	dia_file.close();
	return 0;
}

#include "gauravlib.h"

void Write ( const double om, const double t)
{
ofstream myfile("output/omega.txt",std::ofstream::app);
if (!myfile) Error("Can't open output file","Omega.txt");
myfile<<std::setprecision(12)<<t<<" "<<om<< endl;
myfile.close();
}
void Write (const vector<CV>& w,const int count, string file)
{
file=file+string("/field_")+to_string(count)+string(".csv");
ofstream myfile(file);
if (!myfile) Error("Can't open output file field_%f",file);
for (int i=2;i<res-2;i++)
		myfile<<std::setprecision(12)<<w[i].x<<","<<w[i].h<<","<<w[i].b<<","<< w[i].u<<","<<w[i].v<<","<<w[i].psi<<"\n";
myfile.close();
}

void Read ( vector<double>& v,string file)
{
ifstream myfile(file);
if (!myfile) Error("Can't open input file",file);
double inp;
while(myfile>>inp)
		v.push_back(inp);
myfile.close();
}

void Error (string s1, string s2)
{
	cout<< s1<<" "<<s2<<endl;
}

void deleteDirectoryContents(const std::string& dir_path)
{
    for (const auto& entry : std::filesystem::directory_iterator(dir_path)) 
        std::filesystem::remove_all(entry.path());
}
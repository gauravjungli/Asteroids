#include "gauravlib.h"

void write ( const double Om, const double t)
{
ofstream myfile("output/Omega.txt");
if (!myfile) error("Can't open output file","Omega.txt");
		myfile<<t<<" "<<Om<< endl;
}
void write (const vector<CV>& w,const double t)
{
string file=string("output/field_")+to_string(t)+string(".txt");
ofstream myfile(file);
if (!myfile) error("Can't open output file field_%f",file);
for (int i=2;i<res-2;i++)
		myfile<<w[i].h<<"  "<< w[i].u<<"  "<<w[i].v<<endl;
}

void read ( vector<double>& v,string file)
{
ifstream myfile(file);
if (!myfile) error("Can't open input file",file);
double inp;
while(myfile>>inp)
		v.push_back(inp);
}

void error (string s1, string s2)
{
	cout<< s1<<" "<<s2<<endl;
}
#include "gauravlib.h"

void Write ( const double om, const double t)
{
ofstream myfile("output/Omega.txt");
if (!myfile) Error("Can't open output file","Omega.txt");
		myfile<<t<<" "<<om<< endl;
}
void Write (const vector<CV>& w,const double t)
{
string file=string("output/field_")+to_string(t)+string(".txt");
ofstream myfile(file);
if (!myfile) Error("Can't open output file field_%f",file);
for (int i=2;i<res-2;i++)
		myfile<<w[i].h<<"  "<< w[i].u<<"  "<<w[i].v<<endl;
}

void Read ( vector<double>& v,string file)
{
ifstream myfile(file);
if (!myfile) Error("Can't open input file",file);
double inp;
while(myfile>>inp)
		v.push_back(inp);
}

void Error (string s1, string s2)
{
	cout<< s1<<" "<<s2<<endl;
}
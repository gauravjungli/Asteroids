#include "gauravlib.h"
void write ( const double Om, const double t)
{
ofstream myfile("Omega.txt");
if (!myfile) error("Can't open output file","Omega.txt");
		myfile<<t<<" "<<Om<< endl;
}
void write (const vector<CV>& v,const double t)
{
string file=string("field_")+to_string(t);
ofstream myfile(file);
if (!myfile) error("Can't open output file field_%f",file);
for (int i=0;i<v.size();i++)
		myfile<<v[i].h<<"  "<< v[i].u<<"  "<<v[i].v<<endl;
}

void read ( vector<double>& v,string file)
{
ifstream myfile(file);
if (!myfile) error("Can't open input file",file);
double inp;
while(myfile>>inp)
		v.push_back(inp);
}
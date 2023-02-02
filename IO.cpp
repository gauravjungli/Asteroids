#include "gauravlib.h"

void Write ( const double om, const double t, string file)
{
ofstream myfile(file,std::ofstream::app);
if (!myfile) Error("Can't open output file","Omega.txt");
myfile<<std::setprecision(18)<<t<<" "<<om<< endl;
myfile.close();
}
void Write ( const double om, string file)
{
ofstream myfile(file,std::ofstream::out);
if (!myfile) Error("Can't open output file","Omega.txt");
myfile<<std::setprecision(18)<<" "<<omega<< endl;//change to om for SSAH
myfile.close();
}
void Write (const vector<CV>& w, string file)
{
ofstream myfile(file);
if (!myfile) Error("Can't open output file field_%f",file);
for (int i=2;i<res-2;i++)
		myfile<<std::setprecision(18)<<w[i].x<<","<<w[i].h<<","<<w[i].b<<","<< w[i].u<<","<<w[i].v<<","<<w[i].psi<<"\n";
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

bool Parameters()
{
	ifstream myfile("parameters");
	if (!myfile) Error("Can't open file", "parameters");
	string line;
	while (getline(myfile, line))  
    {
        istringstream iss(line);
        string word1;
        if (!(iss >> word1))
            continue;  // line had no words
        string word2;
        if (!(iss >> word2))
            continue;  // line only had one word

		par[word1]=stod(word2);
    }
    myfile.close();
    myfile.open("output/omega.txt");
    if (!myfile) Error("Can't open file", "omega.txt");
    getline(myfile, line);
    istringstream iss(line);
    string word;
    iss >> word;
    cout<<word;
    par["omega"]=stod(word);
    return true;
}
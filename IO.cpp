#include "gauravlib.h"

void Write ( const double om, const double t, string file)
{
ofstream myfile(file,std::ofstream::app);
if (!myfile) Error("Can't open output file","Omega.txt");
myfile<<std::setprecision(18)<<past_time+t<<" "<<om<< endl;
myfile.close();
}

std::string to_string(double value, int precision) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(precision) << value;
    return oss.str();
}

void Write ( const double om, string file)
{
ofstream myfile(file,std::ofstream::out);
if (!myfile) Error("Can't open output file","Omega.txt");
myfile<<std::setprecision(18)<<" "<<om<< endl;
myfile.close();
}


void Write (const vector<CV>& w, string file)
{
ofstream myfile(file);
if (!myfile) Error("Can't open output file field",file);
for (int i=2;i<res-2;i++)
		myfile<<std::setprecision(18)<<w[i].x<<","<<w[i].h<<","<<w[i].b<<","<< w[i].q<<","<<w[i].r<<","<<w[i].psi<<","<<dia<<"\n";
myfile.close();
}

void Write (const vector<double>& x, const vector<CV>& w, string file)
{
ofstream myfile(file+"/base.txt");
if (!myfile) Error("Can't open output file field",file);
for (int i=0;i<res;i++)
		myfile<<std::setprecision(18)<<x[i]<<","<<(w[i].b+epsilon/gamma*w[i].h)<<"\n";
myfile.close();
}

void Write ( std::map <std::string, string> par, string file)
{
ofstream myfile(file,std::ofstream::out);
if (!myfile) 
{
    Error("Can't open the file","parameters");
    return;
}

for (auto i = par.begin(); i != par.end(); i++)
    {
        myfile<<std::string(50,'-')<<"\n";
		myfile<<left<<std::setw(25)<< i->first << i->second<<endl;
    }
myfile<<std::string(50,'-')<<"\n";
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
        if (line.find("--")!=std::string::npos)
            continue;
        istringstream iss(line);
        string word1;
        if (!(iss >> word1))
            continue;  // line had no words
        string word2;
        if (!(iss >> word2))
            continue;  // line only had one word

		par[word1]= word2;
     
    
    }
    myfile.close();
    return true;
}


// Function to read a 2D array from a file
void Read_grav( vector<Grav>& g, const string& file)
{
    std::ifstream f(file);
    std::vector<double> row;
    if (!f) Error("Can't open file", "grav.txt");
    double num;
    int i=2;
    while (f >> num)
    {
        row.push_back(num);

        // Check if the row is complete
        if (row.size() == 2) 
        {
            g[i].X1=row[0];
            g[i].X2=row[1];
            g[i].X3=0;
            row.clear();
            i++;
        }
    }
}
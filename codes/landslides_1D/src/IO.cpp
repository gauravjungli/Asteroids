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


void Write_data (const vector<CV>& w, string file)
{
ofstream myfile(file);
if (!myfile) Error("Can't open output file field",file);
for (int i=0;i<res;i++)
	myfile<<std::setprecision(18)<<w[i].x<<","<<w[i].b<<","<<w[i].h<<","<< w[i].db
    <<","<< w[i].ddb<<","<< w[i].u<<","<<w[i].v<<","<<w[i].psi<<"\n";
myfile.close();
}

void Write_final_data (const vector<CV>& w, string file)
{
ofstream myfile(file);
if (!myfile) Error("Can't open output file field",file);
for (int i=0;i<res;i++)
	myfile<<std::setprecision(18)<<w[i].x<<","<<w[i].b<<","<<w[i].h<<","<< w[i].db
    <<","<< w[i].ddb<<","<< w[i].u<<","<<w[i].v<<","<<w[i].psi<<","<<w[i].g.X1<<","<<w[i].g.X2<<"\n";
myfile.close();
}

void Write_base ( const vector<CV>& w, string file)
{   
    fs::path base_path=file;
    fs::path file_name= string("base.txt");
	fs::path full_path = base_path / file_name;
	string file2=full_path.string();
    ofstream myfile(file2);
    if (!myfile) 
    {
        Error("Can't open output file field",file2);
        return;
    }
    for (int i=0;i<res;i++)
        myfile<<std::setprecision(18)<<w[i].x<<","<<w[i].b<<","<<w[i].h<<","<<w[i].db<<","<<w[i].ddb<<"\n";
    myfile.close();
}

void Write_par ( std::map <std::string, string> par)
{
ofstream myfile(par_add,std::ofstream::out);
if (!myfile) 
{
    Error("Can't open the file","parameters");
    return;
}

for (auto i = par.begin(); i != par.end(); i++)
    {
        myfile<<std::string(100,'-')<<"\n";
		myfile<<left<<std::setw(50)<< i->first<<"\t" << i->second<<endl;
    }
myfile<<std::string(100,'-')<<"\n";
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
	std::cout<< s1<<" "<<s2<<endl;
}


void deleteDirectoryContents(const std::string& dir_path)
{
    for (const auto& entry : std::filesystem::directory_iterator(dir_path)) 
        std::filesystem::remove_all(entry.path());
}


bool Parameters()
{   ifstream myfile;
    fs::path filePath1 = fs::path(par_add) ;
    fs::path filePath2 = fs::path("/home/g/Asteroids/output/Debug/run1/parameters");
    if (fs::exists(filePath1)) {
         myfile.open(filePath1);
        std::cout << "File found in first directory.\n";
    } 
    
     else if (Debug and fs::exists(filePath2)) 
    {
        myfile.open(filePath2);
        std::cout << "File found in second directory.\n"; 
    }    
    else {
        std::cout << "File not found in either directory.\n";
        return false;  // Exit if the file doesn't exist in either directory
    }
	
	if (!myfile) Error("Can't open file", "parameters");
	string line;
	while (getline(myfile, line))  
    {
        if (line.find("--")!=std::string::npos)
            continue;

        line = regex_replace(line, regex("^\\t+|\\t+$"), ""); 

        // Split the line based on multiple spaces
        istringstream iss(line);
        string key, value;

        // Get the key (potentially with spaces)
        getline(iss, key, '\t'); 
         key = regex_replace(key, regex("^\\s+|\\s+$"), "");
        // Discard multiple spaces
        while (iss.peek() == '\t') {
            iss.get(); 
        }

        // Get the remaining part as the value
        getline(iss, value); 
        par[key] = value;
     
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
    int i=0;
    while (f >> num)
    {
        row.push_back(num);

        // Check if the row is complete
        if (row.size() == 2) 
        {
            g[i].X1 = row[0]; 
            g[i].X2 = row[1];  
            g[i].X3 = 0;
            row.clear();
            i++;
        }
    }
}

void Read_data ( vector<double>& x,vector<double>& b,vector<double>& db,vector<double>& ddb,vector<double>& h,string file)
{
    fs::path base_path = Output_folder;

	fs::path file_name= file;
	fs::path full_path = base_path.parent_path()/ file_name;
	string file1=	full_path.string();
    ifstream myfile(file1);

    std::string line;
    int i=0;
    while(getline(myfile,line))
    {
        istringstream iss(line);
        string word1,word2,word3,word4,word5;
        getline(iss, word1, ',');
        getline(iss, word2, ',');     
        getline(iss, word3, ','); 
        getline(iss, word4, ',');     
        getline(iss, word5, ',');  

        x[i] = stod(word1);
        b[i] = stod(word2); 
        h[i] = stod(word3);   
        db[i] = stod(word4);  
        ddb[i]= stod(word5);

        i++;
    }
    myfile.close();

}
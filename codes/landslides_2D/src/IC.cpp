#include "gauravlib.h"


void Grid(vector<double> & x, vector <double> & y)
{
    for (int i=0;i<rows;i++)
    for (int j=0;j<cols;j++)
    {
        {   int i1 = index(i,j);
            x[i1]= offset+dx * (i+0.5);
            y[i1]= dy * (j+0.5);
        }
    }
}


void Uniform_IC (vector<CV> & w, vector<double> & x, vector<double> & y, vector<Grav>& g)
{   
    vector<double>  b(rows*cols,0);
    vector <double> h(rows*cols,1), u(rows*cols,0),v(rows*cols,0);
    
//Uncomment this one only if you want special initial conditions

    Base(b,h,x,y);

    // for (int j=0;j<res;j++)
    // {
    //     double max_h=1;
    //     if (b[j]<-max_h)
    //     {
    //         h[j]=min_h;
    //         b[j]=b[j]+epsilon/Gamma*uni_h;
    //     }

    //     if (b[j]>max_h)
    //     { 
    //         h[j]=h[j]+Gamma/epsilon*(b[j]-max_h);
    //         b[j]=max_h;
    //     }
        
    // } 
    
     
        for(int i=0;i<rows;i++)
    {
	    for (int j = 0; j < cols; j++)
        {
            int i1 = index(i,j);
            CV temp(h[i1],u[i1],v[i1],b[i1],g[i1],x[i1],y[i1],false);
            w.push_back(temp);
        } 
    }
	    
}


void Base ( vector<double>& b,vector<double>& h,vector<double>& x, vector<double>& y)
{
    fs::path base_path = Output_folder;

	fs::path file_name= "base.txt";
	fs::path full_path = base_path.parent_path()/ file_name;
	string file1=full_path.string();
    ifstream myfile(file1);

    std::string line;
    int i=0;
    while(getline(myfile,line))
    {
        istringstream iss(line);
        string word1,word2,word3,word4;
        getline(iss, word1, ',');
        getline(iss, word2, ','); 
        getline(iss, word3, ',');        
        getline(iss, word4, ',');
        x[i] = stod(word1);
        y[i] = stod(word2);
        b[i] = stod(word3); 
        h[i] = stod(word4); 
        i++;
    }
    myfile.close();

}

#include "gauravlib.h"

int main()
{	

	double om=omega;
    vector<Grav> g(res);
	Grav_sph(g);
	vector<double> x(res);
	Grid(x);
	vector<CV> w;
	deleteDirectoryContents("output");
	for (int i=0;i<=slides;i++)
	{	
		Uniform_IC(w,x,g,om);
		if(i==0){Write(w,i,"output");continue;}
		
		March(w,om,i*finalt,i);
		Write(w,i,"output");
		cout<<i<<endl;
	}
	
	return 0;
	
}
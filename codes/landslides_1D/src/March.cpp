#include "gauravlib.h"

	// predictor step for interior
void Predictor(vector<CV>& w,  vector<CV>& wl, vector<CV>& wr, double dt)
{  	
	BC(w);
	Edge(w,wl,wr);
	vector<CV> wtemp(w);
	
	for (int j = 2; j < res-2; j++)
	{	
		FS hl= Hx(wr[j-1],wl[j]);
		FS hr= Hx(wr[j],wl[j+1]);
		
		FS source=Source( wtemp[j], wr[j-1], wl[j], wr[j], wl[j+1]);
		w[j].Modify(wtemp[j].p - ((hr.p - hl.p) / dx - source.p) * dt 
		, wtemp[j].q - ((hr.q - hl.q) / dx - source.q) * dt
    	, wtemp[j].r - ((hr.r-hl.r) / dx - source.r) * dt);	
	}
}

// corrector step for interior
void Corrector(vector<CV>& w,  vector<CV>& wl, vector<CV>& wr, vector<CV>& w_init, double dt)
{	
	BC(w);
	Edge(w,wl,wr);
    vector<CV> wtemp(w);

	for (int j=2; j < res-2; j++)
	{
		FS hr= Hx(wr[j],wl[j+1]);
		FS hl= Hx(wr[j-1],wl[j]);
		FS source=Source( wtemp[j], wr[j-1], wl[j], wr[j], wl[j+1]);	
		w[j].Modify ( w_init[j].p * weight + (1 - weight) * (wtemp[j].p - ((hr.p - hl.p) / dx - source.p) * dt)
		, w_init[j].q * weight + (1 - weight) * (wtemp[j].q - ((hr.q - hl.q) / dx - source.q) * dt)
		,  w_init[j].r * weight + (1 - weight) * (wtemp[j].r - ((hr.r - hl.r) / dx - source.r) * dt));
	}	
}

void March (vector<CV>& w, double& Ang_Shed)
{	
	double dt = dx / 8;
	static int timesteps=0;
	double sum1=0,sum=0;
	vector<CV> wl(w),wr(w);
	
	string file1=par["verbose_dir"];
	if(!filesystem::exists(file1))
		filesystem::create_directory(file1);
	static double t=0;
	double check_t=1; 
	while(t<finalt)
	{
		if(timesteps%dump==0 && (verbose=="Yes"|| verbose=="yes" ))
		{	
			fs::path base_path = verbose_dir;
			fs::path file_name= string("field_")+to_string(int(timesteps/dump))+string(".csv");
			fs::path full_path = base_path / file_name;
			string file2=	full_path.string();
			Write(w, file2);
		}
		
		vector<CV> w_init(w);
		if (fric_type!="constant" && fric_type!="Constant")
		{
			if (t<seismic_time)
				delta= 0;
			else
				delta = std::min(25.0,Delta);
		}
		
		Shed(w, Ang_Shed);
		Predictor(w,wl,wr,dt);	
		Shed(w, Ang_Shed);
		Corrector(w, wl, wr, w_init, dt); 
		
		Time_step(wl,wr,dt,t,timesteps);

		sum=0;
		for (int i=2;i<res-2;i++)
			sum+=PI/2*(w[i].v*(pow(1+Gamma*w[i].b+epsilon*w[i].h,4)-pow(1+Gamma*w[i].b,4)))*dx;
	/*	if (t>check_t)
		{	
			if (abs(sum)<epsilon*epsilon )
			{	
				break;
			}
			sum1=sum;
			check_t+=0.1;
		}
		//std::cout<<std::setprecision(18)<<t<<"  "<<sum<<"  "<<delta<<endl;
	*/
	}
	fs::path base_path = verbose_dir;
	fs::path file_name= string("log.txt");
	fs::path full_path = base_path / file_name;
	string file2=	full_path.string();
	
	ofstream myfile(file2,std::ofstream::app);
	Shed(w, Ang_Shed);
	myfile<<"Simulation ran for time --> " <<t<<endl;
	myfile<<"Residual Angular Momentum --> " <<sum<<endl;
	myfile<<"Rate of change of Angular Momentum -->"<<sum1-sum<<endl;
	Ang_Shed+=sum;
	
}


void Time_step(vector <CV>& wl, vector <CV>& wr, double & dt, double & t, int & timesteps)
{	
	CFL(wl,wr,dt);
	// time updation
	t = t + dt;
	timesteps = timesteps + 1;
}

void CFL(vector<CV>& wl,vector<CV>& wr, double & dt)
{

	double maxspeed = 0.00000001;
	// to evaluate dt from maximum speeds (CFL condition)
	for (int j = 2; j < res - 2; j++)
	{
		double eig = max(abs(Ax(wr[j-1],wl[j],"max")),abs(Ax(wr[j-1],wl[j],"min")));
		if (maxspeed < eig)
			maxspeed = eig;
	}
	//std::cout<<maxspeed<<endl;
	dt = min(dx/4,dx/4/maxspeed);
}
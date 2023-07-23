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

void March (vector<CV>& w, double final_t)
{	
	double dt = dx / 4;
	static int timesteps=0;
	double sum1=0;
	vector<CV> wl(w),wr(w);
	
	string file=string("output/files_")+to_string(delta)+string("_")+to_string(omega_initial)+string("/data");
	if(!filesystem::exists(file))
		filesystem::create_directory(file);
	static double t=0;
	int check_t=1; 
	while(t<final_t)
	{
		if(timesteps%dump==0)
		{
			string file1=file+string("/field_")+to_string(int((slides)*1000+timesteps/dump))+string(".csv");
			Write(w, file1);
		}
		
		vector<CV> w_init(w);
		
		
		Predictor(w,wl,wr,dt);
		Corrector(w, wl, wr, w_init, dt);
		Shed(w);
		Time_step(wl,wr,dt,t,timesteps);

		double sum=0;
		for (int i=2;i<res-2;i++)
			sum+=w[i].r*dx;
		if (t>check_t)
		{	
			if (sum<0.1*epsilon)
				return;
			check_t++;
		}
		cout<<std::setprecision(18)<<t<<"  "<<sum<<endl;
		sum1=sum;
		
	}
		par["time"]=to_string(past_time+t);
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
	//cout<<maxspeed<<endl;
	dt = min(dx/4,dx/4/maxspeed);
}
#include "gauravlib.h"

	// predictor step for interior
void Predictor(vector<CV>& w,  vector<CV>& wl, vector<CV>& wr, double& om, double dt, AMB amb)
{  	
	vector<CV> wtemp(w);
	BC(w);
	Edge(w,wl,wr);
	double integral1=0, integral2=0;
	Integrals(integral1,integral2,w,wl,wr);
	Ang_mom(amb,integral1,integral2,true);
	double inertia=AMB_corrector(w,wl,wr);
	double alpha=Alpha( om , dt, amb, inertia );
	double sum1=0,sum2=0;
	for (int j = 2; j < res-2; j++)
	{	FS hr= Hx(wr[j],wl[j+1]);
		FS hl= Hx(wr[j-1],wl[j]);
		FS source=Source( wtemp[j], wr[j-1], wl[j], wr[j], wl[j+1], om, alpha );
		w[j].Modify(wtemp[j].p - ((hr.p - hl.p) / dx - source.p) * dt
		, wtemp[j].q - ((hr.q - hl.q) / dx - source.q) * dt
    	, wtemp[j].r - ((hr.r-hl.r) / dx - source.r) * dt);
		//sum1+=om*(hr.p-hl.p)*(sin(w[j].x-dx/2)*sin(w[j].x-dx/2)+sin(w[j].x+dx/2)*sin(w[j].x+dx/2))/2*dt;
		//sum2+=source.r*dt*dx;

	}
	
	om=om+alpha*dt;
	//cout<<Alpha( om , dt, amb )<<endl;
}

// corrector step for interior
void Corrector(vector<CV>& w,  vector<CV>& wl, vector<CV>& wr, vector<CV>& w_init, double& om, double om_init, double dt, AMB amb)
{	
	BC(w);
	Edge(w,wl,wr);
    vector<CV> wtemp(w);
	double integral1=0, integral2=0;
	Integrals(integral1,integral2,w,wl,wr);
	Ang_mom(amb,integral1,integral2,false);
	double inertia=AMB_corrector(w,wl,wr);
	double alpha=Alpha( om , dt, amb, inertia );
	for (int j=2; j < res-2; j++)
	{
	FS hr= Hx(wr[j],wl[j+1]);
	FS hl= Hx(wr[j-1],wl[j]);
	FS source=Source( wtemp[j], wr[j-1], wl[j], wr[j], wl[j+1], om, alpha);	
	w[j].Modify ( w_init[j].p * weight + (1 - weight) * (wtemp[j].p - ((hr.p - hl.p) / dx - source.p) * dt)
	, w_init[j].q * weight + (1 - weight) * (wtemp[j].q - ((hr.q - hl.q) / dx - source.q) * dt)
	,  w_init[j].r * weight + (1 - weight) * (wtemp[j].r - ((hr.r - hl.r) / dx - source.r) * dt));
	}
	om=om_init*weight+(1-weight)*(om+ alpha*dt);

}

double March (vector<CV>& w, double final_t)
{	
	double dt = dx / 4;
	static int timesteps=0;
	vector<CV> wl(w),wr(w);
	double om=0,inertia=0,integral1=0,integral2=0;
	Integrals(integral1,integral2,w,wl,wr);
	AMB amb(integral1,integral2);
	
	
	string file=string("output/files_")+to_string(delta)+string("_")+to_string(omega)+string("/data");
	filesystem::create_directory(file);
	static double t=0;
	while(t<final_t)
	{
		if(timesteps%dump==0)
		{
			string file1=file+string("/field_")+to_string(timesteps/dump)+string(".csv");
			Write(w, file1);
		    file1=string("output/files_")+to_string(delta)+string("_")+to_string(omega)+string("/omega.txt");
			Write(omega+epsilon*om,t,file1);
		}
		vector<CV> w_init(w);
		double om_init=om;
		
		Predictor(w,wl,wr,om,dt,amb);

		Corrector(w,wl,wr,w_init,om,om_init,dt,amb);
	
		Shed(w);

		Time_step(wl,wr,dt,t,timesteps);

		double sum=0;
		for (int i=2;i<res-2;i++)
		{
			sum+=w[i].h*(1+0*2*epsilon*w[i].lambda)*sin(w[i].x);
		}
		cout<<std::setprecision(12)<<t<<"  "<<sum<<endl;
		
	}
	
	return omega+epsilon*om;
}


void Time_step(vector <CV>& wl, vector <CV>& wr, double & dt, double & t, int & timesteps)
{	
	
	CFL(wl,wr,dt);

	// time updation
	t = t + dt;
	//timep = timep + dt;
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
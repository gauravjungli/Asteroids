#include "gauravlib.h"

	// predictor step for interior
void Predictor(vector<CV>& w,  vector<CV>& wl, vector<CV>& wr, double& om, double alpha, double dt)
{  	
	vector<CV> wtemp(w);
	


	for (int j = 2; j < res-2; j++)
	{	
		FS hl= Hx(wr[j-1],wl[j]);
		FS hr= Hx(wr[j],wl[j+1]);
		
		FS source=Source( wtemp[j], wr[j-1], wl[j], wr[j], wl[j+1], om, alpha );
		w[j].Modify(wtemp[j].p - ((hr.p - hl.p) / dx - source.p) * dt
		, wtemp[j].q - ((hr.q - hl.q) / dx - source.q) * dt
    	, wtemp[j].r - ((hr.r-hl.r) / dx - source.r) * dt);
		
	}

	om=om+alpha*dt;
}

// corrector step for interior
void Corrector(vector<CV>& w,  vector<CV>& wl, vector<CV>& wr, vector<CV>& w_init, double& om, double alpha, double om_init, double dt)
{	
	BC(w);
	Edge(w,wl,wr);
    vector<CV> wtemp(w);

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
	double alpha=0,alpha_sys=0, om=0;
	static int timesteps=0;
	double sum1=0;
	vector<CV> wl(w),wr(w);
	
	double inertia=0,integral1=0,integral2=0;
	
	AMB amb(integral1,integral2);
	
	
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
		    file1=string("output/files_")+to_string(delta)+string("_")+to_string(omega_initial)+string("/omega.txt");
			Write(omega+epsilon*om,t,file1);
		}
		if (t>check_t)
		{	
			double max_u=1e-5;
			for(int j=0;j<res;j++)
			{
				if (w[j].u>max_u)
					max_u=w[j].u;	
			}
			if (max_u<5e-4)
			{
				par["time"]=past_time+t;
				return omega+epsilon*om;
			}
			check_t++;
		}
		vector<CV> w_init(w);
		double om_init=om;
		
		BC(w);
		Edge(w,wl,wr);
		if (t<1e-12)
		{
			Integrals(integral1,integral2,w,wl,wr);
			amb=AMB(integral1,integral2);
			
		}
		
		Alpha_fric(w,wl,wr,om,alpha);
		Predictor(w,wl,wr,om,alpha,dt);

		Integrals(integral1,integral2,w,wl,wr);
		Ang_mom(amb,integral1,integral2,true);
		inertia=AMB_corrector(w,wl,wr,true);
		alpha_sys=Alpha_sys( om , dt, amb, inertia );
	

		double alpha_init=alpha;
		Alpha_fric(w,wl,wr,om,alpha);	
		Corrector(w, wl, wr, w_init, om, alpha, om_init, dt);

		Integrals(integral1,integral2,w,wl,wr);
		Ang_mom(amb,integral1,integral2,true);
		inertia=AMB_corrector(w,wl,wr,true);
		alpha_sys=Alpha_sys( om , dt, amb, inertia );
		cout<<"angular acc "<<alpha_sys-((1-weight)*alpha-weight*(alpha_init))<<"  "<<alpha<<endl;
		Shed(w);
		
		Time_step(wl,wr,dt,t,timesteps);

		double sum=0;
		for (int i=2;i<res-2;i++)
			sum+=w[i].h*(1+2*epsilon*w[i].lambda)*sin(w[i].x);


		cout<<std::setprecision(18)<<t<<"  "<<sum-sum1<<endl;
		sum1=sum;
		
	}
	par["time"]=past_time+t;
	return omega+epsilon*om;
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
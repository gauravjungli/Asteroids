#include "gauravlib.h"

	// predictor step for interior
void Predictor(vector<CV>& w,  vector<CV>& wl, vector<CV>& wr, double om, double dt, AMB amb)
{  	vector<CV> wtemp(w);
	BC(wtemp,om);

	for (int j = 2; j < res-2; j++)
	{	FS hr= Hx(wr[j],wl[j+1]);
		FS hl= Hx(wr[j-1],wl[j]);
		FS source=Source( wtemp[j],om);
		w[j].Modify(wtemp[j].p - ((hr.p - hl.p) / dx - source.p) * dt
		, wtemp[j].q - ((hr.q - hl.q) / dx - source.q) * dt
    	, wtemp[j].r - ((hr.r-hl.r) / dx - source.r) * dt,om);
		om=om+Alpha( om , dt, amb )*dt;
	}
}

// corrector step for interior
void Corrector(vector<CV>& w,  vector<CV>& wl, vector<CV>& wr, vector<CV>& w_init, double om, double om_init, double dt, AMB amb)
{	
   vector<CV> wtemp(w);
	BC(wtemp,om);
	for (int j=2; j < res-2; j++)
	{
	FS hr= Hx(wr[j],wl[j+1]);
	FS hl= Hx(wr[j-1],wl[j]);
	FS source=Source( wtemp[j],om);	
	w[j].Modify ( w_init[j].p * weight + (1 - weight) * (wtemp[j].p - ((hr.p - hl.p) / dx - source.p) * dt)
	, w_init[j].q * weight + (1 - weight) * (wtemp[j].q - ((hr.q - hl.q) / dx - source.q) * dt)
	,  w_init[j].r * weight + (1 - weight) * (wtemp[j].r - ((hr.r - hl.r) / dx - source.r) * dt),om);
	}
	om=om_init*weight+(1-weight)*(om+ Alpha( om , dt, amb )*dt);
}

void March (vector<CV>& w, double & om, double finalt)
{	
	double dt = dx / 4;
	static int timesteps=0;
	double alpha=0;
	
	vector<CV> wl(w),wr(w);
	Edge(w,wl,wr,om);
	double integral1=0,integral2=0;
	Integrals(integral1,integral2,w,wl,wr);
	AMB amb(integral1,integral2);
	
	double t=0;
	while(t<finalt)
	{
		if(timesteps%dump==0)
		{
			Write(w,t);
			Write(om,t);
		}
		vector<CV> w_init(w);
		double om_init=om;
		Edge(w,wl,wr,om);
		Predictor(w,wl,wr,om,dt,amb);
		Corrector(w,wl,wr,w_init,om,om_init,dt,amb);
		Integrals(integral1,integral2,w,wl,wr);
		Ang_mom(amb,integral1,integral2);
		Shed(w,om);
		Time_step(wl,wr,dt,t,timesteps);
	}
	
	
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
		double eig = Ax(wr[j-1],wl[j]);
		if (maxspeed < abs(eig))
		maxspeed = abs(eig);
	}

	dt = min(dx/8,dx/8/maxspeed);

}
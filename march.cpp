#include "gauravlib.h"

	// predictor step for interior
void predictor(vector<CV>& w,  vector<CV>& wl, vector<CV>& wr, double Om, double dt)
{  	vector<CV> wtemp(w);
	bc(wtemp);

	for (int j = 0; j < res; j++)
	{	FS hr= Hx(wl[j+1],wr[j+1]);
		FS hl= Hx(wl[j],wr[j]);
		FS source=Source( wtemp[j]);
		w[j].modify(wtemp[j].p - ((hr.p - hl.p) / dx - source.p) * dt
		, wtemp[j].q - ((hr.q - hl.q) / dx - source.q) * dt
    	, wtemp[j].r - ((hr.r-hl.r) / dx - source.r) * dt);
	}
}

// corrector step for interior
void corrector(vector<CV>& w,  vector<CV>& wl, vector<CV>& wr, vector<CV>& wtemphat, double Om, double dt)
{	
   vector<CV> wtemp(w);
	bc(wtemp);
	for (int j=0; j < res; j++)
	{
	FS hr= Hx(wl[j+1],wr[j+1]);
	FS hl= Hx(wl[j],wr[j]);
	FS source=Source( wtemp[j]);	
	w[j].modify ( wtemphat[j].r * weight + (1 - weight) * (wtemp[j].p - ((hr.p - hl.p) / dx - source.p) * dt)
	, wtemphat[j].r * weight + (1 - weight) * (wtemp[j].q - ((hr.q - hl.q) / dx - source.q) * dt)
	,  wtemphat[j].r * weight + (1 - weight) * (wtemp[j].r - ((hr.r - hl.r) / dx - source.r) * dt));
	}
}

void march (vector<CV>& w, double & Om, double finalt)
{	static int timesteps=0;
	double dt = dx / 4;
	double momincb=0;
	double t=0;
	while(t<finalt)
	{
		if(timesteps%dump==0)
		{
			write(w,t);
			write(Om,t);
		}
		vector<CV> wtemp(w);
		vector<CV> wtemphat(wtemp);
		vector<CV> wl(w),wr(w);
		edge(w,wl,wr);
		predictor(w,wl,wr,Om,dt);
		corrector(w,wl,wr,wtemphat,Om,dt);
		Ang_mom(w,wtemphat,Om,dt,momincb);
		shed(w);
		time_step(wl,wr,dt,t,timesteps);
	}
	
	
}

void time_step(vector <CV>& wl, vector <CV>& wr, double & dt, double & t, int & timesteps)
{	
	
	cfl(wl,wr,dt);

	// time updation
	t = t + dt;
	//timep = timep + dt;
	timesteps = timesteps + 1;

}

void cfl(vector<CV>& wl,vector<CV>& wr, double & dt)
{

	double maxspeed = 0.00000001;
	// to evaluate dt from maximum speeds (CFL condition)
	for (int j = 1; j < res - 1; j++)
	{
		double eig = ax(wl[j],wr[j]);
		if (maxspeed < abs(eig))
		maxspeed = abs(eig);
	}

	dt = min(dx/8,dx/8/maxspeed);

}
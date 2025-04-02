#include "gauravlib.h"

	// predictor step for interior
void Predictor(vector<CV>& w,  vector<CV>& wl, vector<CV>& wr, vector<CV>& wb, vector<CV>& wt, double dt)
{  	
	BC(w);
	Edge(w,wl,wr,wb,wt);
	vector<CV> wtemp(w);
	
	for(int i=2;i<rows-2;i++)
	{
	for (int j = 0; j < cols; j++)
	{	
		int i0 = index(i,j);
		int i_1 = index(i-1,j);
		int i1 = index(i+1,j);

		int j0 = index(i,j);
		int j_1 = index(i,j-1);
		int j1 = index(i,j+1);

		FS hl= Hx(wr[i_1],wl[i0]);
		FS hr= Hx(wr[i0],wl[i1]);
		FS hb= Hy(wt[j_1],wb[j0]);
		FS ht= Hy(wt[j0],wb[j1]);
		FS source=Source( wtemp[i0], wr[i_1], wl[i0], wr[i0], wl[i1]);
		w[i0].Modify(wtemp[i0].p - ((hr.p - hl.p) / dx + (ht.p - hb.p) / dy - source.p) * dt
		, wtemp[i0].q - ((hr.q - hl.q) / dx + (ht.q - hb.q) / dy - source.q) * dt
		, wtemp[i0].r - ((hr.r - hl.r) / dx + (ht.r - hb.r) / dy - source.r) * dt);
	}
	}
}

// corrector step for interior
void Corrector(vector<CV>& w,  vector<CV>& wl, vector<CV>& wr, vector<CV>& wb, vector<CV>& wt, vector<CV>& w_init, double dt)
{	
	BC(w);
	Edge(w,wl,wr,wb,wt);
    vector<CV> wtemp(w);

	for(int i=2;i<rows-2;i++)
	{
	for (int j = 0; j < cols; j++)
	{	
		int i0 = index(i,j);
		int i_1 = index(i-1,j);
		int i1 = index(i+1,j);

		int j0 = index(i,j);
		int j_1 = index(i,j-1);
		int j1 = index(i,j+1);

		FS hl= Hx(wr[i_1],wl[i0]);
		FS hr= Hx(wr[i0],wl[i1]);
		FS hb= Hy(wt[j_1],wb[j0]);
		FS ht= Hy(wt[j0],wb[j1]);
		FS source=Source( wtemp[i0], wr[i_1], wl[i0], wr[i0], wl[i1]);
		w[j].Modify ( w_init[i0].p * weight + (1 - weight) * (wtemp[i0].p - ((hr.p - hl.p) / dx + (ht.p - hb.p) / dy - source.p) * dt)
		, w_init[i0].q * weight + (1 - weight) * (wtemp[i0].q - ((hr.q - hl.q) / dx + (ht.q - hb.q) / dy - source.q) * dt)
		,  w_init[i0].r * weight + (1 - weight) * (wtemp[i0].r - ((hr.r - hl.r) / dx + (ht.r - hb.r) / dy - source.r) * dt));
	}
	}
}

void March (vector<CV>& w, double& Ang_Shed)
{	
	double dt = min(dx / 8, dy/8);
	double maximum=0;
	static int timesteps=0;
	vector<CV> wl(w),wr(w),wb(w),wt(w);
	
	string file1 = par["verbose_dir"];
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
			std::cout<<std::setprecision(18)<<int(timesteps/dump)<<"  "<<t<<"  "<<maximum<<"  "<<delta<<endl;
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
		Predictor(w,wl,wr,wb,wt,dt);	
		Shed(w, Ang_Shed);
		Corrector(w, wl, wr,wb,wt, w_init, dt); 
		
		Time_step(wl,wr,wb,wt,dt,t,timesteps);


		for (int i=0;i<rows*cols;i++)
			maximum =max({maximum,w[i].u*w[i].h,w[i].v*w[i].h});
		if (t>check_t)
		{	
			if (abs(maximum)<epsilon*epsilon )
			{	
				break;
			}
			check_t+=1;
		}
		
	
	}
	fs::path base_path = verbose_dir;
	fs::path file_name= string("log.txt");
	fs::path full_path = base_path / file_name;
	string file2=	full_path.string();
	
	ofstream myfile(file2,std::ofstream::app);
	Shed(w, Ang_Shed);
	myfile<<"Simulation ran for time --> " <<t<<endl;
	myfile<<"Maximum momentum --> " <<maximum<<endl;

	//Ang_Shed+=sum;
	
}


void Time_step(vector <CV>& wl, vector <CV>& wr, vector<CV>& wb,vector<CV>& wt, double & dt, double & t, int & timesteps)
{	
	CFL(wl,wr,wb,wt,dt);
	// time updation
	t = t + dt;
	timesteps = timesteps + 1;
}

void CFL(vector<CV>& wl,vector<CV>& wr, vector<CV>& wb,vector<CV>& wt, double & dt)
{

	double maxspeedx = 0.00000001;
	double maxspeedy = 0.00000001;
	// to evaluate dt from maximum speeds (CFL condition)
	for(int i=2;i<rows-2;i++)
	{
		for (int j = 0; j < cols; j++)
		{	
			int i0 = index(i,j);
			double eigx = max(abs(Ax(wl[i0],wr[i0],"max")),abs(Ax(wl[i0],wr[i0],"min")));
			double eigy = max(abs(Ay(wb[i0],wt[i0],"max")),abs(Ay(wb[i0],wt[i0],"min")));
			if (maxspeedx < eigx)
				maxspeedx = eigx;
			if (maxspeedy < eigy)
				maxspeedy = eigy;
			{
				/* code */
			}
			
		}
	}
	//std::cout<<maxspeed<<endl;
	dt = min({1e-2,dx/8/maxspeedx,dy/8/maxspeedy});
}
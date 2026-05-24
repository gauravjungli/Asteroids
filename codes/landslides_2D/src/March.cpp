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
		FS source=Source( wtemp[i0], wr[i_1], wl[i0], wr[i0], wl[i1], wt[j_1], wb[j0], wt[j0], wb[j1]);
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
		FS source=Source( wtemp[i0], wr[i_1], wl[i0], wr[i0], wl[i1], wt[j_1], wb[j0], wt[j0], wb[j1]);
		w[i0].Modify ( w_init[i0].p * weight + (1 - weight) * (wtemp[i0].p - ((hr.p - hl.p) / dx + (ht.p - hb.p) / dy - source.p) * dt)
		, w_init[i0].q * weight + (1 - weight) * (wtemp[i0].q - ((hr.q - hl.q) / dx + (ht.q - hb.q) / dy - source.q) * dt)
		,  w_init[i0].r * weight + (1 - weight) * (wtemp[i0].r - ((hr.r - hl.r) / dx + (ht.r - hb.r) / dy - source.r) * dt));
	}
	}
}

void RK3(vector<CV>& w,  vector<CV>& wl, vector<CV>& wr, vector<CV>& wb, vector<CV>& wt, double dt)
{  	
	BC(w);
	Edge(w,wl,wr,wb,wt);
	vector<CV> wtemp(w);
	vector<CV> w1(w);
	vector<CV> w2(w);

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

		FS source=Source( wtemp[i0], wr[i_1], wl[i0], wr[i0], wl[i1], wt[j_1], wb[j0], wt[j0], wb[j1]);

		w1[i0].p = w[i0].p - ((hr.p - hl.p) / dx + (ht.p - hb.p) / dy - source.p) * dt ;
		w1[i0].q = w[i0].q - ((hr.q - hl.q) / dx + (ht.q - hb.q) / dy- source.q) * dt;
		w1[i0].r = w[i0].r - ((hr.r - hl.r) / dx + (ht.r - hb.r) / dy - source.r) * dt;

		w1[i0].Modify( w1[i0].p,  w1[i0].q , w1[i0].r );

		
		FS friction = Friction(w[i0]);

		if (w[i0].psi>0 and mu>epsilon)
		{
		 	CV w_new (w1[i0]);
			
			w_new.Modify(w_new.p,w_new.q-friction.q * dt, w_new.r - friction.r * dt);

 		 	  if(w_new.u*w1[i0].u<0)
				w_new = CV(w_new.h,min_u*sign(w1[i0].u),w_new.v,w[i0].b,w[i0].Ra,w[i0].dR,w[i0].ddR,w[i0].g,w[i0].x,w[i0].y);

 			 if(w_new.v*w1[i0].v<0)
				w_new = CV(w_new.h,w_new.u,min_u*sign(w1[i0].v),w[i0].b,w[i0].Ra,w[i0].dR,w[i0].ddR,w[i0].g,w[i0].x,w[i0].y);  
 
			w1[i0] = w_new; 

		}
		}
	}

	BC(w1);
	Edge(w1,wl,wr,wb,wt);

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

		FS source = Source( w1[i0], wr[i_1], wl[i0], wr[i0], wl[i1], wt[j_1], wb[j0], wt[j0], wb[j1]);

		w2[i0].p = 0.75*w[i0].p + 0.25*(w1[i0].p - ((hr.p - hl.p) / dx  + (ht.p - hb.p) / dy - source.p) * dt);
		w2[i0].q = 0.75*w[i0].q + 0.25*(w1[i0].q - ((hr.q - hl.q) / dx  + (ht.q - hb.q) / dy - source.q) * dt);
		w2[i0].r = 0.75*w[i0].r + 0.25*(w1[i0].r - ((hr.r - hl.r) / dx + + (ht.r - hb.r) /dy - source.r) * dt);

		w2[i0].Modify( w2[i0].p,  w2[i0].q , w2[i0].r );

		
		FS friction = Friction(w1[i0]);

		if (w1[i0].psi>0 and mu>epsilon)
		{
		 	CV w_new (w2[i0]);
			
			w_new.Modify(w_new.p,w_new.q-0.25*friction.q * dt, w_new.r - 0.25*friction.r * dt);

 		 	  if(w_new.u*w2[i0].u<0)
				w_new = CV(w_new.h,min_u*sign(w2[i0].u),w_new.v,w[i0].b,w[i0].Ra,w[i0].dR,w[i0].ddR,w[i0].g,w[i0].x,w[i0].y);

 			 if(w_new.v*w2[i0].v<0)
				w_new = CV(w_new.h,w_new.u,min_u*sign(w2[i0].v),w[i0].b,w[i0].Ra,w[i0].dR,w[i0].ddR,w[i0].g,w[i0].x,w[i0].y);  
 
			w2[i0]=w_new; 

		}
		}
	}


	BC(w2);
	Edge(w2,wl,wr,wb,wt);

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

		FS source = Source( w2[i0], wr[i_1], wl[i0], wr[i0], wl[i1], wt[j_1], wb[j0], wt[j0], wb[j1]);

		w[i0].p = (1.0/3.0)*w[i0].p + (2.0/3.0)*(w2[i0].p - ((hr.p - hl.p) / dx + (ht.p - hb.p) / dy - source.p) * dt);
		w[i0].q = (1.0/3.0)*w[i0].q + (2.0/3.0)*(w2[i0].q - ((hr.q - hl.q) / dx + (ht.q - hb.q) / dy - source.q) * dt);
		w[i0].r = (1.0/3.0)*w[i0].r + (2.0/3.0)*(w2[i0].r - ((hr.r - hl.r) / dx + (ht.r - hb.r) / dy - source.r) * dt);

		w[i0].Modify( w[i0].p,  w[i0].q , w[i0].r );

		
		FS friction = Friction(w2[i0]);

		if (w2[i0].psi > 0 and mu > epsilon)
		{
		 	CV w_new (w[i0]);
			
			w_new.Modify(w_new.p, w_new.q-(2.0/3.0)*friction.q * dt, w_new.r - (2.0/3.0)*friction.r * dt);

 		 	  if(w_new.u*w[i0].u<0)
				w_new = CV(w_new.h,min_u*sign(w[i0].u),w_new.v,w[i0].b,w[i0].Ra,w[i0].dR,w[i0].ddR,w[i0].g,w[i0].x,w[i0].y);

 			 if(w_new.v*w[i0].v<0)
				w_new = CV(w_new.h,w_new.u,min_u*sign(w[i0].v),w[i0].b,w[i0].Ra,w[i0].dR,w[i0].ddR,w[i0].g,w[i0].x,w[i0].y);  
 
			w[i0]=w_new; 

		}
		}
	}
}

void March (vector<CV>& w, vector<CV>& wl, vector<CV>& wr, vector<CV>& wb, vector<CV>& wt, double& Ang_Shed)
{	
	double dt = min(dx / 8, dy/8);
	static int timesteps=0;
	double sum1=0,sum=0,ang_sum=0,KE=0;

	if(!filesystem::exists(verbose_dir))
	filesystem::create_directory(verbose_dir);
	
	static double t=0;
	double check_t=1; 
	fs::path base_path = Output_folder;
	fs::path file_name= string("log.txt");
	fs::path full_path = base_path / file_name;
	string file1=	full_path.string();
	base_path = verbose_dir;
	double print_time = 0;
	while(t<finalt)
	{
		if(t>=print_time)//(timesteps%dump==0 && (verbose=="Yes"|| verbose=="yes" ))//change
		{	
			fs::path base_path = verbose_dir;
			fs::path file_name= string("field_")+to_string(int(timesteps))+string(".csv");//change divide it by dump
			fs::path full_path = base_path / file_name;
			string file2=	full_path.string();
			Write_data(w, file2);
			print_time += 0.1;
		}
		
		vector<CV> w_init(w);
		if (fric_type=="Variable")
		{	
			double min_mu =  0.1 ;
			double alpha =2.0/3 ;

			if (t<100)
				mu = max(min_mu,Mu*(1-alpha*max(0.0,Gamma_max*exp(-k_d*t/2)-0.25)));
			else
				mu = tan(PI/4); 
		}

		else if (fric_type=="Constant")
			mu = Mu;
		else
		{
			cout<<"No friction type selected. Using the constant friction angle"<<endl;
		}


		 if (solver=="Predictor Corrector")
		{
			Predictor(w, wl, wr,wb,wt, dt);	
			Corrector(w, wl, wr,wb,wt, w_init, dt); 
		}

		else if (solver == "Runge-Kutta 3")
	    	RK3(w, wl, wr,wb,wt, dt);

		else
		{
			cout<<"No valid solver selected. Ending simulation"<<endl;
			return ;
		}
		
		Shed(w, Ang_Shed);
		Edge(w,wl,wr,wb,wt);
		Time_step(wl,wr,wb,wt,dt,t,timesteps);

		if (restart)
			break;
		sum1=ang_sum;
		ang_sum=0;
		sum=0;
		sum1=0;
		for(int i=2;i<rows-2;i++)
		{
		for (int j = 0; j < cols; j++)
		{	
			int i0 = index(i,j);
			sum += w[i0].J*w[i0].h*(w[i0].u*w[i0].u+w[i0].v*w[i0].v)*dx*dy;//+Jinertia1_reg(w[i])*omega)*dx;
			ang_sum +=  w[i0].J*w[i0].R_phi*w[i0].V*w[i0].h*dx*dy;//Ang_mom_reg(w[i])*dx;
			sum1 += w[i0].J*w[i0].h*dx*dy;
		}
		}
		

		if (t>check_t)
		{	
			
			if (abs(sum)<1E-4)
			{	
				break;
			}	
			check_t+=0.1;
			std::cout<<std::setprecision(18)<<t<<"  "<<sum<<"  "<<sum1<<"  "<<dt<<"  "<<endl;
		}
		
	
	}
	file_name= string("field_")+to_string(int(timesteps/dump))+string(".csv");
	 full_path = base_path / file_name;
	 string file2=	full_path.string();
	Write_data(w, file2);

	
	ofstream myfile(file1,std::ofstream::app);
	Shed(w, Ang_Shed);
	myfile<<"Simulation ran for time --> " <<t<<endl;
	myfile<<"Residual Angular Momentum --> " <<ang_sum<<endl;
	myfile<<"Residual Energy -->"<<sum<<endl;

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
		}
	}
	//std::cout<<maxspeed<<endl;
	dt = min({1e-2,dx/8/maxspeedx,dy/8/maxspeedy});
}
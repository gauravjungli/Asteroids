#include "gauravlib.h"

void LXF(vector<CV>& w,  vector<CV>& wl, vector<CV>& wr, double dt)//Needs derivation gain and fixing
{  	
	BC(w,w);
	Edge(w,wl,wr);
	vector<CV> wtemp(w);
	for (int j = 2; j < res-2; j++)
	{	
		double el = dx/dt;//max(abs(Ax(wl[j],wr[j-1],"min")),abs(Ax(wl[j],wr[j-1],"max"))); // change to dx/dt 
		double er = dx/dt;//max(abs(Ax(wl[j+1],wr[j],"min")),abs(Ax(wl[j+1],wr[j],"max")));
		FS hr = Flux(wtemp[j+1]); 
	    FS hl = Flux(wtemp[j-1]);
		FS sourcel = Source( wtemp[j-1], wtemp[j-1], wtemp[j-1], wtemp[j], wtemp[j]);
		FS sourcer = Source( wtemp[j+1], wtemp[j], wtemp[j], wtemp[j+1], wtemp[j+1]);
		FS source = (sourcel +sourcer)/2;
		double lambda = dt/dx;
		w[j].p = (wtemp[j+1].h*wr[j].J*er+wtemp[j-1].h*wl[j].J*el)*lambda/2 + 
		wtemp[j].h*(2*w[j].J-lambda*(wr[j].J*er+wl[j].J*el))/2
			- ((hr.p - hl.p) / (2*dx) - source.p) * dt;

		w[j].q = (wtemp[j+1].h*wtemp[j+1].u*wr[j].J*er + wtemp[j-1].h*wtemp[j-1].u*wl[j].J*el)/2*lambda + 
		wtemp[j].h*wtemp[j].u*(2*w[j].J-(wr[j].J*er+wl[j].J*el)*lambda)/2 
		- ((hr.q - hl.q) / (2*dx) - source.q) * dt;

		w[j].r = (wtemp[j+1].h*wtemp[j+1].v*wr[j].J/wr[j].phi*er + wtemp[j-1].h*wtemp[j-1].v*wl[j].J/wl[j].phi*el)/2*lambda + 
		wtemp[j].h*wtemp[j].v*(2*w[j].J/w[j].phi-(wr[j].J/wr[j].phi*er+wl[j].J/wl[j].phi*el)*lambda)/2 
		- ((hr.r - hl.r) / (2*dx) - source.r) * dt
		+ omega*((wtemp[j+1].h*wr[j].J/pow(wr[j].phi,2)*er + wtemp[j-1].h*wl[j].J/pow(wl[j].phi,2)*el)/2*lambda + 
		wtemp[j].h*(2*w[j].J/pow(w[j].phi,2)-(wr[j].J/pow(wr[j].phi,2)*er+wl[j].J/pow(wl[j].phi,2)*el)*lambda)/2 );


		w[j].Modify( w[j].p, w[j].q, w[j].r);

		FS friction = Friction(w[j]);

		if (wtemp[j].psi>0 and mu>epsilon)
		{
		 	CV w_new (w[j]);
			
			w_new.Modify(w_new.p,w_new.q-friction.q * dt, w_new.r - friction.r * dt);

 		 	if(w_new.u*w[j].u<0)
				w_new = CV(w_new.h,min_u*sign(w[j].u),w_new.v,w[j].b,w[j].db,w[j].ddb,w[j].g,w[j].x);

 			if(w_new.v*w[j].v<0)
			{ 
				w_new = CV(w_new.h,w_new.u,min_u*sign(w[j].v),w[j].b,w[j].db,w[j].ddb,w[j].g,w[j].x);  
			}
			w[j]=w_new; 

		}
		
	}

}

void NT(vector<CV>& w,  vector<CV>& wl, vector<CV>& wr, double dt)
{  	
	BC(w,w);
	Edge(w,wl,wr);
	vector<CV> wtemp(w);
	for (int j = 2; j < res-2; j++)
	{	
		FS hr = Flux(wtemp[j+1]); 
	    FS hl = Flux(wtemp[j-1]);
		FS sourcel = Source( wtemp[j-1], wtemp[j-1], wtemp[j-1], wtemp[j], wtemp[j]);
		FS sourcer = Source( wtemp[j+1], wtemp[j], wtemp[j], wtemp[j+1], wtemp[j+1]);
		FS source = (sourcel +sourcer)/2;
		w[j].p = (wtemp[j+1].h*wr[j].J+wtemp[j-1].h*wl[j].J)/2 + wtemp[j].h*(2*w[j].J-wr[j].J-wl[j].J)/2 
			- ((hr.p - hl.p) / (2*dx) - source.p) * dt;

		w[j].q = (wtemp[j+1].h*wtemp[j+1].u*wr[j].J + wtemp[j-1].h*wtemp[j-1].u*wl[j].J)/2 + 
		wtemp[j].h*wtemp[j].u*(2*w[j].J-wr[j].J-wl[j].J)/2 - ((hr.q - hl.q) / (2*dx) - source.q) * dt;

		w[j].r = (wtemp[j+1].h*wtemp[j+1].v*wr[j].J/wr[j].phi + wtemp[j-1].h*wtemp[j-1].v*wl[j].J/wl[j].phi)/2 + 
		wtemp[j].h*wtemp[j].v*(2*w[j].J/w[j].phi-wr[j].J/wr[j].phi-wl[j].J/wl[j].phi)/2 
		- ((hr.r - hl.r) / (2*dx) - source.r) * dt
		+ omega*((wtemp[j+1].h*wr[j].J/pow(wr[j].phi,2) + wtemp[j-1].h*wl[j].J/pow(wl[j].phi,2))/2 + 
		wtemp[j].h*(2*w[j].J/pow(w[j].phi,2)-wr[j].J/pow(wr[j].phi,2)-wl[j].J/pow(wl[j].phi,2))/2 );


		w[j].Modify( w[j].p, w[j].q, w[j].r);

		FS friction = Friction(w[j]);


		if (wtemp[j].psi>0 and mu>epsilon)
		{
		 	CV w_new (w[j]);
			
			w_new.Modify(w_new.p,w_new.q-friction.q * dt, w_new.r - friction.r * dt);

 		 	if(w_new.u*w[j].u<0)
				w_new = CV(w_new.h,min_u*sign(w[j].u),w_new.v,w[j].b,w[j].db,w[j].ddb,w[j].g,w[j].x);

 			if(w_new.v*w[j].v<0)
			{ 
				w_new = CV(w_new.h,w_new.u,min_u*sign(w[j].v),w[j].b,w[j].db,w[j].ddb,w[j].g,w[j].x);  
			}
			w[j]=w_new; 

		}
		
	}

}


void Predictor(vector<CV>& w,  vector<CV>& wl, vector<CV>& wr, double dt)
{  	
	BC(w,w);
	Edge(w,wl,wr);
	vector<CV> wtemp(w);
	for (int j = 2; j < res-2; j++)
	{	
		FS hl = Hx(wr[j-1],wl[j]);
		FS hr = Hx(wr[j],wl[j+1]);

		FS source = Source( wtemp[j], wr[j-1], wl[j], wr[j], wl[j+1]);
		w[j].Modify(wtemp[j].p - ((hr.p - hl.p) / dx - source.p) * dt 
		, wtemp[j].q - ((hr.q - hl.q) / dx - source.q) * dt
    	, wtemp[j].r - ((hr.r-hl.r) / dx - source.r) * dt);

		FS friction = Friction(wtemp[j]);

		if (wtemp[j].psi>0 and mu>epsilon)
		{
		 	CV w_new (w[j]);
			
			w_new.Modify(w_new.p,w_new.q-friction.q * dt, w_new.r - friction.r * dt);

 		 	  if(w_new.u*w[j].u<0)
			  {
				w_new = CV(w_new.h,min_u*sign(w[j].u),w_new.v,w[j].b,w[j].db,w[j].ddb,w[j].g,w[j].x);
			  }
 			 if(w_new.v*w[j].v<0)
				w_new = CV(w_new.h,w_new.u,min_u*sign(w[j].v),w[j].b,w[j].db,w[j].ddb,w[j].g,w[j].x);  
 
			w[j]=w_new; 

		}
		}
	
}

// corrector step for interior
void Corrector(vector<CV>& w,  vector<CV>& wl, vector<CV>& wr, vector<CV>& w_init, double dt)
{	
	BC(w,w);
	Edge(w,wl,wr);
    vector<CV> wtemp(w);

	for (int j=2; j < res-2; j++)
	{	
		FS hr = Hx(wr[j], wl[j+1]);
		FS hl = Hx(wr[j-1], wl[j]);
		FS source=Source( wtemp[j], wr[j-1], wl[j], wr[j], wl[j+1]);	
		w[j].Modify ( w_init[j].p * weight + (1 - weight) * (wtemp[j].p - ((hr.p - hl.p) / dx - source.p) * dt)
		, w_init[j].q * weight + (1 - weight) * (wtemp[j].q - ((hr.q - hl.q) / dx - source.q) * dt)
		,  w_init[j].r * weight + (1 - weight) * (wtemp[j].r - ((hr.r - hl.r) / dx - source.r) * dt));


		FS friction = Friction(wtemp[j]);

		if (wtemp[j].psi>0 and mu>epsilon)
		{
			CV w_new (w[j]);
			
			w_new.Modify(w_new.p, w_new.q-(1-weight)*friction.q*dt, w_new.r - (1 - weight)*friction.r*dt);

		 	 if(w_new.u*w[j].u<0)
				w_new = CV(w_new.h,min_u*sign(w[j].u),w_new.v,w[j].b,w[j].db,w[j].ddb,w[j].g,w[j].x);

 		 	 if(w_new.v*w[j].v<0)
			 { 
				w_new = CV(w_new.h,w_new.u,min_u*sign(w[j].v),w[j].b,w[j].db,w[j].ddb,w[j].g,w[j].x);    
			 }
			w[j] = w_new; 
		
		} 

	}	
}


void RK3(vector<CV>& w,  vector<CV>& wl, vector<CV>& wr, double dt)
{  	
	
	BC(w);
	Edge(w,wl,wr);
	vector<CV> w0(w);
	vector<CV> w1(w);
	vector<CV> w2(w);

	for (int j = 2; j < res-2; j++)
	{	
		
		FS hl = Hx(wr[j-1],wl[j]);
		FS hr = Hx(wr[j],wl[j+1]);

		FS source = Source( w[j], wr[j-1], wl[j], wr[j], wl[j+1]);

		w1[j].p = w[j].p - ((hr.p - hl.p) / dx - source.p) * dt ;
		w1[j].q = w[j].q - ((hr.q - hl.q) / dx - source.q) * dt;
		w1[j].r = w[j].r - ((hr.r - hl.r) / dx - source.r) * dt;

		w1[j].Modify( w1[j].p,  w1[j].q , w1[j].r );
	
	}

	BC(w1);
	Edge(w1,wl,wr);

	for (int j = 2; j < res-2; j++)
	{	
		FS hl = Hx(wr[j-1],wl[j]);
		FS hr = Hx(wr[j],wl[j+1]);

		FS source = Source( w1[j], wr[j-1], wl[j], wr[j], wl[j+1]);

		w2[j].p = 0.75*w[j].p + 0.25*(w1[j].p - ((hr.p - hl.p) / dx - source.p) * dt);
		w2[j].q = 0.75*w[j].q + 0.25*(w1[j].q - ((hr.q - hl.q) / dx - source.q) * dt);
		w2[j].r = 0.75*w[j].r + 0.25*(w1[j].r - ((hr.r - hl.r) / dx - source.r) * dt);

		w2[j].Modify( w2[j].p,  w2[j].q , w2[j].r );


	}

	BC(w2);
	Edge(w2,wl,wr);

	for (int j = 2; j < res-2; j++)
	{	
		FS hl = Hx(wr[j-1],wl[j]);
		FS hr = Hx(wr[j],wl[j+1]);

		FS source = Source( w2[j], wr[j-1], wl[j], wr[j], wl[j+1]);

		w[j].p = (1.0/3.0)*w[j].p + (2.0/3.0)*(w2[j].p - ((hr.p - hl.p) / dx - source.p) * dt);
		w[j].q = (1.0/3.0)*w[j].q + (2.0/3.0)*(w2[j].q - ((hr.q - hl.q) / dx - source.q) * dt);
		w[j].r = (1.0/3.0)*w[j].r + (2.0/3.0)*(w2[j].r - ((hr.r - hl.r) / dx - source.r) * dt);

		w[j].Modify( w[j].p,  w[j].q , w[j].r );


		
		FS friction = Friction(w[j]);

		if (w[j].psi > 0 and mu > epsilon)
		{
		 	CV w_new (w[j]);
			
			w_new.Modify(w_new.p, w_new.q-(3.0/3.0)*friction.q * dt, w_new.r - (3.0/3.0)*friction.r * dt);

 		 	  if(w_new.u*w[j].u<0)
				w_new = CV(w_new.h,min_u*sign(w[j].u),w_new.v,w[j].b,w[j].db,w[j].ddb,w[j].g,w[j].x);

 			 if(w_new.v*w[j].v<0)
				w_new = CV(w_new.h,w_new.u,min_u*sign(w[j].v),w[j].b,w[j].db,w[j].ddb,w[j].g,w[j].x);  
 
			w[j]=w_new; 

		}
	}	
	
}

void March (vector<CV>& w, vector<CV>& wl, vector<CV>& wr, double& Ang_Shed)
{	
	double dt = dx / 4;
    int timesteps=0;
	double sum1=0,sum=0,ang_sum=0,KE=0;
	
	if(!filesystem::exists(verbose_dir))
		filesystem::create_directory(verbose_dir);
    double t=0;
	double check_t= 1; 
	fs::path base_path = Output_folder;
	fs::path file_name = string("log.txt");
	fs::path full_path = base_path / file_name;
	string file1=	full_path.string();

	 file_name = string("time_log.txt");
	 full_path = base_path / file_name;
	 string file0=	full_path.string();
	 ofstream myfile0(file0,std::ofstream::app);

	base_path = verbose_dir;
	
	while(t<finalt)
	{
		if(timesteps%dump==0 && (verbose=="Yes"|| verbose=="yes" ))
		{	
			fs::path base_path = verbose_dir;
			fs::path file_name= string("field_")+to_string(int(timesteps/dump))+string(".csv");
			fs::path full_path = base_path / file_name;
			string file2=	full_path.string();
			Write_data(w, file2);
			
		}
		vector<CV> w_init(w);
		if (fric_type=="Variable")
		{	
			double min_mu =  0.1 ;
			double alpha =2.0/3 ;

			if (t<10)
				mu = max(min_mu,Mu*(1-alpha*max(0.0,Gamma_max*exp(-k_d*t/2)-0.25)));
			else
				mu = tan(PI/4); 
		}
		else if (fric_type=="Constant")
			if (t<10)
				mu = Mu;
			else
				mu = tan(1.5*Static_Delta* PI / 180);
			
		else
		{
			cout<<"No friction type selected. Using the constant friction angle"<<endl;
			mu = Mu;
		}

		if (solver == "LXF first order")
			LXF (w,wl,wr,dt); 

		else if (solver=="Predictor Corrector")
		{
			Predictor(w, wl, wr, dt);	
			Corrector(w, wl, wr, w_init, dt); 
		}

		else if (solver == "Runge-Kutta 3")
	    	RK3(w, wl, wr, dt);

		else
		{
			cout<<"No valid solver selected. Ending simulation"<<endl;
			return ;
		}
		
		Shed(w, Ang_Shed); 
		Edge(w,wl,wr);

		Time_step(wl,wr,dt,t,timesteps);

		if (restart)
			break;
		sum1=ang_sum;
		ang_sum=0;
		sum=0;
		for (int i=2;i<res-2;i++)
		{
			sum += w[i].J*w[i].h*(w[i].u*w[i].u+w[i].v*w[i].v)/2*dx;//+Jinertia1_reg(w[i])*omega)*dx;
			ang_sum += Ang_mom_reg(w[i])*dx;
		}

		if (t>check_t)
		{	
			myfile0<<std::setprecision(18)<<t<<"  "<<sum<<"  "<<dt<<"  "<<mass_shed<<endl;

			if (abs(sum)<1E-7)
			{	
				break;
			}	
			check_t+=1;
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
	myfile0<<"-------------------------------------------------------------------------"<<endl;
	cout << "Simulation ended ----------------------------------------------------------"<<endl;
	myfile.close();
	myfile0.close();
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
	if (dt<1E-8)
	{	
		cout << " The timestep becomes "<<dt<< endl;
		restart = true;
		limiter = false;
	}
}
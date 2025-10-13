#include "gauravlib.h"

void LXF(vector<CV>& w,  vector<CV>& wl, vector<CV>& wr, double dt)
{  	
	BC(w,w);
	Edge(w,wl,wr);
	vector<CV> wtemp(w);
	for (int j = 2; j < res-2; j++)
	{	
		double el=dx/dt;//max(abs(Ax(wl[j],wr[j-1],"min")),abs(Ax(wl[j],wr[j-1],"max")));
		double er=dx/dt;//max(abs(Ax(wl[j+1],wr[j],"min")),abs(Ax(wl[j+1],wr[j],"max")));
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


		if (wtemp[j].psi>0 and delta>epsilon)
		{
		 	CV w_new (w[j]);
			
			w_new.Modify(w_new.p,w_new.q-friction.q * dt, w_new.r - friction.r * dt);

 		 	if(w_new.u*w[j].u<0)
				w_new = CV(w_new.h,min_u*sign(w[j].u),min_u*sign(w[j].u),w_new.v,w[j].b,w[j].db,w[j].ddb,w[j].g,w[j].x);

 			if(w_new.v*w[j].v<0)
			{ 
				w_new = CV(w_new.h,w_new.u,w_new.u,min_u*sign(w[j].v),w[j].b,w[j].db,w[j].ddb,w[j].g,w[j].x);  
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


		if (wtemp[j].psi>0 and delta>epsilon)
		{
		 	CV w_new (w[j]);
			
			w_new.Modify(w_new.p,w_new.q-friction.q * dt, w_new.r - friction.r * dt);

 		 	if(w_new.u*w[j].u<0)
				w_new = CV(w_new.h,min_u*sign(w[j].u),min_u*sign(w[j].u),w_new.v,w[j].b,w[j].db,w[j].ddb,w[j].g,w[j].x);

 			if(w_new.v*w[j].v<0)
			{ 
				w_new = CV(w_new.h,w_new.u,w_new.u,min_u*sign(w[j].v),w[j].b,w[j].db,w[j].ddb,w[j].g,w[j].x);  
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

		
		FS friction = Friction(w[j]);

		if (wtemp[j].psi>0 and delta>epsilon)
		{
		 	CV w_new (w[j]);
			
			w_new.Modify(w_new.p,w_new.q-friction.q * dt, w_new.r - friction.r * dt);

 		 	  if(w_new.u*w[j].u<0)
				w_new = CV(w_new.h,min_u*sign(w[j].u),min_u*sign(w[j].u),w_new.v,w[j].b,w[j].db,w[j].ddb,w[j].g,w[j].x);

 			 if(w_new.v*w[j].v<0)
				w_new = CV(w_new.h,w_new.u,w_new.u,min_u*sign(w[j].v),w[j].b,w[j].db,w[j].ddb,w[j].g,w[j].x);  
 
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

		FS friction = Friction(w[j]);

		if (wtemp[j].psi>0 and delta>epsilon)
		{
			CV w_new (w[j]);
			
			w_new.Modify(w_new.p,w_new.q-(1-weight)*friction.q * dt, w_new.r - (1 - weight) *  friction.r * dt);

		 	if(w_new.u*w[j].u<0)
				w_new = CV(w_new.h,min_u*sign(w[j].u),min_u*sign(w[j].u),w_new.v,w[j].b,w[j].db,w[j].ddb,w[j].g,w[j].x);

 			 if(w_new.v*w[j].v<0)
				w_new = CV(w_new.h,w_new.u,w_new.u,min_u*sign(w[j].v),w[j].b,w[j].db,w[j].ddb,w[j].g,w[j].x);  

			w[j]=w_new; 
		
		} 

	}	
}

void March (vector<CV>& w, vector<CV>& wl, vector<CV>& wr, double& Ang_Shed)
{	
	double dt = dx / 8;
	static int timesteps=0;
	double sum1=0,sum=0,max_sum=0;
	
	
	if(!filesystem::exists(verbose_dir))
		filesystem::create_directory(verbose_dir);
	static double t=0;
	double check_t= 1;
	fs::path base_path = Output_folder;
	fs::path file_name= string("log.txt");
	fs::path full_path = base_path / file_name;
	string file1=	full_path.string();
	base_path = verbose_dir;
	while(t<finalt)
	{
		if(timesteps%dump==0 && (verbose=="Yes"|| verbose=="yes" ))
		{	
			
			fs::path file_name= string("field_")+to_string(int(timesteps/dump))+string(".csv");
			fs::path full_path = base_path / file_name;
			string file2=	full_path.string();
			Write_data(w, file2);
		}
		
		vector<CV> w_init(w);
		if (fric_type!="constant" && fric_type!="Constant")
		{
			if (t<seismic_time)
				delta= 1;
			else if (t<10)
				delta = 24*(1-exp(-k_d*(t-seismic_time)))+1;//change make it 30 again
			else
				delta =60;
		}

		LXF (w,wl,wr,dt);

	//	Predictor(w, wl, wr, dt);	
		
	//	Corrector(w, wl, wr, w_init, dt); 
		
		Shed(w, Ang_Shed); 

		Time_step(wl,wr,dt,t,timesteps);

		sum1=sum;
		sum=0;
		for (int i=2;i<res-2;i++)
			sum+=(Ang_mom_reg(w[i])*dx);//+Jinertia1_reg(w[i])*omega)*dx;
			
		if(abs(sum)>max_sum)
			max_sum =abs(sum);
		if (t>check_t)
		{	
			if (abs(sum)<pow(epsilon,2))
			{	
				break;
			}
			
			check_t+=1;
		}
		
	
		//std::cout<<std::setprecision(18)<<t<<"  "<<sum<<"  "<<delta<<endl;
	}

	 file_name= string("field_")+to_string(int(timesteps/dump))+string(".csv");
	 full_path = base_path / file_name;
	 string file2=	full_path.string();
	Write_data(w, file2);

	
	
	ofstream myfile(file1,std::ofstream::app);
	Shed(w, Ang_Shed);
	myfile<<"Simulation ran for time --> " <<t<<endl;
	myfile<<"Residual Angular Momentum --> " <<sum<<endl;
	myfile<<"Rate of change of Angular Momentum -->"<<(sum-sum1)/dt<<endl;
	//Ang_Shed+=sum;
	
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
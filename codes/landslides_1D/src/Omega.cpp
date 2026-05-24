#include "gauravlib.h"

double Inertia(vector<CV>& w, int no)
{
	double sum1 = Jinertia1_reg(w[0])+Jinertia1_reg(w[res-1])+Jinertia1_ast(w[0])+Jinertia1_ast(w[res-1]); 

	double sum2 = Jinertia2_reg(w[0])+Jinertia2_reg(w[res-1])+Jinertia2_ast(w[0])+Jinertia2_ast(w[res-1]);

	for (int i=1;i<res-1;i++)
	{
		sum1 += 2.0*(Jinertia1_ast(w[i])+Jinertia1_reg(w[i]));
	
		sum2 += 2.0*(Jinertia2_ast(w[i])+Jinertia2_reg(w[i]));
	}

	if (no==1)
		return dx/2*sum1;
	else
		return dx/2*sum2;
}

void Shed(vector<CV>& w, double& Ang_Shed)
{
bool shed = false;
for ( int j = 2; j < res-2; j++)
{
 	if (w[j].psi <= epsilon*epsilon ) 
	{	
		Ang_Shed += (Ang_mom_reg(w[j])+Jinertia1_reg(w[j])*omega)*dx;
		mass_shed +=  2*PI*epsilon*pow(dia/2,3)*w[j].p*dx;
		w[j]=CV(  min_h,sign(w[j].u)*min_u,sign(w[j].v)*min_u,w[j].b,w[j].db,w[j].ddb,w[j].g,w[j].x  );
		//w[j].psi = epsilon*epsilon;
		shed = true;
	} 

}
//if (shed && Ang_Shed>0)
//	cout<<"Total angular momentum shed till now " << Ang_Shed <<"      "<<past_time<<endl;

}

double Reg_Inertia(vector<CV>& w)
{
	double sum1 = Jinertia1_reg(w[0])+Jinertia1_reg(w[res-1]); 


	for (int i=1;i<res-1;i++)
	{
		sum1 += 2.0*(Jinertia1_reg(w[i]));
	
	}

	return dx/2*sum1;

}


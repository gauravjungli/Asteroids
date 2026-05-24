#include "gauravlib.h"

double Inertia(vector<CV>& w)
{
	double sum1 = 0;

for (int i =0;i<rows;i++)
{
	for ( int j = 0; j < cols; j++)
	{
		int i0 = index(i,j);
		sum1 += (Jinertia1_ast(w[i0])+Jinertia1_reg(w[i0]));
	
	}
}
	return dx*dy*sum1;

}

void Shed(vector<CV>& w, double& Ang_Shed)
{
bool shed = false;
for (int i =2;i<rows-2;i++)
{
	for ( int j = 0; j < cols; j++)
	{
	int i0 = index(i,j);

 	if (w[i0].psi <= epsilon*epsilon ) 
	{	
		Ang_Shed += (Ang_mom_reg(w[i0])+Jinertia1_reg(w[i0])*omega)*dx*dy;
		mass_shed +=  epsilon*pow(dia/2,3)*w[i0].p*dx*dy;
		w[i0]=CV(  min_h,sign(w[i0].u)*min_u,sign(w[i0].v)*min_u,w[i0].b,w[i0].Ra,w[i0].dR,w[i0].ddR,w[i0].g,w[i0].x,w[i0].y  );

		shed = true;
	} 

}
}
if (shed && Ang_Shed>0)
	cout<<"Total angular momentum shed till now " << Ang_Shed <<"      "<<past_time<<endl;

}



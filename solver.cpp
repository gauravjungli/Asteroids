#include "gauravlib.h"
void cfl(vector<double>& x,double& dt)
{

	double maxspeed = 0.00000001;
	// to evaluate dt from maximum speeds (CFL condition)
	for (int j = 1; j < res - 1; j++)
	{
		double eig = ax(j);
		if (maxspeed < abs(eig))
		maxspeed = abs(eig);
	}

	dt = min(dx/8,dx/8/maxspeed);
}

double maxeigenx(CV wl,CV wr,int j, int pm)
// pm is 1,2 for plus and minus of q vector
{
	double e1, e2, e3; //eigenvalues
	double root;
	if (pm == 1)
	{
		root = 0;//change for the sphere
		if (root < 0)
		{//cout<<"non-hperbolic";
			root = 0;
		}
		else
		{
			root = 0;//chnage for the sphere
		}
		e1 = abs(wr.q / wr.p);
		e2 = abs(e1 - root);
		e3 = abs(e1 + root);

		return (std::max(e1, std::max(e2, e3)));
	}
	else
	{
		double tau = 0;//change for the sphere
		root = 0;//change for the sphere
		//  cout<<phi<<"  "<<jj<<"  "<<root<<endl;
		if (root < 0)
		{//cout<<"non-hperbolic";
			root = 0;
		}
		else
		{
			root = sqrt(root) / (pow(wl.p, 2) * pow(tau, 3));
		}

		e1 = abs(wl.q / wl.p);
		e2 = abs(e1 - root);
		e3 = abs(e1 + root);
		return (std::max(e1, std::max(e2, e3)));
	}
}

double ax(int jj)
{    	//	cout<<maxeigenx(jj, 1)<<"   "<< maxeigenx(jj, 2)<<"   "<<"max2"<<endl;

	return (std::max(maxeigenx(jj, 1), maxeigenx(jj, 2))); // 1,2 for plus,minus
}

double Hx(int j, int no, CV wl, CV wr) // this xx, yy corresponds to boundaries of cells
// no is 1,2,3 for denoting p,q,r
{

	

	double fp = flux(wr, j); // first 1 is for f and second 1 is for plus
	double fm = flux(wl ,j);

	if (no == 1)
	{  //cout<<Qp_p<<"    "<<Qm_p<<"     "<<"   "<<jj<<endl;
		return ((fp + fm) / 2 - ax(j - 1) * (wr.p - wl.p) / 2);
	}
	else if (no == 2)
	{
		return ((fp + fm) / 2 - ax(j - 1) * (wr.q - wl.q) / 2);
	}
	else
	{
		return ((fp + fm) / 2 - ax(j - 1) * (wr.r - wl.r) / 2);
	}

}

void cell_boundary(vector<CV>& wl, vector<CV>& wr, vector<CV>& wtemp)
{for (int j=0; j<res;j++)
	{
	double Qp_p = wtemp[j].p - dx * derivative(j, 1,wtemp) / 2;
	double Qm_p = wtemp[j - 1].p + dx * derivative(j - 1, 1,wtemp) / 2;
	double Qp_q = wtemp[j].p - dx * derivative(j, 2,wtemp) / 2;
	double Qm_q = wtemp[j - 1].q + dx * derivative(j - 1, 2,wtemp) / 2;
	double Qp_r = wtemp[j].r - dx * derivative(j, 3,wtemp) / 2;
	double Qm_r = wtemp[j - 1].r + dx * derivative(j - 1, 3,wtemp) / 2;
	CV temp1(Qm_p,Qm_q/Qm_p,Qm_r/Qm_p);
	wl.push_back(temp1);
	CV temp2(Qp_p,Qp_q/Qp_p,Qp_r/Qp_p);
	wr.push_back(temp2);
	}
}



#include "gauravlib.h"

double minmod(double a, double b, double c) // calculate minmod
{
	double temps;
	//temps = std::min(abs(a), std::min(abs(b), abs(c)));

	//return temps;

	if (a < 0 && b < 0 && c < 0)
		return std::max(a, std::max(b, c));
	else if (a > 0 && b > 0 && c > 0)
		return std::min(a, std::min(b, c));
	else
		return 0;
}

double derivative( int j, int no, vector<CV>& wtemp)
{
	if (j == res + 3)
		cout << "Daya kuch to gadbad hain" << endl;
	if (no == 1)
	{
		return minmod(theta * (wtemp[j].p - wtemp[j - 1].p) / dx, (wtemp[j + 1].p - wtemp[j - 1].p) / 2 / dx, theta * (wtemp[j + 1].p - wtemp[j].p) / dx);
	}
	else if (no == 2)
	{
		return minmod(theta * (wtemp[j].q - wtemp[j - 1].q) / dx, (wtemp[j + 1].q - wtemp[j - 1].q) / 2 / dx, theta * (wtemp[j + 1].q - wtemp[j].q) / dx);
	}
	else
	{
		return minmod(theta * (wtemp[j].r - wtemp[j - 1].r) / dx, (wtemp[j + 1].r - wtemp[j - 1].r) / 2 / dx, theta * (wtemp[j + 1].r - wtemp[j].r) / dx);
	}

}

void reconstruction(vector<CV>& w, vector<double>& x, double Om)//change for the sphere
{
	for (int j = 0; j < res; j++)
	{
			w[j].h = w[j].p;
			w[j].u = w[j].q / w[j].p;
			w[j].v = w[j].r;
	}
}
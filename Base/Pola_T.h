//====================LIBRARY===================================================
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//====================THE POLARIZABILITY FUNCTION===============================
double step(double q1, double k1)
{
	if (fabs(q1) <= fabs(2. * k1))
		return 0.;
	else
		return 1. / (q1 * sqrt(pow(q1, 2.) - 4. * pow(k1, 2.)));
}
double fermi_E(double k1, double t1)
{
	double hs;
	hs = (k1 * k1 - 1.) / (t1 + 0.001);
	if ((hs <= 100.))
		return 1. / (1. + exp(hs));
	else
		return 0.;
}
double fermi_E_2mu(double k1, double t1)
{
	double hs;
	hs = (k1 * k1 + 1.) / (t1 + 0.001);
	if ((hs <= 100.))
		return 1. / (1. + exp(hs));
	else
		return 0.;
}
double fpolar_k(double k1, double q1, double t1)
{ // Ham fpolar chinh la ham Pi(q,T)/N0
	double kq;
	kq = (1. / pow(k1, 3.)) * (sqrt(4. * pow(k1, 4.) + pow(q1, 4.)) - pow(k1, 2.) - fabs(pow(k1, 2.) - pow(q1, 2)) + (fermi_E(k1, t1) + fermi_E_2mu(k1, t1)) * (2. * pow(k1, 2.) - sqrt(4. * pow(k1, 4.) + pow(q1, 4.)) + pow(2 * pow(k1, 2.) - pow(q1, 2.), 2.) * step(q1, k1)));
	return kq;
}
void gaulegf2(double x1, double x2, double x[], double w[], int n)
{
	int i, j, m;
	double eps = 3.0E-14;
	double p1, p2, p3, pp, xl, xm, z, z1;

	m = (n + 1) / 2;
	xm = 0.5 * (x2 + x1);
	xl = 0.5 * (x2 - x1);
	for (i = 1; i <= m; i++)
	{
		z = cos(3.141592654 * ((double)i - 0.25) / ((double)n + 0.5));
		while (1)
		{
			p1 = 1.0;
			p2 = 0.0;
			for (j = 1; j <= n; j++)
			{
				p3 = p2;
				p2 = p1;
				p1 = ((2.0 * (double)j - 1.0) * z * p2 - ((double)j - 1.0) * p3) /
					 (double)j;
			}
			pp = (double)n * (z * p1 - p2) / (z * z - 1.0);
			z1 = z;
			z = z1 - p1 / pp;
			if (abs(z - z1) <= eps)
				break;
		}
		x[i] = xm - xl * z;
		x[n + 1 - i] = xm + xl * z;
		w[i] = 2.0 * xl / ((1.0 - z * z) * pp * pp);
		w[n + 1 - i] = w[i];
	}
} /* end gaulegf2 */
double polar(double a, double b, double q1, double t1)
{
	int i, j;
	double sum;
	double k1[2000], w[2000];
	i = 5.;
	gaulegf2(a, b, k1, w, i);
	sum = 0.0;
	for (j = 1; j <= i; j++)
	{
		sum = sum + w[j] * fpolar_k(k1[j], q1, t1);
	}
	return sum;
}
double polarinfi(double q1, double t1)
{
	double kt, kq, gt1, gt2, gt, a, b;
	int i;
	if (q1 <= 0.001)
		kq = -1.;
	else
	{
		i = 1;
		a = 0.00001;
		b = .005;
		gt1 = polar(a, b, q1, t1);
		a = b;
		b = b + .005;
		gt2 = polar(a, b, q1, t1);
		gt = gt1 + gt2;
		kt = 1.;
		while (kt > .00000001 && i < 10000.)
		{
			a = b;
			b = b + .005;
			gt2 = polar(a, b, q1, t1);
			gt = gt + gt2;
			kt = fabs(gt2 / gt);
			i += 1;
		}
		kq = -1 * (gt + gt2);
	}
	return kq;
}
//==============================================================================

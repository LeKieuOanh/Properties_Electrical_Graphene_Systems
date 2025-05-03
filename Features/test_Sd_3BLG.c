//====================LIBRARY===================================================
#include <stdio.h>
#include <math.h>
//====================CONSTANTS(CGS)============================================
#define pi					3.141592654
#define kB					1.3806503*1e-16
#define hbar				1.054571628*1e-27
#define msao				0.033*9.10938188*1e-28
#define msao1				3.326573713*1e28 //1 chia msao
#define elec				4.80320427*1e-10
//====================ROUTINE===================================================
double PI(double q1, double t1);
double step(double q1, double k1);
double fermi_E(double k1, double t1);
double fermi_E_2mu(double k1, double t1);
double fpolar_k(double k1, double q1, double t1);
double polar(double a, double b, double q1, double t1);
double polarinfi(double q1, double t1);
double df(double E1, double t1);

double heso(double q1);
double V11(double q1, double d);
double V12(double q1, double d);
double V22(double q1, double d);
double epsilon(double q1, double d, double t1);
double W11(double q1, double d, double t1);
double ftau_rev(double q1, double k1, double d, double t1);
double tau_rev(double d, double t1, double E1);
double sigma1(double d, double t1);
//====================MAIN======================================================
double kappa1 = 1;
double kappa2 = 1;
double kappa3 = 1;

double Xi1 = 1;
double Xi2 = 1;

double gsgv = 4;
double n	= 1e12;		//n mat do hat tai
double ni	= 1e11;		//mat do hat tai tren lop thu i
double s	= 1;		// The hoa hoc /mu
//------------------------------------------------------------------------------
int main(){
	double E1, d, t1, i1, i2, i, Sd; //x la E/EF
	d = 1; //nm
	FILE *output;
	output=fopen ("2BLG_L1_Sd.dat","w+");
		for (t1 = 0.01; t1 <= 0.1; t1 += 0.01){
 			printf("%.2f\n	", t1);						//Hien thi tren terminal
 			fprintf(output, "%f ", t1);
 			i = i1 = i2 = Sd = 0;
 			for (E1 = 1; E1 <= 30; E1 += 0.5){
 				printf("%2.1f ", E1);					//Hien thi tren terminal
 				i = df(E1, t1) * E1/tau_rev(d*1e-7, t1, E1);
				i2 = i2 + i;
				i1 = i1 + i*E1;
			}
			Sd = (i1/i2-s)/t1;							//kB/e, K^-1
			//Sd = kB*(i1/i2-s)/(elec*t1);				//V/K
			printf("%.20e\n", Sd);						//Hien thi tren terminal
 			fprintf(output,"%.20e", Sd);
			fprintf(output, "\n");
        }
        for (t1 = 0.2; t1 <= 0.9; t1 += 0.1){
 			printf("%.2f\n	", t1);						//Hien thi tren terminal
 			fprintf(output, "%f ", t1);
 			i = i1 = i2 = Sd = 0;
 			for (E1 = 1; E1 <= 30; E1 += 0.5){
 				printf("%2.1f ", E1);					//Hien thi tren terminal
 				i = df(E1, t1) * E1/tau_rev(d*1e-7, t1, E1);
				i2 = i2 + i;
				i1 = i1 + i*E1;
			}
			Sd = (i1/i2-s)/t1;							//kB/e, K^-1
			printf("%.20e\n", Sd);						//Hien thi tren terminal
 			fprintf(output,"%.20e", Sd);
			fprintf(output, "\n");
        }
	fclose(output);
	printf("\nfin\a");
	return 0;
}
//====================THE POLARIZABILITY FUNCTION===============================
void gaulegf(double x1, double x2, double x[], double w[], int n){
	int i, j, m;
	double eps = 3.0E-14;
	double p1, p2, p3, pp, xl, xm, z, z1;
	m = (n+1)/2;
	xm = 0.5*(x2+x1);
	xl = 0.5*(x2-x1);
	for (i = 1; i <= m; i++){
    z = cos(3.141592654*((double)i-0.25)/((double)n+0.5));
    while(1){
    	p1 = 1.0;
    	p2 = 0.0;
    	for (j = 1; j <= n; j++){
    		p3 = p2;
    		p2 = p1;
        	p1 = ((2.0*(double)j-1.0)*z*p2-((double)j-1.0)*p3)/
            (double)j;
    	}
    	pp = (double)n*(z*p1-p2)/(z*z-1.0);
    	z1 = z;
    	z = z1 - p1/pp;
    	if (fabs(z-z1) <= eps) break;
    }
    x[i] = xm - xl*z;
    x[n+1-i] = xm + xl*z;
    w[i] = 2.0*xl/((1.0-z*z)*pp*pp);
    w[n+1-i] = w[i];
	}
} //end gaulegf
double step(double q1, double k1){
	if	(fabs(q1) <= fabs(2.*k1))
		return 0.;
	else
		return 1./(q1*sqrt(pow(q1, 2.)-4.*pow(k1, 2.)));
}
double fermi_E(double k1,double t1){
	double hs, kq;
	hs = (k1*k1-1.)/(t1+0.001);
	if (hs <= 100.)
		return 1./(1.+exp(hs));
	else
		return 0.;
}
double fermi_E_2mu(double k1, double t1){
	double hs,kq;
	hs = (k1*k1+1.)/(t1+0.001);
	if (hs<=100.)
		return 1./(1.+exp(hs));
	else
		return 0.;
}
double fpolar_k(double k1, double q1, double t1){  // Ham fpolar chinh la ham Pi(q,T)/N0
	double kq;
	kq = (1./pow(k1, 3.))*(sqrt(4.*pow(k1, 4.)+pow(q1, 4.))-pow(k1, 2.)-fabs(pow(k1, 2.)-pow(q1, 2))+(fermi_E(k1, t1)+fermi_E_2mu(k1, t1))*(2.*pow(k1, 2.)-sqrt(4.*pow(k1, 4.)+pow(q1, 4.))+pow(2*pow(k1, 2.)-pow(q1, 2.),2.)*step(q1, k1)));
	return kq;
}
double polar(double a, double b, double q1, double t1){
	int i, j;
	double sum;
	double k1[2000], w[2000];
	i = 5.;
	gaulegf(a, b, k1, w, i);
	sum = 0.0;
	for (j = 1; j <= i; j++){
    	sum = sum + w[j]*fpolar_k(k1[j], q1, t1);
	}
	return sum;
}
double PI(double q1, double t1){ //polarinfi //Pi/D0
	double kt, kq, gt1, gt2, gt, a, b;
	int i;
	if (q1 <= 0.001)
    	kq = -1.;
	else{
		i = 1;
		a = .00001;
		b = .005;
			gt1 = polar(a, b, q1, t1);
		a = b;
		b = b+.005;
			gt2 = polar(a, b, q1, t1);
		gt = gt1+gt2;
		kt = 1.;
	while (kt>.00000001 && i<10000.){
		a = b;
		b = b+.005;
		gt2 = polar(a, b, q1, t1);
		gt = gt+gt2;
			kt = fabs(gt2/gt);
    	i += 1;
		}
		kq = -1*(gt+gt2);
	}
  return kq;
}
double df(double E1, double t1){// E.(-df/dE)
	double kq, hs;
	hs = (E1-s)/t1;
	if (hs <= 350.) //700
		kq = exp(hs)*E1 / (pow(1+exp(hs), 2)*t1);
	else
		kq = E1/(t1*exp(350.));
  return kq;
}
//==============================================================================
double heso(double q1){
	double kq;
	kq = (4*pi*pow(elec, 2)) / (q1*sqrt(4*pi*n/gsgv));
	return kq;
} //OK
double V11(double q1, double d){
	double kq;
	kq = heso(q1)*((kappa2 + kappa1*tanh(q1*sqrt(4*pi*n/gsgv)*d))/
		(kappa2*(kappa3 + kappa1) + (pow(kappa2, 2) + kappa3*kappa1)*tanh(q1*sqrt(4*pi*n/gsgv)*d))
		);
	return kq;
} //OK
double V22(double q1, double d){
	double kq;
	kq = heso(q1)*((kappa2 + kappa3*tanh(q1*sqrt(4*pi*n/gsgv)*d))/(
		 kappa2*(kappa3 + kappa1) + (pow(kappa2, 2) + kappa3*kappa1)*tanh(q1*sqrt(4*pi*n/gsgv)*d)
		)
	);
	return kq;
} //OK
double V12(double q1, double d){
	double kq;
	kq = heso(q1)*(kappa2/cosh(q1*sqrt(4*pi*n/gsgv)*d))/(
		kappa2*(kappa3 + kappa1) + (pow(kappa2, 2) + kappa3*kappa1)*tanh(q1*sqrt(4*pi*n/gsgv)*d)
		);
	return kq;
} //OK
double epsilon(double q1, double d, double t1){
	double kq, pi1, pi2;
	pi1 = pi2 = (PI(q1, t1) * gsgv*msao) / (2*pi*pow(hbar,2));
	kq = (1 - pi1*V11(q1, d))*(1 - pi2*V22(q1, d)) - pow(V12(q1, d), 2)*pi1*pi2;
	return kq;
} //OK
double W11(double q1, double d, double t1){
	double kq, pi1, pi2;
	pi1 = pi2 = (PI(q1, t1) * gsgv*msao) / (2*pi*pow(hbar,2));
	kq = Xi1*(
			V11(q1, d) + (pow(V12(q1, d), 2) - V11(q1, d)*V22(q1, d))*pi2
		) / epsilon(q1, d, t1);
	return kq;
} //OK
//==============================================================================
double ftau_rev(double q1, double k1, double d, double t1){ //ham duoi dau tich phan cua tau
	double kq;
	kq = (q1*q1*	pow((1-2*pow(q1/(2*k1), 2)), 2)	/	(sqrt(4*k1*k1-q1*q1)))	*	pow(W11(q1, d, t1), 2); //tau cho T khac 0: k = k1.kF = sqrt(E1)*kF //chua co he so kF: kF.ftau_rev
	return kq;
} //OK
double tau_rev(double d, double t1, double E1){ // ham tich phan chua co he so phia truoc
	int i, j, ka;
	double sum;
	double z[1000], w[1000];
	i = 100.;
	gaulegf(0, 2*sqrt(E1), z, w, i);
	sum = 0.;
	for (j = 1; j <= i; j++){
    	sum = sum + w[j] * ftau_rev(z[j], sqrt(E1), d, t1); //k1=sqrt(E1) dq=kF.dq1; kF^2.ftau_rev.dq1
    }
	return sum;
}

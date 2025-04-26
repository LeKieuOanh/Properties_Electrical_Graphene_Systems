

//=============THE LIBRARIES DECLARATION==========
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>


//===========CONSTANT============================
// CGS - units
#define pi        3.141592653589793
#define hbar      1.054571628*1e-27 //cm2g/s
#define msao      0.033*9.10938215*1e-28 //the bare electron mass : gams
#define elec      4.80320427*1e-10 // charge electron 
#define n2        1e12   // the number electron at second layers

// epsilon1 = 1



//============FUNCTIONS=========================
double step(double q1);
double f0(double q1);
double g0(double q1);
double PI(double q1);
double V11(double q1, double epsilon2, double epsilon3, double d, double n1);
double V22(double q1, double epsilon2, double epsilon3, double d, double n1);
double V12(double q1, double epsilon2, double epsilon3, double d, double n1);
double epsilon(double q1, double d, double epsilon2, double epsilon3, double n1);
double W11(double q1, double d, double epsilon2, double epsilon3, double n1);
double S(double q1,double ni1, double n1, double r0);
double ftau(double q1,double ni1, double d, double epsilon2, double epsilon3,double n1, double r0);
void gaulegf(double x1, double x2, double x[], double w[], int n);
double tau(double ni1,double d, double epsilon2, double epsilon3, double n1, double r0);
void Legendre(int n, double x, double *d,double *f);
double J(double x);
//================EXECUTION======================

//q= q1*kF , k = k1 *kF, T = 0K => k1 = 1,  kF = sqrt(pi * n )

double f0(double q1)
{
    return ((2. + q1*q1) * sqrt(q1*q1 - 4.) / (2.*q1) + log((q1 - sqrt(q1*q1 - 4.)) / (q1 + sqrt(q1*q1 - 4.)) ) );
}

double g0(double q1)
{
    return (( sqrt(4. + pow(q1,4)) / 2.) - log((1. + sqrt( 1. + pow(q1,4)/4.)) / 2.)) ;
}

  
  

double PI(double q1)
{
    double N0  = ( 2. * msao ) / (pi * pow(hbar,2));
    if (q1 <= 2.)
        return -N0 * ( g0(q1));
    else 
       return -N0 * ( g0(q1) - f0(q1));
    
 }


double V11(double q1, double epsilon2, double epsilon3, double d, double n1)
{
    double coeff = 4*pi*pow(elec,2)/ (q1 * sqrt(pi * n1)) ;
    double deno = (epsilon2*(1+epsilon3) + ( pow(epsilon2,2) + epsilon3)*tanh(q1 * sqrt(pi * n1) * d));
    return coeff * (epsilon2 + epsilon3 * tanh(q1 * sqrt(pi * n1)*d)) / deno;
    

}

double V22(double q1, double epsilon2, double epsilon3, double d, double n1)
{
    double coeff = 4*pi*pow(elec,2)/ (q1 * sqrt(pi * n1)) ;
    double deno = (epsilon2*(1+epsilon3) + (pow(epsilon2,2) + epsilon3)*tanh(q1 * sqrt(pi * n1)*d));
    return  coeff * (epsilon2 + 1.0*tanh(q1*sqrt(pi * n1)*d)) / deno ;
    
}

double V12(double q1, double epsilon2, double epsilon3, double d, double n1)
{
    double coeff = (4*pi*pow(elec,2))/ (q1 * sqrt(pi * n1)) ;
    double deno = (epsilon2*(1+epsilon3) + (pow(epsilon2,2) + epsilon3)*tanh(q1 * sqrt(pi * n1)*d));
    return coeff * epsilon2 / ( cosh(q1 * sqrt(pi * n1)*d) * deno ) ;
    
}

double epsilon(double q1, double d, double epsilon2, double epsilon3, double n1)
{
    return ((1 - PI(q1)* V11(q1,epsilon2,epsilon3,d,n1) ) *  ( 1 - PI(q1*sqrt(n1/n2))*V22(q1,epsilon2,epsilon3,d,n1)) - (V12(q1,epsilon2,epsilon3,d,n1)*V12(q1,epsilon2,epsilon3,d,n1) * PI(q1) * PI(q1*sqrt(n1/n2))));
}

double W11(double q1, double d, double epsilon2, double epsilon3, double n1)
{
    
    return ( V11(q1,epsilon2,epsilon3,d,n1) + ( pow(V12(q1,epsilon2,epsilon3,d,n1),2) - V11(q1,epsilon2,epsilon3,d,n1) * V22(q1,epsilon2,epsilon3,d,n1))*PI(q1*sqrt(n1/n2))  / (epsilon(q1,d,epsilon2,epsilon3,n1)));
    
}


double prod(int n)
{
    double Pi=1;
    for (int i = 1; i <= n; i++)
    {
        Pi = Pi * ((double)i) ;
    }
    return Pi;
}

double gam(int n)
{
    int i;
    double prod = 1;
    for (i = n ; i > 1 ; i--)
    {
        prod = prod * ((double)i - 1) ;
    }
    return prod;
}

double J(double x)
{
    double s = 0, s1 = 1, eps = 1e-14;
    int i = 0;
    
        while (fabs(s1) > fabs(s)*eps)
        {
            s1 = (double)pow(-1,i) * pow(x/2,2*i +1) / ( gam(i+2)* prod(i) )  ;
            s = s + s1 ;
            i = i + 1;
        }
    return s ;
}


double S(double q1,double ni1, double n1, double r0)
{
    return 1. - 2 * pi * ni1 * r0 * J(q1*r0*sqrt(pi*n1)) / (q1*sqrt(pi * n1))  ;
}

void gaulegf(double x1, double x2, double x[], double w[], int n)
{
  int i, j, m;
  double eps = 3.0E-14;
  double p1, p2, p3, pp, xl, xm, z, z1;
  
  m = (n+1)/2;
  xm = 0.5*(x2+x1);
  xl = 0.5*(x2-x1);
  for(i=1; i<=m; i++)
  {
    z = cos(3.141592654*((double)i-0.25)/((double)n+0.5));
    while(1)
    {
      p1 = 1.0;
      p2 = 0.0;
      for(j=1; j<=n; j++)
      {
        p3 = p2;
        p2 = p1;
        p1 = ((2.0*(double)j-1.0)*z*p2-((double)j-1.0)*p3)/
              (double)j;
      }
      pp = (double)n*(z*p1-p2)/(z*z-1.0);
      z1 = z;
      z = z1 - p1/pp;
      if(fabs(z-z1) <= eps) break;
    }
    x[i] = xm - xl*z;
    x[n+1-i] = xm + xl*z;
    w[i] = 2.0*xl/((1.0-z*z)*pp*pp);
    w[n+1-i] = w[i];
  }
}

double ftau(double q1,double ni1,double d, double epsilon2, double epsilon3, double n1, double r0)
{
   
    double EF = pi*n1*hbar*hbar/(2*msao);
    return ( pi * n1 * ni1 * q1*q1 * pow(1. - q1*q1/2.,2)* pow(fabs(W11(q1,d,epsilon2,epsilon3,n1)),2) * S(q1, ni1, n1,r0)) / (sqrt(4 - q1*q1) * 2 * pi * hbar * EF) ;
}

double tau(double ni1,double d, double epsilon2, double epsilon3, double n1, double r0)
{
  int i, j;
  double result;
  double x[9000], w[9000];
    i=100.; 
    gaulegf(0., 2., x, w, i);
    result = 0.;
    double EF = pi*n1*hbar*hbar/(2*msao);
    for(j=1; j<=i; j++)
    {
        result +=  w[j] *  ftau(x[j],ni1,d,epsilon2,epsilon3,n1,r0); 
    }
    
  return result; 
}
double sigma(double ni1,double d,  double epsilon2, double epsilon3, double n1, double r0)
{
    double N0  = ( 2. * msao ) / (pi * pow(hbar,2));
    double EF = pi*n1*hbar*hbar / (2*msao);
    return (N0 * EF * hbar * 2. * pi ) / ( msao * tau(ni1,d,epsilon2,epsilon3,n1,r0)) ;            
}
int main()
{
 
    double  a0, r0,l;
    double ni1, n1;
    
    
    double N0  = ( 2. * msao ) / (pi * pow(hbar,2));
    a0 = 4.92 * 1e-8 ;
    l = 1e-7;
    double r02_list[5] = {0., 5.*a0, 7*a0, 8*a0, 10.*a0};
    double d_list[3] = {1*l, 3*l, 100*l};
    double d, epsilon2, epsilon3;
    epsilon2 = 4.0; 
    epsilon3 = 12.53;   
    FILE *output;
    for (int k =0; k < 3; k++)
    {
        d = d_list[k];
       for (int i= 0; i<=4; i++)
    {
        ni1 = 0.95 * 1e12;
        r0 = r02_list[i];
        char s1[] = "./results/BLG_systems_0K/test/d=";
        char s3[10];
        gcvt(d/l,10,s3);
        char s4[] = "nm/0.95e12/sigmar0=";
        char s2[20];
        gcvt(r0/a0,20,s2);
        char s[100]="";
        strcat(s,s1);
        strcat(s,s3);
        strcat(s,s4);
        strcat(s,s2);
        strcat(s,".txt");
    output = fopen(s,"w+");

         for (n1 = 0; n1 <= 5e12; n1 += 0.02e12)
    {
        fprintf(output,"%f \t %f \t %f\n ",n1,sigma(ni1,d,epsilon2,epsilon3,n1,r0),tau(ni1,d,epsilon2,epsilon3,n1,r0));
    }
    }
    }
    
    
    fclose(output);
    return 0;
}
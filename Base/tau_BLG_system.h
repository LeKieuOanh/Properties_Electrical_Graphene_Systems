#include <stdio.h>
#include <math.h>

#include"W_2layer_2BLG.h"
#include"W_2layer_BLG_MLG.h"
#include"W_2layer_BLG_2DEG.h"
#include"S_func.h"
#include"gaulegf.h"


// 2BLG 

double ftau_2BLG(double q1,double ni1,double d, double e1, double e2, double e3, double n1, double n2, double r0)
{
   
    double EF = pi*n1*hbar*hbar / (2*msao);
    return ( pi * n1 * ni1 * q1*q1 * pow(1. - q1*q1/2.,2)* pow(fabs(W11_2BLG(q1,d,e1,e2,e3,n1,n2)),2) * S(q1, ni1, n1,r0)) / (sqrt(4 - q1*q1) * 2 * pi * hbar * EF) ;
}

double tau_2BLG(double ni1,double d, double e1, double e2, double e3, double n1, double n2, double r0)
{
  int i, j;
  double result;
  double x[9000], w[9000];
    i=100.; 
    gaulegf(0., 2., x, w, i);
    result = 0.;
    double EF = pi*n1*hbar*hbar / (2*msao);
    for(j=1; j<=i; j++)
    {
        result +=  w[j] *  ftau_2BLG(x[j],ni1,d,e1,e2,e3,n1,n2,r0); 
    }
    
  return result; 
}

// BLG-MLG

double ftau_BLG_MLG(double q1,double ni1,double d, double e1, double e2, double e3, double n1, double n2, double r0)
{
   
    double EF = pi*n1*hbar*hbar / (2*msao);
    return ( pi * n1 * ni1 * q1*q1 * pow(1. - q1*q1/2.,2)* pow(fabs(W11_BLG_MLG(q1,d,e1,e2,e3,n1,n2)),2) * S(q1, ni1, n1,r0)) / (sqrt(4 - q1*q1) * 2 * pi * hbar * EF) ;
}

double tau_BLG_MLG(double ni1,double d, double e1, double e2, double e3, double n1, double n2, double r0)
{
  int i, j;
  double result;
  double x[9000], w[9000];
    i=100.; 
    gaulegf(0., 2., x, w, i);
    result = 0.;
    double EF = pi*n1*hbar*hbar / (2*msao);
    for(j=1; j<=i; j++)
    {
        result +=  w[j] *  ftau_BLG_MLG(x[j],ni1,d,e1,e2,e3,n1,n2,r0); 
    }
    
  return result; 
}

// BLG-2DEG

double ftau_BLG_2DEG(double q1,double ni1,double d, double e1, double e2, double e3, double n1, double n2, double r0, double m2DEG)
{
   
    double EF = pi*n1*hbar*hbar / (2*msao);
    return ( pi * n1 * ni1 * q1*q1 * pow(1. - q1*q1/2.,2)* pow(fabs(W11_BLG_2DEG(q1,d,e1,e2,e3,n1,n2,m2DEG)),2) * S(q1, ni1, n1,r0)) / (sqrt(4 - q1*q1) * 2 * pi * hbar * EF) ;
}

double tau_BLG_2DEG(double ni1,double d, double e1, double e2, double e3, double n1, double n2, double r0, double m2DEG)
{
  int i, j;
  double result;
  double x[9000], w[9000];
    i=100.; 
    gaulegf(0., 2., x, w, i);
    result = 0.;
    double EF = pi*n1*hbar*hbar / (2*msao);
    for(j=1; j<=i; j++)
    {
        result +=  w[j] *  ftau_BLG_2DEG(x[j],ni1,d,e1,e2,e3,n1,n2,r0,m2DEG); 
    }
    
  return result; 
}
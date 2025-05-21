#include <stdio.h>
#include <math.h>

#include"W_layers.h"

#include"S_func.h"
#include"gaulegf.h"


// 2BLG - 1spacer

double ftau_2BLG_1spacer(double q1,double ni1,double d, double e1, double e2, double e3, double n1, double n2, double r0)
{
   
    double EF = pi*n1*hbar*hbar / (2*msao);
    return ( pi * n1 * ni1 * q1*q1 * pow(1. - q1*q1/2.,2)* pow(fabs(W11_2BLG(q1,d,e1,e2,e3,n1,n2)),2) * S(q1, ni1, n1,r0)) / (sqrt(4 - q1*q1) * 2 * pi * hbar * EF) ;
}

double tau_2BLG_1spacer(double ni1,double d, double e1, double e2, double e3, double n1, double n2, double r0)
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
        result +=  w[j] *  ftau_2BLG_1spacer(x[j],ni1,d,e1,e2,e3,n1,n2,r0); 
    }
    
  return result; 
}

// BLG-MLG - 1spacer

double ftau_BLG_MLG_1spacer(double q1,double ni1,double d, double e1, double e2, double e3, double n1, double n2, double r0)
{
   
    double EF = pi*n1*hbar*hbar / (2*msao);
    return ( pi * n1 * ni1 * q1*q1 * pow(1. - q1*q1/2.,2)* pow(fabs(W11_BLG_MLG_1spacer(q1,d,e1,e2,e3,n1,n2)),2) * S(q1, ni1, n1,r0)) / (sqrt(4 - q1*q1) * 2 * pi * hbar * EF) ;
}

double tau_BLG_MLG_1spacer(double ni1,double d, double e1, double e2, double e3, double n1, double n2, double r0)
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
        result +=  w[j] *  ftau_BLG_MLG_1spacer(x[j],ni1,d,e1,e2,e3,n1,n2,r0); 
    }
    
  return result; 
}

// BLG-2DEG - 1spacer

double ftau_BLG_2DEG_1spacer(double q1,double ni1,double d, double e1, double e2, double e3, double n1, double n2, double r0, double m2DEG)
{
   
    double EF = pi*n1*hbar*hbar / (2*msao);
    return ( pi * n1 * ni1 * q1*q1 * pow(1. - q1*q1/2.,2)* pow(fabs(W11_BLG_2DEG_1spacer(q1,d,e1,e2,e3,n1,n2,m2DEG)),2) * S(q1, ni1, n1,r0)) / (sqrt(4 - q1*q1) * 2 * pi * hbar * EF) ;
}

double tau_BLG_2DEG_1spacer(double ni1,double d, double e1, double e2, double e3, double n1, double n2, double r0, double m2DEG)
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
        result +=  w[j] *  ftau_BLG_2DEG_1spacer(x[j],ni1,d,e1,e2,e3,n1,n2,r0,m2DEG); 
    }
    
  return result; 
}


// 2BLG - 2spacer

double ftau_2BLG_2spacer(double q1, double w, double n1, double n2,double e1, double e21 , double e22, double e3, double d, double r0, double ni1)
{
   
    double EF = pi*n1*hbar*hbar / (2*msao);
    return ( pi * n1 * ni1 * q1*q1 * pow(1. - q1*q1/2.,2)* pow(fabs(W11_2layer_2spacer_2BLG(q1,w,n1,n2,e1,e21,e22,e3,d)),2) * S(q1, ni1, n1,r0)) / (sqrt(4 - q1*q1) * 2 * pi * hbar * EF) ;
}

double tau_2BLG_2spacer(double w, double n1, double n2,double e1, double e21 , double e22, double e3, double d, double r0, double ni1)
{
  int i, j;
  double result;
  double x[9000], wi[9000];
    i=100.; 
    gaulegf(0., 2., x, wi, i);
    result = 0.;
    double EF = pi*n1*hbar*hbar / (2*msao);
    for(j=1; j<=i; j++)
    {
        result +=  wi[j] *  ftau_2BLG_2spacer(x[j],w,n1,n2,e1,e21,e22,e3,d,r0,ni1); 
    }
    
  return result; 
}

// BLG_MLG - 2spacer

double ftau_BLG_MLG_2spacer(double q1, double w, double n1, double n2,double e1, double e21 , double e22, double e3, double d, double r0, double ni1)
{
   
    double EF = pi*n1*hbar*hbar / (2*msao);
    return ( pi * n1 * ni1 * q1*q1 * pow(1. - q1*q1/2.,2)* pow(fabs(W11_2layer_2spacer_BLG_MLG(q1,w,n1,n2,e1,e21,e22,e3,d)),2) * S(q1, ni1, n1,r0)) / (sqrt(4 - q1*q1) * 2 * pi * hbar * EF) ;
}

double tau_BLG_MLG_2spacer(double w, double n1, double n2,double e1, double e21 , double e22, double e3, double d, double r0, double ni1)
{
  int i, j;
  double result;
  double x[9000], wi[9000];
    i=100.; 
    gaulegf(0., 2., x, wi, i);
    result = 0.;
    double EF = pi*n1*hbar*hbar / (2*msao);
    for(j=1; j<=i; j++)
    {
        result +=  wi[j] *  ftau_BLG_MLG_2spacer(x[j],w,n1,n2,e1,e21,e22,e3,d,r0,ni1); 
    }
    
  return result; 
}


//******************************************************************************************************************************************* */
// 2BLG - 3spacer

double ftau11_2BLG_3spacer(double q1, double w1, double w2, double n1, double n2,double e1, double e21 , double e22, double e23, double e3, double d, double r01, double ni1)
{
   
    double EF = pi*n1*hbar*hbar / (2*msao);
    return ( pi * n1 * ni1 * q1*q1 * pow(1. - q1*q1/2.,2)* pow(fabs(W11_2layer_3spacer_2BLG(q1,w1,w2,n1,n2,e1,e21,e22,e23,e3,d)),2) * S(q1, ni1, n1,r01)) / (sqrt(4 - q1*q1) * 2 * pi * hbar * EF) ;
}

double tau11_2BLG_3spacer(double w1, double w2, double n1, double n2,double e1, double e21 , double e22, double e23, double e3, double d, double r01, double ni1)
{
  int i, j;
  double result;
  double x[9000], wi[9000];
    i=100.; 
    gaulegf(0., 2., x, wi, i);
    result = 0.;
    double EF = pi*n1*hbar*hbar / (2*msao);
    for(j=1; j<=i; j++)
    {
        result +=  wi[j] *  ftau11_2BLG_3spacer(x[j],w1,w2,n1,n2,e1,e21,e22,e23,e3,d,r01,ni1); 
    }
    
  return result; 
}

double ftau12_2BLG_3spacer(double q1, double w1, double w2, double n1, double n2,double e1, double e21 , double e22, double e23, double e3, double d, double r02, double ni2)
{
   
    double EF = pi*n1*hbar*hbar / (2*msao);
    return ( pi * n1 * ni2 * q1*q1 * pow(1. - q1*q1/2.,2)* pow(fabs(W11_2layer_3spacer_2BLG(q1,w1,w2,n1,n2,e1,e21,e22,e23,e3,d)),2) * S(q1, ni2, n1,r02)) / (sqrt(4 - q1*q1) * 2 * pi * hbar * EF) ;
}

double tau12_2BLG_3spacer(double w1, double w2, double n1, double n2,double e1, double e21 , double e22, double e23, double e3, double d, double r02, double ni2)
{
  int i, j;
  double result;
  double x[9000], wi[9000];
    i=100.; 
    gaulegf(0., 2., x, wi, i);
    result = 0.;
    double EF = pi*n1*hbar*hbar / (2*msao);
    for(j=1; j<=i; j++)
    {
        result +=  wi[j] *  ftau12_2BLG_3spacer(x[j],w1,w2,n1,n2,e1,e21,e22,e23,e3,d,r02,ni2); 
    }
    
  return result; 
}

double tauEf1_2BLG_3spacer(double w1, double w2, double n1, double n2, double e1, double e21, double e22, double e23, double e3, double d, double ni2, double r02, double ni1, double r01)
{
    return tau11_2BLG_3spacer(w1,w2,n1,n2,e1,e21,e22,e23,e3,d,r01,ni1) + tau12_2BLG_3spacer(w1,w2,n1,n2,e1,e21,e22,e23,e3,d,r02,ni2) ;

}

//*********************************************************************************************************************************/
// BLG_MLG - 3spacer

double ftau11_BLG_MLG_3spacer(double q1, double w1, double w2, double n1, double n2,double e1, double e21 , double e22, double e23, double e3, double d, double r01, double ni1)
{
   
    double EF = pi*n1*hbar*hbar / (2*msao);
    return ( pi * n1 * ni1 * q1*q1 * pow(1. - q1*q1/2.,2)* pow(fabs(W11_2layer_3spacer_BLG_MLG(q1,w1,w2,n1,n2,e1,e21,e22,e23,e3,d)),2) * S(q1, ni1, n1,r01)) / (sqrt(4 - q1*q1) * 2 * pi * hbar * EF) ;
}

double tau11_BLG_MLG_2spacer(double w1, double w2, double n1, double n2,double e1, double e21 , double e22, double e23, double e3, double d, double r01, double ni1)
{
  int i, j;
  double result;
  double x[9000], wi[9000];
    i=100.; 
    gaulegf(0., 2., x, wi, i);
    result = 0.;
    double EF = pi*n1*hbar*hbar / (2*msao);
    for(j=1; j<=i; j++)
    {
        result +=  wi[j] *  ftau11_BLG_MLG_3spacer(x[j],w1,w2,n1,n2,e1,e21,e22,e23,e3,d,r01,ni1); 
    }
    
  return result; 
}

double ftau12_BLG_MLG_3spacer(double q1, double w1, double w2, double n1, double n2,double e1, double e21 , double e22, double e23, double e3, double d, double r02, double ni2)
{
   
    double EF = pi*n1*hbar*hbar / (2*msao);
    return ( pi * n1 * ni2 * q1*q1 * pow(1. - q1*q1/2.,2)* pow(fabs(W12_2layer_3spacer_BLG_MLG(q1,w1,w2,n1,n2,e1,e21,e22,e23,e3,d)),2) * S(q1, ni2, n1,r02)) / (sqrt(4 - q1*q1) * 2 * pi * hbar * EF) ;
}

double tau12_BLG_MLG_3spacer(double w1, double w2, double n1, double n2,double e1, double e21 , double e22, double e23, double e3, double d, double r02, double ni2)
{
  int i, j;
  double result;
  double x[9000], wi[9000];
    i=100.; 
    gaulegf(0., 2., x, wi, i);
    result = 0.;
    double EF = pi*n1*hbar*hbar / (2*msao);
    for(j=1; j<=i; j++)
    {
        result +=  wi[j] *  ftau12_BLG_MLG_3spacer(x[j],w1,w2,n1,n2,e1,e21,e22,e23,e3,d,r02,ni2); 
    }
    
  return result; 
}

double tauEf1_BLG_MLG_3spacer(double w1, double w2, double n1, double n2, double e1, double e21, double e22, double e23, double e3, double d, double ni2, double r02, double ni1, double r01)
{
    return tau11_BLG_MLG_2spacer(w1,w2,n1,n2,e1,e21,e22,e23,e3,d,r01,ni1) + tau12_BLG_MLG_3spacer(w1,w2,n1,n2,e1,e21,e22,e23,e3,d,r02,ni2) ;

}

//**************************************************************************************************************************/
// 3BLG - T - 1spacer
double ftau11_3BLG_T_1spacer(double q, double ni, double d, double e1, double e2, double e3, double e4, double n, double r0, double t)
{
    double EF = pi * n * hbar * hbar / (2 * msao);
    return (pi * n * ni * q * q * pow(1. - q * q / 2., 2) * pow(fabs(W11_3BLG_T_1spacer(q, n, d, e1, e2, e3, e4, t)), 2) * S(q, ni, n, r0)) / (sqrt(4 - q * q) * 2 * pi * hbar * EF);
}

double tau11_3BLG_T_1spacer(double ni, double d, double e1, double e2, double e3, double e4, double n, double r0, double t)
{
    int i, j;
    double result;
    double x[9000], w[9000];
    i = 100.;
    gaulegf(0., 2., x, w, i);
    result = 0.;
    for (j = 1; j <= i; j++)
    {
        result += w[j] * ftau11_3BLG_T_1spacer(x[j], ni, d, e1, e2, e3, e4, n, r0, t);
    }

    return result;
}
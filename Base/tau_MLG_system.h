#include <stdio.h>
#include <math.h>

#include"W_2layer_2MLG.h"
#include"W_2layer_MLG_BLG.h"
#include"W_2layer_MLG_2DEG.h"
#include"W_2layer_3spacer_2MLG.h"
#include"S_func.h"
#include"gaulegf.h"


// 2MLG 

double ftau_2MLG(double q1,double ni1,double d, double e1, double e2, double e3, double n1, double n2, double r0)
{
   return (ni1  * pow(W11_2MLG(q1, d, e1, e2, e3, n1, n2),2) * S(q1,ni1,n1,r0) * q1 * q1 * sqrt(4 - q1 * q1) * sqrt(n1)) / (4 * sqrt(pi) * hbar * hbar * vF ) ;

}

double tau_2MLG(double ni1,double d, double e1, double e2, double e3, double n1, double n2, double r0)
{
  int i, j;
  double result;
  double x[9000], w[9000];
    i=100.; 
    gaulegf(0., 2., x, w, i);
    result = 0.;
    
    for(j=1; j<=i; j++)
    {
        result +=  w[j] *  ftau_2MLG(x[j],ni1,d,e1,e2,e3,n1,n2,r0); 
    }
    
  return result; 
}

// MLG-BLG

double ftau_MLG_BLG(double q1,double ni1,double d, double e1, double e2, double e3, double n1, double n2, double r0)
{
   return (ni1  * pow(W11_MLG_BLG(q1, d, e1, e2, e3, n1, n2),2) * S(q1,ni1,n1,r0) * q1 * q1 * sqrt(4 - q1 * q1) * sqrt(n1)) / (4 * sqrt(pi) * hbar * hbar * vF ) ;

}

double tau_MLG_BLG(double ni1,double d, double e1, double e2, double e3, double n1, double n2, double r0)
{
  int i, j;
  double result;
  double x[9000], w[9000];
    i=100.; 
    gaulegf(0., 2., x, w, i);
    result = 0.;
    
    for(j=1; j<=i; j++)
    {
        result +=  w[j] *  ftau_MLG_BLG(x[j],ni1,d,e1,e2,e3,n1,n2,r0); 
    }
    
  return result; 
}


// MLG-2DEG

double ftau_MLG_2DEG(double q1,double ni1,double d, double e1, double e2, double e3, double n1, double n2, double r0, double m2DEG)
{
   return (ni1  * pow(W11_MLG_2DEG(q1, d, e1, e2, e3, n1, n2, m2DEG),2) * S(q1,ni1,n1,r0) * q1 * q1 * sqrt(4 - q1 * q1) * sqrt(n1)) / (4 * sqrt(pi) * hbar * hbar * vF ) ;

}

double tau_MLG_2DEG(double ni1,double d, double e1, double e2, double e3, double n1, double n2, double r0, double m2DEG)
{
  int i, j;
  double result;
  double x[9000], w[9000];
    i=100.; 
    gaulegf(0., 2., x, w, i);
    result = 0.;
    
    for(j=1; j<=i; j++)
    {
        result +=  w[j] *  ftau_MLG_2DEG(x[j],ni1,d,e1,e2,e3,n1,n2,r0,m2DEG); 
    }
    
  return result; 
}


// 2MLG - 3 spacer
double ftau11_2layer_3spacer_2MLG (double q1, double w1, double w2, double n1, double n2, double e1, double e21 , double e22, double e23, double e3, double d, double ni1, double r01)
{
    return (ni1 * 2 * sqrt(n1) * q1 * q1 * sqrt(1 - q1*q1/4) * pow(W11_2layer_3spacer_2MLG(q1,w1,w2,n1,n2,e1,e21,e22,e23,e3,d),2) * S(q1,ni1,n1,r01)) / (4 * hbar * hbar * vF * sqrt(pi)) ;  
}
	
double tau11_2layer_3spacer_2MLG(double w1, double w2, double n1, double n2, double e1, double e21 , double e22, double e23, double e3, double d, double ni1, double r01)
{
  int i, j;
  double result;
  double x[9000], w[9000];
    i=100.;
    gaulegf(0., 2., x, w, i);
    result = 0.;
    
    for(j=1; j<=i; j++)
    {
        result +=  w[j] *  ftau11_2layer_3spacer_2MLG(x[j],w1,w2,n1,n2,e1,e21,e22,e23,e3,d,ni1,r01); 
    }
    
  return result;
}

double ftau12_2layer_3spacer_2MLG(double q1, double w1, double w2, double n1, double n2,double e1, double e21 , double e22, double e23, double e3, double d, double ni2, double r02)
{
    return (ni2 * 2 * sqrt(n1) * q1 * q1 * sqrt(1 - q1*q1/4) * pow(W12_2layer_3spacer_2MLG(q1,w1,w2,n1,n2,e1,e21,e22,e23,e3,d),2) * S(q1,ni2,n1,r02)) / (4 * hbar * hbar * vF * sqrt(pi)) ;  
}

double tau12_2layer_3spacer_2MLG(double w1, double w2, double n1, double n2,double e1, double e21 , double e22, double e23, double e3, double d, double ni2, double r02)
{
  int i, j;
  double result;
  double x[9000], w[9000];
    i=100.; 
    gaulegf(0., 2., x, w, i);
    result = 0.;
    
    for(j=1; j<=i; j++)
    {
        result +=  w[j] *  ftau12_2layer_3spacer_2MLG(x[j],w1,w2,n1,n2,e1,e21,e22,e23,e3,d,ni2,r02); 
    }
    
  return result; 
}


double tauEf1_2layer_3spacer_2MLG(double w1, double w2, double n1, double n2, double e1, double e21, double e22, double e23, double e3, double d, double ni2, double r02, double ni1, double r01)
{
    return tau11_2layer_3spacer_2MLG(w1,w2,n1,n2,e1,e21,e22,e23,e3,d,ni1,r01) + tau12_2layer_3spacer_2MLG(w1,w2,n1,n2,e1,e21,e22,e23,e3,d,ni2,r02) ;

}

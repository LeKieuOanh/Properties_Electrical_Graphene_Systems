#include <stdio.h>
#include <math.h>

#include"W_2layer_2MLG.h"
#include"W_2layer_MLG_BLG.h"
#include"W_2layer_MLG_2DEG.h"
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

#ifndef W_2layer_2MLG_h
#define W_2layer_2MLG_h

#include <stdio.h>
#include <math.h>
#include"V_2layer.h"
#include"Pola_MLG_0K.h"

#define elec      4.80320427*1e-10 // charge electron 
#define pi        3.141592653589793


double epsilon_2MLG(double q1, double d, double e1, double e2, double e3, double n1, double n2)
{
    return ((1 - PIMLG(q1,n1)* V11(q1,e1,e2,e3,d,n1) ) *  ( 1 - PIMLG(q1*sqrt(n1/1e12),n2)*V22(q1,e1,e2,e3,d,n1)) - (V12(q1,e1,e2,e3,d,n1)*V12(q1,e1,e2,e3,d,n1) * PIMLG(q1,n1) * PIMLG(q1*sqrt(n1/1e12),n2)));
}

double W11_2MLG(double q1, double d, double e1, double e2, double e3, double n1, double n2)
{
  return (  V11(q1,e1,e2,e3,d,n1) + ( pow(V12(q1,e1,e2,e3,d,n1),2) - V11(q1,e1,e2,e3,d,n1) * V22(q1,e1,e2,e3,d,n1))*PIMLG(q1*sqrt(n1/n2),n2)) / (epsilon_2MLG(q1,d,e1,e2,e3,n1,n2));
    
}

double W22_2MLG(double q, double d, double e1, double e2, double e3, double n1, double n2)
{
    return (V22(q, e1, e2, e3, d, n1) + ( pow(V12(q, e1, e2, e3, d, n1),2) - V11(q, e1, e2, e3, d, n1) * V22(q, e1, e2, e3, d, n1))*PIMLG(q,n1)) / (epsilon_2MLG(q, d, e1, e2, e3, n1, n2));
}

double W12_2MLG(double q, double d, double e1, double e2, double e3, double n1, double n2)
{
    return V12(q, e1, e2, e3, d, n1)/ (epsilon_2MLG(q, d, e1, e2, e3, n1, n2));
}

#endif
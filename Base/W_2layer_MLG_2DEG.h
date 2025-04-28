#ifndef W_2layer_MLG_2DEG_h
#define W_2layer_MLG_2DEG_h

#include <stdio.h>
#include <math.h>
#include"V_2layer.h"
#include"Pola_MLG_0K.h"
#include"Pola_2DEG_0K.h"

#define elec      4.80320427*1e-10 // charge electron 
#define pi        3.141592653589793


double epsilon_MLG_2DEG(double q1, double d, double e1, double e2, double e3, double n1, double n2, double m2DEG)
{
    return ((1 - PIMLG(q1,n1)* V11(q1,e1,e2,e3,d,n1) ) *  ( 1 - PI2DEG(q1*sqrt(n1/n2),m2DEG)*V22(q1,e1,e2,e3,d,n1)) - (V12(q1,e1,e2,e3,d,n1)*V12(q1,e1,e2,e3,d,n1) * PIMLG(q1,n1) * PI2DEG(q1*sqrt(n1/n2),m2DEG)));
}

double W11_MLG_2DEG(double q1, double d, double e1, double e2, double e3, double n1, double n2, double m2DEG)
{
  return (  V11(q1,e1,e2,e3,d,n1) + ( pow(V12(q1,e1,e2,e3,d,n1),2) - V11(q1,e1,e2,e3,d,n1) * V22(q1,e1,e2,e3,d,n1))*PI2DEG(q1*sqrt(n1/n2),m2DEG)) / (epsilon_MLG_2DEG(q1,d,e1,e2,e3,n1,n2,m2DEG));
    
}

double W22_MLG_2DEG(double q, double d, double e1, double e2, double e3, double n1, double n2, double m2DEG)
{
    return (V22(q, e1, e2, e3, d, n1) + ( pow(V12(q, e1, e2, e3, d, n1),2) - V11(q, e1, e2, e3, d, n1) * V22(q, e1, e2, e3, d, n1))*PIMLG(q,n1)) / (epsilon_MLG_2DEG(q, d, e1, e2, e3, n1, n2, m2DEG));
}

double W12_MLG_2DEG(double q, double d, double e1, double e2, double e3, double n1, double n2, double m2DEG)
{
    return V12(q, e1, e2, e3, d, n1)/ (epsilon_MLG_2DEG(q, d, e1, e2, e3, n1, n2, m2DEG));
}

#endif
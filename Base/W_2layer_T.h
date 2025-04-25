#include <stdio.h>
#include <math.h>
#include"V_2layer.h"
#include"Pola_BLG_T.h"

#define elec      4.80320427*1e-10 // charge electron 
#define pi        3.141592653589793


double epsilon(double q, double d, double e1, double e2, double e3, double n, double n2, double t)
{
    return ((1 - polarinfi(q, t)* V11(q, e1, e2, e3, d, n) ) *  ( 1 - polarinfi(q*sqrt(n/n2), t)*V22(q, e1, e2, e3, d, n)) - pow(V12(q, e1, e2, e3, d, n),2) * polarinfi(q, t) * polarinfi(q*sqrt(n/n2), t));
}

double W11(double q, double d, double e1, double e2, double e3, double n, double n2, double t)
{
    return (V11(q, e1, e2, e3, d, n)) + ( pow(V12(q, e1, e2, e3, d, n),2) - V11(q, e1, e2, e3, d, n) * V22(q, e1, e2, e3, d, n))*polarinfi(q*sqrt(n/n2), t) / (epsilon(q, d, e1, e2, e3, n, n2, t));
}

double W22(double q, double d, double e1, double e2, double e3, double n, double n2, double t)
{
    return (V22(q, e1, e2, e3, d, n)) + ( pow(V12(q, e1, e2, e3, d, n),2) - V11(q, e1, e2, e3, d, n) * V22(q, e1, e2, e3, d, n))*polarinfi(q, t) / (epsilon(q, d, e1, e2, e3, n, n2, t));
}

double W12(double q, double d, double e1, double e2, double e3, double n, double n2, double t)
{
    return V12(q, e1, e2, e3, d, n)/ (epsilon(q, d, e1, e2, e3, n, n2, t));
}
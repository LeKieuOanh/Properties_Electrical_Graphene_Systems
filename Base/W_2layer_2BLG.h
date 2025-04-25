#ifndef W_2layer_2BLG_h
#define W_2layer_2BLG_h

#include <stdio.h>
#include <math.h>
#include"V_2layer.h"
#include"Pola_BLG_0K.h"

#define elec      4.80320427*1e-10 // charge electron 
#define pi        3.141592653589793


double epsilon_2BLG(double q, double d, double e1, double e2, double e3, double n, double n2)
{
    return ((1 - PIBLG(q)* V11(q, e1, e2, e3, d, n) ) *  ( 1 - PIBLG(q*sqrt(n/n2))*V22(q, e1, e2, e3, d, n)) - pow(V12(q, e1, e2, e3, d, n),2) * PIBLG(q) * PIBLG(q*sqrt(n/n2)));
}

double W11_2BLG(double q, double d, double e1, double e2, double e3, double n, double n2)
{
    return (V11(q, e1, e2, e3, d, n)) + ( pow(V12(q, e1, e2, e3, d, n),2) - V11(q, e1, e2, e3, d, n) * V22(q, e1, e2, e3, d, n))*PIBLG(q*sqrt(n/n2)) / (epsilon_2BLG(q, d, e1, e2, e3, n, n2));
}

double W22_2BLG(double q, double d, double e1, double e2, double e3, double n, double n2)
{
    return (V22(q, e1, e2, e3, d, n)) + ( pow(V12(q, e1, e2, e3, d, n),2) - V11(q, e1, e2, e3, d, n) * V22(q, e1, e2, e3, d, n))*PIBLG(q) / (epsilon_2BLG(q, d, e1, e2, e3, n, n2));
}

double W12_2BLG(double q, double d, double e1, double e2, double e3, double n, double n2)
{
    return V12(q, e1, e2, e3, d, n)/ (epsilon_2BLG(q, d, e1, e2, e3, n, n2));
}

#endif
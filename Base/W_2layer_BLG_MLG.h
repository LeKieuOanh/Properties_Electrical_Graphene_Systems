#ifndef W_2layer_BLG_MLG_h
#define W_2layer_BLG_MLG_h

#include <stdio.h>
#include <math.h>
#include"V_2layer.h"
#include"Pola_BLG_0K.h"
#include"Pola_MLG_0K.h"

#define elec      4.80320427*1e-10 // charge electron 
#define pi        3.141592653589793

// 1: BLG ----- 2: MLG

double epsilon_BLG_MLG(double q, double d, double e1, double e2, double e3, double n, double n2)
{
    return ((1 - PIBLG(q)* V11(q, e1, e2, e3, d, n) ) *  ( 1 - PIMLG(q*sqrt(n/n2), n2)*V22(q, e1, e2, e3, d, n)) - pow(V12(q, e1, e2, e3, d, n),2) * PIBLG(q) * PIMLG(q*sqrt(n/n2), n2));
}

double W11_BLG_MLG(double q, double d, double e1, double e2, double e3, double n, double n2)
{
    return (V11(q, e1, e2, e3, d, n) + ( pow(V12(q, e1, e2, e3, d, n),2) - V11(q, e1, e2, e3, d, n) * V22(q, e1, e2, e3, d, n))*PIMLG(q*sqrt(n/n2), n2)) / (epsilon_BLG_MLG(q, d, e1, e2, e3, n, n2));
}

double W22_BLG_MLG(double q, double d, double e1, double e2, double e3, double n, double n2)
{
    return (V22(q, e1, e2, e3, d, n) + ( pow(V12(q, e1, e2, e3, d, n),2) - V11(q, e1, e2, e3, d, n) * V22(q, e1, e2, e3, d, n))*PIBLG(q)) / (epsilon_BLG_MLG(q, d, e1, e2, e3, n, n2));
}

double W12_BLG_MLG(double q, double d, double e1, double e2, double e3, double n, double n2)
{
    return V12(q, e1, e2, e3, d, n)/ (epsilon_BLG_MLG(q, d, e1, e2, e3, n, n2));
}


#endif
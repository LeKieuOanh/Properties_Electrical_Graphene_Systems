#ifndef W_2layer_BLG_2DEG_h
#define W_2layer_BLG_2DEG_h

#include <stdio.h>
#include <math.h>
#include"V_2layer.h"
#include"Pola_BLG_0K.h"
#include"Pola_2DEG_0K.h"

#define elec      4.80320427*1e-10 // charge electron 
#define pi        3.141592653589793

// 1: BLG ----- 2: GaAs

double epsilon_BLG_2DEG(double q, double d, double e1, double e2, double e3, double n, double n2, double m2DEG)
{
    return ((1 - PIBLG(q)* V11(q, e1, e2, e3, d, n) ) *  ( 1 - PI2DEG(q*sqrt(n/n2), m2DEG)*V22(q, e1, e2, e3, d, n)) - pow(V12(q, e1, e2, e3, d, n),2) * PIBLG(q) * PI2DEG(q*sqrt(n/n2), m2DEG));
}

double W11_BLG_2DEG(double q, double d, double e1, double e2, double e3, double n, double n2, double m2DEG)
{
    return (V11(q, e1, e2, e3, d, n)) + ( pow(V12(q, e1, e2, e3, d, n),2) - V11(q, e1, e2, e3, d, n) * V22(q, e1, e2, e3, d, n))*PI2DEG(q*sqrt(n/n2), m2DEG) / (epsilon_BLG_2DEG(q, d, e1, e2, e3, n, n2, m2DEG));
}

double W22_BLG_2DEG(double q, double d, double e1, double e2, double e3, double n, double n2, double m2DEG)
{
    return (V22(q, e1, e2, e3, d, n)) + ( pow(V12(q, e1, e2, e3, d, n),2) - V11(q, e1, e2, e3, d, n) * V22(q, e1, e2, e3, d, n))*PIBLG(q) / (epsilon_BLG_2DEG(q, d, e1, e2, e3, n, n2, m2DEG));
}

double W12_BLG_2DEG(double q, double d, double e1, double e2, double e3, double n, double n2, double m2DEG)
{
    return V12(q, e1, e2, e3, d, n)/ (epsilon_BLG_2DEG(q, d, e1, e2, e3, n, n2, m2DEG));
}

#endif
#ifndef W_layers_h
#define W_layers_h

#include <stdio.h>
#include <math.h>
#include"V_layers.h"
#include"Pola_0K.h"

// 2layer - 1spacer - 2MLG

double epsilon_2MLG_1spacer(double q1, double d, double e1, double e2, double e3, double n1, double n2)
{
    return ((1 - PIMLG(q1,n1)* V11_2layer_1spacer(q1,e1,e2,e3,d,n1) ) *  ( 1 - PIMLG(q1*sqrt(n1/1e12),n2)*V22_2layer_1spacer(q1,e1,e2,e3,d,n1)) - (V12_2layer_1spacer(q1,e1,e2,e3,d,n1)*V12_2layer_1spacer(q1,e1,e2,e3,d,n1) * PIMLG(q1,n1) * PIMLG(q1*sqrt(n1/1e12),n2)));
}

double W11_2MLG_1spacer(double q1, double d, double e1, double e2, double e3, double n1, double n2)
{
  return (  V11_2layer_1spacer(q1,e1,e2,e3,d,n1) + ( pow(V12_2layer_1spacer(q1,e1,e2,e3,d,n1),2) - V11_2layer_1spacer(q1,e1,e2,e3,d,n1) * V22_2layer_1spacer(q1,e1,e2,e3,d,n1))*PIMLG(q1*sqrt(n1/n2),n2)) / (epsilon_2MLG_1spacer(q1,d,e1,e2,e3,n1,n2));
    
}

double W22_2MLG_1spacer(double q, double d, double e1, double e2, double e3, double n1, double n2)
{
    return (V22_2layer_1spacer(q, e1, e2, e3, d, n1) + ( pow(V12_2layer_1spacer(q, e1, e2, e3, d, n1),2) - V11_2layer_1spacer(q, e1, e2, e3, d, n1) * V22_2layer_1spacer(q, e1, e2, e3, d, n1))*PIMLG(q,n1)) / (epsilon_2MLG_1spacer(q, d, e1, e2, e3, n1, n2));
}

double W12_2MLG_1spacer(double q, double d, double e1, double e2, double e3, double n1, double n2)
{
    return V12_2layer_1spacer(q, e1, e2, e3, d, n1)/ (epsilon_2MLG_1spacer(q, d, e1, e2, e3, n1, n2));
}

// 2layer - 1 spacer - MLG_BLG

double epsilon_MLG_BLG_1spacer(double q1, double d, double e1, double e2, double e3, double n1, double n2)
{
    return ((1 - PIMLG(q1,n1)* V11_2layer_1spacer(q1,e1,e2,e3,d,n1) ) *  ( 1 - PIBLG(q1*sqrt(n1/n2))*V22_2layer_1spacer(q1,e1,e2,e3,d,n1)) - (V12_2layer_1spacer(q1,e1,e2,e3,d,n1)*V12_2layer_1spacer(q1,e1,e2,e3,d,n1) * PIMLG(q1,n1) * PIBLG(q1*sqrt(n1/n2))));
}

double W11_MLG_BLG_1spacer(double q1, double d, double e1, double e2, double e3, double n1, double n2)
{
  return (  V11_2layer_1spacer(q1,e1,e2,e3,d,n1) + ( pow(V12_2layer_1spacer(q1,e1,e2,e3,d,n1),2) - V11_2layer_1spacer(q1,e1,e2,e3,d,n1) * V22_2layer_1spacer(q1,e1,e2,e3,d,n1))*PIBLG(q1*sqrt(n1/n2))) / (epsilon_MLG_BLG_1spacer(q1,d,e1,e2,e3,n1,n2));
    
}

double W22_MLG_BLG_1spacer(double q, double d, double e1, double e2, double e3, double n1, double n2)
{
    return (V22_2layer_1spacer(q, e1, e2, e3, d, n1) + ( pow(V12_2layer_1spacer(q, e1, e2, e3, d, n1),2) - V11_2layer_1spacer(q, e1, e2, e3, d, n1) * V22_2layer_1spacer(q, e1, e2, e3, d, n1))*PIMLG(q,n1)) / (epsilon_MLG_BLG_1spacer(q, d, e1, e2, e3, n1, n2));
}

double W12_MLG_BLG_1spacer(double q, double d, double e1, double e2, double e3, double n1, double n2)
{
    return V12_2layer_1spacer(q, e1, e2, e3, d, n1)/ (epsilon_MLG_BLG_1spacer(q, d, e1, e2, e3, n1, n2));
}

// 2layer - 1spacer - MLG_2DEG

double epsilon_MLG_2DEG_1spacer(double q1, double d, double e1, double e2, double e3, double n1, double n2, double m2DEG)
{
    return ((1 - PIMLG(q1,n1)* V11_2layer_1spacer(q1,e1,e2,e3,d,n1) ) *  ( 1 - PI2DEG(q1*sqrt(n1/n2),m2DEG)*V22_2layer_1spacer(q1,e1,e2,e3,d,n1)) - (V12_2layer_1spacer(q1,e1,e2,e3,d,n1)*V12_2layer_1spacer(q1,e1,e2,e3,d,n1) * PIMLG(q1,n1) * PI2DEG(q1*sqrt(n1/n2),m2DEG)));
}

double W11_MLG_2DEG_1spacer(double q1, double d, double e1, double e2, double e3, double n1, double n2, double m2DEG)
{
  return (  V11_2layer_1spacer(q1,e1,e2,e3,d,n1) + ( pow(V12_2layer_1spacer(q1,e1,e2,e3,d,n1),2) - V11_2layer_1spacer(q1,e1,e2,e3,d,n1) * V22_2layer_1spacer(q1,e1,e2,e3,d,n1))*PI2DEG(q1*sqrt(n1/n2),m2DEG)) / (epsilon_MLG_2DEG_1spacer(q1,d,e1,e2,e3,n1,n2,m2DEG));
    
}

double W22_MLG_2DEG_1spacer(double q, double d, double e1, double e2, double e3, double n1, double n2, double m2DEG)
{
    return (V22_2layer_1spacer(q, e1, e2, e3, d, n1) + ( pow(V12_2layer_1spacer(q, e1, e2, e3, d, n1),2) - V11_2layer_1spacer(q, e1, e2, e3, d, n1) * V22_2layer_1spacer(q, e1, e2, e3, d, n1))*PIMLG(q,n1)) / (epsilon_MLG_2DEG_1spacer(q, d, e1, e2, e3, n1, n2, m2DEG));
}

double W12_MLG_2DEG_1spacer(double q, double d, double e1, double e2, double e3, double n1, double n2, double m2DEG)
{
    return V12_2layer_1spacer(q, e1, e2, e3, d, n1)/ (epsilon_MLG_2DEG_1spacer(q, d, e1, e2, e3, n1, n2, m2DEG));
}

// 2layer - 1spacer - 2BLG

double epsilon_2BLG_1spacer(double q, double d, double e1, double e2, double e3, double n, double n2)
{
    return ((1 - PIBLG(q)* V11_2layer_1spacer(q, e1, e2, e3, d, n) ) *  ( 1 - PIBLG(q*sqrt(n/n2))*V22_2layer_1spacer(q, e1, e2, e3, d, n)) - pow(V12_2layer_1spacer(q, e1, e2, e3, d, n),2) * PIBLG(q) * PIBLG(q*sqrt(n/n2)));
}

double W11_2BLG_1spacer(double q, double d, double e1, double e2, double e3, double n, double n2)
{
    return ( V11_2layer_1spacer(q, e1, e2, e3, d, n) + ( pow(V12_2layer_1spacer(q, e1, e2, e3, d, n),2) - V11_2layer_1spacer(q, e1, e2, e3, d, n) * V22_2layer_1spacer(q, e1, e2, e3, d, n))*PIBLG(q*sqrt(n/n2)) ) / (epsilon_2BLG_1spacer(q, d, e1, e2, e3, n, n2));
}

double W22_2BLG_1spacer(double q, double d, double e1, double e2, double e3, double n, double n2)
{
    return (V22_2layer_1spacer(q, e1, e2, e3, d, n) + ( pow(V12_2layer_1spacer(q, e1, e2, e3, d, n),2) - V11_2layer_1spacer(q, e1, e2, e3, d, n) * V22_2layer_1spacer(q, e1, e2, e3, d, n))*PIBLG(q)) / (epsilon_2BLG_1spacer(q, d, e1, e2, e3, n, n2));
}

double W12_2BLG_1spacer(double q, double d, double e1, double e2, double e3, double n, double n2)
{
    return V12_2layer_1spacer(q, e1, e2, e3, d, n)/ (epsilon_2BLG_1spacer(q, d, e1, e2, e3, n, n2));
}


// 2layer - 1spacer - BLG_MLG
double epsilon_BLG_MLG_1spacer(double q, double d, double e1, double e2, double e3, double n, double n2)
{
    return ((1 - PIBLG(q)* V11_2layer_1spacer(q, e1, e2, e3, d, n) ) *  ( 1 - PIMLG(q*sqrt(n/n2), n2)*V22_2layer_1spacer(q, e1, e2, e3, d, n)) - pow(V12_2layer_1spacer(q, e1, e2, e3, d, n),2) * PIBLG(q) * PIMLG(q*sqrt(n/n2), n2));
}

double W11_BLG_MLG_1spacer(double q, double d, double e1, double e2, double e3, double n, double n2)
{
    return (V11_2layer_1spacer(q, e1, e2, e3, d, n) + ( pow(V12_2layer_1spacer(q, e1, e2, e3, d, n),2) - V11_2layer_1spacer(q, e1, e2, e3, d, n) * V22_2layer_1spacer(q, e1, e2, e3, d, n))*PIMLG(q*sqrt(n/n2), n2)) / (epsilon_BLG_MLG_1spacer(q, d, e1, e2, e3, n, n2));
}

double W22_BLG_MLG_1spacer(double q, double d, double e1, double e2, double e3, double n, double n2)
{
    return (V22_2layer_1spacer(q, e1, e2, e3, d, n) + ( pow(V12_2layer_1spacer(q, e1, e2, e3, d, n),2) - V11_2layer_1spacer(q, e1, e2, e3, d, n) * V22_2layer_1spacer(q, e1, e2, e3, d, n))*PIBLG(q)) / (epsilon_BLG_MLG_1spacer(q, d, e1, e2, e3, n, n2));
}

double W12_BLG_MLG_1spacer(double q, double d, double e1, double e2, double e3, double n, double n2)
{
    return V12_2layer_1spacer(q, e1, e2, e3, d, n)/ (epsilon_BLG_MLG_1spacer(q, d, e1, e2, e3, n, n2));
}


// 2layer - 1spacer - BLG_2DEG

double epsilon_BLG_2DEG_1spacer(double q, double d, double e1, double e2, double e3, double n, double n2, double m2DEG)
{
    return ((1 - PIBLG(q)* V11_2layer_1spacer(q, e1, e2, e3, d, n) ) *  ( 1 - PI2DEG(q*sqrt(n/n2), m2DEG)*V22_2layer_1spacer(q, e1, e2, e3, d, n)) - pow(V12_2layer_1spacer(q, e1, e2, e3, d, n),2) * PIBLG(q) * PI2DEG(q*sqrt(n/n2), m2DEG));
}

double W11_BLG_2DEG_1spacer(double q, double d, double e1, double e2, double e3, double n, double n2, double m2DEG)
{
    return (V11_2layer_1spacer(q, e1, e2, e3, d, n) + ( pow(V12_2layer_1spacer(q, e1, e2, e3, d, n),2) - V11_2layer_1spacer(q, e1, e2, e3, d, n) * V22_2layer_1spacer(q, e1, e2, e3, d, n))*PI2DEG(q*sqrt(n/n2), m2DEG)) / (epsilon_BLG_2DEG_1spacer(q, d, e1, e2, e3, n, n2, m2DEG));
}

double W22_BLG_2DEG_1spacer(double q, double d, double e1, double e2, double e3, double n, double n2, double m2DEG)
{
    return (V22_2layer_1spacer(q, e1, e2, e3, d, n) + ( pow(V12_2layer_1spacer(q, e1, e2, e3, d, n),2) - V11_2layer_1spacer(q, e1, e2, e3, d, n) * V22_2layer_1spacer(q, e1, e2, e3, d, n))*PIBLG(q)) / (epsilon_BLG_2DEG_1spacer(q, d, e1, e2, e3, n, n2, m2DEG));
}

double W12_BLG_2DEG_1spacer(double q, double d, double e1, double e2, double e3, double n, double n2, double m2DEG)
{
    return V12_2layer_1spacer(q, e1, e2, e3, d, n)/ (epsilon_BLG_2DEG_1spacer(q, d, e1, e2, e3, n, n2, m2DEG));
}


// 2layer - 3spacer - 2MLG

double epsilon_2layer_3spacer_2MLG(double q1, double w1, double w2, double n1, double n2, double e1, double e21, double e22, double e23, double e3, double d)
{
    return (1 - PIMLG(q1, n1)*V11_2layer_3spacer(q1*sqrt(pi*n1),w1, w2 ,e1,e21,e22,e23,e3,d)) * (1 - PIMLG(q1*sqrt(n1/n2),n2)*V22_2layer_3spacer(q1*sqrt(pi*n1),w1, w2 , e1, e21, e22, e23, e3, d)) - V12_2layer_3spacer(q1*sqrt(pi*n1),w1, w2,e1,e21,e22,e23,e3,d)*V12_2layer_3spacer(q1*sqrt(pi*n1),w1,w2,e1,e21,e22,e23,e3,d)*PIMLG(q1,n1)*PIMLG(q1*sqrt(n1/n2),n2);
}

double W11_2layer_3spacer_2MLG(double q1, double w1, double w2, double n1, double n2, double e1, double e21, double e22, double e23, double e3, double d)
{
    return (V11_2layer_3spacer(q1*sqrt(pi*n1),w1,w2,e1,e21,e22,e23,e3,d) + (V12_2layer_3spacer(q1*sqrt(pi*n1),w1, w2,e1,e21,e22,e23, e3,d)*V12_2layer_3spacer(q1*sqrt(pi*n1),w1, w2,e1,e21,e22,e23,e3,d) - V11_2layer_3spacer(q1*sqrt(pi*n1),w1,w2,e1,e21,e22,e23,e3,d)*V22_2layer_3spacer(q1*sqrt(pi*n1),w1,w2,e1,e21,e22,e23,e3,d))*PIMLG(q1*sqrt(n1/n2),n2))/(epsilon_2layer_3spacer_2MLG(q1,w1,w2,n1,n2,e1,e21,e22,e23,e3,d));
}

double W12_2layer_3spacer_2MLG(double q1, double w1, double w2, double n1, double n2, double e1, double e21, double e22, double e23, double e3, double d)
{
    return V12_2layer_3spacer(q1*sqrt(pi*n1),w1,w2,e1,e21,e22,e23,e3,d) / epsilon_2layer_3spacer_2MLG(q1,w1,w2,n1,n2,e1,e21,e22,e23,e3,d);
}

#endif
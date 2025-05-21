#ifndef V_layers_h
#define V_layers_h

#include <stdio.h>
#include <math.h>

#define elec      4.80320427*1e-10 // charge electron 
#define pi        3.141592653589793


// 2layer - 1spacer
double heso_2layer_1spacer(double n, double q){
    return 4*pi*pow(elec,2)/ (q * sqrt(pi * n)) ;
}

double MS_2layer_1spacer(double q,double e1, double e2, double e3, double d, double n){
    return (e2*(e1+e3) + ( pow(e2,2) + e1*e3)*tanh(q * sqrt(pi * n) * d));
}

double V11_2layer_1spacer(double q,double e1, double e2, double e3, double d, double n){
    return heso_2layer_1spacer(n, q) * (e2 + e3 * tanh(q * sqrt(pi * n)*d)) / MS_2layer_1spacer(q, e1, e2, e3, d, n);
}

double V22_2layer_1spacer(double q,double e1, double e2, double e3, double d, double n){
    return heso_2layer_1spacer(n, q) * (e2 + e1*tanh(q*sqrt(pi * n)*d)) / MS_2layer_1spacer(q, e1, e2, e3, d, n);
}

double V12_2layer_1spacer(double q,double e1, double e2, double e3, double d, double n){
    return heso_2layer_1spacer(n, q) * e2 / (cosh(q * sqrt(pi * n)*d) *MS_2layer_1spacer(q, e1, e2, e3, d, n));
}


// 2layer - 2spacer
double V22_2layer_2spacer(double q, double w, double e1, double e21, double e22, double e3, double d) // V_bb
{
    double coeff = 2*pi*pow(elec,2)/ q;
    double MS = exp(2*q*(d+w))*(e1 + e21)*(e21 + e22)*(e3 + e22) + exp(2*q*d)*(e1 + e21)*(e21 - e22)*(e22 - e3) + exp(4*q*w)*(e1 - e21)*(e21 - e22)*(e3 + e22) + exp(2*q*w)*(e1 - e21)*(e21 + e22)*(e22 - e3);
    double e_bb = 2*(exp(2*q*(d+w))*(e1 + e21)*(e21 + e22) + exp(2*q*d)*(e1 + e21)*(e22 - e21) + exp(4*q*w)*(e21 - e1)*(e22 - e21) + exp(2*q*w)*(e21 - e1)*(e22 + e22));
    return coeff * e_bb / MS;
}

double V11_2layer_2spacer(double q, double w, double e1, double e21, double e22, double e3, double d)  // V_tt
{
    double coeff = 2*pi*pow(elec,2)/ q;
    double MS = exp(2*q*(d+w))*(e1 + e21)*(e21 + e22)*(e3 + e22) + exp(2*q*d)*(e1 + e21)*(e21 - e22)*(e22 - e3) + exp(4*q*w)*(e1 - e21)*(e21 - e22)*(e3 + e22) + exp(2*q*w)*(e1 - e21)*(e21 + e22)*(e22 - e3);
    double e_tt = 2*(exp(2*q*(d+w))*(e3 + e22)*(e21 + e22) + exp(2*q*d)*(e22 - e3)*(e21 - e22) + exp(4*q*w)*(e21 - e22)*(e22 + e3) + exp(2*q*w)*(e22 - e3)*(e21 + e22));
    return coeff * e_tt / MS;
}

double V12_2layer_2spacer(double q, double w, double e1, double e21, double e22, double e3, double d)
{
    double coeff = 2*pi*pow(elec,2)/ (q );
    double MS = exp(2*q*(d+w))*(e1 + e21)*(e21 + e22)*(e3 + e22) + exp(2*q*d)*(e1 + e21)*(e21 - e22)*(e22 - e3) + exp(4*q*w)*(e1 - e21)*(e21 - e22)*(e3 + e22) + exp(2*q*w)*(e1 - e21)*(e21 + e22)*(e22 - e3);
    double e_tb = 8*(exp(q*(d+2*w))*e21*e22);
    return coeff * e_tb / MS;
}


// 2layer - 3spacer
double V22_2layer_3spacer(double q, double w1, double w2, double e1, double e21, double e22, double e23, double e3, double d) // V_bb
{
    double coeff = 4*pi*pow(elec,2)/ q;
    double MS = (exp(2*q*(2*w1+w2))*(e3 + e23)*(e21+e22)*(e1 - e21)*(e22-e23) + exp(2*q*(w1+2*w2))*(e3 + e23)*(e22 + e23)*(e1 - e21)*(e21 - e22) + exp(2*q*(d+2*w1))*(e1+e21)*(e3+e23)*(e21-e22)*(e22-e23) + exp(2*q*(d+w1+w2))*(e1+e21)*(e3+e23)*(e21+e22)*(e22+e23) - exp(2*q*(d+w1))*(e1+e21)*(e22+e23)*(e3-e23)*(e21-e22) - exp(2*q*(d+w2))*(e1+e21)*(e21+e22)*(e3-e23)*(e22-e23) - exp(2*q*(w1+w2))*(e21+e22)*(e22+e23)*(e1-e21)*(e3-e23) - exp(4*q*w2)*(e1-e21)*(e3-e23)*(e21-e22)*(e22-e23));
    double e_bb = exp(2*q*(d+w1+w2))*(e1+e21)*(e21+e22)*(e22+e23) - exp(2*q*(w1+w2))*(e21+e22)*(e22+e23)*(e1-e21) + exp(2*q*(d+2*w1))*(e1+e21)*(e21-e22)*(e22-e23) - exp(4*q*w2)*(e1-e21)*(e21-e22)*(e22-e23) - exp(2*q*(d+w2))*(e1+e21)*(e21+e22)*(e22-e23) - exp(2*q*(d+w1))*(e1+e21)*(e22+e23)*(e21-e22) + exp(2*q*(2*w1+w2))*(e21+e22)*(e1-e21)*(e22-e23) + exp(2*q*(w1+2*w2))*(e22+e23)*(e1-e21)*(e21-e22);
    return coeff * e_bb / MS;
}

double V11_2layer_3spacer(double q, double w1, double w2, double e1, double e21, double e22, double e23, double e3, double d)  // V_tt
{
    double coeff = 4*pi*pow(elec,2)/ q;
    double MS = (exp(2*q*(2*w1+w2))*(e3 + e23)*(e21+e22)*(e1 - e21)*(e22-e23) + exp(2*q*(w1+2*w2))*(e3 + e23)*(e22 + e23)*(e1 - e21)*(e21 - e22) + exp(2*q*(d+2*w1))*(e1+e21)*(e3+e23)*(e21-e22)*(e22-e23) + exp(2*q*(d+w1+w2))*(e1+e21)*(e3+e23)*(e21+e22)*(e22+e23) - exp(2*q*(d+w1))*(e1+e21)*(e22+e23)*(e3-e23)*(e21-e22) - exp(2*q*(d+w2))*(e1+e21)*(e21+e22)*(e3-e23)*(e22-e23) - exp(2*q*(w1+w2))*(e21+e22)*(e22+e23)*(e1-e21)*(e3-e23) - exp(4*q*w2)*(e1-e21)*(e3-e23)*(e21-e22)*(e22-e23));
    double e_tt = exp(2*q*(d+w1+w2))*(e3+e23)*(e21+e22)*(e22+e23) - exp(2*q*(w1+w2))*(e21+e22)*(e22+e23)*(e3-e23) + exp(2*q*(d+2*w1))*(e3+e23)*(e21-e22)*(e22-e23) - exp(4*q*w2)*(e3-e23)*(e21-e22)*(e22-e23) - exp(2*q*(d+w2))*(e21+e22)*(e3-e23)*(e22-e23) - exp(2*q*(d+w1))*(e22+e23)*(e3-e23)*(e21-e22) + exp(2*q*(2*w1+w2))*(e3+e23)*(e21+e22)*(e22-e23) + exp(2*q*(w1+2*w2))*(e3+e23)*(e22+e23)*(e21-e22);
    return coeff * e_tt / MS;
}

double V12_2layer_3spacer(double q, double w1, double w2, double e1, double e21, double e22, double e23, double e3, double d)
{
    double coeff = 4*pi*pow(elec,2)/q ;
    double MS = (exp(2*q*(2*w1+w2))*(e3 + e23)*(e21+e22)*(e1 - e21)*(e22-e23) + exp(2*q*(w1+2*w2))*(e3 + e23)*(e22 + e23)*(e1 - e21)*(e21 - e22) + exp(2*q*(d+2*w1))*(e1+e21)*(e3+e23)*(e21-e22)*(e22-e23) + exp(2*q*(d+w1+w2))*(e1+e21)*(e3+e23)*(e21+e22)*(e22+e23) - exp(2*q*(d+w1))*(e1+e21)*(e22+e23)*(e3-e23)*(e21-e22) - exp(2*q*(d+w2))*(e1+e21)*(e21+e22)*(e3-e23)*(e22-e23) - exp(2*q*(w1+w2))*(e21+e22)*(e22+e23)*(e1-e21)*(e3-e23) - exp(4*q*w2)*(e1-e21)*(e3-e23)*(e21-e22)*(e22-e23));
    double e_tb = 8*e21*e22*e23*exp(q*(d+2*w1+2*w2));
    return coeff * e_tb / MS;
}


// 3layer - 1spacer

double f11(double q, double d, double k2, double k3, double k4)
{
    return 2 * ((k2 + k3) * (k3 - k4) + 2 * k3 * (k2 - k3) * exp(2 * q * d) + (k2 + k3) * (k3 + k4) * exp(4 * q * d));
}

double f22(double q, double d, double k1, double k2, double k3, double k4)
{
    return 8 * exp(2 * q * d) * (k1 * cosh(q * d) + k2 * sinh(q * d)) * (k3 * cosh(q * d) + k4 * sinh(q * d));
}

double f33(double q, double d, double k1, double k2, double k3)
{
    return 2 * ((k2 + k3) * (k2 - k1) + 2 * k2 * (k3 - k2) * exp(2 * q * d) + (k1 + k2) * (k2 + k3) * exp(4 * q * d));
}

double f12(double q, double d, double k2, double k3, double k4)
{
    return 8 * k2 * exp(2 * q * d) * (k3 * cosh(q * d) + k4 * sinh(q * d));
}

double f13(double q, double d, double k2, double k3)
{
    return 8 * k2 * (k3 * exp(2 * q * d));
}

double f23(double q, double d, double k1, double k2, double k3)
{
    return 8 * k3 * exp(q * d) * (k2 * cosh(q * d) + k1 * sinh(q * d));
}

double heso_3layer(double q)
{
    return 2 * pi * pow(elec, 2) / q;
}

double MS_V_3layer(double q, double d, double k1, double k2, double k3, double k4)
{
    return (k1 - k2) * (k2 + k3) * (k3 - k4) + 2 * exp(2 * q * d) * (k2 - k3) * (k1 * k3 - k2 * k4) + exp(4 * q * d) * (k1 + k2) * (k2 + k3) * (k3 + k4);
}

double V11_3layer(double q, double d, double k1, double k2, double k3, double k4)
{
    return heso_3layer(q) * f11(q, d, k2, k3, k4) / MS_V_3layer(q, d, k1, k2, k3, k4);
}

double V22_3layer(double q, double d, double k1, double k2, double k3, double k4)
{
    return heso_3layer(q) * f22(q, d, k1, k2, k3, k4) / MS_V_3layer(q, d, k1, k2, k3, k4);
}

double V33_3layer(double q, double d, double k1, double k2, double k3, double k4)
{
    return heso_3layer(q) * f33(q, d, k1, k2, k3) / MS_V_3layer(q, d, k1, k2, k3, k4);
}

double V12_3layer(double q, double d, double k1, double k2, double k3, double k4)
{
    return heso_3layer(q) * f12(q, d, k2, k3, k4) / MS_V_3layer(q, d, k1, k2, k3, k4);
}

double V13_3layer(double q, double d, double k1, double k2, double k3, double k4)
{
    return heso_3layer(q) * f13(q, d, k2, k3) / MS_V_3layer(q, d, k1, k2, k3, k4);
}

double V23_3layer(double q, double d, double k1, double k2, double k3, double k4)
{
    return heso_3layer(q) * f23(q, d, k1, k2, k3) / MS_V_3layer(q, d, k1, k2, k3, k4);
}


#endif
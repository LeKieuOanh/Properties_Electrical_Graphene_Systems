#ifndef V_2layer_3spacer_h
#define V_2layer_3spacer_h

#include <stdio.h>
#include <math.h>

#define elec      4.80320427*1e-10 // charge electron 
#define pi        3.141592653589793


// layer 1: MLG, layer2: MLG/BLG/2DEG
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


#endif

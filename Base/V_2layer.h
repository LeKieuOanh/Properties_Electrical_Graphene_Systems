#ifndef V_2layer_h
#define V_2layer_h

#include <stdio.h>
#include <math.h>

#define elec      4.80320427*1e-10 // charge electron 
#define pi        3.141592653589793

double heso(double n, double q){
    return 4*pi*pow(elec,2)/ (q * sqrt(pi * n)) ;
}

double MS(double q,double e1, double e2, double e3, double d, double n){
    return (e2*(e1+e3) + ( pow(e2,2) + e1*e3)*tanh(q * sqrt(pi * n) * d));
}

double V11(double q,double e1, double e2, double e3, double d, double n){
    return heso(n, q) * (e2 + e3 * tanh(q * sqrt(pi * n)*d)) / MS(q, e1, e2, e3, d, n);
}

double V22(double q,double e1, double e2, double e3, double d, double n){
    return heso(n, q) * (e2 + e1*tanh(q*sqrt(pi * n)*d)) / MS(q, e1, e2, e3, d, n);
}

double V12(double q,double e1, double e2, double e3, double d, double n){
    return heso(n, q) * e2 / (cosh(q * sqrt(pi * n)*d) *MS(q, e1, e2, e3, d, n));
}


#endif
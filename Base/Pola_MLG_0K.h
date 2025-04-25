#ifndef Pola_MLG_0K_h
#define Pola_MLG_0K_h

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define pi        3.141592653589793
#define hbar      1.054571628*1e-27 //cm2g/s
#define msao      0.033*9.10938215*1e-28 //the bare electron mass : gams
#define vF        1e8 


double PIMLG(double q1, double n1)     
{
    double N01 = (2*sqrt(n1)) / (sqrt(pi) * hbar *vF);
    if (q1 <= 2.)
    return -N01;
    else
    return -N01 * (1 + pi*q1 /8. - sqrt(q1*q1 - 4) / (2. * q1) - (q1 * asin(2./q1) / 4.)) ; 
    
}

#endif
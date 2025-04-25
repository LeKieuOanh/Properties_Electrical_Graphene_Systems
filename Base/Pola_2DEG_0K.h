#ifndef Pola_2DEG_0K_h
#define Pola_2DEG_0K_h

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define pi        3.141592653589793
#define hbar      1.054571628*1e-27 //cm2g/s


double PI2DEG(double q1, double m2DEG)
{
    
    double N02 = (m2DEG) / (pi * hbar * hbar);
    if (q1 <= 2.)
    return -N02;
    else
    return -N02 * (1 - sqrt(q1*q1 - 4) / q1) ; 
    
}

#endif
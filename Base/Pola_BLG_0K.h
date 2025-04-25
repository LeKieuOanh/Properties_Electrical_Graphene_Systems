#ifndef Pola_BLG_0K_h
#define Pola_BLG_0K_h

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define pi        3.141592653589793
#define hbar      1.054571628*1e-27 //cm2g/s
#define msao      0.033*9.10938215*1e-28 //the bare electron mass : gams


// Calculate polarization for BLG
double f0(double q1)
{
    return ((2. + q1*q1) * sqrt(q1*q1 - 4.) / (2.*q1) + log((q1 - sqrt(q1*q1 - 4.)) / (q1 + sqrt(q1*q1 - 4.)) ) );
}

double g0(double q1)
{
    return (( sqrt(4. + pow(q1,4)) / 2.) - log((1. + sqrt( 1. + pow(q1,4)/4.)) / 2.)) ;
}


double PIBLG(double q1)
{
    double N0  = ( 2. * msao ) / (pi * pow(hbar,2));
    if (q1 <= 2.)
        return -N0 * ( g0(q1));
    else 
       return -N0 * ( g0(q1) - f0(q1));
    
}

#endif
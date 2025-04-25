#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include"tau_BLG_system.h"

double sigma_2BLG(double ni1, double d, double e1, double e2, double e3, double n1, double n2, double r0)
{
 
    double N01  = ( 2. * msao ) / (pi * pow(hbar,2));
    double EF = pi*n1*hbar*hbar / (2*msao);
    return (N01 * EF * hbar * 2. * pi ) / ( msao * tau_2BLG(ni1,d,e1,e2,e3,n1,n2,r0)) ;            
}

double sigma_BLG_MLG(double ni1, double d, double e1, double e2, double e3, double n1, double n2, double r0)
{
 
    double N01  = ( 2. * msao ) / (pi * pow(hbar,2));
    double EF = pi*n1*hbar*hbar / (2*msao);
    return (N01 * EF * hbar * 2. * pi ) / ( msao * tau_BLG_MLG(ni1,d,e1,e2,e3,n1,n2,r0)) ;            
}

double sigma_BLG_2DEG(double ni1, double d, double e1, double e2, double e3, double n1, double n2, double r0, double m2DEG)
{
 
    double N01  = ( 2. * msao ) / (pi * pow(hbar,2));
    double EF = pi*n1*hbar*hbar / (2*msao);
    return (N01 * EF * hbar * 2. * pi ) / ( msao * tau_BLG_2DEG(ni1,d,e1,e2,e3,n1,n2,r0,m2DEG)) ;            
}
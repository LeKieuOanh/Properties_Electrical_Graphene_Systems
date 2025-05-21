#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include"tau_MLG_system.h"

double sigma_2MLG(double ni1, double d, double e1, double e2, double e3, double n1, double n2, double r0)
{
    double EF = hbar * vF * sqrt(pi*n1);
    return (2 * EF) / ( hbar * tau_2MLG(ni1, d, e1, e2, e3, n1, n2, r0));             
}

double sigma_MLG_BLG(double ni1, double d, double e1, double e2, double e3, double n1, double n2, double r0)
{
    double EF = hbar * vF * sqrt(pi*n1);
    return (2 * EF) / ( hbar * tau_MLG_BLG(ni1, d, e1, e2, e3, n1, n2, r0));         
}

double sigma_MLG_2DEG(double ni1, double d, double e1, double e2, double e3, double n1, double n2, double r0, double m2DEG)
{
    double EF = hbar * vF * sqrt(pi*n1);
    return (2 * EF) / ( hbar * tau_MLG_2DEG(ni1, d, e1, e2, e3, n1, n2, r0, m2DEG));           
}

double sigma_3spacer_2MLG(double w1, double w2, double n1, double n2, double e1, double e21, double e22, double e23, double e3, double d, double ni2, double r02, double ni1, double r01)
{
    double EF = hbar * vF * sqrt(pi*n1);
    return (2 * EF) / (hbar * tauEf1_2layer_3spacer_2MLG(w1,w2,n1,n2,e1,e21,e22,e23,e3,d,ni2,r02,ni1,r01));
}

double muy_3spacer_2MLG(double w1, double w2, double n1, double n2, double e1, double e21, double e22, double e23, double e3, double d, double ni2, double r02, double ni1, double r01)
{
    return sigma_3spacer_2MLG(w1,w2,n1,n2,e1,e21,e22,e23,e3,d,ni2,r02,ni1,r01) / (n1 *elec) ;
}
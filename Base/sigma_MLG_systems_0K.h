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
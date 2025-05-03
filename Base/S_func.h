#ifndef S_func_h
#define S_func_h

#include <stdio.h>
#include <math.h>

// #define elec      4.80320427*1e-10 // charge electron
#define pi 3.141592653589793

double prod(int n)
{
    double Pi = 1;
    for (int i = 1; i <= n; i++)
    {
        Pi = Pi * ((double)i);
    }
    return Pi;
}

double gam(int n)
{
    int i;
    double prod = 1;
    for (i = n; i > 1; i--)
    {
        prod = prod * ((double)i - 1);
    }
    return prod;
}

double J(double x)
{
    double s = 0, s1 = 1, eps = 1e-14;
    int i = 0;

    while (fabs(s1) > fabs(s) * eps)
    {
        s1 = (double)pow(-1, i) * pow(x / 2, 2 * i + 1) / (gam(i + 2) * prod(i));
        s = s + s1;
        i = i + 1;
    }
    return s;
}

double S(double q, double ni, double n, double r0)
{
    return 1. - 2 * pi * ni * r0 * J(q * r0 * sqrt(pi * n)) / (q * sqrt(pi * n));
}

#endif
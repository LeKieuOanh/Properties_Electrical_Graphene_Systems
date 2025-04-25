#include <stdio.h>
#include <math.h>

#define elec 4.80320427 * 1e-10
#define pi 3.141592653589793

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

double heso(double q)
{
    return 2 * pi * pow(elec, 2) / q;
}

double MS_V(double q, double d, double k1, double k2, double k3, double k4)
{
    return (k1 - k2) * (k2 + k3) * (k3 - k4) + 2 * exp(2 * q * d) * (k2 - k3) * (k1 * k3 - k2 * k4) + exp(4 * q * d) * (k1 + k2) * (k2 + k3) * (k3 + k4);
}

double V11(double q, double d, double k1, double k2, double k3, double k4)
{
    return heso(q) * f11(q, d, k2, k3, k4) / MS_V(q, d, k1, k2, k3, k4);
}

double V22(double q, double d, double k1, double k2, double k3, double k4)
{
    return heso(q) * f22(q, d, k1, k2, k3, k4) / MS_V(q, d, k1, k2, k3, k4);
}

double V33(double q, double d, double k1, double k2, double k3, double k4)
{
    return heso(q) * f33(q, d, k1, k2, k3) / MS_V(q, d, k1, k2, k3, k4);
}

double V12(double q, double d, double k1, double k2, double k3, double k4)
{
    return heso(q) * f12(q, d, k2, k3, k4) / MS_V(q, d, k1, k2, k3, k4);
}

double V13(double q, double d, double k1, double k2, double k3, double k4)
{
    return heso(q) * f13(q, d, k2, k3) / MS_V(q, d, k1, k2, k3, k4);
}

double V23(double q, double d, double k1, double k2, double k3, double k4)
{
    return heso(q) * f23(q, d, k1, k2, k3) / MS_V(q, d, k1, k2, k3, k4);
}
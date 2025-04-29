#ifndef W_3layer_3BLG_T_h
#define W_3layer_3BLG_T_h

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "Pola_BLG_T.h"
#include "V_3layer.h"

#define hbar 1.054571628 * 1e-27        // cm2g/s
#define msao 0.033 * 9.10938215 * 1e-28 // the bare electron mass : gams
#define pi 3.141592653589793
#define gsgv 4

double MS_3BLG_T(double q, double n1, double d, double e1, double e2, double e3, double e4, double t)
{
    double N0 = (2. * msao) / (pi * pow(hbar, 2));
    return (-1 + V33(q * sqrt(n1), d, e1, e2, e3, e4) * polarinfi(q, t) * N0 + V22(q * sqrt(n1), d, e1, e2, e3, e4) * polarinfi(q, t) * N0 + V11(q * sqrt(n1), d, e1, e2, e3, e4) * polarinfi(q, t) * N0 - V22(q * sqrt(n1), d, e1, e2, e3, e4) * V33(q * sqrt(n1), d, e1, e2, e3, e4) * polarinfi(q, t) * N0 * polarinfi(q, t) * N0 - V11(q * sqrt(n1), d, e1, e2, e3, e4) * V33(q * sqrt(n1), d, e1, e2, e3, e4) * polarinfi(q, t) * N0 * polarinfi(q, t) * N0 - V11(q * sqrt(n1), d, e1, e2, e3, e4) * V22(q * sqrt(n1), d, e1, e2, e3, e4) * polarinfi(q, t) * N0 * polarinfi(q, t) * N0 + V11(q * sqrt(n1), d, e1, e2, e3, e4) * V22(q * sqrt(n1), d, e1, e2, e3, e4) * V33(q * sqrt(n1), d, e1, e2, e3, e4) * polarinfi(q, t) * N0 * polarinfi(q, t) * N0 * polarinfi(q, t) * N0 + V12(q * sqrt(n1), d, e1, e2, e3, e4) * V12(q * sqrt(n1), d, e1, e2, e3, e4) * polarinfi(q, t) * N0 * polarinfi(q, t) * N0 - V12(q * sqrt(n1), d, e1, e2, e3, e4) * V12(q * sqrt(n1), d, e1, e2, e3, e4) * V33(q * sqrt(n1), d, e1, e2, e3, e4) * polarinfi(q, t) * N0 * polarinfi(q, t) * N0 * polarinfi(q, t) * N0 + 2 * V12(q * sqrt(n1), d, e1, e2, e3, e4) * V23(q * sqrt(n1), d, e1, e2, e3, e4) * V13(q * sqrt(n1), d, e1, e2, e3, e4) * polarinfi(q, t) * N0 * polarinfi(q, t) * N0 * polarinfi(q, t) * N0 + V13(q * sqrt(n1), d, e1, e2, e3, e4) * V13(q * sqrt(n1), d, e1, e2, e3, e4) * polarinfi(q, t) * N0 * polarinfi(q, t) * N0 - V13(q * sqrt(n1), d, e1, e2, e3, e4) * V13(q * sqrt(n1), d, e1, e2, e3, e4) * V22(q * sqrt(n1), d, e1, e2, e3, e4) * polarinfi(q, t) * N0 * polarinfi(q, t) * N0 * polarinfi(q, t) * N0 + V23(q * sqrt(n1), d, e1, e2, e3, e4) * V23(q * sqrt(n1), d, e1, e2, e3, e4) * polarinfi(q, t) * N0 * polarinfi(q, t) * N0 - V11(q * sqrt(n1), d, e1, e2, e3, e4) * V23(q * sqrt(n1), d, e1, e2, e3, e4) * V23(q * sqrt(n1), d, e1, e2, e3, e4) * polarinfi(q, t) * N0 * polarinfi(q, t) * N0 * polarinfi(q, t) * N0);
}

double W11_3BLG_T(double q, double n1, double d, double e1, double e2, double e3, double e4, double t)
{
    double N0 = (2. * msao) / (pi * pow(hbar, 2));
    double C11 = (polarinfi(q, t) * N0 * V22(q * sqrt(n1), d, e1, e2, e3, e4) + polarinfi(q, t) * N0 * V33(q * sqrt(n1), d, e1, e2, e3, e4) - polarinfi(q, t) * N0 * polarinfi(q, t) * N0 * V22(q * sqrt(n1), d, e1, e2, e3, e4) * V33(q * sqrt(n1), d, e1, e2, e3, e4) + polarinfi(q, t) * N0 * polarinfi(q, t) * N0 * V23(q * sqrt(n1), d, e1, e2, e3, e4) * V23(q * sqrt(n1), d, e1, e2, e3, e4) - 1);
    double C12 = -(polarinfi(q, t) * N0 * V12(q * sqrt(n1), d, e1, e2, e3, e4) - polarinfi(q, t) * N0 * polarinfi(q, t) * N0 * V12(q * sqrt(n1), d, e1, e2, e3, e4) * V33(q * sqrt(n1), d, e1, e2, e3, e4) + polarinfi(q, t) * N0 * polarinfi(q, t) * N0 * V13(q * sqrt(n1), d, e1, e2, e3, e4) * V23(q * sqrt(n1), d, e1, e2, e3, e4));
    double C13 = -(polarinfi(q, t) * N0 * V13(q * sqrt(n1), d, e1, e2, e3, e4) + polarinfi(q, t) * N0 * polarinfi(q, t) * N0 * V12(q * sqrt(n1), d, e1, e2, e3, e4) * V23(q * sqrt(n1), d, e1, e2, e3, e4) - polarinfi(q, t) * N0 * polarinfi(q, t) * N0 * V13(q * sqrt(n1), d, e1, e2, e3, e4) * V22(q * sqrt(n1), d, e1, e2, e3, e4));
    double V11SC = (C11 * V11(q * sqrt(n1), d, e1, e2, e3, e4) + C12 * V12(q * sqrt(n1), d, e1, e2, e3, e4) + C13 * V13(q * sqrt(n1), d, e1, e2, e3, e4));
    return V11SC / MS_3BLG_T(q, n1, d, e1, e2, e3, e4, t);
}

double W22_3BLG_T(double q, double n1, double d, double e1, double e2, double e3, double e4, double t)
{
    double N0 = (2. * msao) / (pi * pow(hbar, 2));
    double C21 = (polarinfi(q, t) * N0 * V11(q * sqrt(n1), d, e1, e2, e3, e4) + polarinfi(q, t) * N0 * V33(q * sqrt(n1), d, e1, e2, e3, e4) - polarinfi(q, t) * N0 * polarinfi(q, t) * N0 * V11(q * sqrt(n1), d, e1, e2, e3, e4) * V33(q * sqrt(n1), d, e1, e2, e3, e4) + polarinfi(q, t) * N0 * polarinfi(q, t) * N0 * V13(q * sqrt(n1), d, e1, e2, e3, e4) * V13(q * sqrt(n1), d, e1, e2, e3, e4) - 1);
    double C22 = -(polarinfi(q, t) * N0 * V12(q * sqrt(n1), d, e1, e2, e3, e4) - polarinfi(q, t) * N0 * polarinfi(q, t) * N0 * V12(q * sqrt(n1), d, e1, e2, e3, e4) * V33(q * sqrt(n1), d, e1, e2, e3, e4) + polarinfi(q, t) * N0 * polarinfi(q, t) * N0 * V23(q * sqrt(n1), d, e1, e2, e3, e4) * V13(q * sqrt(n1), d, e1, e2, e3, e4));
    double C23 = -(polarinfi(q, t) * N0 * V23(q * sqrt(n1), d, e1, e2, e3, e4) - polarinfi(q, t) * N0 * polarinfi(q, t) * N0 * V11(q * sqrt(n1), d, e1, e2, e3, e4) * V23(q * sqrt(n1), d, e1, e2, e3, e4) + polarinfi(q, t) * N0 * polarinfi(q, t) * N0 * V13(q * sqrt(n1), d, e1, e2, e3, e4) * V12(q * sqrt(n1), d, e1, e2, e3, e4));
    double V22SC = (C21 * V22(q * sqrt(n1), d, e1, e2, e3, e4) + C22 * V12(q * sqrt(n1), d, e1, e2, e3, e4) + C23 * V23(q * sqrt(n1), d, e1, e2, e3, e4));
    return V22SC / MS_3BLG_T(q, n1, d, e1, e2, e3, e4, t);
}

double W33_3BLG_T(double q, double n1, double d, double e1, double e2, double e3, double e4, double t)
{
    double N0 = (2. * msao) / (pi * pow(hbar, 2));
    double C13 = (V22(q * sqrt(n1), d, e1, e2, e3, e4) * V13(q * sqrt(n1), d, e1, e2, e3, e4) * polarinfi(q, t) * N0 * polarinfi(q, t) * N0 - V13(q * sqrt(n1), d, e1, e2, e3, e4) * polarinfi(q, t) * N0 - V12(q * sqrt(n1), d, e1, e2, e3, e4) * V23(q * sqrt(n1), d, e1, e2, e3, e4) * polarinfi(q, t) * N0 * polarinfi(q, t) * N0);
    double C23 = (V11(q * sqrt(n1), d, e1, e2, e3, e4) * V23(q * sqrt(n1), d, e1, e2, e3, e4) * polarinfi(q, t) * N0 * polarinfi(q, t) * N0 - V23(q * sqrt(n1), d, e1, e2, e3, e4) * polarinfi(q, t) * N0 - V12(q * sqrt(n1), d, e1, e2, e3, e4) * V13(q * sqrt(n1), d, e1, e2, e3, e4) * polarinfi(q, t) * N0 * polarinfi(q, t) * N0);
    double C33 = (V12(q * sqrt(n1), d, e1, e2, e3, e4) * V12(q * sqrt(n1), d, e1, e2, e3, e4) * polarinfi(q, t) * N0 * polarinfi(q, t) * N0 + V22(q * sqrt(n1), d, e1, e2, e3, e4) * polarinfi(q, t) * N0 - V11(q * sqrt(n1), d, e1, e2, e3, e4) * V22(q * sqrt(n1), d, e1, e2, e3, e4) * polarinfi(q, t) * N0 * polarinfi(q, t) * N0 + V11(q * sqrt(n1), d, e1, e2, e3, e4) * polarinfi(q, t) * N0 - 1);
    double V33SC = (C13 * V13(q * sqrt(n1), d, e1, e2, e3, e4) + C23 * V23(q * sqrt(n1), d, e1, e2, e3, e4) + C33 * V33(q * sqrt(n1), d, e1, e2, e3, e4));
    return V33SC / MS_3BLG_T(q, n1, d, e1, e2, e3, e4, t);
}


#endif
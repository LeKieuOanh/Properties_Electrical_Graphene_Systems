#ifndef df0dE_h
#define df0dE_h

#include <stdio.h>
#include <math.h>

// E * df0 / dE
double df1(double x, double y) // x = E/E0, y = T/TF;  E*(-df0/dE) (MS)
{
     double s = 1;
     double kq, hs;
     hs = (x - s) / y;
     if (hs <= 700.)
          kq = exp(hs) * x / (pow(1 + exp(hs), 2) * y);
     else
          kq = x / (y * exp(700.));
     return kq;
}

// E^2 * df0 / dE
double df2(double x, double y) // x = E/E0, y = T/TF ; E^2*(-df0/dE) (TS)
{
     double s = 1;
     double kq, hs;
     hs = (x - s) / y;
     if (hs <= 700.)
          kq = exp(hs) * x * x / (pow(1 + exp(hs), 2) * y);
     else
          kq = x / (y * exp(700.));
     return kq;
}

#endif
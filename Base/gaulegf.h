#ifndef gaulegf_h
#define gaulegf_h

#include<stdio.h>
#include<math.h>

void gaulegf(double x1, double x2, double x[], double w[], int n)
{
  int i, j, m;
  double eps = 3.0E-14;
  double p1, p2, p3, pp, xl, xm, z, z1;
  
  m = (n+1)/2;
  xm = 0.5*(x2+x1);
  xl = 0.5*(x2-x1);
  for(i=1; i<=m; i++)
  {
    z = cos(3.141592654*((double)i-0.25)/((double)n+0.5));
    while(1)
    {
      p1 = 1.0;
      p2 = 0.0;
      for(j=1; j<=n; j++)
      {
        p3 = p2;
        p2 = p1;
        p1 = ((2.0*(double)j-1.0)*z*p2-((double)j-1.0)*p3)/
              (double)j;
      }
      pp = (double)n*(z*p1-p2)/(z*z-1.0);
      z1 = z;
      z = z1 - p1/pp;
      if(fabs(z-z1) <= eps) break;
    }
    x[i] = xm - xl*z;
    x[n+1-i] = xm + xl*z;
    w[i] = 2.0*xl/((1.0-z*z)*pp*pp);
    w[n+1-i] = w[i];
  }
}

#endif
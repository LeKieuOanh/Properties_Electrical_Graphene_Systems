#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "df0dE.h"
#include "../Base/tau_BLG_system.h"

int main()
{
    double E, kq, ni, d, e1, e2, e3, e4, n, r0, t;
    E = 85;
    r0 = 0;
    n = 1e12;
    // n2 = 1e12;
    // t = 0.02;
    ni = 1e11;
    d = 1e-7;
    e1 = e2 = e3 = e4 = 1;
    FILE *output;
    output = fopen("./result/data/3BLG-sigma-T-0.01-0.1.txt", "a+");
    for (t = 0.03; t <= 0.03; t += 0.01)
    {
        kq = 0;
        for (int i = 0.; i <= E; i += 1)
        {
            kq += df1(i, t) * tau11_3BLG_T(ni, d, e1, e2, e3, e4, n, r0, t);
        }
        fprintf(output, "%f %.20f\n", t, kq);
    }

    fclose(output);
    return 0;
}

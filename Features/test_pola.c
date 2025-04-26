

//=============THE LIBRARIES DECLARATION==========
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "../Base/gaulegf.h"
#include "../Base/Pola_BLG_0K.h"


int main()
{
 
    double a0;
    double n1;
    double d;
    a0 = 4.92 * 1e-8 ;
    n1 = 3e12;
    d = 100e-7;   
     
    FILE *output;
    output = fopen("test_Pola_BLG.txt","w+");
    double x[1000], w[1000];
    gaulegf(0,4, x, w, 100.);
    double N0  = ( 2. * msao ) / (pi * pow(hbar,2));
    for ( int i = 0 ; i <= 100; i++)
    {
        fprintf(output,"%f\t %f\n",x[i],-PIBLG(x[i])/N0);
    }
     fclose(output);

    
    return 0;
}
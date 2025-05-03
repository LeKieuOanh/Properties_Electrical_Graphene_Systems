#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include"../Base/electrical_BLG_systems_0K.h"

int main()
{
    double  a0, r0,l, ni1, n1, n2;
    // double N0  = ( 2. * msao ) / (pi * pow(hbar,2));
    a0 = 4.92 * 1e-8 ;
    l = 1e-7;
    n2 = 1e12;
    double r0_list[5] = {0., 5.*a0, 7*a0, 8*a0, 10.*a0};
    double d_list[3] = {1*l, 3*l, 100*l};
    double d, e1, e2, e3;
    e1 = 1;
    e2 = 4.0; 
    e3 = 12.53;   
    FILE *output;
    for (int k = 0; k < 3; k++)
    {
        d = d_list[k];
        for (int i= 0; i<=4; i++)
    {
        ni1 = 0.95*1e12;
        r0 = r0_list[i];
        char s1[] = "./results/BLG_systems_0K/BLG_MLG/d=";
        char s3[10];
        gcvt(d/l,10,s3);
        char s4[] = "nm/0.95e12/sigmar0=";
        char s2[20];
        gcvt(r0/a0,20,s2);
        char s[100]="";
        strcat(s,s1);
        strcat(s,s3);
        strcat(s,s4);
        strcat(s,s2);
        strcat(s,".txt");
    output = fopen(s,"w+");

         for (n1 = 0; n1 <= 5e12; n1 += 0.02e12)
    {
        fprintf(output,"%f \t %.5e \n",n1,sigma_BLG_MLG(ni1,d,e1,e2,e3,n1,n2,r0));
    }
    }
       for (int i= 0; i<=4; i++)
    {
        ni1 = 0.5 * 1e12;
        r0 = r0_list[i];
        char s1[] = "./results/BLG_systems_0K/BLG_MLG/d=";
        char s3[10];
        gcvt(d/l,10,s3);
        char s4[] = "nm/0.5e12/sigmar0=";
        char s2[20];
        gcvt(r0/a0,20,s2);
        char s[100]="";
        strcat(s,s1);
        strcat(s,s3);
        strcat(s,s4);
        strcat(s,s2);
        strcat(s,".txt");
    output = fopen(s,"w+");

         for (n1 = 0; n1 <= 5e12; n1 += 0.02e12)
    {
        fprintf(output,"%f \t %.5e\n",n1,sigma_BLG_MLG(ni1,d,e1,e2,e3,n1,n2,r0));
    }
    }
    }
    
    
    fclose(output);
    return 0;
}
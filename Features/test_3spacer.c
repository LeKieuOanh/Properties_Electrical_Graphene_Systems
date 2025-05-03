#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include"../Base/electrical_MLG_systems_0K.h"

int main()
{
    // 2layer - 3spacer - 2MLG
    double a0, e1, e21, e22, e23, e3, n1, n2, ni1, ni2, r01, r02, ni, d, w1, w2;
    a0 = 4.92 * 1e-8;
    e1 = 1.0;
    e3 = 22;
    n1 = n2 = 1.e12;
    ni = 0.5e11;
    r01 = r02 = 10* a0;
    ni1 = ni2 = 0.5e12;
    e21 = 4;
    e22 = 10;
    d = 3e-7;
    double e23_list[3] = {10, 14, 16};
    w1 = 1e-7;
    w2 = 2e-7;
    FILE *output;	

    // Calculate for d, w
        // Calculate for e21, e22
        for (int k =0; k<=2; k++)
        {
            e23 = e23_list[k];
            // Calcalate for r0 different;
            for (int j=0; j<=4; j++)
            {
                ni = 0.5e12;
                ni1 = ni2 = ni; 
                char s1[] = "./results/MLG_systems_0K/MLG_MLG/3spacer/d-r01=";
                char s2[20];
                gcvt(r01/a0, 20, s2);
                char s3[] = "-r02=";
                char s4[20];
                gcvt(r02/a0, 20, s4);
                char s5[] = "-e21=";
                char s6[20];
                gcvt(e21, 20, s6);
                char s7[] = "-d=";
                char s8[20];
                gcvt(d/1e-7, 20, s8);
                char s9[] = "-e22=";
                char s10[20];
                gcvt(e22, 20, s10);
                char s11[] = "-e23=";
                char s12[20];
                gcvt(e23, 20, s12);

                char s[100]="";
                strcat(s,s1);
                strcat(s,s2);
                strcat(s,s3);
                strcat(s,s4);                        
                strcat(s,s5);
                strcat(s,s6);
                strcat(s,s7);
                strcat(s,s8);
                strcat(s,s9);
                strcat(s,s10);
                strcat(s,s11);
                strcat(s,s12);
                strcat(s,".txt");
                output = fopen(s,"w+");

                for (n1 = 0; n1 <= 5e12; n1 += 0.5e11)
                {
                    n2 = n1;
                    fprintf(output,"%f \t %f \t %f \n",n1,muy_2layer_3spacer_2MLG(w1,w2,n1,n2,e1,e21,e22,e23,e3,d,ni2,r02,ni1,r01), sigma_2layer_3spacer_2MLG(w1,w2,n1,n2,e1,e21,e22,e23,e3,d,ni2,r02,ni1,r01));
                }
            }
        }    

    fclose(output);
    return 0;
}
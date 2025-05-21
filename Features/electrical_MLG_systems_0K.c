#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include"../Base/electrical_MLG_systems_0K.h"

int main()
{
    // 2 layer - 1spacer (**********uncomment went use*********)
    // double  a0, r0,l, ni1, n1, n2;
    // a0 = 4.92 * 1e-8 ;
    // l = 1e-7;
    // n2 = 1e12;
    // double r0_list[5] = {0., 5.*a0, 7*a0, 8*a0, 10.*a0};
    // double d_list[3] = {1*l, 3*l, 100*l};
    // double d, e1, e2, e3;
    // e1 = 1;
    // e2 = 4.0; 
    // e3 = 12.53;   
    // FILE *output;

    // // MLG-MLG
    // for (int k = 0; k < 3; k++)
    // {
    //     d = d_list[k];
    //     for (int i= 0; i<=4; i++)
    // {
    //     ni1 = 0.95*1e12;
    //     r0 = r0_list[i];
    //     char s1[] = "./results/MLG_systems_0K/MLG_MLG/d=";
    //     char s3[10];
    //     gcvt(d/l,10,s3);
    //     char s4[] = "nm/0.95e12/sigmar0=";
    //     char s2[20];
    //     gcvt(r0/a0,20,s2);
    //     char s[100]="";
    //     strcat(s,s1);
    //     strcat(s,s3);
    //     strcat(s,s4);
    //     strcat(s,s2);
    //     strcat(s,".txt");
    // output = fopen(s,"w+");

    //      for (n1 = 0; n1 <= 5e12; n1 += 0.02e12)
    // {
    //     fprintf(output,"%f \t %.5e \n",n1,sigma_2MLG(ni1,d,e1,e2,e3,n1,n2,r0));
    // }
    // }
    //    for (int i= 0; i<=4; i++)
    // {
    //     ni1 = 0.5 * 1e12;
    //     r0 = r0_list[i];
    //     char s1[] = "./results/MLG_systems_0K/MLG_MLG/d=";
    //     char s3[10];
    //     gcvt(d/l,10,s3);
    //     char s4[] = "nm/0.5e12/sigmar0=";
    //     char s2[20];
    //     gcvt(r0/a0,20,s2);
    //     char s[100]="";
    //     strcat(s,s1);
    //     strcat(s,s3);
    //     strcat(s,s4);
    //     strcat(s,s2);
    //     strcat(s,".txt");
    // output = fopen(s,"w+");

    //      for (n1 = 0; n1 <= 5e12; n1 += 0.02e12)
    // {
    //     fprintf(output,"%f \t %.5e\n",n1,sigma_2MLG(ni1,d,e1,e2,e3,n1,n2,r0));
    // }
    // }
    // }

    // // MLG-BLG
    // for (int k = 0; k < 3; k++)
    // {
    //     d = d_list[k];
    //     for (int i= 0; i<=4; i++)
    // {
    //     ni1 = 0.95*1e12;
    //     r0 = r0_list[i];
    //     char s1[] = "./results/MLG_systems_0K/MLG_BLG/d=";
    //     char s3[10];
    //     gcvt(d/l,10,s3);
    //     char s4[] = "nm/0.95e12/sigmar0=";
    //     char s2[20];
    //     gcvt(r0/a0,20,s2);
    //     char s[100]="";
    //     strcat(s,s1);
    //     strcat(s,s3);
    //     strcat(s,s4);
    //     strcat(s,s2);
    //     strcat(s,".txt");
    // output = fopen(s,"w+");

    //      for (n1 = 0; n1 <= 5e12; n1 += 0.02e12)
    // {
    //     fprintf(output,"%f \t %.5e \n",n1,sigma_MLG_BLG(ni1,d,e1,e2,e3,n1,n2,r0));
    // }
    // }
    //    for (int i= 0; i<=4; i++)
    // {
    //     ni1 = 0.5 * 1e12;
    //     r0 = r0_list[i];
    //     char s1[] = "./results/MLG_systems_0K/MLG_BLG/d=";
    //     char s3[10];
    //     gcvt(d/l,10,s3);
    //     char s4[] = "nm/0.5e12/sigmar0=";
    //     char s2[20];
    //     gcvt(r0/a0,20,s2);
    //     char s[100]="";
    //     strcat(s,s1);
    //     strcat(s,s3);
    //     strcat(s,s4);
    //     strcat(s,s2);
    //     strcat(s,".txt");
    // output = fopen(s,"w+");

    //      for (n1 = 0; n1 <= 5e12; n1 += 0.02e12)
    // {
    //     fprintf(output,"%f \t %.5e\n",n1,sigma_MLG_BLG(ni1,d,e1,e2,e3,n1,n2,r0));
    // }
    // }
    // }
    
    // // 2layer - 3spacer - 2MLG
    // double a0, e1, e21, e22, e23, e3, n1, n2, ni1, ni2, r01, r02, ni, d, w1, w2;
    // a0 = 4.92 * 1e-8;
    // e1 = 1.0;
    // e3 = 22;
    // n1 = n2 = 1.e12;
    // ni = 0.5e11;
    // d = 3e-7;
    // double r1_list[5] = {5*a0, 10*a0, 5*a0, 10*a0, 0};
    // double r2_list[5] = {5*a0, 5*a0, 10*a0, 10*a0, 0};
    // double ni1_list[5] = {ni, ni, 5*ni, ni, 10*ni};
    // double ni2_list[5] = {ni, 5*ni, ni, 10*ni, ni};
    // e21 = 4;
    // e22 = 10;
    // double e23_list[4] = {10, 12, 14, 16};
    // w1 = 1e-7;
    // w2 = 2e-7;
    // FILE *output;	

    // // Calculate for e21, e22
    // for (int k =0; k<=3; k++)
    // {
    //     e23 = e23_list[k];
    //     // Calcalate for r0 different;
    //     for (int j=0; j<=4; j++)
    //     {
    //         ni = 0.5e12;
    //         ni1 = ni2 = ni; 
    //         char s1[] = "./results/MLG_systems_0K/MLG_MLG/3spacer/d-ni1=ni2=0.5e12-r01=";
    //         char s2[20];
    //         gcvt(r1_list[j]/a0, 20, s2);
    //         char s3[] = "-r02=";
    //         char s4[20];
    //         gcvt(r2_list[j]/a0, 20, s4);
    //         char s5[] = "-e21=";
    //         char s6[20];
    //         gcvt(e21, 20, s6);
    //         char s9[] = "-e22=";
    //         char s10[20];
    //         gcvt(e22, 20, s10);
    //         char s11[] = "-e23=";
    //         char s12[20];
    //         gcvt(e23, 20, s12);

    //         char s[100]="";
    //         strcat(s,s1);
    //         strcat(s,s2);
    //         strcat(s,s3);
    //         strcat(s,s4);                        
    //         strcat(s,s5);
    //         strcat(s,s6);
    //         strcat(s,s9);
    //         strcat(s,s10);
    //         strcat(s,s11);
    //         strcat(s,s12);
    //         strcat(s,".txt");
    //         output = fopen(s,"w+");

    //         for (n1 = 0; n1 <= 5e12; n1 += 0.5e11)
    //         {
    //             n2 = n1;
    //             fprintf(output,"%f \t %f \t %f \n",n1,muy_2layer_3spacer_2MLG(w1,w2,n1,n2,e1,e21,e22,e23,e3,d,ni2,r2_list[j],ni1,r1_list[j]), sigma_2layer_3spacer_2MLG(w1,w2,n1,n2,e1,e21,e22,e23,e3,d,ni2,r2_list[j],ni1,r1_list[j]));
    //         }
    //     }

    //     // Calcalate for ni different;

    //     for (int j=0; j<= 4; j++)
    //     {
    //         ni = 0.5e11;
    //         r01 = r02 = 10*a0;    
    //         char s1[] = "./results/MLG_systems_0K/MLG_MLG/3spacer/d-r01=r02=10a0-ni1=";
    //         char s2[20];
    //         gcvt(ni1_list[j]/ni, 20, s2);
    //         char s3[] = "-ni2=";
    //         char s4[20];
    //         gcvt(ni2_list[j]/ni, 20, s4);
    //         char s5[] = "-e21=";
    //         char s6[20];
    //         gcvt(e21, 20, s6);
    //         char s7[] = "-d=";
    //         char s8[20];
    //         gcvt(d/1e-7, 20, s8);
    //         char s9[] = "-e22=";
    //         char s10[20];
    //         gcvt(e22, 20, s10);
    //         char s11[] = "-e23=";
    //         char s12[20];
    //         gcvt(e23, 20, s12);

    //         char s[100]="";
    //         strcat(s,s1);
    //         strcat(s,s2);
    //         strcat(s,s3);
    //         strcat(s,s4);                        
    //         strcat(s,s5);
    //         strcat(s,s6);
    //         strcat(s,s7);
    //         strcat(s,s8);
    //         strcat(s,s9);
    //         strcat(s,s10);
    //         strcat(s,s11);
    //         strcat(s,s12);
    //         strcat(s,".txt");
    //         output = fopen(s,"w+");

    //         for (n1 = 0; n1 <= 5e12; n1 += 0.5e11)
    //         {
    //             n2 = n1;
    //             fprintf(output,"%f \t %f \t %f\n",n1,muy_2layer_3spacer_2MLG(w1,w2,n1,n2,e1,e21,e22,e23,e3,d,ni2_list[j],r02,ni1_list[j],r01), sigma_2layer_3spacer_2MLG(w1,w2,n1,n2,e1,e21,e22,e23,e3,d,ni2_list[j],r02,ni1_list[j],r01));
    //         }
    //     }
    // } 
    
    // 2layer - 2spacer - 2MLG
    double a0, e1, e21, e22, e3, n1, n2, ni1, r01, r02, ni2, d1, w1;
    a0 = 4.92 * 1e-8;
    e1 = 1.0;
    e3 = 12.53;
    ni1 = 0.45e12;
    ni2 = 0;
    n2 = 1e12;
    d1 = 6e-7;
    double r1_list[3] = {0, 5*a0, 10*a0};
    r02 = 0;
    e21 = 4;
    double e22_list[2] = {7.5, 22};
    w1 = 2e-7;
    FILE *output;	

        // Calcalate for r0 different;
        for (int j=0; j<=1; j++)
        {
            e22 = e22_list[j];
            char s1[] = "./results/MLG_systems_0K/MLG_MLG/2spacer/d=6-ni1=0.45e12";
            char s5[] = "-e21=";
            char s6[20];
            gcvt(e21, 20, s6);
            char s9[] = "-e22=";
            char s10[20];
            gcvt(e22, 20, s10);

            char s[100]="";
            strcat(s,s1);
            strcat(s,s5);
            strcat(s,s6);
            strcat(s,s9);
            strcat(s,s10);
            strcat(s,".txt");
            output = fopen(s,"w+");

            for (n1 = 1e12; n1 <= 5e12; n1 += 0.5e11)
            {
                fprintf(output,"%f \t %.5e \t %.5e  \t %.5e \n",n1/1e12,sigma_2spacer_2MLG(w1,n1,n2,e1,e21,e22,e3,d1,ni1,r1_list[0],ni2,r02), sigma_2spacer_2MLG(w1,n1,n2,e1,e21,e22,e3,d1,ni1,r1_list[1],ni2,r02), sigma_2spacer_2MLG(w1,n1,n2,e1,e21,e22,e3,d1,ni1,r1_list[2],ni2,r02));
            }
        }



    fclose(output);
    return 0;
}
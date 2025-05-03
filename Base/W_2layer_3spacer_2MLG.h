#ifndef W_2layer_3spacer_2MLG_h
#define W_2layer_3spacer_2MLG_h

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include <string.h>


#include"Pola_MLG_0K.h"
#include"S_func.h"
#include"V_2layer_3spacer.h"
#include"gaulegf.h"


// unscreen-potential: q -> q*sqrt(pi*n1)
double epsilon_2layer_3spacer_2MLG(double q1, double w1, double w2, double n1, double n2, double e1, double e21, double e22, double e23, double e3, double d)
{
    return (1 - PIMLG(q1, n1)*V11_2layer_3spacer(q1*sqrt(pi*n1),w1, w2 ,e1,e21,e22,e23,e3,d)) * (1 - PIMLG(q1*sqrt(n1/n2),n2)*V22_2layer_3spacer(q1*sqrt(pi*n1),w1, w2 , e1, e21, e22, e23, e3, d)) - V12_2layer_3spacer(q1*sqrt(pi*n1),w1, w2,e1,e21,e22,e23,e3,d)*V12_2layer_3spacer(q1*sqrt(pi*n1),w1,w2,e1,e21,e22,e23,e3,d)*PIMLG(q1,n1)*PIMLG(q1*sqrt(n1/n2),n2);
}

double W11_2layer_3spacer_2MLG(double q1, double w1, double w2, double n1, double n2, double e1, double e21, double e22, double e23, double e3, double d)
{
    return (V11_2layer_3spacer(q1*sqrt(pi*n1),w1,w2,e1,e21,e22,e23,e3,d) + (V12_2layer_3spacer(q1*sqrt(pi*n1),w1, w2,e1,e21,e22,e23, e3,d)*V12_2layer_3spacer(q1*sqrt(pi*n1),w1, w2,e1,e21,e22,e23,e3,d) - V11_2layer_3spacer(q1*sqrt(pi*n1),w1,w2,e1,e21,e22,e23,e3,d)*V22_2layer_3spacer(q1*sqrt(pi*n1),w1,w2,e1,e21,e22,e23,e3,d))*PIMLG(q1*sqrt(n1/n2),n2))/(epsilon_2layer_3spacer_2MLG(q1,w1,w2,n1,n2,e1,e21,e22,e23,e3,d));
}

double W12_2layer_3spacer_2MLG(double q1, double w1, double w2, double n1, double n2, double e1, double e21, double e22, double e23, double e3, double d)
{
    return V12_2layer_3spacer(q1*sqrt(pi*n1),w1,w2,e1,e21,e22,e23,e3,d) / epsilon_2layer_3spacer_2MLG(q1,w1,w2,n1,n2,e1,e21,e22,e23,e3,d);
}

#endif

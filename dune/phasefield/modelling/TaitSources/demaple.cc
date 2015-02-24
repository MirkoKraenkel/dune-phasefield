double evalRho(double x)
{
  {
    return(0.4958034541*x+0.1065766549);
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t16;
  double t17;
  double t18;
  double t20;
  double t22;
  double t29;
  double t3;
  double t4;
  double t7;
  {
    t3 = phi*phi;
    t4 = t3*t3;
    t7 = t3*phi;
    t16 = 6.0*t4*phi;
    t17 = 15.0*t4;
    t18 = 10.0*t7;
    t20 = rho*rho;
    t22 = t20*t20;
    t29 = log(rho);
    return(2.0/delta*A*(t4+(2.0*alpha-2.0)*t7+(-3.0*alpha+1.0)*t3+alpha)+(t16-
t17+t18)*(0.35*t22*t20*rho-0.1170549602*rho+0.6043849679E-1)+(1.0-t16+t17-t18)*
(0.1753186565*rho*t29+0.2172006685*rho+0.186848759E-1));
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t17;
  double t21;
  double t22;
  double t24;
  double t3;
  double t31;
  double t4;
  {
    t3 = phi*phi;
    t4 = t3*phi;
    t17 = t3*t3;
    t21 = 30.0*t17-60.0*t4+30.0*t3;
    t22 = rho*rho;
    t24 = t22*t22;
    t31 = log(rho);
    return(2.0/delta*A*(4.0*t4+3.0*(2.0*alpha-2.0)*t3+2.0*(-3.0*alpha+1.0)*phi)
+t21*(0.35*t24*t22*rho-0.1170549602*rho+0.6043849679E-1)-t21*(0.1753186565*rho*
t31+0.2172006685*rho+0.186848759E-1));
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  double t17;
  double t18;
  double t20;
  double t27;
  double t3;
  {
    t3 = phi*phi;
    t17 = 120.0*t3*phi-180.0*t3+60.0*phi;
    t18 = rho*rho;
    t20 = t18*t18;
    t27 = log(rho);
    return(2.0/delta*A*(12.0*t3+6.0*(2.0*alpha-2.0)*phi-6.0*alpha+2.0)+t17*(
0.35*t20*t18*rho-0.1170549602*rho+0.6043849679E-1)-t17*(0.1753186565*rho*t27+
0.2172006685*rho+0.186848759E-1));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t13;
  double t17;
  double t2;
  double t3;
  double t4;
  double t5;
  double t6;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = t2*phi;
    t4 = rho*rho;
    t5 = t4*t4;
    t6 = t5*t4;
    t13 = t1*phi;
    t17 = log(rho);
    return(0.147E2*t3*t6-0.3057445711E1*t3-0.3675E2*t2*t6+0.7643614278E1*t2+
0.245E2*t13*t6-0.5095742852E1*t13+0.1753186565*t17+0.392519325-0.1051911939E1*
t3*t17+0.2629779848E1*t2*t17-0.1753186565E1*t13*t17);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  double t1;
  double t11;
  double t2;
  double t3;
  double t4;
  double t5;
  double t6;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = t2*phi;
    t4 = rho*rho;
    t5 = t4*t4;
    t6 = t5*t4;
    t11 = t1*phi;
    return(0.5E-9*(0.1764E12*t3*t6-0.441E12*t2*t6+0.294E12*t11*t6+350637313.0
-2103823878.0*t3+5259559695.0*t2-3506373130.0*t11)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t16;
  double t2;
  double t3;
  double t4;
  double t5;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = rho*rho;
    t4 = t3*t3;
    t5 = t4*t3;
    t9 = t1*phi;
    t16 = log(rho);
    return(0.735E2*t2*t5-0.1528722856E2*t2-147.0*t9*t5+0.3057445711E2*t9+
0.735E2*t1*t5-0.1528722856E2*t1-0.5259559695E1*t2*t16+0.1051911939E2*t9*t16
-0.5259559695E1*t1*t16);
  }
}

#include <math.h>
double pressure(double rho,double phi)
{
  double t1;
  double t11;
  double t2;
  double t22;
  double t4;
  double t5;
  double t6;
  double t7;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t6 = t1*phi;
    t7 = 10.0*t6;
    t9 = rho*rho;
    t11 = t9*t9;
    t22 = log(rho);
    return((t4-t5+t7)*(-0.35*t11*t9*rho+0.1170549602*rho-0.6043849679E-1+rho*(
0.245E1*t11*t9-0.1170549602))+(1.0-t4+t5-t7)*(-0.1753186565*rho*t22
-0.2172006685*rho-0.186848759E-1+rho*(0.1753186565*t22+0.392519325))-2.0/delta*
A*(t2+(2.0*alpha-2.0)*t6+(-3.0*alpha+1.0)*t1+alpha));
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t1;
  double t11;
  double t18;
  double t2;
  double t3;
  double t4;
  double t5;
  double t6;
  {
    t1 = rho*rho;
    t2 = t1*t1;
    t3 = t2*t1;
    t4 = phi*phi;
    t5 = t4*t4;
    t6 = t5*phi;
    t11 = t4*phi;
    t18 = sqrt(0.882E12*t3*t6-0.2205E13*t5*t3+0.147E13*t11*t3+1753186565.0
-0.1051911939E11*t6+0.2629779848E11*t5-0.1753186565E11*t11);
    return(0.1E-4*t18);
  }
}

double evalRho(double x)
{
  {
    return(0.4958034541*x+0.1065766549);
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t16;
  double t17;
  double t18;
  double t20;
  double t22;
  double t29;
  double t3;
  double t4;
  double t7;
  {
    t3 = phi*phi;
    t4 = t3*t3;
    t7 = t3*phi;
    t16 = 6.0*t4*phi;
    t17 = 15.0*t4;
    t18 = 10.0*t7;
    t20 = rho*rho;
    t22 = t20*t20;
    t29 = log(rho);
    return(2.0/delta*A*(t4+(2.0*alpha-2.0)*t7+(-3.0*alpha+1.0)*t3+alpha)+(t16-
t17+t18)*(0.8575855211E-1*t22*t20*rho-0.2868132559E-1*rho+0.1480890852E-1)+(1.0
-t16+t17-t18)*(0.1015E1+0.9523661641E1*rho*t29+0.1179877668E2*rho));
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t17;
  double t21;
  double t22;
  double t24;
  double t3;
  double t31;
  double t4;
  {
    t3 = phi*phi;
    t4 = t3*phi;
    t17 = t3*t3;
    t21 = 30.0*t17-60.0*t4+30.0*t3;
    t22 = rho*rho;
    t24 = t22*t22;
    t31 = log(rho);
    return(2.0/delta*A*(4.0*t4+3.0*(2.0*alpha-2.0)*t3+2.0*(-3.0*alpha+1.0)*phi)
+t21*(0.8575855211E-1*t24*t22*rho-0.2868132559E-1*rho+0.1480890852E-1)-t21*(
0.1015E1+0.9523661641E1*rho*t31+0.1179877668E2*rho));
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  double t17;
  double t18;
  double t20;
  double t27;
  double t3;
  {
    t3 = phi*phi;
    t17 = 120.0*t3*phi-180.0*t3+60.0*phi;
    t18 = rho*rho;
    t20 = t18*t18;
    t27 = log(rho);
    return(2.0/delta*A*(12.0*t3+6.0*(2.0*alpha-2.0)*phi-6.0*alpha+2.0)+t17*(
0.8575855211E-1*t20*t18*rho-0.2868132559E-1*rho+0.1480890852E-1)-t17*(0.1015E1+
0.9523661641E1*rho*t27+0.1179877668E2*rho));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t13;
  double t17;
  double t2;
  double t3;
  double t4;
  double t5;
  double t6;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = t2*phi;
    t4 = rho*rho;
    t5 = t4*t4;
    t6 = t5*t4;
    t13 = t1*phi;
    t17 = log(rho);
    return(0.3601859189E1*t3*t6-0.1281067179E3*t3-0.9004647972E1*t2*t6+
0.3202667947E3*t2+0.6003098648E1*t13*t6-0.2135111965E3*t13+0.2132243832E2+
0.9523661641E1*t17-0.5714196985E2*t3*t17+0.1428549246E3*t2*t17-0.9523661641E2*
t13*t17);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  double t1;
  double t11;
  double t2;
  double t3;
  double t4;
  double t5;
  double t6;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = t2*phi;
    t4 = rho*rho;
    t5 = t4*t4;
    t6 = t5*t4;
    t11 = t1*phi;
    return(0.1E-8*(0.2161115513E11*t3*t6-0.5402788784E11*t2*t6+0.3601859189E11*
t11*t6+9523661641.0-0.5714196985E11*t3+0.1428549246E12*t2-0.9523661641E11*t11)/
rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t16;
  double t2;
  double t3;
  double t4;
  double t5;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = rho*rho;
    t4 = t3*t3;
    t5 = t4*t3;
    t9 = t1*phi;
    t16 = log(rho);
    return(0.1800929594E2*t2*t5-0.6405335894E3*t2-0.3601859189E2*t9*t5+
0.1281067179E4*t9+0.1800929594E2*t1*t5-0.6405335894E3*t1-0.2857098492E3*t2*t16+
0.5714196985E3*t9*t16-0.2857098492E3*t1*t16);
  }
}

#include <math.h>
double pressure(double rho,double phi)
{
  double t1;
  double t11;
  double t2;
  double t22;
  double t4;
  double t5;
  double t6;
  double t7;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t6 = t1*phi;
    t7 = 10.0*t6;
    t9 = rho*rho;
    t11 = t9*t9;
    t22 = log(rho);
    return((t4-t5+t7)*(-0.8575855211E-1*t11*t9*rho+0.2868132559E-1*rho
-0.1480890852E-1+rho*(0.6003098648*t11*t9-0.2868132559E-1))+(1.0-t4+t5-t7)*(
-0.1015E1-0.9523661641E1*rho*t22-0.1179877668E2*rho+rho*(0.2132243832E2+
0.9523661641E1*t22))-2.0/delta*A*(t2+(2.0*alpha-2.0)*t6+(-3.0*alpha+1.0)*t1+
alpha));
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t1;
  double t11;
  double t18;
  double t2;
  double t3;
  double t4;
  double t5;
  double t6;
  {
    t1 = rho*rho;
    t2 = t1*t1;
    t3 = t2*t1;
    t4 = phi*phi;
    t5 = t4*t4;
    t6 = t5*phi;
    t11 = t4*phi;
    t18 = sqrt(2161115513.0*t3*t6-5402788784.0*t5*t3+3601859189.0*t11*t3+
952366164.0-5714196984.0*t6+0.1428549246E11*t5-9523661640.0*t11);
    return(0.1E-3*t18);
  }
}

double evalRho(double x)
{
  {
    return(0.4958034541*x+0.1065766549);
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t16;
  double t17;
  double t18;
  double t20;
  double t22;
  double t29;
  double t3;
  double t4;
  double t7;
  {
    t3 = phi*phi;
    t4 = t3*t3;
    t7 = t3*phi;
    t16 = 6.0*t4*phi;
    t17 = 15.0*t4;
    t18 = 10.0*t7;
    t20 = rho*rho;
    t22 = t20*t20;
    t29 = log(rho);
    return(2.0/delta*A*(t4+(2.0*alpha-2.0)*t7+(-3.0*alpha+1.0)*t3+alpha)+(t16-
t17+t18)*(0.8575855211E-1*t22*t20*rho-0.2868132559E-1*rho+0.1480890852E-1)+(1.0
-t16+t17-t18)*(0.15001E-1+0.1407531515*rho*t29+0.1743777827*rho));
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t17;
  double t21;
  double t22;
  double t24;
  double t3;
  double t31;
  double t4;
  {
    t3 = phi*phi;
    t4 = t3*phi;
    t17 = t3*t3;
    t21 = 30.0*t17-60.0*t4+30.0*t3;
    t22 = rho*rho;
    t24 = t22*t22;
    t31 = log(rho);
    return(2.0/delta*A*(4.0*t4+3.0*(2.0*alpha-2.0)*t3+2.0*(-3.0*alpha+1.0)*phi)
+t21*(0.8575855211E-1*t24*t22*rho-0.2868132559E-1*rho+0.1480890852E-1)-t21*(
0.15001E-1+0.1407531515*rho*t31+0.1743777827*rho));
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  double t17;
  double t18;
  double t20;
  double t27;
  double t3;
  {
    t3 = phi*phi;
    t17 = 120.0*t3*phi-180.0*t3+60.0*phi;
    t18 = rho*rho;
    t20 = t18*t18;
    t27 = log(rho);
    return(2.0/delta*A*(12.0*t3+6.0*(2.0*alpha-2.0)*phi-6.0*alpha+2.0)+t17*(
0.8575855211E-1*t20*t18*rho-0.2868132559E-1*rho+0.1480890852E-1)-t17*(
0.15001E-1+0.1407531515*rho*t27+0.1743777827*rho));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t13;
  double t17;
  double t2;
  double t3;
  double t4;
  double t5;
  double t6;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = t2*phi;
    t4 = rho*rho;
    t5 = t4*t4;
    t6 = t5*t4;
    t13 = t1*phi;
    t17 = log(rho);
    return(0.3601859189E1*t3*t6-0.2062873559E1*t3-0.9004647972E1*t2*t6+
0.5157183897E1*t2+0.6003098648E1*t13*t6-0.3438122598E1*t13+0.3151309342+
0.1407531515*t17-0.844518909*t3*t17+0.2111297272E1*t2*t17-0.1407531515E1*t13*
t17);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  double t1;
  double t11;
  double t2;
  double t3;
  double t4;
  double t5;
  double t6;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = t2*phi;
    t4 = rho*rho;
    t5 = t4*t4;
    t6 = t5*t4;
    t11 = t1*phi;
    return(0.5E-9*(0.4322231027E11*t3*t6-0.1080557757E12*t2*t6+0.7203718378E11*
t11*t6+281506303.0-1689037818.0*t3+4222594545.0*t2-2815063030.0*t11)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t16;
  double t2;
  double t3;
  double t4;
  double t5;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = rho*rho;
    t4 = t3*t3;
    t5 = t4*t3;
    t9 = t1*phi;
    t16 = log(rho);
    return(0.1800929594E2*t2*t5-0.1031436779E2*t2-0.3601859189E2*t9*t5+
0.2062873559E2*t9+0.1800929594E2*t1*t5-0.1031436779E2*t1-0.4222594545E1*t16*t2+
0.844518909E1*t16*t9-0.4222594545E1*t1*t16);
  }
}

#include <math.h>
double pressure(double rho,double phi)
{
  double t1;
  double t11;
  double t2;
  double t22;
  double t4;
  double t5;
  double t6;
  double t7;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t6 = t1*phi;
    t7 = 10.0*t6;
    t9 = rho*rho;
    t11 = t9*t9;
    t22 = log(rho);
    return((t4-t5+t7)*(-0.8575855211E-1*t11*t9*rho+0.2868132559E-1*rho
-0.1480890852E-1+rho*(0.6003098648*t11*t9-0.2868132559E-1))+(1.0-t4+t5-t7)*(
-0.15001E-1-0.1407531515*rho*t22-0.1743777827*rho+rho*(0.3151309342+
0.1407531515*t22))-2.0/delta*A*(t2+(2.0*alpha-2.0)*t6+(-3.0*alpha+1.0)*t1+alpha
));
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t1;
  double t11;
  double t18;
  double t2;
  double t3;
  double t4;
  double t5;
  double t6;
  {
    t1 = rho*rho;
    t2 = t1*t1;
    t3 = t2*t1;
    t4 = phi*phi;
    t5 = t4*t4;
    t6 = t5*phi;
    t11 = t4*phi;
    t18 = sqrt(0.2161115513E12*t3*t6-0.5402788784E12*t3*t5+0.3601859189E12*t3*
t11+1407531515.0-8445189090.0*t6+0.2111297272E11*t5-0.1407531515E11*t11);
    return(0.1E-4*t18);
  }
}

double evalRho(double x)
{
  {
    return(0.4958034541*x+0.1065766549);
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t16;
  double t17;
  double t18;
  double t20;
  double t22;
  double t29;
  double t3;
  double t4;
  double t7;
  {
    t3 = phi*phi;
    t4 = t3*t3;
    t7 = t3*phi;
    t16 = 6.0*t4*phi;
    t17 = 15.0*t4;
    t18 = 10.0*t7;
    t20 = rho*rho;
    t22 = t20*t20;
    t29 = log(rho);
    return(2.0/delta*A*(t4+(2.0*alpha-2.0)*t7+(-3.0*alpha+1.0)*t3+alpha)+(t16-
t17+t18)*(0.8575855211E-1*t22*t20*rho-0.2868132559E-1*rho+0.1480890852E-1)+(1.0
-t16+t17-t18)*(0.1500001E-1+0.1407438624*rho*t29+0.1743662745*rho));
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t17;
  double t21;
  double t22;
  double t24;
  double t3;
  double t31;
  double t4;
  {
    t3 = phi*phi;
    t4 = t3*phi;
    t17 = t3*t3;
    t21 = 30.0*t17-60.0*t4+30.0*t3;
    t22 = rho*rho;
    t24 = t22*t22;
    t31 = log(rho);
    return(2.0/delta*A*(4.0*t4+3.0*(2.0*alpha-2.0)*t3+2.0*(-3.0*alpha+1.0)*phi)
+t21*(0.8575855211E-1*t24*t22*rho-0.2868132559E-1*rho+0.1480890852E-1)-t21*(
0.1500001E-1+0.1407438624*rho*t31+0.1743662745*rho));
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  double t17;
  double t18;
  double t20;
  double t27;
  double t3;
  {
    t3 = phi*phi;
    t17 = 120.0*t3*phi-180.0*t3+60.0*phi;
    t18 = rho*rho;
    t20 = t18*t18;
    t27 = log(rho);
    return(2.0/delta*A*(12.0*t3+6.0*(2.0*alpha-2.0)*phi-6.0*alpha+2.0)+t17*(
0.8575855211E-1*t20*t18*rho-0.2868132559E-1*rho+0.1480890852E-1)-t17*(
0.1500001E-1+0.1407438624*rho*t27+0.1743662745*rho));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t13;
  double t17;
  double t2;
  double t3;
  double t4;
  double t5;
  double t6;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = t2*phi;
    t4 = rho*rho;
    t5 = t4*t4;
    t6 = t5*t4;
    t13 = t1*phi;
    t17 = log(rho);
    return(0.3601859189E1*t3*t6-0.2062748775E1*t3-0.9004647972E1*t2*t6+
0.5156871937E1*t2+0.6003098648E1*t13*t6-0.3437914625E1*t13+0.3151101369+
0.1407438624*t17-0.8444631744*t3*t17+0.2111157936E1*t2*t17-0.1407438624E1*t13*
t17);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  double t1;
  double t11;
  double t2;
  double t3;
  double t4;
  double t5;
  double t6;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = t2*phi;
    t4 = rho*rho;
    t5 = t4*t4;
    t6 = t5*t4;
    t11 = t1*phi;
    return(0.2E-9*(0.1080557757E12*t3*t6-0.2701394392E12*t2*t6+0.1800929594E12*
t11*t6+703719312.0-4222315872.0*t3+0.1055578968E11*t2-7037193120.0*t11)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t16;
  double t2;
  double t3;
  double t4;
  double t5;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = rho*rho;
    t4 = t3*t3;
    t5 = t4*t3;
    t9 = t1*phi;
    t16 = log(rho);
    return(0.1800929594E2*t2*t5-0.1031374387E2*t2-0.3601859189E2*t9*t5+
0.2062748775E2*t9+0.1800929594E2*t1*t5-0.1031374387E2*t1-0.4222315872E1*t2*t16+
0.8444631744E1*t9*t16-0.4222315872E1*t1*t16);
  }
}

#include <math.h>
double pressure(double rho,double phi)
{
  double t1;
  double t11;
  double t2;
  double t22;
  double t4;
  double t5;
  double t6;
  double t7;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t6 = t1*phi;
    t7 = 10.0*t6;
    t9 = rho*rho;
    t11 = t9*t9;
    t22 = log(rho);
    return((t4-t5+t7)*(-0.1480890852E-1-0.8575855211E-1*t11*t9*rho+
0.2868132559E-1*rho+rho*(-0.2868132559E-1+0.6003098648*t11*t9))+(1.0-t4+t5-t7)*
(-0.1500001E-1-0.1407438624*rho*t22-0.1743662745*rho+rho*(0.3151101369+
0.1407438624*t22))-2.0/delta*A*(t2+(2.0*alpha-2.0)*t6+(-3.0*alpha+1.0)*t1+alpha
));
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t1;
  double t11;
  double t18;
  double t2;
  double t3;
  double t4;
  double t5;
  double t6;
  {
    t1 = rho*rho;
    t2 = t1*t1;
    t3 = t2*t1;
    t4 = phi*phi;
    t5 = t4*t4;
    t6 = t5*phi;
    t11 = t4*phi;
    t18 = sqrt(0.1350697196E11*t3*t6-0.337674299E11*t3*t5+0.2251161993E11*t3*
t11+87964914.0-527789484.0*t6+1319473710.0*t5-879649140.0*t11);
    return(0.4E-4*t18);
  }
}

double evalRho(double x)
{
  {
    return(0.4958034541*x+0.1065766549);
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t16;
  double t17;
  double t18;
  double t20;
  double t22;
  double t29;
  double t3;
  double t4;
  double t7;
  {
    t3 = phi*phi;
    t4 = t3*t3;
    t7 = t3*phi;
    t16 = 6.0*t4*phi;
    t17 = 15.0*t4;
    t18 = 10.0*t7;
    t20 = rho*rho;
    t22 = t20*t20;
    t29 = log(rho);
    return(2.0/delta*A*(t4+(2.0*alpha-2.0)*t7+(-3.0*alpha+1.0)*t3+alpha)+(t16-
t17+t18)*(0.8575855211E-1*t22*t20*rho-0.2868132559E-1*rho+0.1480890852E-1)+(1.0
-t16+t17-t18)*(0.1E-7+0.9429832461E-7*rho*t29+0.1168253257E-6*rho));
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t17;
  double t21;
  double t22;
  double t24;
  double t3;
  double t31;
  double t4;
  {
    t3 = phi*phi;
    t4 = t3*phi;
    t17 = t3*t3;
    t21 = 30.0*t17-60.0*t4+30.0*t3;
    t22 = rho*rho;
    t24 = t22*t22;
    t31 = log(rho);
    return(2.0/delta*A*(4.0*t4+3.0*(2.0*alpha-2.0)*t3+2.0*(-3.0*alpha+1.0)*phi)
+t21*(0.8575855211E-1*t24*t22*rho-0.2868132559E-1*rho+0.1480890852E-1)-t21*(
0.1E-7+0.9429832461E-7*rho*t31+0.1168253257E-6*rho));
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  double t17;
  double t18;
  double t20;
  double t27;
  double t3;
  {
    t3 = phi*phi;
    t17 = 120.0*t3*phi-180.0*t3+60.0*phi;
    t18 = rho*rho;
    t20 = t18*t18;
    t27 = log(rho);
    return(2.0/delta*A*(12.0*t3+6.0*(2.0*alpha-2.0)*phi-6.0*alpha+2.0)+t17*(
0.8575855211E-1*t20*t18*rho-0.2868132559E-1*rho+0.1480890852E-1)-t17*(0.1E-7+
0.9429832461E-7*rho*t27+0.1168253257E-6*rho));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t13;
  double t17;
  double t2;
  double t3;
  double t4;
  double t5;
  double t6;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = t2*phi;
    t4 = rho*rho;
    t5 = t4*t4;
    t6 = t5*t4;
    t13 = t1*phi;
    t17 = log(rho);
    return(0.3601859189E1*t3*t6-0.1720892203*t3-0.9004647972E1*t2*t6+
0.4302230507*t2+0.6003098648E1*t13*t6-0.2868153671*t13+0.2111236503E-6+
0.9429832461E-7*t17-0.5657899477E-6*t3*t17+0.1414474869E-5*t2*t17
-0.9429832461E-6*t13*t17);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  double t1;
  double t11;
  double t2;
  double t3;
  double t4;
  double t5;
  double t6;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = t2*phi;
    t4 = rho*rho;
    t5 = t4*t4;
    t6 = t5*t4;
    t11 = t1*phi;
    return(0.1E-16*(0.2161115513E19*t3*t6-0.5402788784E19*t2*t6+0.3601859189E19
*t11*t6+9429832461.0-0.5657899477E11*t3+0.1414474869E12*t2-0.9429832461E11*t11)
/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t16;
  double t2;
  double t3;
  double t4;
  double t5;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = rho*rho;
    t4 = t3*t3;
    t5 = t4*t3;
    t9 = t1*phi;
    t16 = log(rho);
    return(0.1800929594E2*t2*t5-0.8604461014*t2-0.3601859189E2*t9*t5+
0.1720892203E1*t9+0.1800929594E2*t1*t5-0.8604461014*t1-0.2828949738E-5*t2*t16+
0.5657899477E-5*t9*t16-0.2828949738E-5*t1*t16);
  }
}

#include <math.h>
double pressure(double rho,double phi)
{
  double t1;
  double t11;
  double t2;
  double t22;
  double t4;
  double t5;
  double t6;
  double t7;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t6 = t1*phi;
    t7 = 10.0*t6;
    t9 = rho*rho;
    t11 = t9*t9;
    t22 = log(rho);
    return((t4-t5+t7)*(-0.8575855211E-1*t11*t9*rho+0.2868132559E-1*rho
-0.1480890852E-1+rho*(0.6003098648*t11*t9-0.2868132559E-1))+(1.0-t4+t5-t7)*(
-0.1E-7-0.9429832461E-7*rho*t22-0.1168253257E-6*rho+rho*(0.2111236503E-6+
0.9429832461E-7*t22))-2.0/delta*A*(t2+(2.0*alpha-2.0)*t6+(-3.0*alpha+1.0)*t1+
alpha));
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t1;
  double t11;
  double t18;
  double t2;
  double t3;
  double t4;
  double t5;
  double t6;
  {
    t1 = rho*rho;
    t2 = t1*t1;
    t3 = t2*t1;
    t4 = phi*phi;
    t5 = t4*t4;
    t6 = t5*phi;
    t11 = t4*phi;
    t18 = sqrt(0.2161115513E18*t3*t6-0.5402788784E18*t3*t5+0.3601859189E18*t3*
t11+942983246.0-5657899476.0*t6+0.1414474869E11*t5-9429832460.0*t11);
    return(0.1E-7*t18);
  }
}

double evalRho(double x)
{
  {
    return(0.4958034541*x+0.1065766549);
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t16;
  double t17;
  double t18;
  double t20;
  double t22;
  double t29;
  double t3;
  double t4;
  double t7;
  {
    t3 = phi*phi;
    t4 = t3*t3;
    t7 = t3*phi;
    t16 = 6.0*t4*phi;
    t17 = 15.0*t4;
    t18 = 10.0*t7;
    t20 = rho*rho;
    t22 = t20*t20;
    t29 = log(rho);
    return(2.0/delta*A*(t4+(2.0*alpha-2.0)*t7+(-3.0*alpha+1.0)*t3+alpha)+(t16-
t17+t18)*(0.8575855211E-1*t22*t20*rho-0.2868132559E-1*rho+0.1480890852E-1)+(1.0
-t16+t17-t18)*(0.1E-1+0.9382917919E-1*rho*t29+0.1162441057*rho));
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t17;
  double t21;
  double t22;
  double t24;
  double t3;
  double t31;
  double t4;
  {
    t3 = phi*phi;
    t4 = t3*phi;
    t17 = t3*t3;
    t21 = 30.0*t17-60.0*t4+30.0*t3;
    t22 = rho*rho;
    t24 = t22*t22;
    t31 = log(rho);
    return(2.0/delta*A*(4.0*t4+3.0*(2.0*alpha-2.0)*t3+2.0*(-3.0*alpha+1.0)*phi)
+t21*(0.8575855211E-1*t24*t22*rho-0.2868132559E-1*rho+0.1480890852E-1)-t21*(
0.1E-1+0.9382917919E-1*rho*t31+0.1162441057*rho));
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  double t17;
  double t18;
  double t20;
  double t27;
  double t3;
  {
    t3 = phi*phi;
    t17 = 120.0*t3*phi-180.0*t3+60.0*phi;
    t18 = rho*rho;
    t20 = t18*t18;
    t27 = log(rho);
    return(2.0/delta*A*(12.0*t3+6.0*(2.0*alpha-2.0)*phi-6.0*alpha+2.0)+t17*(
0.8575855211E-1*t20*t18*rho-0.2868132559E-1*rho+0.1480890852E-1)-t17*(0.1E-1+
0.9382917919E-1*rho*t27+0.1162441057*rho));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t13;
  double t17;
  double t2;
  double t3;
  double t4;
  double t5;
  double t6;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = t2*phi;
    t4 = rho*rho;
    t5 = t4*t4;
    t6 = t5*t4;
    t13 = t1*phi;
    t17 = log(rho);
    return(0.3601859189E1*t3*t6-0.1432527663E1*t3-0.9004647972E1*t2*t6+
0.3581319157E1*t2+0.6003098648E1*t13*t6-0.2387546105E1*t13+0.2100732849+
0.9382917919E-1*t17-0.5629750751*t3*t17+0.1407437688E1*t2*t17-0.9382917919*t13*
t17);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  double t1;
  double t11;
  double t2;
  double t3;
  double t4;
  double t5;
  double t6;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = t2*phi;
    t4 = rho*rho;
    t5 = t4*t4;
    t6 = t5*t4;
    t11 = t1*phi;
    return(0.1E-10*(0.2161115513E13*t3*t6-0.5402788784E13*t2*t6+0.3601859189E13
*t11*t6+9382917919.0-0.5629750751E11*t3+0.1407437688E12*t2-0.9382917919E11*t11)
/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t16;
  double t2;
  double t3;
  double t4;
  double t5;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = rho*rho;
    t4 = t3*t3;
    t5 = t4*t3;
    t9 = t1*phi;
    t16 = log(rho);
    return(0.1800929594E2*t2*t5-0.7162638315E1*t2-0.3601859189E2*t9*t5+
0.1432527663E2*t9+0.1800929594E2*t1*t5-0.7162638315E1*t1-0.2814875376E1*t2*t16+
0.5629750751E1*t9*t16-0.2814875376E1*t1*t16);
  }
}

#include <math.h>
double pressure(double rho,double phi)
{
  double t1;
  double t11;
  double t2;
  double t22;
  double t4;
  double t5;
  double t6;
  double t7;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t6 = t1*phi;
    t7 = 10.0*t6;
    t9 = rho*rho;
    t11 = t9*t9;
    t22 = log(rho);
    return((t4-t5+t7)*(-0.8575855211E-1*t11*t9*rho+0.2868132559E-1*rho
-0.1480890852E-1+rho*(0.6003098648*t11*t9-0.2868132559E-1))+(1.0-t4+t5-t7)*(
-0.1E-1-0.9382917919E-1*rho*t22-0.1162441057*rho+rho*(0.2100732849+
0.9382917919E-1*t22))-2.0/delta*A*(t2+(2.0*alpha-2.0)*t6+(-3.0*alpha+1.0)*t1+
alpha));
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t1;
  double t11;
  double t18;
  double t2;
  double t3;
  double t4;
  double t5;
  double t6;
  {
    t1 = rho*rho;
    t2 = t1*t1;
    t3 = t2*t1;
    t4 = phi*phi;
    t5 = t4*t4;
    t6 = t5*phi;
    t11 = t4*phi;
    t18 = sqrt(0.1350697196E11*t3*t6-0.337674299E11*t3*t5+0.2251161993E11*t3*
t11+58643237.0-351859422.0*t6+879648555.0*t5-586432370.0*t11);
    return(0.4E-4*t18);
  }
}

double evalRho(double x)
{
  {
    return(0.4958034541*x+0.1065766549);
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t16;
  double t17;
  double t18;
  double t20;
  double t22;
  double t3;
  double t30;
  double t31;
  double t33;
  double t4;
  double t7;
  {
    t3 = phi*phi;
    t4 = t3*t3;
    t7 = t3*phi;
    t16 = 6.0*t4*phi;
    t17 = 15.0*t4;
    t18 = 10.0*t7;
    t20 = rho*rho;
    t22 = t20*t20;
    t30 = dv(rho-0.1065766549);
    t31 = t30*t30;
    t33 = log(rho);
    return(2.0/delta*A*(t4+(2.0*alpha-2.0)*t7+(-3.0*alpha+1.0)*t3+alpha)+(t16-
t17+t18)*(0.8575855211E-1*t22*t20*rho-0.2868132559E-1*rho+0.1480890852E-1)+(1.0
-t16+t17-t18)*(t31+av*rho*t33+(bv-av)*rho+0.15-0.1135858337E-1*dv));
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t17;
  double t21;
  double t22;
  double t24;
  double t3;
  double t32;
  double t33;
  double t35;
  double t4;
  {
    t3 = phi*phi;
    t4 = t3*phi;
    t17 = t3*t3;
    t21 = 30.0*t17-60.0*t4+30.0*t3;
    t22 = rho*rho;
    t24 = t22*t22;
    t32 = dv(rho-0.1065766549);
    t33 = t32*t32;
    t35 = log(rho);
    return(2.0/delta*A*(4.0*t4+3.0*(2.0*alpha-2.0)*t3+2.0*(-3.0*alpha+1.0)*phi)
+t21*(0.8575855211E-1*t24*t22*rho-0.2868132559E-1*rho+0.1480890852E-1)-t21*(t33
+av*rho*t35+(bv-av)*rho+0.15-0.1135858337E-1*dv));
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  double t17;
  double t18;
  double t20;
  double t28;
  double t29;
  double t3;
  double t31;
  {
    t3 = phi*phi;
    t17 = 120.0*t3*phi-180.0*t3+60.0*phi;
    t18 = rho*rho;
    t20 = t18*t18;
    t28 = dv(rho-0.1065766549);
    t29 = t28*t28;
    t31 = log(rho);
    return(2.0/delta*A*(12.0*t3+6.0*(2.0*alpha-2.0)*phi-6.0*alpha+2.0)+t17*(
0.8575855211E-1*t20*t18*rho-0.2868132559E-1*rho+0.1480890852E-1)-t17*(t29+av*
rho*t31+(bv-av)*rho+0.15-0.1135858337E-1*dv));
  }
}






double evalRho(double x)
{
  {
    return(0.4958034541*x+0.1065766549);
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t16;
  double t17;
  double t18;
  double t20;
  double t22;
  double t29;
  double t3;
  double t4;
  double t7;
  {
    t3 = phi*phi;
    t4 = t3*t3;
    t7 = t3*phi;
    t16 = 6.0*t4*phi;
    t17 = 15.0*t4;
    t18 = 10.0*t7;
    t20 = rho*rho;
    t22 = t20*t20;
    t29 = log(rho);
    return(2.0/delta*A*(t4+(2.0*alpha-2.0)*t7+(-3.0*alpha+1.0)*t3+alpha)+(t16-
t17+t18)*(0.8575855211E-1*t22*t20*rho-0.2868132559E-1*rho+0.1480890852E-1)+(1.0
-t16+t17-t18)*(0.2386414166E-1+0.2239152818*rho*t29+0.2774065798*rho));
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t17;
  double t21;
  double t22;
  double t24;
  double t3;
  double t31;
  double t4;
  {
    t3 = phi*phi;
    t4 = t3*phi;
    t17 = t3*t3;
    t21 = 30.0*t17-60.0*t4+30.0*t3;
    t22 = rho*rho;
    t24 = t22*t22;
    t31 = log(rho);
    return(2.0/delta*A*(4.0*t4+3.0*(2.0*alpha-2.0)*t3+2.0*(-3.0*alpha+1.0)*phi)
+t21*(0.8575855211E-1*t24*t22*rho-0.2868132559E-1*rho+0.1480890852E-1)-t21*(
0.2386414166E-1+0.2239152818*rho*t31+0.2774065798*rho));
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  double t17;
  double t18;
  double t20;
  double t27;
  double t3;
  {
    t3 = phi*phi;
    t17 = 120.0*t3*phi-180.0*t3+60.0*phi;
    t18 = rho*rho;
    t20 = t18*t18;
    t27 = log(rho);
    return(2.0/delta*A*(12.0*t3+6.0*(2.0*alpha-2.0)*phi-6.0*alpha+2.0)+t17*(
0.8575855211E-1*t20*t18*rho-0.2868132559E-1*rho+0.1480890852E-1)-t17*(
0.2386414166E-1+0.2239152818*rho*t27+0.2774065798*rho));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t13;
  double t17;
  double t2;
  double t3;
  double t4;
  double t5;
  double t6;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = t2*phi;
    t4 = rho*rho;
    t5 = t4*t4;
    t6 = t5*t4;
    t13 = t1*phi;
    t17 = log(rho);
    return(0.3601859189E1*t3*t6-0.3180019123E1*t3-0.9004647972E1*t2*t6+
0.7950047808E1*t2+0.6003098648E1*t13*t6-0.5300031872E1*t13+0.5013218616+
0.2239152818*t17-0.1343491691E1*t3*t17+0.3358729227E1*t2*t17-0.2239152818E1*t13
*t17);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  double t1;
  double t11;
  double t2;
  double t3;
  double t4;
  double t5;
  double t6;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = t2*phi;
    t4 = rho*rho;
    t5 = t4*t4;
    t6 = t5*t4;
    t11 = t1*phi;
    return(0.2E-9*(0.1080557757E12*t3*t6-0.2701394392E12*t2*t6+0.1800929594E12*
t11*t6+1119576409.0-6717458454.0*t3+0.1679364614E11*t2-0.1119576409E11*t11)/rho
);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t16;
  double t2;
  double t3;
  double t4;
  double t5;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = rho*rho;
    t4 = t3*t3;
    t5 = t4*t3;
    t9 = t1*phi;
    t16 = log(rho);
    return(0.1800929594E2*t2*t5-0.1590009562E2*t2-0.3601859189E2*t9*t5+
0.3180019123E2*t9+0.1800929594E2*t1*t5-0.1590009562E2*t1-0.6717458454E1*t2*t16+
0.1343491691E2*t9*t16-0.6717458454E1*t1*t16);
  }
}

#include <math.h>
double pressure(double rho,double phi)
{
  double t1;
  double t11;
  double t2;
  double t22;
  double t4;
  double t5;
  double t6;
  double t7;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t6 = t1*phi;
    t7 = 10.0*t6;
    t9 = rho*rho;
    t11 = t9*t9;
    t22 = log(rho);
    return((t4-t5+t7)*(-0.8575855211E-1*t11*t9*rho+0.2868132559E-1*rho
-0.1480890852E-1+rho*(0.6003098648*t11*t9-0.2868132559E-1))+(1.0-t4+t5-t7)*(
-0.2386414166E-1-0.2239152818*rho*t22-0.2774065798*rho+rho*(0.5013218616+
0.2239152818*t22))-2.0/delta*A*(t2+(2.0*alpha-2.0)*t6+(-3.0*alpha+1.0)*t1+alpha
));
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t1;
  double t11;
  double t18;
  double t2;
  double t3;
  double t4;
  double t5;
  double t6;
  {
    t1 = rho*rho;
    t2 = t1*t1;
    t3 = t2*t1;
    t4 = phi*phi;
    t5 = t4*t4;
    t6 = t5*phi;
    t11 = t4*phi;
    t18 = sqrt(0.2161115513E12*t3*t6-0.5402788784E12*t3*t5+0.3601859189E12*t3*
t11+2239152818.0-0.1343491691E11*t6+0.3358729227E11*t5-0.2239152818E11*t11);
    return(0.1E-4*t18);
  }
}

double evalRho(double x)
{
  {
    return(0.4958034541*x+0.1065766549);
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t16;
  double t17;
  double t18;
  double t20;
  double t22;
  double t29;
  double t3;
  double t4;
  double t7;
  {
    t3 = phi*phi;
    t4 = t3*t3;
    t7 = t3*phi;
    t16 = 6.0*t4*phi;
    t17 = 15.0*t4;
    t18 = 10.0*t7;
    t20 = rho*rho;
    t22 = t20*t20;
    t29 = log(rho);
    return(2.0/delta*A*(t4+(2.0*alpha-2.0)*t7+(-3.0*alpha+1.0)*t3+alpha)+(t16-
t17+t18)*(0.8575855211E-1*t22*t20*rho-0.2868132559E-1*rho+0.1480890852E-1)+(1.0
-t16+t17-t18)*(0.1498641417E-1+0.1406162938*rho*t29+0.1742082311*rho));
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t17;
  double t21;
  double t22;
  double t24;
  double t3;
  double t31;
  double t4;
  {
    t3 = phi*phi;
    t4 = t3*phi;
    t17 = t3*t3;
    t21 = 30.0*t17-60.0*t4+30.0*t3;
    t22 = rho*rho;
    t24 = t22*t22;
    t31 = log(rho);
    return(2.0/delta*A*(4.0*t4+3.0*(2.0*alpha-2.0)*t3+2.0*(-3.0*alpha+1.0)*phi)
+t21*(0.8575855211E-1*t24*t22*rho-0.2868132559E-1*rho+0.1480890852E-1)-t21*(
0.1498641417E-1+0.1406162938*rho*t31+0.1742082311*rho));
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  double t17;
  double t18;
  double t20;
  double t27;
  double t3;
  {
    t3 = phi*phi;
    t17 = 120.0*t3*phi-180.0*t3+60.0*phi;
    t18 = rho*rho;
    t20 = t18*t18;
    t27 = log(rho);
    return(2.0/delta*A*(12.0*t3+6.0*(2.0*alpha-2.0)*phi-6.0*alpha+2.0)+t17*(
0.8575855211E-1*t20*t18*rho-0.2868132559E-1*rho+0.1480890852E-1)-t17*(
0.1498641417E-1+0.1406162938*rho*t27+0.1742082311*rho));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t13;
  double t17;
  double t2;
  double t3;
  double t4;
  double t5;
  double t6;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = t2*phi;
    t4 = rho*rho;
    t5 = t4*t4;
    t6 = t5*t4;
    t13 = t1*phi;
    t17 = log(rho);
    return(0.3601859189E1*t3*t6-0.2061035103E1*t3-0.9004647972E1*t2*t6+
0.5152587757E1*t2+0.6003098648E1*t13*t6-0.3435058505E1*t13+0.3148245249+
0.1406162938*t17-0.8436977628*t3*t17+0.2109244407E1*t2*t17-0.1406162938E1*t13*
t17);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  double t1;
  double t11;
  double t2;
  double t3;
  double t4;
  double t5;
  double t6;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = t2*phi;
    t4 = rho*rho;
    t5 = t4*t4;
    t6 = t5*t4;
    t11 = t1*phi;
    return(0.2E-9*(0.1080557757E12*t3*t6-0.2701394392E12*t2*t6+0.1800929594E12*
t11*t6+703081469.0-4218488814.0*t3+0.1054622204E11*t2-7030814690.0*t11)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t16;
  double t2;
  double t3;
  double t4;
  double t5;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = rho*rho;
    t4 = t3*t3;
    t5 = t4*t3;
    t9 = t1*phi;
    t16 = log(rho);
    return(0.1800929594E2*t2*t5-0.1030517551E2*t2-0.3601859189E2*t9*t5+
0.2061035103E2*t9+0.1800929594E2*t1*t5-0.1030517551E2*t1-0.4218488814E1*t2*t16+
0.8436977628E1*t9*t16-0.4218488814E1*t1*t16);
  }
}

#include <math.h>
double pressure(double rho,double phi)
{
  double t1;
  double t11;
  double t2;
  double t22;
  double t4;
  double t5;
  double t6;
  double t7;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t6 = t1*phi;
    t7 = 10.0*t6;
    t9 = rho*rho;
    t11 = t9*t9;
    t22 = log(rho);
    return((t4-t5+t7)*(-0.8575855211E-1*t11*t9*rho+0.2868132559E-1*rho
-0.1480890852E-1+rho*(0.6003098648*t11*t9-0.2868132559E-1))+(1.0-t4+t5-t7)*(
-0.1498641417E-1-0.1406162938*rho*t22-0.1742082311*rho+rho*(0.3148245249+
0.1406162938*t22))-2.0/delta*A*(t2+(2.0*alpha-2.0)*t6+(-3.0*alpha+1.0)*t1+alpha
));
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t1;
  double t11;
  double t18;
  double t2;
  double t3;
  double t4;
  double t5;
  double t6;
  {
    t1 = rho*rho;
    t2 = t1*t1;
    t3 = t2*t1;
    t4 = phi*phi;
    t5 = t4*t4;
    t6 = t5*phi;
    t11 = t4*phi;
    t18 = sqrt(0.2161115513E12*t3*t6-0.5402788784E12*t3*t5+0.3601859189E12*t3*
t11+1406162938.0-8436977628.0*t6+0.2109244407E11*t5-0.1406162938E11*t11);
    return(0.1E-4*t18);
  }
}

double evalRho(double x)
{
  {
    return(0.4958034541*x+0.1065766549);
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t16;
  double t17;
  double t18;
  double t20;
  double t22;
  double t29;
  double t3;
  double t4;
  double t7;
  {
    t3 = phi*phi;
    t4 = t3*t3;
    t7 = t3*phi;
    t16 = 6.0*t4*phi;
    t17 = 15.0*t4;
    t18 = 10.0*t7;
    t20 = rho*rho;
    t22 = t20*t20;
    t29 = log(rho);
    return(2.0/delta*A*(t4+(2.0*alpha-2.0)*t7+(-3.0*alpha+1.0)*t3+alpha)+(t16-
t17+t18)*(0.8575855211E-1*t22*t20*rho-0.2868132559E-1*rho+0.1480890852E-1)+(1.0
-t16+t17-t18)*(0.1499887414E-1+0.1407332047*rho*t29+0.1743530709*rho));
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t17;
  double t21;
  double t22;
  double t24;
  double t3;
  double t31;
  double t4;
  {
    t3 = phi*phi;
    t4 = t3*phi;
    t17 = t3*t3;
    t21 = 30.0*t17-60.0*t4+30.0*t3;
    t22 = rho*rho;
    t24 = t22*t22;
    t31 = log(rho);
    return(2.0/delta*A*(4.0*t4+3.0*(2.0*alpha-2.0)*t3+2.0*(-3.0*alpha+1.0)*phi)
+t21*(0.8575855211E-1*t24*t22*rho-0.2868132559E-1*rho+0.1480890852E-1)-t21*(
0.1499887414E-1+0.1407332047*rho*t31+0.1743530709*rho));
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  double t17;
  double t18;
  double t20;
  double t27;
  double t3;
  {
    t3 = phi*phi;
    t17 = 120.0*t3*phi-180.0*t3+60.0*phi;
    t18 = rho*rho;
    t20 = t18*t18;
    t27 = log(rho);
    return(2.0/delta*A*(12.0*t3+6.0*(2.0*alpha-2.0)*phi-6.0*alpha+2.0)+t17*(
0.8575855211E-1*t20*t18*rho-0.2868132559E-1*rho+0.1480890852E-1)-t17*(
0.1499887414E-1+0.1407332047*rho*t27+0.1743530709*rho));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t13;
  double t17;
  double t2;
  double t3;
  double t4;
  double t5;
  double t6;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = t2*phi;
    t4 = rho*rho;
    t5 = t4*t4;
    t6 = t5*t4;
    t13 = t1*phi;
    t17 = log(rho);
    return(0.3601859189E1*t3*t6-0.2062605607E1*t3-0.9004647972E1*t2*t6+
0.5156514018E1*t2+0.6003098648E1*t13*t6-0.3437676012E1*t13+0.3150862756+
0.1407332047*t17-0.8443992282*t3*t17+0.211099807E1*t2*t17-0.1407332047E1*t13*
t17);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  double t1;
  double t11;
  double t2;
  double t3;
  double t4;
  double t5;
  double t6;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = t2*phi;
    t4 = rho*rho;
    t5 = t4*t4;
    t6 = t5*t4;
    t11 = t1*phi;
    return(0.1E-9*(0.2161115513E12*t3*t6-0.5402788784E12*t2*t6+0.3601859189E12*
t11*t6+1407332047.0-8443992282.0*t3+0.211099807E11*t2-0.1407332047E11*t11)/rho
);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t16;
  double t2;
  double t3;
  double t4;
  double t5;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = rho*rho;
    t4 = t3*t3;
    t5 = t4*t3;
    t9 = t1*phi;
    t16 = log(rho);
    return(0.1800929594E2*t2*t5-0.1031302804E2*t2-0.3601859189E2*t9*t5+
0.2062605607E2*t9+0.1800929594E2*t1*t5-0.1031302804E2*t1-0.4221996141E1*t2*t16+
0.8443992282E1*t9*t16-0.4221996141E1*t1*t16);
  }
}

#include <math.h>
double pressure(double rho,double phi)
{
  double t1;
  double t11;
  double t2;
  double t22;
  double t4;
  double t5;
  double t6;
  double t7;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t6 = t1*phi;
    t7 = 10.0*t6;
    t9 = rho*rho;
    t11 = t9*t9;
    t22 = log(rho);
    return((t4-t5+t7)*(-0.8575855211E-1*t11*t9*rho+0.2868132559E-1*rho
-0.1480890852E-1+rho*(0.6003098648*t11*t9-0.2868132559E-1))+(1.0-t4+t5-t7)*(
-0.1499887414E-1-0.1407332047*rho*t22-0.1743530709*rho+rho*(0.3150862756+
0.1407332047*t22))-2.0/delta*A*(t2+(2.0*alpha-2.0)*t6+(-3.0*alpha+1.0)*t1+alpha
));
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t1;
  double t11;
  double t18;
  double t2;
  double t3;
  double t4;
  double t5;
  double t6;
  {
    t1 = rho*rho;
    t2 = t1*t1;
    t3 = t2*t1;
    t4 = phi*phi;
    t5 = t4*t4;
    t6 = t5*phi;
    t11 = t4*phi;
    t18 = sqrt(0.2161115513E12*t3*t6-0.5402788784E12*t3*t5+0.3601859189E12*t3*
t11+1407332047.0-8443992282.0*t6+0.211099807E11*t5-0.1407332047E11*t11);
    return(0.1E-4*t18);
  }
}

double evalRho(double x)
{
  {
    return(0.4958034541*x+0.1065766549);
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t16;
  double t17;
  double t18;
  double t20;
  double t22;
  double t3;
  double t30;
  double t32;
  double t34;
  double t4;
  double t7;
  {
    t3 = phi*phi;
    t4 = t3*t3;
    t7 = t3*phi;
    t16 = 6.0*t4*phi;
    t17 = 15.0*t4;
    t18 = 10.0*t7;
    t20 = rho*rho;
    t22 = t20*t20;
    t30 = e(rho-0.1065766549);
    t32 = pow(t30-5.0,2.0);
    t34 = log(rho);
    return(2.0/delta*A*(t4+(2.0*alpha-2.0)*t7+(-3.0*alpha+1.0)*t3+alpha)+(t16-
t17+t18)*(0.8575855211E-1*t22*t20*rho-0.2868132559E-1*rho+0.1480890852E-1)+(1.0
-t16+t17-t18)*(t32+av*rho*t34+(bv-av)*rho+0.7179291685E-1-0.1135858337E-1*e));
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t17;
  double t21;
  double t22;
  double t24;
  double t3;
  double t32;
  double t34;
  double t36;
  double t4;
  {
    t3 = phi*phi;
    t4 = t3*phi;
    t17 = t3*t3;
    t21 = 30.0*t17-60.0*t4+30.0*t3;
    t22 = rho*rho;
    t24 = t22*t22;
    t32 = e(rho-0.1065766549);
    t34 = pow(t32-5.0,2.0);
    t36 = log(rho);
    return(2.0/delta*A*(4.0*t4+3.0*(2.0*alpha-2.0)*t3+2.0*(-3.0*alpha+1.0)*phi)
+t21*(0.8575855211E-1*t24*t22*rho-0.2868132559E-1*rho+0.1480890852E-1)-t21*(t34
+av*rho*t36+(bv-av)*rho+0.7179291685E-1-0.1135858337E-1*e));
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  double t17;
  double t18;
  double t20;
  double t28;
  double t3;
  double t30;
  double t32;
  {
    t3 = phi*phi;
    t17 = 120.0*t3*phi-180.0*t3+60.0*phi;
    t18 = rho*rho;
    t20 = t18*t18;
    t28 = e(rho-0.1065766549);
    t30 = pow(t28-5.0,2.0);
    t32 = log(rho);
    return(2.0/delta*A*(12.0*t3+6.0*(2.0*alpha-2.0)*phi-6.0*alpha+2.0)+t17*(
0.8575855211E-1*t20*t18*rho-0.2868132559E-1*rho+0.1480890852E-1)-t17*(t30+av*
rho*t32+(bv-av)*rho+0.7179291685E-1-0.1135858337E-1*e));
  }
}






double evalRho(double x)
{
  {
    return(0.4958034541*x+0.1065766549);
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t16;
  double t17;
  double t18;
  double t20;
  double t22;
  double t29;
  double t3;
  double t4;
  double t7;
  {
    t3 = phi*phi;
    t4 = t3*t3;
    t7 = t3*phi;
    t16 = 6.0*t4*phi;
    t17 = 15.0*t4;
    t18 = 10.0*t7;
    t20 = rho*rho;
    t22 = t20*t20;
    t29 = log(rho);
    return(2.0/delta*A*(t4+(2.0*alpha-2.0)*t7+(-3.0*alpha+1.0)*t3+alpha)+(t16-
t17+t18)*(0.8575855211E-1*t22*t20*rho-0.2868132559E-1*rho+0.1480890852E-1)+(1.0
-t16+t17-t18)*(0.15E-1+0.1407437685*rho*t29+0.1743661584*rho));
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t17;
  double t21;
  double t22;
  double t24;
  double t3;
  double t31;
  double t4;
  {
    t3 = phi*phi;
    t4 = t3*phi;
    t17 = t3*t3;
    t21 = 30.0*t17-60.0*t4+30.0*t3;
    t22 = rho*rho;
    t24 = t22*t22;
    t31 = log(rho);
    return(2.0/delta*A*(4.0*t4+3.0*(2.0*alpha-2.0)*t3+2.0*(-3.0*alpha+1.0)*phi)
+t21*(0.8575855211E-1*t24*t22*rho-0.2868132559E-1*rho+0.1480890852E-1)-t21*(
0.15E-1+0.1407437685*rho*t31+0.1743661584*rho));
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  double t17;
  double t18;
  double t20;
  double t27;
  double t3;
  {
    t3 = phi*phi;
    t17 = 120.0*t3*phi-180.0*t3+60.0*phi;
    t18 = rho*rho;
    t20 = t18*t18;
    t27 = log(rho);
    return(2.0/delta*A*(12.0*t3+6.0*(2.0*alpha-2.0)*phi-6.0*alpha+2.0)+t17*(
0.8575855211E-1*t20*t18*rho-0.2868132559E-1*rho+0.1480890852E-1)-t17*(0.15E-1+
0.1407437685*rho*t27+0.1743661584*rho));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t13;
  double t17;
  double t2;
  double t3;
  double t4;
  double t5;
  double t6;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = t2*phi;
    t4 = rho*rho;
    t5 = t4*t4;
    t6 = t5*t4;
    t13 = t1*phi;
    t17 = log(rho);
    return(0.3601859189E1*t3*t6-0.2062747515E1*t3-0.9004647972E1*t2*t6+
0.5156868787E1*t2+0.6003098648E1*t13*t6-0.3437912525E1*t13+0.3151099269+
0.1407437685*t17-0.844462611*t3*t17+0.2111156528E1*t2*t17-0.1407437685E1*t13*
t17);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  double t1;
  double t11;
  double t2;
  double t3;
  double t4;
  double t5;
  double t6;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = t2*phi;
    t4 = rho*rho;
    t5 = t4*t4;
    t6 = t5*t4;
    t11 = t1*phi;
    return(0.5E-9*(0.4322231027E11*t3*t6-0.1080557757E12*t2*t6+0.7203718378E11*
t11*t6+281487537.0-1688925222.0*t3+4222313055.0*t2-2814875370.0*t11)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t16;
  double t2;
  double t3;
  double t4;
  double t5;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = rho*rho;
    t4 = t3*t3;
    t5 = t4*t3;
    t9 = t1*phi;
    t16 = log(rho);
    return(0.1800929594E2*t2*t5-0.1031373757E2*t2-0.3601859189E2*t9*t5+
0.2062747515E2*t9+0.1800929594E2*t1*t5-0.1031373757E2*t1-0.4222313055E1*t16*t2+
0.844462611E1*t16*t9-0.4222313055E1*t1*t16);
  }
}

#include <math.h>
double pressure(double rho,double phi)
{
  double t1;
  double t11;
  double t2;
  double t22;
  double t4;
  double t5;
  double t6;
  double t7;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t6 = t1*phi;
    t7 = 10.0*t6;
    t9 = rho*rho;
    t11 = t9*t9;
    t22 = log(rho);
    return((t4-t5+t7)*(-0.8575855211E-1*t11*t9*rho+0.2868132559E-1*rho
-0.1480890852E-1+rho*(0.6003098648*t11*t9-0.2868132559E-1))+(1.0-t4+t5-t7)*(
-0.15E-1-0.1407437685*rho*t22-0.1743661584*rho+rho*(0.3151099269+0.1407437685*
t22))-2.0/delta*A*(t2+(2.0*alpha-2.0)*t6+(-3.0*alpha+1.0)*t1+alpha));
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t1;
  double t11;
  double t18;
  double t2;
  double t3;
  double t4;
  double t5;
  double t6;
  {
    t1 = rho*rho;
    t2 = t1*t1;
    t3 = t2*t1;
    t4 = phi*phi;
    t5 = t4*t4;
    t6 = t5*phi;
    t11 = t4*phi;
    t18 = sqrt(0.2161115513E12*t3*t6-0.5402788784E12*t3*t5+0.3601859189E12*t3*
t11+1407437685.0-8444626110.0*t6+0.2111156528E11*t5-0.1407437685E11*t11);
    return(0.1E-4*t18);
  }
}

double evalRho(double x)
{
  {
    return(0.4958034541*x+0.1065766549);
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t16;
  double t17;
  double t18;
  double t20;
  double t22;
  double t3;
  double t30;
  double t4;
  double t7;
  {
    t3 = phi*phi;
    t4 = t3*t3;
    t7 = t3*phi;
    t16 = 6.0*t4*phi;
    t17 = 15.0*t4;
    t18 = 10.0*t7;
    t20 = rho*rho;
    t22 = t20*t20;
    t30 = log(rho);
    return(2.0/delta*A*(t4+(2.0*alpha-2.0)*t7+(-3.0*alpha+1.0)*t3+alpha)+(t16-
t17+t18)*(0.8575855211E-1*t22*t20*rho-0.2868132559E-1*rho+0.1480890852E-1)+(1.0
-t16+t17-t18)*(av*rho*t30+(bv-av)*rho+cv));
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t17;
  double t21;
  double t22;
  double t24;
  double t3;
  double t32;
  double t4;
  {
    t3 = phi*phi;
    t4 = t3*phi;
    t17 = t3*t3;
    t21 = 30.0*t17-60.0*t4+30.0*t3;
    t22 = rho*rho;
    t24 = t22*t22;
    t32 = log(rho);
    return(2.0/delta*A*(4.0*t4+3.0*(2.0*alpha-2.0)*t3+2.0*(-3.0*alpha+1.0)*phi)
+t21*(0.8575855211E-1*t24*t22*rho-0.2868132559E-1*rho+0.1480890852E-1)-t21*(av*
rho*t32+(bv-av)*rho+cv));
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  double t17;
  double t18;
  double t20;
  double t28;
  double t3;
  {
    t3 = phi*phi;
    t17 = 120.0*t3*phi-180.0*t3+60.0*phi;
    t18 = rho*rho;
    t20 = t18*t18;
    t28 = log(rho);
    return(2.0/delta*A*(12.0*t3+6.0*(2.0*alpha-2.0)*phi-6.0*alpha+2.0)+t17*(
0.8575855211E-1*t20*t18*rho-0.2868132559E-1*rho+0.1480890852E-1)-t17*(av*rho*
t28+(bv-av)*rho+cv));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t13;
  double t17;
  double t2;
  double t3;
  double t4;
  double t5;
  double t6;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = t2*phi;
    t4 = rho*rho;
    t5 = t4*t4;
    t6 = t5*t4;
    t13 = t1*phi;
    t17 = log(rho);
    return(0.3601859189E1*t3*t6-0.1720879535*t3-0.9004647972E1*t2*t6+
0.4302198838*t2+0.6003098648E1*t13*t6-0.2868132559*t13+av*t17+bv-6.0*t3*av*t17
-6.0*t3*bv+15.0*t2*av*t17+15.0*t2*bv-10.0*t13*av*t17-10.0*t13*bv);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  double t1;
  double t11;
  double t2;
  double t3;
  double t4;
  double t5;
  double t6;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = t2*phi;
    t4 = rho*rho;
    t5 = t4*t4;
    t6 = t5*t4;
    t11 = t1*phi;
    return(-0.1E-8*(-0.2161115513E11*t3*t6+0.5402788784E11*t2*t6
-0.3601859189E11*t11*t6-1000000000.0*av+6000000000.0*av*t3-0.15E11*av*t2+0.1E11
*av*t11)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t17;
  double t2;
  double t3;
  double t4;
  double t5;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = rho*rho;
    t4 = t3*t3;
    t5 = t4*t3;
    t9 = t1*phi;
    t17 = log(rho);
    return(0.1800929594E2*t2*t5-0.8604397677*t2-0.3601859189E2*t9*t5+
0.1720879535E1*t9+0.1800929594E2*t1*t5-0.8604397677*t1-30.0*t2*av*t17-30.0*t2*
bv+60.0*t9*av*t17+60.0*t9*bv-30.0*t1*av*t17-30.0*t1*bv);
  }
}

#include <math.h>
double pressure(double rho,double phi)
{
  double t1;
  double t11;
  double t2;
  double t23;
  double t4;
  double t5;
  double t6;
  double t7;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t6 = t1*phi;
    t7 = 10.0*t6;
    t9 = rho*rho;
    t11 = t9*t9;
    t23 = log(rho);
    return((t4-t5+t7)*(-0.8575855211E-1*t11*t9*rho+0.2868132559E-1*rho
-0.1480890852E-1+rho*(0.6003098648*t11*t9-0.2868132559E-1))+(1.0-t4+t5-t7)*(-av
*rho*t23-(bv-av)*rho-cv+rho*(av*t23+bv))-2.0/delta*A*(t2+(2.0*alpha-2.0)*t6+(
-3.0*alpha+1.0)*t1+alpha));
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t1;
  double t11;
  double t2;
  double t22;
  double t3;
  double t4;
  double t5;
  double t6;
  {
    t1 = rho*rho;
    t2 = t1*t1;
    t3 = t2*t1;
    t4 = phi*phi;
    t5 = t4*t4;
    t6 = t5*phi;
    t11 = t4*phi;
    t22 = sqrt(2161115513.0*t3*t6-5402788784.0*t5*t3+3601859189.0*t11*t3+
100000000.0*av-600000000.0*av*t6+1500000000.0*av*t5-1000000000.0*av*t11);
    return(0.1E-3*t22);
  }
}

double evalRho(double x)
{
  {
    return(0.4958034541*x+0.1065766549);
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t16;
  double t17;
  double t18;
  double t20;
  double t22;
  double t29;
  double t3;
  double t4;
  double t7;
  {
    t3 = phi*phi;
    t4 = t3*t3;
    t7 = t3*phi;
    t16 = 6.0*t4*phi;
    t17 = 15.0*t4;
    t18 = 10.0*t7;
    t20 = rho*rho;
    t22 = t20*t20;
    t29 = log(rho);
    return(2.0/delta*A*(t4+(2.0*alpha-2.0)*t7+(-3.0*alpha+1.0)*t3+alpha)+(t16-
t17+t18)*(0.8575855211E-1*t22*t20*rho-0.2868132559E-1*rho+0.1480890852E-1)+(1.0
-t16+t17-t18)*(0.3198902536E-2+0.3001504024E-1*rho*t29+0.3718535686E-1*rho));
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t17;
  double t21;
  double t22;
  double t24;
  double t3;
  double t31;
  double t4;
  {
    t3 = phi*phi;
    t4 = t3*phi;
    t17 = t3*t3;
    t21 = 30.0*t17-60.0*t4+30.0*t3;
    t22 = rho*rho;
    t24 = t22*t22;
    t31 = log(rho);
    return(2.0/delta*A*(4.0*t4+3.0*(2.0*alpha-2.0)*t3+2.0*(-3.0*alpha+1.0)*phi)
+t21*(0.8575855211E-1*t24*t22*rho-0.2868132559E-1*rho+0.1480890852E-1)-t21*(
0.3198902536E-2+0.3001504024E-1*rho*t31+0.3718535686E-1*rho));
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  double t17;
  double t18;
  double t20;
  double t27;
  double t3;
  {
    t3 = phi*phi;
    t17 = 120.0*t3*phi-180.0*t3+60.0*phi;
    t18 = rho*rho;
    t20 = t18*t18;
    t27 = log(rho);
    return(2.0/delta*A*(12.0*t3+6.0*(2.0*alpha-2.0)*phi-6.0*alpha+2.0)+t17*(
0.8575855211E-1*t20*t18*rho-0.2868132559E-1*rho+0.1480890852E-1)-t17*(
0.3198902536E-2+0.3001504024E-1*rho*t27+0.3718535686E-1*rho));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t13;
  double t17;
  double t2;
  double t3;
  double t4;
  double t5;
  double t6;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = t2*phi;
    t4 = rho*rho;
    t5 = t4*t4;
    t6 = t5*t4;
    t13 = t1*phi;
    t17 = log(rho);
    return(0.3601859189E1*t3*t6-0.5752903361*t3-0.9004647972E1*t2*t6+
0.143822584E1*t2+0.6003098648E1*t13*t6-0.9588172269*t13+0.672003971E-1+
0.3001504024E-1*t17-0.1800902414*t3*t17+0.4502256036*t2*t17-0.3001504024*t13*
t17);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  double t1;
  double t11;
  double t2;
  double t3;
  double t4;
  double t5;
  double t6;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = t2*phi;
    t4 = rho*rho;
    t5 = t4*t4;
    t6 = t5*t4;
    t11 = t1*phi;
    return(0.4E-10*(0.5402788784E12*t3*t6-0.1350697196E13*t2*t6+0.9004647972E12
*t11*t6+750376006.0-4502256036.0*t3+0.1125564009E11*t2-7503760060.0*t11)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t16;
  double t2;
  double t3;
  double t4;
  double t5;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = rho*rho;
    t4 = t3*t3;
    t5 = t4*t3;
    t9 = t1*phi;
    t16 = log(rho);
    return(0.1800929594E2*t2*t5-0.2876451681E1*t2-0.3601859189E2*t9*t5+
0.5752903361E1*t9+0.1800929594E2*t1*t5-0.2876451681E1*t1-0.9004512072*t2*t16+
0.1800902414E1*t9*t16-0.9004512072*t1*t16);
  }
}

#include <math.h>
double pressure(double rho,double phi)
{
  double t1;
  double t11;
  double t2;
  double t22;
  double t4;
  double t5;
  double t6;
  double t7;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t6 = t1*phi;
    t7 = 10.0*t6;
    t9 = rho*rho;
    t11 = t9*t9;
    t22 = log(rho);
    return((t4-t5+t7)*(-0.8575855211E-1*t11*t9*rho+0.2868132559E-1*rho
-0.1480890852E-1+rho*(0.6003098648*t11*t9-0.2868132559E-1))+(1.0-t4+t5-t7)*(
-0.3198902536E-2-0.3001504024E-1*rho*t22-0.3718535686E-1*rho+rho*(
0.672003971E-1+0.3001504024E-1*t22))-2.0/delta*A*(t2+(2.0*alpha-2.0)*t6+(-3.0*
alpha+1.0)*t1+alpha));
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t1;
  double t11;
  double t18;
  double t2;
  double t3;
  double t4;
  double t5;
  double t6;
  {
    t1 = rho*rho;
    t2 = t1*t1;
    t3 = t2*t1;
    t4 = phi*phi;
    t5 = t4*t4;
    t6 = t5*phi;
    t11 = t4*phi;
    t18 = sqrt(0.5402788782E13*t3*t6-0.1350697196E14*t3*t5+0.9004647972E13*t3*
t11+7503760060.0-0.4502256035E11*t6+0.1125564009E12*t5-0.750376006E11*t11);
    return(0.2E-5*t18);
  }
}

double evalRho(double x)
{
  {
    return(0.4958034541*x+0.1065766549);
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t16;
  double t17;
  double t18;
  double t20;
  double t22;
  double t29;
  double t3;
  double t4;
  double t7;
  {
    t3 = phi*phi;
    t4 = t3*t3;
    t7 = t3*phi;
    t16 = 6.0*t4*phi;
    t17 = 15.0*t4;
    t18 = 10.0*t7;
    t20 = rho*rho;
    t22 = t20*t20;
    t29 = log(rho);
    return(2.0/delta*A*(t4+(2.0*alpha-2.0)*t7+(-3.0*alpha+1.0)*t3+alpha)+(t16-
t17+t18)*(0.8575855211E-1*t22*t20*rho-0.2868132559E-1*rho+0.1480890852E-1)+(1.0
-t16+t17-t18)*(0.3001504024E-1*rho*t29+0.3718535686E-1*rho+0.3198902536E-2));
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t17;
  double t21;
  double t22;
  double t24;
  double t3;
  double t31;
  double t4;
  {
    t3 = phi*phi;
    t4 = t3*phi;
    t17 = t3*t3;
    t21 = 30.0*t17-60.0*t4+30.0*t3;
    t22 = rho*rho;
    t24 = t22*t22;
    t31 = log(rho);
    return(2.0/delta*A*(4.0*t4+3.0*(2.0*alpha-2.0)*t3+2.0*(-3.0*alpha+1.0)*phi)
+t21*(0.8575855211E-1*t24*t22*rho-0.2868132559E-1*rho+0.1480890852E-1)-t21*(
0.3001504024E-1*rho*t31+0.3718535686E-1*rho+0.3198902536E-2));
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  double t17;
  double t18;
  double t20;
  double t27;
  double t3;
  {
    t3 = phi*phi;
    t17 = 120.0*t3*phi-180.0*t3+60.0*phi;
    t18 = rho*rho;
    t20 = t18*t18;
    t27 = log(rho);
    return(2.0/delta*A*(12.0*t3+6.0*(2.0*alpha-2.0)*phi-6.0*alpha+2.0)+t17*(
0.8575855211E-1*t20*t18*rho-0.2868132559E-1*rho+0.1480890852E-1)-t17*(
0.3001504024E-1*rho*t27+0.3718535686E-1*rho+0.3198902536E-2));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t13;
  double t17;
  double t2;
  double t3;
  double t4;
  double t5;
  double t6;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = t2*phi;
    t4 = rho*rho;
    t5 = t4*t4;
    t6 = t5*t4;
    t13 = t1*phi;
    t17 = log(rho);
    return(0.3601859189E1*t3*t6-0.5752903361*t3-0.9004647972E1*t2*t6+
0.143822584E1*t2+0.6003098648E1*t13*t6-0.9588172269*t13+0.3001504024E-1*t17+
0.672003971E-1-0.1800902414*t3*t17+0.4502256036*t2*t17-0.3001504024*t13*t17);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  double t1;
  double t11;
  double t2;
  double t3;
  double t4;
  double t5;
  double t6;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = t2*phi;
    t4 = rho*rho;
    t5 = t4*t4;
    t6 = t5*t4;
    t11 = t1*phi;
    return(0.4E-10*(0.5402788784E12*t3*t6-0.1350697196E13*t2*t6+0.9004647972E12
*t11*t6+750376006.0-4502256036.0*t3+0.1125564009E11*t2-7503760060.0*t11)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t16;
  double t2;
  double t3;
  double t4;
  double t5;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = rho*rho;
    t4 = t3*t3;
    t5 = t4*t3;
    t9 = t1*phi;
    t16 = log(rho);
    return(0.1800929594E2*t2*t5-0.2876451681E1*t2-0.3601859189E2*t9*t5+
0.5752903361E1*t9+0.1800929594E2*t1*t5-0.2876451681E1*t1-0.9004512072*t2*t16+
0.1800902414E1*t9*t16-0.9004512072*t1*t16);
  }
}

#include <math.h>
double pressure(double rho,double phi)
{
  double t1;
  double t11;
  double t2;
  double t22;
  double t4;
  double t5;
  double t6;
  double t7;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t6 = t1*phi;
    t7 = 10.0*t6;
    t9 = rho*rho;
    t11 = t9*t9;
    t22 = log(rho);
    return((t4-t5+t7)*(-0.8575855211E-1*t11*t9*rho+0.2868132559E-1*rho
-0.1480890852E-1+rho*(0.6003098648*t11*t9-0.2868132559E-1))+(1.0-t4+t5-t7)*(
-0.3001504024E-1*rho*t22-0.3718535686E-1*rho-0.3198902536E-2+rho*(
0.3001504024E-1*t22+0.672003971E-1))-2.0/delta*A*(t2+(2.0*alpha-2.0)*t6+(-3.0*
alpha+1.0)*t1+alpha));
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t1;
  double t11;
  double t18;
  double t2;
  double t3;
  double t4;
  double t5;
  double t6;
  {
    t1 = rho*rho;
    t2 = t1*t1;
    t3 = t2*t1;
    t4 = phi*phi;
    t5 = t4*t4;
    t6 = t5*phi;
    t11 = t4*phi;
    t18 = sqrt(0.5402788782E13*t3*t6-0.1350697196E14*t3*t5+0.9004647972E13*t3*
t11+7503760060.0-0.4502256035E11*t6+0.1125564009E12*t5-0.750376006E11*t11);
    return(0.2E-5*t18);
  }
}


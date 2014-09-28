#include <math.h>
double rhosol(double x)
{
  {
    return((0.113810037E2-0.5075641578E2*pow(0.5*tanh(1/delta*x)+0.5,5.0)+
0.1268910394E3*pow(0.5*tanh(1/delta*x)+0.5,4.0)-0.845940263E2*pow(0.5*tanh(1/
delta*x)+0.5,3.0))/(-0.9E1*pow(0.5*tanh(1/delta*x)+0.5,5.0)+0.225E2*pow(0.5*
tanh(1/delta*x)+0.5,4.0)-0.15E2*pow(0.5*tanh(1/delta*x)+0.5,3.0)+3.0));
  }
}

#include <math.h>
double gradrho(double x)
{
  {
    return((-0.1268910394E3*pow(0.5*tanh(1/delta*x)+0.5,4.0)*(1.0-pow(tanh(1/
delta*x),2.0))/delta+0.2537820788E3*pow(0.5*tanh(1/delta*x)+0.5,3.0)*(1.0-pow(
tanh(1/delta*x),2.0))/delta-0.1268910394E3*pow(0.5*tanh(1/delta*x)+0.5,2.0)*(
1.0-pow(tanh(1/delta*x),2.0))/delta)/(-0.9E1*pow(0.5*tanh(1/delta*x)+0.5,5.0)+
0.225E2*pow(0.5*tanh(1/delta*x)+0.5,4.0)-0.15E2*pow(0.5*tanh(1/delta*x)+0.5,3.0
)+3.0)-(0.113810037E2-0.5075641578E2*pow(0.5*tanh(1/delta*x)+0.5,5.0)+
0.1268910394E3*pow(0.5*tanh(1/delta*x)+0.5,4.0)-0.845940263E2*pow(0.5*tanh(1/
delta*x)+0.5,3.0))/pow(-0.9E1*pow(0.5*tanh(1/delta*x)+0.5,5.0)+0.225E2*pow(0.5*
tanh(1/delta*x)+0.5,4.0)-0.15E2*pow(0.5*tanh(1/delta*x)+0.5,3.0)+3.0,2.0)*(
-0.225E2*pow(0.5*tanh(1/delta*x)+0.5,4.0)*(1.0-pow(tanh(1/delta*x),2.0))/delta+
0.45E2*pow(0.5*tanh(1/delta*x)+0.5,3.0)*(1.0-pow(tanh(1/delta*x),2.0))/delta
-0.225E2*pow(0.5*tanh(1/delta*x)+0.5,2.0)*(1.0-pow(tanh(1/delta*x),2.0))/delta)
);
  }
}

#include <math.h>
double gradrho(double x)
{
  {
    return((-0.1268910394E3*pow(0.5*tanh(1/delta*x)+0.5,4.0)*(1.0-pow(tanh(1/
delta*x),2.0))/delta+0.2537820788E3*pow(0.5*tanh(1/delta*x)+0.5,3.0)*(1.0-pow(
tanh(1/delta*x),2.0))/delta-0.1268910394E3*pow(0.5*tanh(1/delta*x)+0.5,2.0)*(
1.0-pow(tanh(1/delta*x),2.0))/delta)/(-0.9E1*pow(0.5*tanh(1/delta*x)+0.5,5.0)+
0.225E2*pow(0.5*tanh(1/delta*x)+0.5,4.0)-0.15E2*pow(0.5*tanh(1/delta*x)+0.5,3.0
)+3.0)-(0.113810037E2-0.5075641578E2*pow(0.5*tanh(1/delta*x)+0.5,5.0)+
0.1268910394E3*pow(0.5*tanh(1/delta*x)+0.5,4.0)-0.845940263E2*pow(0.5*tanh(1/
delta*x)+0.5,3.0))/pow(-0.9E1*pow(0.5*tanh(1/delta*x)+0.5,5.0)+0.225E2*pow(0.5*
tanh(1/delta*x)+0.5,4.0)-0.15E2*pow(0.5*tanh(1/delta*x)+0.5,3.0)+3.0,2.0)*(
-0.225E2*pow(0.5*tanh(1/delta*x)+0.5,4.0)*(1.0-pow(tanh(1/delta*x),2.0))/delta+
0.45E2*pow(0.5*tanh(1/delta*x)+0.5,3.0)*(1.0-pow(tanh(1/delta*x),2.0))/delta
-0.225E2*pow(0.5*tanh(1/delta*x)+0.5,2.0)*(1.0-pow(tanh(1/delta*x),2.0))/delta)
);
  }
}

#include <math.h>
double veloRhs(double x)
{
  double t1;
  double t104;
  double t11;
  double t116;
  double t118;
  double t13;
  double t17;
  double t18;
  double t19;
  double t20;
  double t21;
  double t23;
  double t24;
  double t25;
  double t26;
  double t3;
  double t35;
  double t38;
  double t4;
  double t41;
  double t43;
  double t44;
  double t48;
  double t49;
  double t5;
  double t50;
  double t55;
  double t57;
  double t6;
  double t63;
  double t65;
  double t67;
  double t68;
  double t7;
  double t78;
  double t79;
  double t8;
  double t80;
  double t81;
  double t91;
  {
    t1 = 1/delta;
    t3 = tanh(t1*x);
    t4 = 0.5*t3;
    t5 = t4+0.5;
    t6 = t5*t5;
    t7 = t6*t6;
    t8 = t7*t5;
    t11 = t6*t5;
    t13 = 0.113810037E2-0.5075641578E2*t8+0.1268910394E3*t7-0.845940263E2*t11;
    t17 = -0.9E1*t8+0.225E2*t7-0.15E2*t11+3.0;
    t18 = 1/t17;
    t19 = t13*t18;
    t20 = delta*delta;
    t21 = 1/t20;
    t23 = 0.5-t4;
    t24 = t23*t23;
    t25 = t3*t3;
    t26 = 1.0-t25;
    t35 = t7*t26*t1;
    t38 = t11*t26*t1;
    t41 = t6*t26*t1;
    t43 = 0.15E2*t35-0.3E2*t38+0.15E2*t41;
    t44 = log(t19);
    t48 = 6.0*t8;
    t49 = 15.0*t7;
    t50 = 10.0*t11;
    t55 = -0.1268910394E3*t35+0.2537820788E3*t38-0.1268910394E3*t41;
    t57 = t17*t17;
    t63 = -0.225E2*t35+0.45E2*t38-0.225E2*t41;
    t65 = t55*t18-t13/t57*t63;
    t67 = 1/t13;
    t68 = t67*t17;
    t78 = t26*t1;
    t79 = t13*t13;
    t80 = 1/t79;
    t81 = t80*t57;
    t91 = x*t26*t1;
    t104 = t1*t13;
    t116 = 30.0*t7-60.0*t11+30.0*t6;
    t118 = t19*t44;
    return(t19*(0.2E1*t21*t5*t24*t26-0.2E1*t21*t6*t23*t26+t43*(-1.0+0.15E1*t44)
+0.15E1*(t48-t49+t50)*t65*t68-t43*(-4.0+3.0*t44)+3.0*(1.0-t48+t49-t50)*t65*t68
-0.5*t78*t81+0.1E1*x*t3*t26*t21*t80*t57+0.1E1*t91/t79/t13*t57*t55-0.1E1*t91*t80
*t17*t63)-0.5*t78*(4.0*t104*t18*t5*t24-4.0*t104*t18*t6*t23+t116*(-0.25E1*t19+
0.15E1*t118+1.0)-t116*(-7.0*t19+3.0*t118+0.945940263E1)+0.1E1*t1*t67*t17*t3*t26
+0.5*t65*t26*t81));
  }
}

#include <math.h>
double gradphi(double x)
{
  {
    return(0.1E2-0.1E2*pow(tanh(0.2E2*x),2.0));
  }
}

      t0 = phiSourcel;
#include <math.h>
double musol(double x)
{
  {
    return(0.4E2*pow(0.5*tanh(0.2E2*x)+0.5,2.0)*pow(0.5-0.5*tanh(0.2E2*x),2.0)+
(6.0*pow(0.5*tanh(0.2E2*x)+0.5,5.0)-15.0*pow(0.5*tanh(0.2E2*x)+0.5,4.0)+10.0*
pow(0.5*tanh(0.2E2*x)+0.5,3.0))*(-1.0+0.15E1*log((0.113810037E2-0.5075641578E2*
pow(0.5*tanh(0.2E2*x)+0.5,5.0)+0.1268910394E3*pow(0.5*tanh(0.2E2*x)+0.5,4.0)
-0.845940263E2*pow(0.5*tanh(0.2E2*x)+0.5,3.0))/(-0.9E1*pow(0.5*tanh(0.2E2*x)+
0.5,5.0)+0.225E2*pow(0.5*tanh(0.2E2*x)+0.5,4.0)-0.15E2*pow(0.5*tanh(0.2E2*x)+
0.5,3.0)+3.0)))+(1.0-6.0*pow(0.5*tanh(0.2E2*x)+0.5,5.0)+15.0*pow(0.5*tanh(0.2E2
*x)+0.5,4.0)-10.0*pow(0.5*tanh(0.2E2*x)+0.5,3.0))*(-4.0+3.0*log((0.113810037E2
-0.5075641578E2*pow(0.5*tanh(0.2E2*x)+0.5,5.0)+0.1268910394E3*pow(0.5*tanh(
0.2E2*x)+0.5,4.0)-0.845940263E2*pow(0.5*tanh(0.2E2*x)+0.5,3.0))/(-0.9E1*pow(0.5
*tanh(0.2E2*x)+0.5,5.0)+0.225E2*pow(0.5*tanh(0.2E2*x)+0.5,4.0)-0.15E2*pow(0.5*
tanh(0.2E2*x)+0.5,3.0)+3.0)))-x*(0.1E2-0.1E2*pow(tanh(0.2E2*x),2.0))/pow(
0.113810037E2-0.5075641578E2*pow(0.5*tanh(0.2E2*x)+0.5,5.0)+0.1268910394E3*pow(
0.5*tanh(0.2E2*x)+0.5,4.0)-0.845940263E2*pow(0.5*tanh(0.2E2*x)+0.5,3.0),2.0)*
pow(-0.9E1*pow(0.5*tanh(0.2E2*x)+0.5,5.0)+0.225E2*pow(0.5*tanh(0.2E2*x)+0.5,4.0
)-0.15E2*pow(0.5*tanh(0.2E2*x)+0.5,3.0)+3.0,2.0));
  }
}

#include <math.h>
double rhosol(double x)
{
  {
    return((0.113810037E2-0.5075641578E2*pow(0.5*tanh(1/delta*x)+0.5,5.0)+
0.1268910394E3*pow(0.5*tanh(1/delta*x)+0.5,4.0)-0.845940263E2*pow(0.5*tanh(1/
delta*x)+0.5,3.0))/(-0.9E1*pow(0.5*tanh(1/delta*x)+0.5,5.0)+0.225E2*pow(0.5*
tanh(1/delta*x)+0.5,4.0)-0.15E2*pow(0.5*tanh(1/delta*x)+0.5,3.0)+3.0));
  }
}

#include <math.h>
double gradrho(double x)
{
  {
    return((-0.1268910394E3*pow(0.5*tanh(1/delta*x)+0.5,4.0)*(1.0-pow(tanh(1/
delta*x),2.0))/delta+0.2537820788E3*pow(0.5*tanh(1/delta*x)+0.5,3.0)*(1.0-pow(
tanh(1/delta*x),2.0))/delta-0.1268910394E3*pow(0.5*tanh(1/delta*x)+0.5,2.0)*(
1.0-pow(tanh(1/delta*x),2.0))/delta)/(-0.9E1*pow(0.5*tanh(1/delta*x)+0.5,5.0)+
0.225E2*pow(0.5*tanh(1/delta*x)+0.5,4.0)-0.15E2*pow(0.5*tanh(1/delta*x)+0.5,3.0
)+3.0)-(0.113810037E2-0.5075641578E2*pow(0.5*tanh(1/delta*x)+0.5,5.0)+
0.1268910394E3*pow(0.5*tanh(1/delta*x)+0.5,4.0)-0.845940263E2*pow(0.5*tanh(1/
delta*x)+0.5,3.0))/pow(-0.9E1*pow(0.5*tanh(1/delta*x)+0.5,5.0)+0.225E2*pow(0.5*
tanh(1/delta*x)+0.5,4.0)-0.15E2*pow(0.5*tanh(1/delta*x)+0.5,3.0)+3.0,2.0)*(
-0.225E2*pow(0.5*tanh(1/delta*x)+0.5,4.0)*(1.0-pow(tanh(1/delta*x),2.0))/delta+
0.45E2*pow(0.5*tanh(1/delta*x)+0.5,3.0)*(1.0-pow(tanh(1/delta*x),2.0))/delta
-0.225E2*pow(0.5*tanh(1/delta*x)+0.5,2.0)*(1.0-pow(tanh(1/delta*x),2.0))/delta)
);
  }
}

#include <math.h>
double gradphi(double x)
{
  {
    return(0.5*(1.0-pow(tanh(1/delta*x),2.0))/delta);
  }
}

#include <math.h>
double thetasol(double x)
{
  double t1;
  double t11;
  double t13;
  double t14;
  double t18;
  double t19;
  double t21;
  double t22;
  double t3;
  double t33;
  double t34;
  double t36;
  double t37;
  double t4;
  double t48;
  double t49;
  double t5;
  double t54;
  double t57;
  double t6;
  double t60;
  double t64;
  double t7;
  double t74;
  double t8;
  {
    t1 = 1/delta;
    t3 = tanh(t1*x);
    t4 = 0.5*t3;
    t5 = t4+0.5;
    t6 = t5*t5;
    t7 = t6*t6;
    t8 = t7*t5;
    t11 = t6*t5;
    t13 = 0.113810037E2-0.5075641578E2*t8+0.1268910394E3*t7-0.845940263E2*t11;
    t14 = t1*t13;
    t18 = -0.9E1*t8+0.225E2*t7-0.15E2*t11+3.0;
    t19 = 1/t18;
    t21 = 0.5-t4;
    t22 = t21*t21;
    t33 = 30.0*t7-60.0*t11+30.0*t6;
    t34 = t13*t19;
    t36 = log(t34);
    t37 = t34*t36;
    t48 = t3*t3;
    t49 = 1.0-t48;
    t54 = t7*t49*t1;
    t57 = t11*t49*t1;
    t60 = t6*t49*t1;
    t64 = t18*t18;
    t74 = t13*t13;
    return(4.0*t14*t19*t5*t22-4.0*t14*t19*t6*t21+t33*(-0.25E1*t34+0.15E1*t37+
1.0)-t33*(-7.0*t34+3.0*t37+0.945940263E1)+0.1E1*t1/t13*t18*t3*t49+0.5*((
-0.1268910394E3*t54+0.2537820788E3*t57-0.1268910394E3*t60)*t19-t13/t64*(
-0.225E2*t54+0.45E2*t57-0.225E2*t60))*t49/t74*t64);
  }
}

      t0 = phiSourcel;
#include <math.h>
double musol(double x)
{
  {
    return(2.0/delta*pow(0.5*tanh(1/delta*x)+0.5,2.0)*pow(0.5-0.5*tanh(1/delta*
x),2.0)+(6.0*pow(0.5*tanh(1/delta*x)+0.5,5.0)-15.0*pow(0.5*tanh(1/delta*x)+0.5,
4.0)+10.0*pow(0.5*tanh(1/delta*x)+0.5,3.0))*(-1.0+0.15E1*log((0.113810037E2
-0.5075641578E2*pow(0.5*tanh(1/delta*x)+0.5,5.0)+0.1268910394E3*pow(0.5*tanh(1/
delta*x)+0.5,4.0)-0.845940263E2*pow(0.5*tanh(1/delta*x)+0.5,3.0))/(-0.9E1*pow(
0.5*tanh(1/delta*x)+0.5,5.0)+0.225E2*pow(0.5*tanh(1/delta*x)+0.5,4.0)-0.15E2*
pow(0.5*tanh(1/delta*x)+0.5,3.0)+3.0)))+(1.0-6.0*pow(0.5*tanh(1/delta*x)+0.5,
5.0)+15.0*pow(0.5*tanh(1/delta*x)+0.5,4.0)-10.0*pow(0.5*tanh(1/delta*x)+0.5,3.0
))*(-4.0+3.0*log((0.113810037E2-0.5075641578E2*pow(0.5*tanh(1/delta*x)+0.5,5.0)
+0.1268910394E3*pow(0.5*tanh(1/delta*x)+0.5,4.0)-0.845940263E2*pow(0.5*tanh(1/
delta*x)+0.5,3.0))/(-0.9E1*pow(0.5*tanh(1/delta*x)+0.5,5.0)+0.225E2*pow(0.5*
tanh(1/delta*x)+0.5,4.0)-0.15E2*pow(0.5*tanh(1/delta*x)+0.5,3.0)+3.0)))-0.5*x*(
1.0-pow(tanh(1/delta*x),2.0))/delta/pow(0.113810037E2-0.5075641578E2*pow(0.5*
tanh(1/delta*x)+0.5,5.0)+0.1268910394E3*pow(0.5*tanh(1/delta*x)+0.5,4.0)
-0.845940263E2*pow(0.5*tanh(1/delta*x)+0.5,3.0),2.0)*pow(-0.9E1*pow(0.5*tanh(1/
delta*x)+0.5,5.0)+0.225E2*pow(0.5*tanh(1/delta*x)+0.5,4.0)-0.15E2*pow(0.5*tanh(
1/delta*x)+0.5,3.0)+3.0,2.0));
  }
}

#include <math.h>
double veloRhs(double x)
{
  double t1;
  double t11;
  double t13;
  double t17;
  double t18;
  double t19;
  double t20;
  double t22;
  double t23;
  double t25;
  double t27;
  double t28;
  double t3;
  double t39;
  double t4;
  double t41;
  double t42;
  double t5;
  double t57;
  double t6;
  double t60;
  double t63;
  double t67;
  double t7;
  double t77;
  double t8;
  {
    t1 = 1/delta;
    t3 = tanh(t1*x);
    t4 = 0.5*t3;
    t5 = t4+0.5;
    t6 = t5*t5;
    t7 = t6*t6;
    t8 = t7*t5;
    t11 = t6*t5;
    t13 = 0.113810037E2-0.5075641578E2*t8+0.1268910394E3*t7-0.845940263E2*t11;
    t17 = -0.9E1*t8+0.225E2*t7-0.15E2*t11+3.0;
    t18 = 1/t17;
    t19 = t13*t18;
    t20 = diff(solmu(x),x);
    t22 = t3*t3;
    t23 = 1.0-t22;
    t25 = t1*t13;
    t27 = 0.5-t4;
    t28 = t27*t27;
    t39 = 30.0*t7-60.0*t11+30.0*t6;
    t41 = log(t19);
    t42 = t19*t41;
    t57 = t7*t23*t1;
    t60 = t11*t23*t1;
    t63 = t6*t23*t1;
    t67 = t17*t17;
    t77 = t13*t13;
    return(t19*t20-0.5*t23*t1*(4.0*t25*t18*t5*t28-4.0*t25*t18*t6*t27+t39*(
-0.25E1*t19+0.15E1*t42+1.0)-t39*(-7.0*t19+3.0*t42+0.945940263E1)+0.1E1*t1/t13*
t17*t3*t23+0.5*((-0.1268910394E3*t57+0.2537820788E3*t60-0.1268910394E3*t63)*t18
-t13/t67*(-0.225E2*t57+0.45E2*t60-0.225E2*t63))*t23/t77*t67));
  }
}

#include <math.h>
double rhosol(double x)
{
  {
    return((0.113810037E2-0.5075641578E2*pow(0.5*tanh(1/delta*x)+0.5,5.0)+
0.1268910394E3*pow(0.5*tanh(1/delta*x)+0.5,4.0)-0.845940263E2*pow(0.5*tanh(1/
delta*x)+0.5,3.0))/(-0.9E1*pow(0.5*tanh(1/delta*x)+0.5,5.0)+0.225E2*pow(0.5*
tanh(1/delta*x)+0.5,4.0)-0.15E2*pow(0.5*tanh(1/delta*x)+0.5,3.0)+3.0));
  }
}

#include <math.h>
double gradrho(double x)
{
  {
    return((-0.1268910394E3*pow(0.5*tanh(1/delta*x)+0.5,4.0)*(1.0-pow(tanh(1/
delta*x),2.0))/delta+0.2537820788E3*pow(0.5*tanh(1/delta*x)+0.5,3.0)*(1.0-pow(
tanh(1/delta*x),2.0))/delta-0.1268910394E3*pow(0.5*tanh(1/delta*x)+0.5,2.0)*(
1.0-pow(tanh(1/delta*x),2.0))/delta)/(-0.9E1*pow(0.5*tanh(1/delta*x)+0.5,5.0)+
0.225E2*pow(0.5*tanh(1/delta*x)+0.5,4.0)-0.15E2*pow(0.5*tanh(1/delta*x)+0.5,3.0
)+3.0)-(0.113810037E2-0.5075641578E2*pow(0.5*tanh(1/delta*x)+0.5,5.0)+
0.1268910394E3*pow(0.5*tanh(1/delta*x)+0.5,4.0)-0.845940263E2*pow(0.5*tanh(1/
delta*x)+0.5,3.0))/pow(-0.9E1*pow(0.5*tanh(1/delta*x)+0.5,5.0)+0.225E2*pow(0.5*
tanh(1/delta*x)+0.5,4.0)-0.15E2*pow(0.5*tanh(1/delta*x)+0.5,3.0)+3.0,2.0)*(
-0.225E2*pow(0.5*tanh(1/delta*x)+0.5,4.0)*(1.0-pow(tanh(1/delta*x),2.0))/delta+
0.45E2*pow(0.5*tanh(1/delta*x)+0.5,3.0)*(1.0-pow(tanh(1/delta*x),2.0))/delta
-0.225E2*pow(0.5*tanh(1/delta*x)+0.5,2.0)*(1.0-pow(tanh(1/delta*x),2.0))/delta)
);
  }
}

#include <math.h>
double gradphi(double x)
{
  {
    return(0.5*(1.0-pow(tanh(1/delta*x),2.0))/delta);
  }
}

#include <math.h>
double thetasol(double x)
{
  double t1;
  double t11;
  double t13;
  double t14;
  double t18;
  double t19;
  double t21;
  double t22;
  double t3;
  double t33;
  double t34;
  double t36;
  double t37;
  double t4;
  double t48;
  double t49;
  double t5;
  double t54;
  double t57;
  double t6;
  double t60;
  double t64;
  double t7;
  double t74;
  double t8;
  {
    t1 = 1/delta;
    t3 = tanh(t1*x);
    t4 = 0.5*t3;
    t5 = t4+0.5;
    t6 = t5*t5;
    t7 = t6*t6;
    t8 = t7*t5;
    t11 = t6*t5;
    t13 = 0.113810037E2-0.5075641578E2*t8+0.1268910394E3*t7-0.845940263E2*t11;
    t14 = t1*t13;
    t18 = -0.9E1*t8+0.225E2*t7-0.15E2*t11+3.0;
    t19 = 1/t18;
    t21 = 0.5-t4;
    t22 = t21*t21;
    t33 = 30.0*t7-60.0*t11+30.0*t6;
    t34 = t13*t19;
    t36 = log(t34);
    t37 = t34*t36;
    t48 = t3*t3;
    t49 = 1.0-t48;
    t54 = t7*t49*t1;
    t57 = t11*t49*t1;
    t60 = t6*t49*t1;
    t64 = t18*t18;
    t74 = t13*t13;
    return(4.0*t14*t19*t5*t22-4.0*t14*t19*t6*t21+t33*(-0.25E1*t34+0.15E1*t37+
1.0)-t33*(-7.0*t34+3.0*t37+0.945940263E1)+0.1E1*t1/t13*t18*t3*t49+0.5*((
-0.1268910394E3*t54+0.2537820788E3*t57-0.1268910394E3*t60)*t19-t13/t64*(
-0.225E2*t54+0.45E2*t57-0.225E2*t60))*t49/t74*t64);
  }
}

      t0 = phiSourcel;
#include <math.h>
double musol(double x)
{
  {
    return(2.0/delta*pow(0.5*tanh(1/delta*x)+0.5,2.0)*pow(0.5-0.5*tanh(1/delta*
x),2.0)+(6.0*pow(0.5*tanh(1/delta*x)+0.5,5.0)-15.0*pow(0.5*tanh(1/delta*x)+0.5,
4.0)+10.0*pow(0.5*tanh(1/delta*x)+0.5,3.0))*(-1.0+0.15E1*log((0.113810037E2
-0.5075641578E2*pow(0.5*tanh(1/delta*x)+0.5,5.0)+0.1268910394E3*pow(0.5*tanh(1/
delta*x)+0.5,4.0)-0.845940263E2*pow(0.5*tanh(1/delta*x)+0.5,3.0))/(-0.9E1*pow(
0.5*tanh(1/delta*x)+0.5,5.0)+0.225E2*pow(0.5*tanh(1/delta*x)+0.5,4.0)-0.15E2*
pow(0.5*tanh(1/delta*x)+0.5,3.0)+3.0)))+(1.0-6.0*pow(0.5*tanh(1/delta*x)+0.5,
5.0)+15.0*pow(0.5*tanh(1/delta*x)+0.5,4.0)-10.0*pow(0.5*tanh(1/delta*x)+0.5,3.0
))*(-4.0+3.0*log((0.113810037E2-0.5075641578E2*pow(0.5*tanh(1/delta*x)+0.5,5.0)
+0.1268910394E3*pow(0.5*tanh(1/delta*x)+0.5,4.0)-0.845940263E2*pow(0.5*tanh(1/
delta*x)+0.5,3.0))/(-0.9E1*pow(0.5*tanh(1/delta*x)+0.5,5.0)+0.225E2*pow(0.5*
tanh(1/delta*x)+0.5,4.0)-0.15E2*pow(0.5*tanh(1/delta*x)+0.5,3.0)+3.0)))-0.5*x*(
1.0-pow(tanh(1/delta*x),2.0))/delta/pow(0.113810037E2-0.5075641578E2*pow(0.5*
tanh(1/delta*x)+0.5,5.0)+0.1268910394E3*pow(0.5*tanh(1/delta*x)+0.5,4.0)
-0.845940263E2*pow(0.5*tanh(1/delta*x)+0.5,3.0),2.0)*pow(-0.9E1*pow(0.5*tanh(1/
delta*x)+0.5,5.0)+0.225E2*pow(0.5*tanh(1/delta*x)+0.5,4.0)-0.15E2*pow(0.5*tanh(
1/delta*x)+0.5,3.0)+3.0,2.0));
  }
}

#include <math.h>
double veloRhs(double x)
{
  double t1;
  double t104;
  double t11;
  double t116;
  double t118;
  double t13;
  double t17;
  double t18;
  double t19;
  double t20;
  double t21;
  double t23;
  double t24;
  double t25;
  double t26;
  double t3;
  double t35;
  double t38;
  double t4;
  double t41;
  double t43;
  double t44;
  double t48;
  double t49;
  double t5;
  double t50;
  double t55;
  double t57;
  double t6;
  double t63;
  double t65;
  double t67;
  double t68;
  double t7;
  double t78;
  double t79;
  double t8;
  double t80;
  double t81;
  double t91;
  {
    t1 = 1/delta;
    t3 = tanh(t1*x);
    t4 = 0.5*t3;
    t5 = t4+0.5;
    t6 = t5*t5;
    t7 = t6*t6;
    t8 = t7*t5;
    t11 = t6*t5;
    t13 = 0.113810037E2-0.5075641578E2*t8+0.1268910394E3*t7-0.845940263E2*t11;
    t17 = -0.9E1*t8+0.225E2*t7-0.15E2*t11+3.0;
    t18 = 1/t17;
    t19 = t13*t18;
    t20 = delta*delta;
    t21 = 1/t20;
    t23 = 0.5-t4;
    t24 = t23*t23;
    t25 = t3*t3;
    t26 = 1.0-t25;
    t35 = t7*t26*t1;
    t38 = t11*t26*t1;
    t41 = t6*t26*t1;
    t43 = 0.15E2*t35-0.3E2*t38+0.15E2*t41;
    t44 = log(t19);
    t48 = 6.0*t8;
    t49 = 15.0*t7;
    t50 = 10.0*t11;
    t55 = -0.1268910394E3*t35+0.2537820788E3*t38-0.1268910394E3*t41;
    t57 = t17*t17;
    t63 = -0.225E2*t35+0.45E2*t38-0.225E2*t41;
    t65 = t55*t18-t13/t57*t63;
    t67 = 1/t13;
    t68 = t67*t17;
    t78 = t26*t1;
    t79 = t13*t13;
    t80 = 1/t79;
    t81 = t80*t57;
    t91 = x*t26*t1;
    t104 = t1*t13;
    t116 = 30.0*t7-60.0*t11+30.0*t6;
    t118 = t19*t44;
    return(t19*(0.2E1*t21*t5*t24*t26-0.2E1*t21*t6*t23*t26+t43*(-1.0+0.15E1*t44)
+0.15E1*(t48-t49+t50)*t65*t68-t43*(-4.0+3.0*t44)+3.0*(1.0-t48+t49-t50)*t65*t68
-0.5*t78*t81+0.1E1*x*t3*t26*t21*t80*t57+0.1E1*t91/t79/t13*t57*t55-0.1E1*t91*t80
*t17*t63)-0.5*t78*(4.0*t104*t18*t5*t24-4.0*t104*t18*t6*t23+t116*(-0.25E1*t19+
0.15E1*t118+1.0)-t116*(-7.0*t19+3.0*t118+0.945940263E1)+0.1E1*t1*t67*t17*t3*t26
+0.5*t65*t26*t81));
  }
}

#include <math.h>
double rhosol(double x)
{
  {
    return((0.113810037E2-0.5075641578E2*pow(0.5*tanh(1/delta*x)+0.5,5.0)+
0.1268910394E3*pow(0.5*tanh(1/delta*x)+0.5,4.0)-0.845940263E2*pow(0.5*tanh(1/
delta*x)+0.5,3.0))/(-0.9E1*pow(0.5*tanh(1/delta*x)+0.5,5.0)+0.225E2*pow(0.5*
tanh(1/delta*x)+0.5,4.0)-0.15E2*pow(0.5*tanh(1/delta*x)+0.5,3.0)+3.0));
  }
}

#include <math.h>
double gradrho(double x)
{
  {
    return((-0.1268910394E3*pow(0.5*tanh(1/delta*x)+0.5,4.0)*(1.0-pow(tanh(1/
delta*x),2.0))/delta+0.2537820788E3*pow(0.5*tanh(1/delta*x)+0.5,3.0)*(1.0-pow(
tanh(1/delta*x),2.0))/delta-0.1268910394E3*pow(0.5*tanh(1/delta*x)+0.5,2.0)*(
1.0-pow(tanh(1/delta*x),2.0))/delta)/(-0.9E1*pow(0.5*tanh(1/delta*x)+0.5,5.0)+
0.225E2*pow(0.5*tanh(1/delta*x)+0.5,4.0)-0.15E2*pow(0.5*tanh(1/delta*x)+0.5,3.0
)+3.0)-(0.113810037E2-0.5075641578E2*pow(0.5*tanh(1/delta*x)+0.5,5.0)+
0.1268910394E3*pow(0.5*tanh(1/delta*x)+0.5,4.0)-0.845940263E2*pow(0.5*tanh(1/
delta*x)+0.5,3.0))/pow(-0.9E1*pow(0.5*tanh(1/delta*x)+0.5,5.0)+0.225E2*pow(0.5*
tanh(1/delta*x)+0.5,4.0)-0.15E2*pow(0.5*tanh(1/delta*x)+0.5,3.0)+3.0,2.0)*(
-0.225E2*pow(0.5*tanh(1/delta*x)+0.5,4.0)*(1.0-pow(tanh(1/delta*x),2.0))/delta+
0.45E2*pow(0.5*tanh(1/delta*x)+0.5,3.0)*(1.0-pow(tanh(1/delta*x),2.0))/delta
-0.225E2*pow(0.5*tanh(1/delta*x)+0.5,2.0)*(1.0-pow(tanh(1/delta*x),2.0))/delta)
);
  }
}

#include <math.h>
double gradphi(double x)
{
  {
    return(0.5*(1.0-pow(tanh(1/delta*x),2.0))/delta);
  }
}

#include <math.h>
double thetasol(double x)
{
  double t1;
  double t11;
  double t13;
  double t14;
  double t18;
  double t19;
  double t21;
  double t22;
  double t3;
  double t33;
  double t34;
  double t36;
  double t37;
  double t4;
  double t48;
  double t49;
  double t5;
  double t54;
  double t57;
  double t6;
  double t60;
  double t64;
  double t7;
  double t74;
  double t8;
  {
    t1 = 1/delta;
    t3 = tanh(t1*x);
    t4 = 0.5*t3;
    t5 = t4+0.5;
    t6 = t5*t5;
    t7 = t6*t6;
    t8 = t7*t5;
    t11 = t6*t5;
    t13 = 0.113810037E2-0.5075641578E2*t8+0.1268910394E3*t7-0.845940263E2*t11;
    t14 = t1*t13;
    t18 = -0.9E1*t8+0.225E2*t7-0.15E2*t11+3.0;
    t19 = 1/t18;
    t21 = 0.5-t4;
    t22 = t21*t21;
    t33 = 30.0*t7-60.0*t11+30.0*t6;
    t34 = t13*t19;
    t36 = log(t34);
    t37 = t34*t36;
    t48 = t3*t3;
    t49 = 1.0-t48;
    t54 = t7*t49*t1;
    t57 = t11*t49*t1;
    t60 = t6*t49*t1;
    t64 = t18*t18;
    t74 = t13*t13;
    return(4.0*t14*t19*t5*t22-4.0*t14*t19*t6*t21+t33*(-0.25E1*t34+0.15E1*t37+
1.0)-t33*(-7.0*t34+3.0*t37+0.945940263E1)+0.1E1*t1/t13*t18*t3*t49+0.5*((
-0.1268910394E3*t54+0.2537820788E3*t57-0.1268910394E3*t60)*t19-t13/t64*(
-0.225E2*t54+0.45E2*t57-0.225E2*t60))*t49/t74*t64);
  }
}

      t0 = phiSourcel;
#include <math.h>
double musol(double x)
{
  {
    return(2.0/delta*pow(0.5*tanh(1/delta*x)+0.5,2.0)*pow(0.5-0.5*tanh(1/delta*
x),2.0)+(6.0*pow(0.5*tanh(1/delta*x)+0.5,5.0)-15.0*pow(0.5*tanh(1/delta*x)+0.5,
4.0)+10.0*pow(0.5*tanh(1/delta*x)+0.5,3.0))*(-1.0+0.15E1*log((0.113810037E2
-0.5075641578E2*pow(0.5*tanh(1/delta*x)+0.5,5.0)+0.1268910394E3*pow(0.5*tanh(1/
delta*x)+0.5,4.0)-0.845940263E2*pow(0.5*tanh(1/delta*x)+0.5,3.0))/(-0.9E1*pow(
0.5*tanh(1/delta*x)+0.5,5.0)+0.225E2*pow(0.5*tanh(1/delta*x)+0.5,4.0)-0.15E2*
pow(0.5*tanh(1/delta*x)+0.5,3.0)+3.0)))+(1.0-6.0*pow(0.5*tanh(1/delta*x)+0.5,
5.0)+15.0*pow(0.5*tanh(1/delta*x)+0.5,4.0)-10.0*pow(0.5*tanh(1/delta*x)+0.5,3.0
))*(-4.0+3.0*log((0.113810037E2-0.5075641578E2*pow(0.5*tanh(1/delta*x)+0.5,5.0)
+0.1268910394E3*pow(0.5*tanh(1/delta*x)+0.5,4.0)-0.845940263E2*pow(0.5*tanh(1/
delta*x)+0.5,3.0))/(-0.9E1*pow(0.5*tanh(1/delta*x)+0.5,5.0)+0.225E2*pow(0.5*
tanh(1/delta*x)+0.5,4.0)-0.15E2*pow(0.5*tanh(1/delta*x)+0.5,3.0)+3.0)))-0.5*x*(
1.0-pow(tanh(1/delta*x),2.0))/delta/pow(0.113810037E2-0.5075641578E2*pow(0.5*
tanh(1/delta*x)+0.5,5.0)+0.1268910394E3*pow(0.5*tanh(1/delta*x)+0.5,4.0)
-0.845940263E2*pow(0.5*tanh(1/delta*x)+0.5,3.0),2.0)*pow(-0.9E1*pow(0.5*tanh(1/
delta*x)+0.5,5.0)+0.225E2*pow(0.5*tanh(1/delta*x)+0.5,4.0)-0.15E2*pow(0.5*tanh(
1/delta*x)+0.5,3.0)+3.0,2.0));
  }
}

      t0 = veloRhs;
#include <math.h>
double rhosol(double x)
{
  {
    return((0.113810037E2-0.5075641578E2*pow(0.5*tanh(1/delta*x)+0.5,5.0)+
0.1268910394E3*pow(0.5*tanh(1/delta*x)+0.5,4.0)-0.845940263E2*pow(0.5*tanh(1/
delta*x)+0.5,3.0))/(-0.9E1*pow(0.5*tanh(1/delta*x)+0.5,5.0)+0.225E2*pow(0.5*
tanh(1/delta*x)+0.5,4.0)-0.15E2*pow(0.5*tanh(1/delta*x)+0.5,3.0)+3.0));
  }
}

#include <math.h>
double gradrho(double x)
{
  {
    return((-0.1268910394E3*pow(0.5*tanh(1/delta*x)+0.5,4.0)*(1.0-pow(tanh(1/
delta*x),2.0))/delta+0.2537820788E3*pow(0.5*tanh(1/delta*x)+0.5,3.0)*(1.0-pow(
tanh(1/delta*x),2.0))/delta-0.1268910394E3*pow(0.5*tanh(1/delta*x)+0.5,2.0)*(
1.0-pow(tanh(1/delta*x),2.0))/delta)/(-0.9E1*pow(0.5*tanh(1/delta*x)+0.5,5.0)+
0.225E2*pow(0.5*tanh(1/delta*x)+0.5,4.0)-0.15E2*pow(0.5*tanh(1/delta*x)+0.5,3.0
)+3.0)-(0.113810037E2-0.5075641578E2*pow(0.5*tanh(1/delta*x)+0.5,5.0)+
0.1268910394E3*pow(0.5*tanh(1/delta*x)+0.5,4.0)-0.845940263E2*pow(0.5*tanh(1/
delta*x)+0.5,3.0))/pow(-0.9E1*pow(0.5*tanh(1/delta*x)+0.5,5.0)+0.225E2*pow(0.5*
tanh(1/delta*x)+0.5,4.0)-0.15E2*pow(0.5*tanh(1/delta*x)+0.5,3.0)+3.0,2.0)*(
-0.225E2*pow(0.5*tanh(1/delta*x)+0.5,4.0)*(1.0-pow(tanh(1/delta*x),2.0))/delta+
0.45E2*pow(0.5*tanh(1/delta*x)+0.5,3.0)*(1.0-pow(tanh(1/delta*x),2.0))/delta
-0.225E2*pow(0.5*tanh(1/delta*x)+0.5,2.0)*(1.0-pow(tanh(1/delta*x),2.0))/delta)
);
  }
}

#include <math.h>
double gradphi(double x)
{
  {
    return(0.5*(1.0-pow(tanh(1/delta*x),2.0))/delta);
  }
}

#include <math.h>
double thetasol(double x)
{
  double t1;
  double t11;
  double t13;
  double t14;
  double t18;
  double t19;
  double t21;
  double t22;
  double t3;
  double t33;
  double t34;
  double t36;
  double t37;
  double t4;
  double t48;
  double t49;
  double t5;
  double t54;
  double t57;
  double t6;
  double t60;
  double t64;
  double t7;
  double t74;
  double t8;
  {
    t1 = 1/delta;
    t3 = tanh(t1*x);
    t4 = 0.5*t3;
    t5 = t4+0.5;
    t6 = t5*t5;
    t7 = t6*t6;
    t8 = t7*t5;
    t11 = t6*t5;
    t13 = 0.113810037E2-0.5075641578E2*t8+0.1268910394E3*t7-0.845940263E2*t11;
    t14 = t1*t13;
    t18 = -0.9E1*t8+0.225E2*t7-0.15E2*t11+3.0;
    t19 = 1/t18;
    t21 = 0.5-t4;
    t22 = t21*t21;
    t33 = 30.0*t7-60.0*t11+30.0*t6;
    t34 = t13*t19;
    t36 = log(t34);
    t37 = t34*t36;
    t48 = t3*t3;
    t49 = 1.0-t48;
    t54 = t7*t49*t1;
    t57 = t11*t49*t1;
    t60 = t6*t49*t1;
    t64 = t18*t18;
    t74 = t13*t13;
    return(4.0*t14*t19*t5*t22-4.0*t14*t19*t6*t21+t33*(-0.25E1*t34+0.15E1*t37+
1.0)-t33*(-7.0*t34+3.0*t37+0.945940263E1)+0.1E1*t1/t13*t18*t3*t49+0.5*((
-0.1268910394E3*t54+0.2537820788E3*t57-0.1268910394E3*t60)*t19-t13/t64*(
-0.225E2*t54+0.45E2*t57-0.225E2*t60))*t49/t74*t64);
  }
}

      t0 = phiSourcel;
#include <math.h>
double musol(double x)
{
  {
    return(2.0/delta*pow(0.5*tanh(1/delta*x)+0.5,2.0)*pow(0.5-0.5*tanh(1/delta*
x),2.0)+(6.0*pow(0.5*tanh(1/delta*x)+0.5,5.0)-15.0*pow(0.5*tanh(1/delta*x)+0.5,
4.0)+10.0*pow(0.5*tanh(1/delta*x)+0.5,3.0))*(-1.0+0.15E1*log((0.113810037E2
-0.5075641578E2*pow(0.5*tanh(1/delta*x)+0.5,5.0)+0.1268910394E3*pow(0.5*tanh(1/
delta*x)+0.5,4.0)-0.845940263E2*pow(0.5*tanh(1/delta*x)+0.5,3.0))/(-0.9E1*pow(
0.5*tanh(1/delta*x)+0.5,5.0)+0.225E2*pow(0.5*tanh(1/delta*x)+0.5,4.0)-0.15E2*
pow(0.5*tanh(1/delta*x)+0.5,3.0)+3.0)))+(1.0-6.0*pow(0.5*tanh(1/delta*x)+0.5,
5.0)+15.0*pow(0.5*tanh(1/delta*x)+0.5,4.0)-10.0*pow(0.5*tanh(1/delta*x)+0.5,3.0
))*(-4.0+3.0*log((0.113810037E2-0.5075641578E2*pow(0.5*tanh(1/delta*x)+0.5,5.0)
+0.1268910394E3*pow(0.5*tanh(1/delta*x)+0.5,4.0)-0.845940263E2*pow(0.5*tanh(1/
delta*x)+0.5,3.0))/(-0.9E1*pow(0.5*tanh(1/delta*x)+0.5,5.0)+0.225E2*pow(0.5*
tanh(1/delta*x)+0.5,4.0)-0.15E2*pow(0.5*tanh(1/delta*x)+0.5,3.0)+3.0)))-0.5*x*(
1.0-pow(tanh(1/delta*x),2.0))/delta/pow(0.113810037E2-0.5075641578E2*pow(0.5*
tanh(1/delta*x)+0.5,5.0)+0.1268910394E3*pow(0.5*tanh(1/delta*x)+0.5,4.0)
-0.845940263E2*pow(0.5*tanh(1/delta*x)+0.5,3.0),2.0)*pow(-0.9E1*pow(0.5*tanh(1/
delta*x)+0.5,5.0)+0.225E2*pow(0.5*tanh(1/delta*x)+0.5,4.0)-0.15E2*pow(0.5*tanh(
1/delta*x)+0.5,3.0)+3.0,2.0));
  }
}

      t0 = veloRhs;
#include <math.h>
double rhosol(double x)
{
  {
    return((0.113810037E2-0.5075641578E2*pow(0.5*tanh(1/delta*x)+0.5,5.0)+
0.1268910394E3*pow(0.5*tanh(1/delta*x)+0.5,4.0)-0.845940263E2*pow(0.5*tanh(1/
delta*x)+0.5,3.0))/(-0.9E1*pow(0.5*tanh(1/delta*x)+0.5,5.0)+0.225E2*pow(0.5*
tanh(1/delta*x)+0.5,4.0)-0.15E2*pow(0.5*tanh(1/delta*x)+0.5,3.0)+3.0));
  }
}

#include <math.h>
double gradrho(double x)
{
  {
    return((-0.1268910394E3*pow(0.5*tanh(1/delta*x)+0.5,4.0)*(1.0-pow(tanh(1/
delta*x),2.0))/delta+0.2537820788E3*pow(0.5*tanh(1/delta*x)+0.5,3.0)*(1.0-pow(
tanh(1/delta*x),2.0))/delta-0.1268910394E3*pow(0.5*tanh(1/delta*x)+0.5,2.0)*(
1.0-pow(tanh(1/delta*x),2.0))/delta)/(-0.9E1*pow(0.5*tanh(1/delta*x)+0.5,5.0)+
0.225E2*pow(0.5*tanh(1/delta*x)+0.5,4.0)-0.15E2*pow(0.5*tanh(1/delta*x)+0.5,3.0
)+3.0)-(0.113810037E2-0.5075641578E2*pow(0.5*tanh(1/delta*x)+0.5,5.0)+
0.1268910394E3*pow(0.5*tanh(1/delta*x)+0.5,4.0)-0.845940263E2*pow(0.5*tanh(1/
delta*x)+0.5,3.0))/pow(-0.9E1*pow(0.5*tanh(1/delta*x)+0.5,5.0)+0.225E2*pow(0.5*
tanh(1/delta*x)+0.5,4.0)-0.15E2*pow(0.5*tanh(1/delta*x)+0.5,3.0)+3.0,2.0)*(
-0.225E2*pow(0.5*tanh(1/delta*x)+0.5,4.0)*(1.0-pow(tanh(1/delta*x),2.0))/delta+
0.45E2*pow(0.5*tanh(1/delta*x)+0.5,3.0)*(1.0-pow(tanh(1/delta*x),2.0))/delta
-0.225E2*pow(0.5*tanh(1/delta*x)+0.5,2.0)*(1.0-pow(tanh(1/delta*x),2.0))/delta)
);
  }
}

#include <math.h>
double gradphi(double x)
{
  {
    return(0.5*(1.0-pow(tanh(1/delta*x),2.0))/delta);
  }
}

#include <math.h>
double thetasol(double x)
{
  double t1;
  double t11;
  double t13;
  double t14;
  double t18;
  double t19;
  double t21;
  double t22;
  double t3;
  double t33;
  double t34;
  double t36;
  double t37;
  double t4;
  double t48;
  double t49;
  double t5;
  double t54;
  double t57;
  double t6;
  double t60;
  double t64;
  double t7;
  double t74;
  double t8;
  {
    t1 = 1/delta;
    t3 = tanh(t1*x);
    t4 = 0.5*t3;
    t5 = t4+0.5;
    t6 = t5*t5;
    t7 = t6*t6;
    t8 = t7*t5;
    t11 = t6*t5;
    t13 = 0.113810037E2-0.5075641578E2*t8+0.1268910394E3*t7-0.845940263E2*t11;
    t14 = t1*t13;
    t18 = -0.9E1*t8+0.225E2*t7-0.15E2*t11+3.0;
    t19 = 1/t18;
    t21 = 0.5-t4;
    t22 = t21*t21;
    t33 = 30.0*t7-60.0*t11+30.0*t6;
    t34 = t13*t19;
    t36 = log(t34);
    t37 = t34*t36;
    t48 = t3*t3;
    t49 = 1.0-t48;
    t54 = t7*t49*t1;
    t57 = t11*t49*t1;
    t60 = t6*t49*t1;
    t64 = t18*t18;
    t74 = t13*t13;
    return(4.0*t14*t19*t5*t22-4.0*t14*t19*t6*t21+t33*(-0.25E1*t34+0.15E1*t37+
1.0)-t33*(-7.0*t34+3.0*t37+0.945940263E1)+0.1E1*t1/t13*t18*t3*t49+0.5*((
-0.1268910394E3*t54+0.2537820788E3*t57-0.1268910394E3*t60)*t19-t13/t64*(
-0.225E2*t54+0.45E2*t57-0.225E2*t60))*t49/t74*t64);
  }
}

      t0 = phiSourcel;
#include <math.h>
double musol(double x)
{
  {
    return(2.0/delta*pow(0.5*tanh(1/delta*x)+0.5,2.0)*pow(0.5-0.5*tanh(1/delta*
x),2.0)+(6.0*pow(0.5*tanh(1/delta*x)+0.5,5.0)-15.0*pow(0.5*tanh(1/delta*x)+0.5,
4.0)+10.0*pow(0.5*tanh(1/delta*x)+0.5,3.0))*(-1.0+0.15E1*log((0.113810037E2
-0.5075641578E2*pow(0.5*tanh(1/delta*x)+0.5,5.0)+0.1268910394E3*pow(0.5*tanh(1/
delta*x)+0.5,4.0)-0.845940263E2*pow(0.5*tanh(1/delta*x)+0.5,3.0))/(-0.9E1*pow(
0.5*tanh(1/delta*x)+0.5,5.0)+0.225E2*pow(0.5*tanh(1/delta*x)+0.5,4.0)-0.15E2*
pow(0.5*tanh(1/delta*x)+0.5,3.0)+3.0)))+(1.0-6.0*pow(0.5*tanh(1/delta*x)+0.5,
5.0)+15.0*pow(0.5*tanh(1/delta*x)+0.5,4.0)-10.0*pow(0.5*tanh(1/delta*x)+0.5,3.0
))*(-4.0+3.0*log((0.113810037E2-0.5075641578E2*pow(0.5*tanh(1/delta*x)+0.5,5.0)
+0.1268910394E3*pow(0.5*tanh(1/delta*x)+0.5,4.0)-0.845940263E2*pow(0.5*tanh(1/
delta*x)+0.5,3.0))/(-0.9E1*pow(0.5*tanh(1/delta*x)+0.5,5.0)+0.225E2*pow(0.5*
tanh(1/delta*x)+0.5,4.0)-0.15E2*pow(0.5*tanh(1/delta*x)+0.5,3.0)+3.0)))-0.5*x*(
1.0-pow(tanh(1/delta*x),2.0))/delta/pow(0.113810037E2-0.5075641578E2*pow(0.5*
tanh(1/delta*x)+0.5,5.0)+0.1268910394E3*pow(0.5*tanh(1/delta*x)+0.5,4.0)
-0.845940263E2*pow(0.5*tanh(1/delta*x)+0.5,3.0),2.0)*pow(-0.9E1*pow(0.5*tanh(1/
delta*x)+0.5,5.0)+0.225E2*pow(0.5*tanh(1/delta*x)+0.5,4.0)-0.15E2*pow(0.5*tanh(
1/delta*x)+0.5,3.0)+3.0,2.0));
  }
}

      t0 = veloRhs;
#include <math.h>
double rhosol(double x)
{
  {
    return((0.113810037E2-0.5075641578E2*pow(0.5+0.5*tanh(1/delta*x),5.0)+
0.1268910394E3*pow(0.5+0.5*tanh(1/delta*x),4.0)-0.845940263E2*pow(0.5+0.5*tanh(
1/delta*x),3.0))/(3.0-0.9E1*pow(0.5+0.5*tanh(1/delta*x),5.0)+0.225E2*pow(0.5+
0.5*tanh(1/delta*x),4.0)-0.15E2*pow(0.5+0.5*tanh(1/delta*x),3.0)));
  }
}

#include <math.h>
double gradrho(double x)
{
  {
    return((-0.1268910394E3*pow(0.5+0.5*tanh(1/delta*x),4.0)*(1.0-pow(tanh(1/
delta*x),2.0))/delta+0.2537820788E3*pow(0.5+0.5*tanh(1/delta*x),3.0)*(1.0-pow(
tanh(1/delta*x),2.0))/delta-0.1268910394E3*pow(0.5+0.5*tanh(1/delta*x),2.0)*(
1.0-pow(tanh(1/delta*x),2.0))/delta)/(3.0-0.9E1*pow(0.5+0.5*tanh(1/delta*x),5.0
)+0.225E2*pow(0.5+0.5*tanh(1/delta*x),4.0)-0.15E2*pow(0.5+0.5*tanh(1/delta*x),
3.0))-(0.113810037E2-0.5075641578E2*pow(0.5+0.5*tanh(1/delta*x),5.0)+
0.1268910394E3*pow(0.5+0.5*tanh(1/delta*x),4.0)-0.845940263E2*pow(0.5+0.5*tanh(
1/delta*x),3.0))/pow(3.0-0.9E1*pow(0.5+0.5*tanh(1/delta*x),5.0)+0.225E2*pow(0.5
+0.5*tanh(1/delta*x),4.0)-0.15E2*pow(0.5+0.5*tanh(1/delta*x),3.0),2.0)*(
-0.225E2*pow(0.5+0.5*tanh(1/delta*x),4.0)*(1.0-pow(tanh(1/delta*x),2.0))/delta+
0.45E2*pow(0.5+0.5*tanh(1/delta*x),3.0)*(1.0-pow(tanh(1/delta*x),2.0))/delta
-0.225E2*pow(0.5+0.5*tanh(1/delta*x),2.0)*(1.0-pow(tanh(1/delta*x),2.0))/delta)
);
  }
}

#include <math.h>
double gradphi(double x)
{
  {
    return(0.5*(1.0-pow(tanh(1/delta*x),2.0))/delta);
  }
}

#include <math.h>
double thetasol(double x)
{
  double t1;
  double t11;
  double t13;
  double t14;
  double t18;
  double t19;
  double t21;
  double t22;
  double t3;
  double t33;
  double t34;
  double t36;
  double t37;
  double t4;
  double t48;
  double t49;
  double t5;
  double t54;
  double t57;
  double t6;
  double t60;
  double t64;
  double t7;
  double t74;
  double t8;
  {
    t1 = 1/delta;
    t3 = tanh(t1*x);
    t4 = 0.5*t3;
    t5 = 0.5+t4;
    t6 = t5*t5;
    t7 = t6*t6;
    t8 = t7*t5;
    t11 = t6*t5;
    t13 = 0.113810037E2-0.5075641578E2*t8+0.1268910394E3*t7-0.845940263E2*t11;
    t14 = t1*t13;
    t18 = 3.0-0.9E1*t8+0.225E2*t7-0.15E2*t11;
    t19 = 1/t18;
    t21 = 0.5-t4;
    t22 = t21*t21;
    t33 = 30.0*t7-60.0*t11+30.0*t6;
    t34 = t13*t19;
    t36 = log(t34);
    t37 = t34*t36;
    t48 = t3*t3;
    t49 = 1.0-t48;
    t54 = t7*t49*t1;
    t57 = t11*t49*t1;
    t60 = t6*t49*t1;
    t64 = t18*t18;
    t74 = t13*t13;
    return(4.0*t14*t19*t5*t22-4.0*t14*t19*t6*t21+t33*(-0.25E1*t34+0.15E1*t37+
1.0)-t33*(-7.0*t34+3.0*t37+0.945940263E1)+0.1E1*t1/t13*t18*t3*t49+0.5*((
-0.1268910394E3*t54+0.2537820788E3*t57-0.1268910394E3*t60)*t19-t13/t64*(
-0.225E2*t54+0.45E2*t57-0.225E2*t60))*t49/t74*t64);
  }
}

      t0 = phiSourcel;
#include <math.h>
double musol(double x)
{
  {
    return(2.0/delta*pow(0.5+0.5*tanh(1/delta*x),2.0)*pow(0.5-0.5*tanh(1/delta*
x),2.0)+(6.0*pow(0.5+0.5*tanh(1/delta*x),5.0)-15.0*pow(0.5+0.5*tanh(1/delta*x),
4.0)+10.0*pow(0.5+0.5*tanh(1/delta*x),3.0))*(-1.0+0.15E1*log((0.113810037E2
-0.5075641578E2*pow(0.5+0.5*tanh(1/delta*x),5.0)+0.1268910394E3*pow(0.5+0.5*
tanh(1/delta*x),4.0)-0.845940263E2*pow(0.5+0.5*tanh(1/delta*x),3.0))/(3.0-0.9E1
*pow(0.5+0.5*tanh(1/delta*x),5.0)+0.225E2*pow(0.5+0.5*tanh(1/delta*x),4.0)
-0.15E2*pow(0.5+0.5*tanh(1/delta*x),3.0))))+(1.0-6.0*pow(0.5+0.5*tanh(1/delta*x
),5.0)+15.0*pow(0.5+0.5*tanh(1/delta*x),4.0)-10.0*pow(0.5+0.5*tanh(1/delta*x),
3.0))*(-4.0+3.0*log((0.113810037E2-0.5075641578E2*pow(0.5+0.5*tanh(1/delta*x),
5.0)+0.1268910394E3*pow(0.5+0.5*tanh(1/delta*x),4.0)-0.845940263E2*pow(0.5+0.5*
tanh(1/delta*x),3.0))/(3.0-0.9E1*pow(0.5+0.5*tanh(1/delta*x),5.0)+0.225E2*pow(
0.5+0.5*tanh(1/delta*x),4.0)-0.15E2*pow(0.5+0.5*tanh(1/delta*x),3.0))))-0.5*x*(
1.0-pow(tanh(1/delta*x),2.0))/delta/pow(0.113810037E2-0.5075641578E2*pow(0.5+
0.5*tanh(1/delta*x),5.0)+0.1268910394E3*pow(0.5+0.5*tanh(1/delta*x),4.0)
-0.845940263E2*pow(0.5+0.5*tanh(1/delta*x),3.0),2.0)*pow(3.0-0.9E1*pow(0.5+0.5*
tanh(1/delta*x),5.0)+0.225E2*pow(0.5+0.5*tanh(1/delta*x),4.0)-0.15E2*pow(0.5+
0.5*tanh(1/delta*x),3.0),2.0));
  }
}

      t0 = veloRhs;

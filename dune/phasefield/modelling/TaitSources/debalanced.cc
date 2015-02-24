#include <math.h>
double rhosol(double x)
{
  double t3;
  {
    t3 = tanh(x/delta);
    return(0.247901727*t3+0.3544783819);
  }
}

#include <math.h>
double gradrho(double x)
{
  double t1;
  double t3;
  double t4;
  {
    t1 = 1/delta;
    t3 = tanh(x*t1);
    t4 = t3*t3;
    return(0.247901727*(1.0-t4)*t1);
  }
}

#include <math.h>
double gradphi(double x)
{
  {
    return(0.5*(1.0-pow(tanh(x/delta),2.0))/delta);
  }
}

#include <math.h>
double thetasol(double x)
{
  double t1;
  double t21;
  double t25;
  double t27;
  double t28;
  double t30;
  double t36;
  double t4;
  double t43;
  double t6;
  double t7;
  double t8;
  {
    t1 = 1/delta;
    t4 = tanh(x*t1);
    t6 = 0.5*t4+0.5;
    t7 = t6*t6;
    t8 = t7*t6;
    t21 = t7*t7;
    t25 = 30.0*t21-60.0*t8+30.0*t7;
    t27 = 0.247901727*t4+0.3544783819;
    t28 = t27*t27;
    t30 = t28*t28;
    t36 = log(t27);
    t43 = t4*t4;
    return(2.0*t1*A*(4.0*t8+3.0*(2.0*alpha-2.0)*t7+2.0*(-3.0*alpha+1.0)*t6)+t25
*(0.35*t30*t28*t27-0.2901812679E-1*t4+0.189450439E-1)-t25*(0.1753186565*t27*t36
+0.5384442083E-1*t4+0.9567781742E-1)+0.1E1*t1*t4*(1.0-t43));
  }
}

#include <math.h>
double phiSource(double x)
{
  double t1;
  double t11;
  double t12;
  double t13;
  double t26;
  double t3;
  double t30;
  double t32;
  double t33;
  double t35;
  double t41;
  double t48;
  double t5;
  double t9;
  {
    t1 = 1/delta;
    t3 = tanh(x*t1);
    t5 = 0.247901727*t3+0.3544783819;
    t9 = tanh(t5*t1);
    t11 = 0.5*t9+0.5;
    t12 = t11*t11;
    t13 = t12*t11;
    t26 = t12*t12;
    t30 = 30.0*t26-60.0*t13+30.0*t12;
    t32 = 0.247901727*t9+0.3544783819;
    t33 = t32*t32;
    t35 = t33*t33;
    t41 = log(t32);
    t48 = t9*t9;
    return(1/t5*(2.0*t1*A*(4.0*t13+3.0*(2.0*alpha-2.0)*t12+2.0*(-3.0*alpha+1.0)
*t11)+t30*(0.35*t35*t33*t32-0.2901812679E-1*t9+0.189450439E-1)-t30*(
0.1753186565*t32*t41+0.5384442083E-1*t9+0.9567781742E-1)+0.1E1*t1*t9*(1.0-t48))
);
  }
}

#include <math.h>
double musol(double x)
{
  double t10;
  double t11;
  double t12;
  double t13;
  double t20;
  double t24;
  double t3;
  double t5;
  double t6;
  double t7;
  double t8;
  {
    t3 = tanh(x/delta);
    t5 = 0.5*t3+0.5;
    t6 = t5*t5;
    t7 = t6*t6;
    t8 = t7*t5;
    t10 = 0.247901727*t3+0.3544783819;
    t11 = t10*t10;
    t12 = t11*t11;
    t13 = t12*t11;
    t20 = t6*t5;
    t24 = log(t10);
    return(0.147E2*t8*t13-0.3057445711E1*t8-0.3675E2*t7*t13+0.7643614278E1*t7+
0.245E2*t20*t13-0.5095742852E1*t20+0.1753186565*t24+0.392519325-0.1051911939E1*
t8*t24+0.2629779848E1*t7*t24-0.1753186565E1*t20*t24);
  }
}

#include <math.h>
double veloSource(double x)
{
  double t1;
  double t10;
  double t101;
  double t11;
  double t24;
  double t28;
  double t3;
  double t30;
  double t31;
  double t33;
  double t39;
  double t4;
  double t5;
  double t51;
  double t55;
  double t56;
  double t6;
  double t60;
  double t69;
  double t81;
  double t88;
  double t9;
  {
    t1 = 1/delta;
    t3 = tanh(x*t1);
    t4 = t3*t3;
    t5 = 1.0-t4;
    t6 = t5*t1;
    t9 = 0.5*t3+0.5;
    t10 = t9*t9;
    t11 = t10*t9;
    t24 = t10*t10;
    t28 = 30.0*t24-60.0*t11+30.0*t10;
    t30 = 0.247901727*t3+0.3544783819;
    t31 = t30*t30;
    t33 = t31*t31;
    t39 = log(t30);
    t51 = t33*t31;
    t55 = t24*t9;
    t56 = t33*t30;
    t60 = t24*t5;
    t69 = t11*t5;
    t81 = 1/t30;
    t88 = t1*t81;
    t101 = 0.3675E2*t24*t51*t6+0.2186493232E2*t55*t56*t6-0.7643614278E1*t60*t1
-0.735E2*t11*t51*t6-0.546623308E2*t24*t56*t6+0.1528722856E2*t69*t1+0.3675E2*t10
*t51*t6+0.3644155387E2*t11*t56*t6-0.7643614278E1*t10*t5*t1+0.4346179772E-1*t6*
t81-0.2629779848E1*t24*t39*t6-0.2607707863*t55*t5*t88+0.5259559696E1*t11*t39*t6
+0.6519269659*t60*t88-0.2629779848E1*t10*t39*t6-0.4346179772*t69*t88;
    return(-0.5*t6*(2.0*t1*A*(4.0*t11+3.0*(2.0*alpha-2.0)*t10+2.0*(-3.0*alpha+
1.0)*t9)+t28*(0.35*t33*t31*t30-0.2901812679E-1*t3+0.189450439E-1)-t28*(
0.1753186565*t30*t39+0.5384442083E-1*t3+0.9567781742E-1)+0.1E1*t1*t3*t5)+t30*
t101);
  }
}

#include <math.h>
double rhosol(double x)
{
  double t3;
  {
    t3 = tanh(x/delta);
    return(0.247901727*t3+0.3544783819);
  }
}

#include <math.h>
double gradrho(double x)
{
  double t1;
  double t3;
  double t4;
  {
    t1 = 1/delta;
    t3 = tanh(x*t1);
    t4 = t3*t3;
    return(0.247901727*(1.0-t4)*t1);
  }
}

#include <math.h>
double gradphi(double x)
{
  {
    return(0.5*(1.0-pow(tanh(x/delta),2.0))/delta);
  }
}

#include <math.h>
double thetasol(double x)
{
  double t1;
  double t21;
  double t25;
  double t27;
  double t28;
  double t30;
  double t36;
  double t4;
  double t43;
  double t6;
  double t7;
  double t8;
  {
    t1 = 1/delta;
    t4 = tanh(x*t1);
    t6 = 0.5*t4+0.5;
    t7 = t6*t6;
    t8 = t7*t6;
    t21 = t7*t7;
    t25 = 30.0*t21-60.0*t8+30.0*t7;
    t27 = 0.3544783819+0.247901727*t4;
    t28 = t27*t27;
    t30 = t28*t28;
    t36 = log(t27);
    t43 = t4*t4;
    return(2.0*t1*A*(4.0*t8+3.0*(2.0*alpha-2.0)*t7+2.0*(-3.0*alpha+1.0)*t6)+t25
*(0.464199863E-2+0.8575855211E-1*t30*t28*t27-0.7110150146E-2*t4)-t25*(
0.5197411266E1+0.9523661641E1*t27*t36+0.2924937115E1*t4)+0.1E1*t1*t4*(1.0-t43)
);
  }
}

#include <math.h>
double phiSource(double x)
{
  double t1;
  double t11;
  double t12;
  double t13;
  double t26;
  double t3;
  double t30;
  double t32;
  double t33;
  double t35;
  double t41;
  double t48;
  double t5;
  double t9;
  {
    t1 = 1/delta;
    t3 = tanh(x*t1);
    t5 = 0.247901727*t3+0.3544783819;
    t9 = tanh(t5*t1);
    t11 = 0.5*t9+0.5;
    t12 = t11*t11;
    t13 = t12*t11;
    t26 = t12*t12;
    t30 = 30.0*t26-60.0*t13+30.0*t12;
    t32 = 0.247901727*t9+0.3544783819;
    t33 = t32*t32;
    t35 = t33*t33;
    t41 = log(t32);
    t48 = t9*t9;
    return(1/t5*(2.0*t1*A*(4.0*t13+3.0*(2.0*alpha-2.0)*t12+2.0*(-3.0*alpha+1.0)
*t11)+t30*(0.8575855211E-1*t35*t33*t32-0.7110150146E-2*t9+0.464199863E-2)-t30*(
0.5197411266E1+0.9523661641E1*t32*t41+0.2924937115E1*t9)+0.1E1*t1*t9*(1.0-t48))
);
  }
}

#include <math.h>
double musol(double x)
{
  double t10;
  double t11;
  double t12;
  double t13;
  double t20;
  double t24;
  double t3;
  double t5;
  double t6;
  double t7;
  double t8;
  {
    t3 = tanh(x/delta);
    t5 = 0.5*t3+0.5;
    t6 = t5*t5;
    t7 = t6*t6;
    t8 = t7*t5;
    t10 = 0.247901727*t3+0.3544783819;
    t11 = t10*t10;
    t12 = t11*t11;
    t13 = t12*t11;
    t20 = t6*t5;
    t24 = log(t10);
    return(0.3601859189E1*t8*t13-0.1281067179E3*t8-0.9004647972E1*t7*t13+
0.3202667947E3*t7+0.6003098648E1*t20*t13-0.2135111965E3*t20+0.2132243832E2+
0.9523661641E1*t24-0.5714196985E2*t8*t24+0.1428549246E3*t7*t24-0.9523661641E2*
t20*t24);
  }
}

#include <math.h>
double veloSource(double x)
{
  double t1;
  double t10;
  double t101;
  double t11;
  double t24;
  double t28;
  double t3;
  double t30;
  double t31;
  double t33;
  double t39;
  double t4;
  double t5;
  double t51;
  double t55;
  double t56;
  double t6;
  double t60;
  double t69;
  double t81;
  double t88;
  double t9;
  {
    t1 = 1/delta;
    t3 = tanh(x*t1);
    t4 = t3*t3;
    t5 = 1.0-t4;
    t6 = t5*t1;
    t9 = 0.5*t3+0.5;
    t10 = t9*t9;
    t11 = t10*t9;
    t24 = t10*t10;
    t28 = 30.0*t24-60.0*t11+30.0*t10;
    t30 = 0.247901727*t3+0.3544783819;
    t31 = t30*t30;
    t33 = t31*t31;
    t39 = log(t30);
    t51 = t33*t31;
    t55 = t24*t9;
    t56 = t33*t30;
    t60 = t24*t5;
    t69 = t11*t5;
    t81 = 1/t30;
    t88 = t1*t81;
    t101 = 0.9004647972E1*t24*t51*t6+0.535744268E1*t55*t56*t6-0.3202667948E3*
t60*t1-0.1800929594E2*t11*t51*t6-0.133936067E2*t24*t56*t6+0.6405335894E3*t69*t1
+0.9004647972E1*t10*t51*t6+0.8929071133E1*t11*t56*t6-0.3202667948E3*t10*t5*t1+
0.2360932168E1*t6*t81-0.1428549246E3*t24*t39*t6-0.1416559301E2*t55*t5*t88+
0.2857098492E3*t11*t39*t6+0.3541398252E2*t60*t88-0.1428549246E3*t10*t39*t6
-0.2360932168E2*t69*t88;
    return(-0.5*t6*(2.0*t1*A*(4.0*t11+3.0*(2.0*alpha-2.0)*t10+2.0*(-3.0*alpha+
1.0)*t9)+t28*(0.8575855211E-1*t33*t31*t30-0.7110150146E-2*t3+0.464199863E-2)-
t28*(0.5197411266E1+0.9523661641E1*t30*t39+0.2924937115E1*t3)+0.1E1*t1*t3*t5)+
t30*t101);
  }
}

#include <math.h>
double rhosol(double x)
{
  double t3;
  {
    t3 = tanh(x/delta);
    return(0.247901727*t3+0.3544783819);
  }
}

#include <math.h>
double gradrho(double x)
{
  double t1;
  double t3;
  double t4;
  {
    t1 = 1/delta;
    t3 = tanh(x*t1);
    t4 = t3*t3;
    return(0.247901727*(1.0-t4)*t1);
  }
}

#include <math.h>
double gradphi(double x)
{
  {
    return(0.5*(1.0-pow(tanh(x/delta),2.0))/delta);
  }
}

#include <math.h>
double thetasol(double x)
{
  double t1;
  double t21;
  double t25;
  double t27;
  double t28;
  double t30;
  double t36;
  double t4;
  double t43;
  double t6;
  double t7;
  double t8;
  {
    t1 = 1/delta;
    t4 = tanh(x*t1);
    t6 = 0.5*t4+0.5;
    t7 = t6*t6;
    t8 = t7*t6;
    t21 = t7*t7;
    t25 = 30.0*t21-60.0*t8+30.0*t7;
    t27 = 0.247901727*t4+0.3544783819;
    t28 = t27*t27;
    t30 = t28*t28;
    t36 = log(t27);
    t43 = t4*t4;
    return(2.0*t1*A*(4.0*t8+3.0*(2.0*alpha-2.0)*t7+2.0*(-3.0*alpha+1.0)*t6)+t25
*(0.8575855211E-1*t30*t28*t27-0.7110150146E-2*t4+0.464199863E-2)-t25*(
0.7681415425E-1+0.1407531515*t27*t36+0.4322855348E-1*t4)+0.1E1*t1*t4*(1.0-t43)
);
  }
}

#include <math.h>
double phiSource(double x)
{
  double t1;
  double t11;
  double t12;
  double t13;
  double t26;
  double t3;
  double t30;
  double t32;
  double t33;
  double t35;
  double t41;
  double t48;
  double t5;
  double t9;
  {
    t1 = 1/delta;
    t3 = tanh(x*t1);
    t5 = 0.247901727*t3+0.3544783819;
    t9 = tanh(t5*t1);
    t11 = 0.5*t9+0.5;
    t12 = t11*t11;
    t13 = t12*t11;
    t26 = t12*t12;
    t30 = 30.0*t26-60.0*t13+30.0*t12;
    t32 = 0.247901727*t9+0.3544783819;
    t33 = t32*t32;
    t35 = t33*t33;
    t41 = log(t32);
    t48 = t9*t9;
    return(1/t5*(2.0*t1*A*(4.0*t13+3.0*(2.0*alpha-2.0)*t12+2.0*(-3.0*alpha+1.0)
*t11)+t30*(0.8575855211E-1*t35*t33*t32-0.7110150146E-2*t9+0.464199863E-2)-t30*(
0.7681415425E-1+0.1407531515*t32*t41+0.4322855348E-1*t9)+0.1E1*t1*t9*(1.0-t48))
);
  }
}

#include <math.h>
double musol(double x)
{
  double t10;
  double t11;
  double t12;
  double t13;
  double t20;
  double t24;
  double t3;
  double t5;
  double t6;
  double t7;
  double t8;
  {
    t3 = tanh(x/delta);
    t5 = 0.5*t3+0.5;
    t6 = t5*t5;
    t7 = t6*t6;
    t8 = t7*t5;
    t10 = 0.247901727*t3+0.3544783819;
    t11 = t10*t10;
    t12 = t11*t11;
    t13 = t12*t11;
    t20 = t6*t5;
    t24 = log(t10);
    return(0.3601859189E1*t8*t13-0.2062873559E1*t8-0.9004647972E1*t7*t13+
0.5157183897E1*t7+0.6003098648E1*t20*t13-0.3438122598E1*t20+0.3151309342+
0.1407531515*t24-0.844518909*t8*t24+0.2111297272E1*t7*t24-0.1407531515E1*t20*
t24);
  }
}

#include <math.h>
double veloSource(double x)
{
  double t1;
  double t10;
  double t101;
  double t11;
  double t24;
  double t28;
  double t3;
  double t30;
  double t31;
  double t33;
  double t39;
  double t4;
  double t5;
  double t51;
  double t55;
  double t56;
  double t6;
  double t60;
  double t69;
  double t81;
  double t88;
  double t9;
  {
    t1 = 1/delta;
    t3 = tanh(x*t1);
    t4 = t3*t3;
    t5 = 1.0-t4;
    t6 = t5*t1;
    t9 = 0.5*t3+0.5;
    t10 = t9*t9;
    t11 = t10*t9;
    t24 = t10*t10;
    t28 = 30.0*t24-60.0*t11+30.0*t10;
    t30 = 0.247901727*t3+0.3544783819;
    t31 = t30*t30;
    t33 = t31*t31;
    t39 = log(t30);
    t51 = t33*t31;
    t55 = t24*t9;
    t56 = t33*t30;
    t60 = t24*t5;
    t69 = t11*t5;
    t81 = 1/t30;
    t88 = t1*t81;
    t101 = 0.9004647972E1*t24*t51*t6+0.535744268E1*t55*t56*t6-0.5157183898E1*
t60*t1-0.1800929594E2*t11*t51*t6-0.133936067E2*t24*t56*t6+0.1031436779E2*t69*t1
+0.9004647972E1*t10*t51*t6+0.8929071133E1*t11*t56*t6-0.5157183897E1*t10*t5*t1+
0.3489294934E-1*t6*t81-0.2111297272E1*t24*t39*t6-0.209357696*t55*t5*t88+
0.4222594544E1*t11*t39*t6+0.5233942399*t60*t88-0.2111297272E1*t10*t39*t6
-0.3489294934*t69*t88;
    return(-0.5*t6*(2.0*t1*A*(4.0*t11+3.0*(2.0*alpha-2.0)*t10+2.0*(-3.0*alpha+
1.0)*t9)+t28*(0.8575855211E-1*t33*t31*t30-0.7110150146E-2*t3+0.464199863E-2)-
t28*(0.7681415425E-1+0.1407531515*t30*t39+0.4322855348E-1*t3)+0.1E1*t1*t3*t5)+
t30*t101);
  }
}

#include <math.h>
double rhosol(double x)
{
  double t3;
  {
    t3 = tanh(x/delta);
    return(0.3544783819+0.247901727*t3);
  }
}

#include <math.h>
double gradrho(double x)
{
  double t1;
  double t3;
  double t4;
  {
    t1 = 1/delta;
    t3 = tanh(x*t1);
    t4 = t3*t3;
    return(0.247901727*(1.0-t4)*t1);
  }
}

#include <math.h>
double gradphi(double x)
{
  {
    return(0.5*(1.0-pow(tanh(x/delta),2.0))/delta);
  }
}

#include <math.h>
double thetasol(double x)
{
  double t1;
  double t21;
  double t25;
  double t27;
  double t28;
  double t30;
  double t36;
  double t4;
  double t43;
  double t6;
  double t7;
  double t8;
  {
    t1 = 1/delta;
    t4 = tanh(x*t1);
    t6 = 0.5*t4+0.5;
    t7 = t6*t6;
    t8 = t7*t6;
    t21 = t7*t7;
    t25 = 30.0*t21-60.0*t8+30.0*t7;
    t27 = 0.247901727*t4+0.3544783819;
    t28 = t27*t27;
    t30 = t28*t28;
    t36 = log(t27);
    t43 = t4*t4;
    return(2.0*t1*A*(4.0*t8+3.0*(2.0*alpha-2.0)*t7+2.0*(-3.0*alpha+1.0)*t6)+t25
*(0.8575855211E-1*t30*t28*t27-0.7110150146E-2*t4+0.464199863E-2)-t25*(
0.7680908484E-1+0.1407438624*t27*t36+0.4322570058E-1*t4)+0.1E1*t1*t4*(1.0-t43)
);
  }
}

#include <math.h>
double phiSource(double x)
{
  double t1;
  double t11;
  double t12;
  double t13;
  double t26;
  double t3;
  double t30;
  double t32;
  double t33;
  double t35;
  double t41;
  double t48;
  double t5;
  double t9;
  {
    t1 = 1/delta;
    t3 = tanh(x*t1);
    t5 = 0.3544783819+0.247901727*t3;
    t9 = tanh(t5*t1);
    t11 = 0.5*t9+0.5;
    t12 = t11*t11;
    t13 = t12*t11;
    t26 = t12*t12;
    t30 = 30.0*t26-60.0*t13+30.0*t12;
    t32 = 0.247901727*t9+0.3544783819;
    t33 = t32*t32;
    t35 = t33*t33;
    t41 = log(t32);
    t48 = t9*t9;
    return(1/t5*(2.0*t1*A*(4.0*t13+3.0*(2.0*alpha-2.0)*t12+2.0*(-3.0*alpha+1.0)
*t11)+t30*(0.8575855211E-1*t35*t33*t32-0.7110150146E-2*t9+0.464199863E-2)-t30*(
0.7680908484E-1+0.1407438624*t32*t41+0.4322570058E-1*t9)+0.1E1*t1*t9*(1.0-t48))
);
  }
}

#include <math.h>
double musol(double x)
{
  double t10;
  double t11;
  double t12;
  double t13;
  double t20;
  double t24;
  double t3;
  double t5;
  double t6;
  double t7;
  double t8;
  {
    t3 = tanh(x/delta);
    t5 = 0.5*t3+0.5;
    t6 = t5*t5;
    t7 = t6*t6;
    t8 = t7*t5;
    t10 = 0.3544783819+0.247901727*t3;
    t11 = t10*t10;
    t12 = t11*t11;
    t13 = t12*t11;
    t20 = t6*t5;
    t24 = log(t10);
    return(0.3601859189E1*t8*t13-0.2062748775E1*t8-0.9004647972E1*t7*t13+
0.5156871937E1*t7+0.6003098648E1*t20*t13-0.3437914625E1*t20+0.3151101369+
0.1407438624*t24-0.8444631744*t8*t24+0.2111157936E1*t7*t24-0.1407438624E1*t20*
t24);
  }
}

#include <math.h>
double veloSource(double x)
{
  double t1;
  double t10;
  double t101;
  double t11;
  double t24;
  double t28;
  double t3;
  double t30;
  double t31;
  double t33;
  double t39;
  double t4;
  double t5;
  double t51;
  double t55;
  double t56;
  double t6;
  double t60;
  double t69;
  double t81;
  double t88;
  double t9;
  {
    t1 = 1/delta;
    t3 = tanh(x*t1);
    t4 = t3*t3;
    t5 = 1.0-t4;
    t6 = t5*t1;
    t9 = 0.5*t3+0.5;
    t10 = t9*t9;
    t11 = t10*t9;
    t24 = t10*t10;
    t28 = 30.0*t24-60.0*t11+30.0*t10;
    t30 = 0.3544783819+0.247901727*t3;
    t31 = t30*t30;
    t33 = t31*t31;
    t39 = log(t30);
    t51 = t33*t31;
    t55 = t24*t9;
    t56 = t33*t30;
    t60 = t24*t5;
    t69 = t11*t5;
    t81 = 1/t30;
    t88 = t1*t81;
    t101 = 0.9004647972E1*t24*t51*t6+0.535744268E1*t55*t56*t6-0.5156871938E1*
t60*t1-0.1800929594E2*t11*t51*t6-0.133936067E2*t24*t56*t6+0.1031374387E2*t69*t1
+0.9004647972E1*t10*t51*t6+0.8929071133E1*t11*t56*t6-0.5156871938E1*t10*t5*t1+
0.3489064655E-1*t6*t81-0.2111157936E1*t24*t39*t6-0.2093438793*t55*t5*t88+
0.4222315872E1*t11*t39*t6+0.5233596983*t60*t88-0.2111157936E1*t10*t39*t6
-0.3489064655*t69*t88;
    return(-0.5*t6*(2.0*t1*A*(4.0*t11+3.0*(2.0*alpha-2.0)*t10+2.0*(-3.0*alpha+
1.0)*t9)+t28*(0.8575855211E-1*t33*t31*t30-0.7110150146E-2*t3+0.464199863E-2)-
t28*(0.7680908484E-1+0.1407438624*t30*t39+0.4322570058E-1*t3)+0.1E1*t1*t3*t5)+
t30*t101);
  }
}

#include <math.h>
double rhosol(double x)
{
  double t3;
  {
    t3 = tanh(x/delta);
    return(0.247901727*t3+0.3544783819);
  }
}

#include <math.h>
double gradrho(double x)
{
  double t1;
  double t3;
  double t4;
  {
    t1 = 1/delta;
    t3 = tanh(x*t1);
    t4 = t3*t3;
    return(0.247901727*(1.0-t4)*t1);
  }
}

#include <math.h>
double gradphi(double x)
{
  {
    return(0.5*(1.0-pow(tanh(x/delta),2.0))/delta);
  }
}

#include <math.h>
double thetasol(double x)
{
  double t1;
  double t21;
  double t25;
  double t27;
  double t28;
  double t30;
  double t36;
  double t4;
  double t43;
  double t6;
  double t7;
  double t8;
  {
    t1 = 1/delta;
    t4 = tanh(x*t1);
    t6 = 0.5*t4+0.5;
    t7 = t6*t6;
    t8 = t7*t6;
    t21 = t7*t7;
    t25 = 30.0*t21-60.0*t8+30.0*t7;
    t27 = 0.247901727*t4+0.3544783819;
    t28 = t27*t27;
    t30 = t28*t28;
    t36 = log(t27);
    t43 = t4*t4;
    return(2.0*t1*A*(4.0*t8+3.0*(2.0*alpha-2.0)*t7+2.0*(-3.0*alpha+1.0)*t6)+t25
*(0.8575855211E-1*t30*t28*t27-0.7110150146E-2*t4+0.464199863E-2)-t25*(
0.5141205242E-7+0.9429832461E-7*t27*t36+0.289612E-7*t4)+0.1E1*t1*t4*(1.0-t43));
  }
}

#include <math.h>
double phiSource(double x)
{
  double t1;
  double t11;
  double t12;
  double t13;
  double t26;
  double t3;
  double t30;
  double t32;
  double t33;
  double t35;
  double t41;
  double t48;
  double t5;
  double t9;
  {
    t1 = 1/delta;
    t3 = tanh(x*t1);
    t5 = 0.247901727*t3+0.3544783819;
    t9 = tanh(t5*t1);
    t11 = 0.5*t9+0.5;
    t12 = t11*t11;
    t13 = t12*t11;
    t26 = t12*t12;
    t30 = 30.0*t26-60.0*t13+30.0*t12;
    t32 = 0.247901727*t9+0.3544783819;
    t33 = t32*t32;
    t35 = t33*t33;
    t41 = log(t32);
    t48 = t9*t9;
    return(1/t5*(2.0*t1*A*(4.0*t13+3.0*(2.0*alpha-2.0)*t12+2.0*(-3.0*alpha+1.0)
*t11)+t30*(0.8575855211E-1*t35*t33*t32-0.7110150146E-2*t9+0.464199863E-2)-t30*(
0.5141205242E-7+0.9429832461E-7*t32*t41+0.289612E-7*t9)+0.1E1*t1*t9*(1.0-t48))
);
  }
}

#include <math.h>
double musol(double x)
{
  double t10;
  double t11;
  double t12;
  double t13;
  double t20;
  double t24;
  double t3;
  double t5;
  double t6;
  double t7;
  double t8;
  {
    t3 = tanh(x/delta);
    t5 = 0.5*t3+0.5;
    t6 = t5*t5;
    t7 = t6*t6;
    t8 = t7*t5;
    t10 = 0.247901727*t3+0.3544783819;
    t11 = t10*t10;
    t12 = t11*t11;
    t13 = t12*t11;
    t20 = t6*t5;
    t24 = log(t10);
    return(0.3601859189E1*t8*t13-0.1720892203*t8-0.9004647972E1*t7*t13+
0.4302230507*t7+0.6003098648E1*t20*t13-0.2868153671*t20+0.2111236503E-6+
0.9429832461E-7*t24-0.5657899477E-6*t8*t24+0.1414474869E-5*t7*t24
-0.9429832461E-6*t20*t24);
  }
}

#include <math.h>
double veloSource(double x)
{
  double t1;
  double t10;
  double t101;
  double t11;
  double t24;
  double t28;
  double t3;
  double t30;
  double t31;
  double t33;
  double t39;
  double t4;
  double t5;
  double t51;
  double t55;
  double t56;
  double t6;
  double t60;
  double t69;
  double t81;
  double t88;
  double t9;
  {
    t1 = 1/delta;
    t3 = tanh(x*t1);
    t4 = t3*t3;
    t5 = 1.0-t4;
    t6 = t5*t1;
    t9 = 0.5*t3+0.5;
    t10 = t9*t9;
    t11 = t10*t9;
    t24 = t10*t10;
    t28 = 30.0*t24-60.0*t11+30.0*t10;
    t30 = 0.247901727*t3+0.3544783819;
    t31 = t30*t30;
    t33 = t31*t31;
    t39 = log(t30);
    t51 = t33*t31;
    t55 = t24*t9;
    t56 = t33*t30;
    t60 = t24*t5;
    t69 = t11*t5;
    t81 = 1/t30;
    t88 = t1*t81;
    t101 = 0.9004647972E1*t24*t51*t6+0.535744268E1*t55*t56*t6-0.4302230508*t60*
t1-0.1800929594E2*t11*t51*t6-0.133936067E2*t24*t56*t6+0.8604461014*t69*t1+
0.9004647972E1*t10*t51*t6+0.8929071133E1*t11*t56*t6-0.4302230506*t10*t5*t1+
0.2337671752E-7*t6*t81-0.1414474869E-5*t24*t39*t6-0.1402603052E-6*t55*t5*t88+
0.2828949738E-5*t11*t39*t6+0.3506507628E-6*t60*t88-0.1414474869E-5*t10*t39*t6
-0.2337671752E-6*t69*t88;
    return(-0.5*t6*(2.0*t1*A*(4.0*t11+3.0*(2.0*alpha-2.0)*t10+2.0*(-3.0*alpha+
1.0)*t9)+t28*(0.8575855211E-1*t33*t31*t30-0.7110150146E-2*t3+0.464199863E-2)-
t28*(0.5141205242E-7+0.9429832461E-7*t30*t39+0.289612E-7*t3)+0.1E1*t1*t3*t5)+
t30*t101);
  }
}

#include <math.h>
double rhosol(double x)
{
  double t3;
  {
    t3 = tanh(x/delta);
    return(0.247901727*t3+0.3544783819);
  }
}

#include <math.h>
double gradrho(double x)
{
  double t1;
  double t3;
  double t4;
  {
    t1 = 1/delta;
    t3 = tanh(x*t1);
    t4 = t3*t3;
    return(0.247901727*(1.0-t4)*t1);
  }
}

#include <math.h>
double gradphi(double x)
{
  {
    return(0.5*(1.0-pow(tanh(x/delta),2.0))/delta);
  }
}

#include <math.h>
double thetasol(double x)
{
  double t1;
  double t21;
  double t25;
  double t27;
  double t28;
  double t30;
  double t36;
  double t4;
  double t43;
  double t6;
  double t7;
  double t8;
  {
    t1 = 1/delta;
    t4 = tanh(x*t1);
    t6 = 0.5*t4+0.5;
    t7 = t6*t6;
    t8 = t7*t6;
    t21 = t7*t7;
    t25 = 30.0*t21-60.0*t8+30.0*t7;
    t27 = 0.247901727*t4+0.3544783819;
    t28 = t27*t27;
    t30 = t28*t28;
    t36 = log(t27);
    t43 = t4*t4;
    return(2.0*t1*A*(4.0*t8+3.0*(2.0*alpha-2.0)*t7+2.0*(-3.0*alpha+1.0)*t6)+t25
*(0.8575855211E-1*t30*t28*t27-0.7110150146E-2*t4+0.464199863E-2)-t25*(
0.5120602249E-1+0.9382917919E-1*t27*t36+0.2881711456E-1*t4)+0.1E1*t1*t4*(1.0-
t43));
  }
}

#include <math.h>
double phiSource(double x)
{
  double t1;
  double t11;
  double t12;
  double t13;
  double t26;
  double t3;
  double t30;
  double t32;
  double t33;
  double t35;
  double t41;
  double t48;
  double t5;
  double t9;
  {
    t1 = 1/delta;
    t3 = tanh(x*t1);
    t5 = 0.247901727*t3+0.3544783819;
    t9 = tanh(t5*t1);
    t11 = 0.5*t9+0.5;
    t12 = t11*t11;
    t13 = t12*t11;
    t26 = t12*t12;
    t30 = 30.0*t26-60.0*t13+30.0*t12;
    t32 = 0.247901727*t9+0.3544783819;
    t33 = t32*t32;
    t35 = t33*t33;
    t41 = log(t32);
    t48 = t9*t9;
    return(1/t5*(2.0*t1*A*(4.0*t13+3.0*(2.0*alpha-2.0)*t12+2.0*(-3.0*alpha+1.0)
*t11)+t30*(0.8575855211E-1*t35*t33*t32-0.7110150146E-2*t9+0.464199863E-2)-t30*(
0.5120602249E-1+0.9382917919E-1*t32*t41+0.2881711456E-1*t9)+0.1E1*t1*t9*(1.0-
t48)));
  }
}

#include <math.h>
double musol(double x)
{
  double t10;
  double t11;
  double t12;
  double t13;
  double t20;
  double t24;
  double t3;
  double t5;
  double t6;
  double t7;
  double t8;
  {
    t3 = tanh(x/delta);
    t5 = 0.5*t3+0.5;
    t6 = t5*t5;
    t7 = t6*t6;
    t8 = t7*t5;
    t10 = 0.247901727*t3+0.3544783819;
    t11 = t10*t10;
    t12 = t11*t11;
    t13 = t12*t11;
    t20 = t6*t5;
    t24 = log(t10);
    return(0.3601859189E1*t8*t13-0.1432527663E1*t8-0.9004647972E1*t7*t13+
0.3581319157E1*t7+0.6003098648E1*t20*t13-0.2387546105E1*t20+0.2100732849+
0.9382917919E-1*t24-0.5629750751*t8*t24+0.1407437688E1*t7*t24-0.9382917919*t20*
t24);
  }
}

#include <math.h>
double veloSource(double x)
{
  double t1;
  double t10;
  double t101;
  double t11;
  double t24;
  double t28;
  double t3;
  double t30;
  double t31;
  double t33;
  double t39;
  double t4;
  double t5;
  double t51;
  double t55;
  double t56;
  double t6;
  double t60;
  double t69;
  double t81;
  double t88;
  double t9;
  {
    t1 = 1/delta;
    t3 = tanh(x*t1);
    t4 = t3*t3;
    t5 = 1.0-t4;
    t6 = t5*t1;
    t9 = 0.5*t3+0.5;
    t10 = t9*t9;
    t11 = t10*t9;
    t24 = t10*t10;
    t28 = 30.0*t24-60.0*t11+30.0*t10;
    t30 = 0.247901727*t3+0.3544783819;
    t31 = t30*t30;
    t33 = t31*t31;
    t39 = log(t30);
    t51 = t33*t31;
    t55 = t24*t9;
    t56 = t33*t30;
    t60 = t24*t5;
    t69 = t11*t5;
    t81 = 1/t30;
    t88 = t1*t81;
    t101 = 0.9004647972E1*t24*t51*t6+0.535744268E1*t55*t56*t6-0.3581319158E1*
t60*t1-0.1800929594E2*t11*t51*t6-0.133936067E2*t24*t56*t6+0.7162638314E1*t69*t1
+0.9004647972E1*t10*t51*t6+0.8929071133E1*t11*t56*t6-0.3581319158E1*t10*t5*t1+
0.2326041556E-1*t6*t81-0.1407437688E1*t24*t39*t6-0.1395624934*t55*t5*t88+
0.2814875376E1*t11*t39*t6+0.3489062335*t60*t88-0.1407437688E1*t10*t39*t6
-0.2326041556*t69*t88;
    return(-0.5*t6*(2.0*t1*A*(4.0*t11+3.0*(2.0*alpha-2.0)*t10+2.0*(-3.0*alpha+
1.0)*t9)+t28*(0.8575855211E-1*t33*t31*t30-0.7110150146E-2*t3+0.464199863E-2)-
t28*(0.5120602249E-1+0.9382917919E-1*t30*t39+0.2881711456E-1*t3)+0.1E1*t1*t3*t5
)+t30*t101);
  }
}

#include <math.h>
double rhosol(double x)
{
  double t3;
  {
    t3 = tanh(x/delta);
    return(0.247901727*t3+0.3544783819);
  }
}

#include <math.h>
double gradrho(double x)
{
  double t1;
  double t3;
  double t4;
  {
    t1 = 1/delta;
    t3 = tanh(x*t1);
    t4 = t3*t3;
    return(0.247901727*(1.0-t4)*t1);
  }
}

#include <math.h>
double gradphi(double x)
{
  {
    return(0.5*(1.0-pow(tanh(x/delta),2.0))/delta);
  }
}

#include <math.h>
double thetasol(double x)
{
  double t1;
  double t21;
  double t25;
  double t26;
  double t27;
  double t28;
  double t30;
  double t37;
  double t38;
  double t4;
  double t40;
  double t48;
  double t6;
  double t7;
  double t8;
  {
    t1 = 1/delta;
    t4 = tanh(x*t1);
    t6 = 0.5*t4+0.5;
    t7 = t6*t6;
    t8 = t7*t6;
    t21 = t7*t7;
    t25 = 30.0*t21-60.0*t8+30.0*t7;
    t26 = 0.247901727*t4;
    t27 = t26+0.3544783819;
    t28 = t27*t27;
    t30 = t28*t28;
    t37 = dv(t26+0.247901727);
    t38 = t37*t37;
    t40 = log(t27);
    t48 = t4*t4;
    return(2.0*t1*A*(4.0*t8+3.0*(2.0*alpha-2.0)*t7+2.0*(-3.0*alpha+1.0)*t6)+t25
*(0.8575855211E-1*t30*t28*t27-0.7110150146E-2*t4+0.464199863E-2)-t25*(t38+av*
t27*t40+(bv-av)*t27+0.15-0.1135858337E-1*dv)+0.1E1*t1*t4*(1.0-t48));
  }
}

#include <math.h>
double phiSource(double x)
{
  double t1;
  double t11;
  double t12;
  double t13;
  double t26;
  double t3;
  double t30;
  double t31;
  double t32;
  double t33;
  double t35;
  double t42;
  double t43;
  double t45;
  double t5;
  double t53;
  double t9;
  {
    t1 = 1/delta;
    t3 = tanh(x*t1);
    t5 = 0.247901727*t3+0.3544783819;
    t9 = tanh(t5*t1);
    t11 = 0.5*t9+0.5;
    t12 = t11*t11;
    t13 = t12*t11;
    t26 = t12*t12;
    t30 = 30.0*t26-60.0*t13+30.0*t12;
    t31 = 0.247901727*t9;
    t32 = t31+0.3544783819;
    t33 = t32*t32;
    t35 = t33*t33;
    t42 = dv(t31+0.247901727);
    t43 = t42*t42;
    t45 = log(t32);
    t53 = t9*t9;
    return(1/t5*(2.0*t1*A*(4.0*t13+3.0*(2.0*alpha-2.0)*t12+2.0*(-3.0*alpha+1.0)
*t11)+t30*(0.8575855211E-1*t35*t33*t32-0.7110150146E-2*t9+0.464199863E-2)-t30*(
t43+av*t32*t45+(bv-av)*t32+0.15-0.1135858337E-1*dv)+0.1E1*t1*t9*(1.0-t53)));
  }
}



#include <math.h>
double rhosol(double x)
{
  double t3;
  {
    t3 = tanh(x/delta);
    return(0.247901727*t3+0.3544783819);
  }
}

#include <math.h>
double gradrho(double x)
{
  double t1;
  double t3;
  double t4;
  {
    t1 = 1/delta;
    t3 = tanh(x*t1);
    t4 = t3*t3;
    return(0.247901727*(1.0-t4)*t1);
  }
}

#include <math.h>
double gradphi(double x)
{
  {
    return(0.5*(1.0-pow(tanh(x/delta),2.0))/delta);
  }
}

#include <math.h>
double thetasol(double x)
{
  double t1;
  double t21;
  double t25;
  double t27;
  double t28;
  double t30;
  double t36;
  double t4;
  double t43;
  double t6;
  double t7;
  double t8;
  {
    t1 = 1/delta;
    t4 = tanh(x*t1);
    t6 = 0.5*t4+0.5;
    t7 = t6*t6;
    t8 = t7*t6;
    t21 = t7*t7;
    t25 = 30.0*t21-60.0*t8+30.0*t7;
    t27 = 0.247901727*t4+0.3544783819;
    t28 = t27*t27;
    t30 = t28*t28;
    t36 = log(t27);
    t43 = t4*t4;
    return(2.0*t1*A*(4.0*t8+3.0*(2.0*alpha-2.0)*t7+2.0*(-3.0*alpha+1.0)*t6)+t25
*(0.8575855211E-1*t30*t28*t27-0.7110150146E-2*t4+0.464199863E-2)-t25*(
0.1221987772+0.2239152818*t27*t36+0.6876957021E-1*t4)+0.1E1*t1*t4*(1.0-t43));
  }
}

#include <math.h>
double phiSource(double x)
{
  double t1;
  double t11;
  double t12;
  double t13;
  double t26;
  double t3;
  double t30;
  double t32;
  double t33;
  double t35;
  double t41;
  double t48;
  double t5;
  double t9;
  {
    t1 = 1/delta;
    t3 = tanh(x*t1);
    t5 = 0.247901727*t3+0.3544783819;
    t9 = tanh(t5*t1);
    t11 = 0.5*t9+0.5;
    t12 = t11*t11;
    t13 = t12*t11;
    t26 = t12*t12;
    t30 = 30.0*t26-60.0*t13+30.0*t12;
    t32 = 0.247901727*t9+0.3544783819;
    t33 = t32*t32;
    t35 = t33*t33;
    t41 = log(t32);
    t48 = t9*t9;
    return(1/t5*(2.0*t1*A*(4.0*t13+3.0*(2.0*alpha-2.0)*t12+2.0*(-3.0*alpha+1.0)
*t11)+t30*(0.8575855211E-1*t35*t33*t32-0.7110150146E-2*t9+0.464199863E-2)-t30*(
0.1221987772+0.2239152818*t32*t41+0.6876957021E-1*t9)+0.1E1*t1*t9*(1.0-t48)));
  }
}

#include <math.h>
double musol(double x)
{
  double t10;
  double t11;
  double t12;
  double t13;
  double t20;
  double t24;
  double t3;
  double t5;
  double t6;
  double t7;
  double t8;
  {
    t3 = tanh(x/delta);
    t5 = 0.5*t3+0.5;
    t6 = t5*t5;
    t7 = t6*t6;
    t8 = t7*t5;
    t10 = 0.247901727*t3+0.3544783819;
    t11 = t10*t10;
    t12 = t11*t11;
    t13 = t12*t11;
    t20 = t6*t5;
    t24 = log(t10);
    return(0.3601859189E1*t8*t13-0.3180019123E1*t8-0.9004647972E1*t7*t13+
0.7950047808E1*t7+0.6003098648E1*t20*t13-0.5300031872E1*t20+0.5013218616+
0.2239152818*t24-0.1343491691E1*t8*t24+0.3358729227E1*t7*t24-0.2239152818E1*t20
*t24);
  }
}

#include <math.h>
double veloSource(double x)
{
  double t1;
  double t10;
  double t101;
  double t11;
  double t24;
  double t28;
  double t3;
  double t30;
  double t31;
  double t33;
  double t39;
  double t4;
  double t5;
  double t51;
  double t55;
  double t56;
  double t6;
  double t60;
  double t69;
  double t81;
  double t88;
  double t9;
  {
    t1 = 1/delta;
    t3 = tanh(x*t1);
    t4 = t3*t3;
    t5 = 1.0-t4;
    t6 = t5*t1;
    t9 = 0.5*t3+0.5;
    t10 = t9*t9;
    t11 = t10*t9;
    t24 = t10*t10;
    t28 = 30.0*t24-60.0*t11+30.0*t10;
    t30 = 0.247901727*t3+0.3544783819;
    t31 = t30*t30;
    t33 = t31*t31;
    t39 = log(t30);
    t51 = t33*t31;
    t55 = t24*t9;
    t56 = t33*t30;
    t60 = t24*t5;
    t69 = t11*t5;
    t81 = 1/t30;
    t88 = t1*t81;
    t101 = 0.9004647972E1*t24*t51*t6+0.535744268E1*t55*t56*t6-0.7950047808E1*
t60*t1-0.1800929594E2*t11*t51*t6-0.133936067E2*t24*t56*t6+0.1590009562E2*t69*t1
+0.9004647972E1*t10*t51*t6+0.8929071133E1*t11*t56*t6-0.7950047808E1*t10*t5*t1+
0.5550898506E-1*t6*t81-0.3358729228E1*t24*t39*t6-0.3330539104*t55*t5*t88+
0.6717458454E1*t11*t39*t6+0.8326347759*t60*t88-0.3358729227E1*t10*t39*t6
-0.5550898506*t69*t88;
    return(-0.5*t6*(2.0*t1*A*(4.0*t11+3.0*(2.0*alpha-2.0)*t10+2.0*(-3.0*alpha+
1.0)*t9)+t28*(0.8575855211E-1*t33*t31*t30-0.7110150146E-2*t3+0.464199863E-2)-
t28*(0.1221987772+0.2239152818*t30*t39+0.6876957021E-1*t3)+0.1E1*t1*t3*t5)+t30*
t101);
  }
}

#include <math.h>
double rhosol(double x)
{
  double t3;
  {
    t3 = tanh(x/delta);
    return(0.247901727*t3+0.3544783819);
  }
}

#include <math.h>
double gradrho(double x)
{
  double t1;
  double t3;
  double t4;
  {
    t1 = 1/delta;
    t3 = tanh(x*t1);
    t4 = t3*t3;
    return(0.247901727*(1.0-t4)*t1);
  }
}

#include <math.h>
double gradphi(double x)
{
  {
    return(0.5*(1.0-pow(tanh(x/delta),2.0))/delta);
  }
}

#include <math.h>
double thetasol(double x)
{
  double t1;
  double t21;
  double t25;
  double t27;
  double t28;
  double t30;
  double t36;
  double t4;
  double t43;
  double t6;
  double t7;
  double t8;
  {
    t1 = 1/delta;
    t4 = tanh(x*t1);
    t6 = 0.5*t4+0.5;
    t7 = t6*t6;
    t8 = t7*t6;
    t21 = t7*t7;
    t25 = 30.0*t21-60.0*t8+30.0*t7;
    t27 = 0.247901727*t4+0.3544783819;
    t28 = t27*t27;
    t30 = t28*t28;
    t36 = log(t27);
    t43 = t4*t4;
    return(2.0*t1*A*(4.0*t8+3.0*(2.0*alpha-2.0)*t7+2.0*(-3.0*alpha+1.0)*t6)+t25
*(0.8575855211E-1*t30*t28*t27-0.7110150146E-2*t4+0.464199863E-2)-t25*(
0.7673946604E-1+0.1406162938*t27*t36+0.4318652135E-1*t4)+0.1E1*t1*t4*(1.0-t43)
);
  }
}

#include <math.h>
double phiSource(double x)
{
  double t1;
  double t11;
  double t12;
  double t13;
  double t26;
  double t3;
  double t30;
  double t32;
  double t33;
  double t35;
  double t41;
  double t48;
  double t5;
  double t9;
  {
    t1 = 1/delta;
    t3 = tanh(x*t1);
    t5 = 0.247901727*t3+0.3544783819;
    t9 = tanh(t5*t1);
    t11 = 0.5*t9+0.5;
    t12 = t11*t11;
    t13 = t12*t11;
    t26 = t12*t12;
    t30 = 30.0*t26-60.0*t13+30.0*t12;
    t32 = 0.247901727*t9+0.3544783819;
    t33 = t32*t32;
    t35 = t33*t33;
    t41 = log(t32);
    t48 = t9*t9;
    return(1/t5*(2.0*t1*A*(4.0*t13+3.0*(2.0*alpha-2.0)*t12+2.0*(-3.0*alpha+1.0)
*t11)+t30*(0.8575855211E-1*t35*t33*t32-0.7110150146E-2*t9+0.464199863E-2)-t30*(
0.7673946604E-1+0.1406162938*t32*t41+0.4318652135E-1*t9)+0.1E1*t1*t9*(1.0-t48))
);
  }
}

#include <math.h>
double musol(double x)
{
  double t10;
  double t11;
  double t12;
  double t13;
  double t20;
  double t24;
  double t3;
  double t5;
  double t6;
  double t7;
  double t8;
  {
    t3 = tanh(x/delta);
    t5 = 0.5*t3+0.5;
    t6 = t5*t5;
    t7 = t6*t6;
    t8 = t7*t5;
    t10 = 0.247901727*t3+0.3544783819;
    t11 = t10*t10;
    t12 = t11*t11;
    t13 = t12*t11;
    t20 = t6*t5;
    t24 = log(t10);
    return(0.3601859189E1*t8*t13-0.2061035103E1*t8-0.9004647972E1*t7*t13+
0.5152587757E1*t7+0.6003098648E1*t20*t13-0.3435058505E1*t20+0.3148245249+
0.1406162938*t24-0.8436977628*t8*t24+0.2109244407E1*t7*t24-0.1406162938E1*t20*
t24);
  }
}

#include <math.h>
double veloSource(double x)
{
  double t1;
  double t10;
  double t101;
  double t11;
  double t24;
  double t28;
  double t3;
  double t30;
  double t31;
  double t33;
  double t39;
  double t4;
  double t5;
  double t51;
  double t55;
  double t56;
  double t6;
  double t60;
  double t69;
  double t81;
  double t88;
  double t9;
  {
    t1 = 1/delta;
    t3 = tanh(x*t1);
    t4 = t3*t3;
    t5 = 1.0-t4;
    t6 = t5*t1;
    t9 = 0.5*t3+0.5;
    t10 = t9*t9;
    t11 = t10*t9;
    t24 = t10*t10;
    t28 = 30.0*t24-60.0*t11+30.0*t10;
    t30 = 0.247901727*t3+0.3544783819;
    t31 = t30*t30;
    t33 = t31*t31;
    t39 = log(t30);
    t51 = t33*t31;
    t55 = t24*t9;
    t56 = t33*t30;
    t60 = t24*t5;
    t69 = t11*t5;
    t81 = 1/t30;
    t88 = t1*t81;
    t101 = 0.9004647972E1*t24*t51*t6+0.535744268E1*t55*t56*t6-0.5152587758E1*
t60*t1-0.1800929594E2*t11*t51*t6-0.133936067E2*t24*t56*t6+0.1030517551E2*t69*t1
+0.9004647972E1*t10*t51*t6+0.8929071133E1*t11*t56*t6-0.5152587758E1*t10*t5*t1+
0.3485902208E-1*t6*t81-0.2109244407E1*t24*t39*t6-0.2091541325*t55*t5*t88+
0.4218488814E1*t11*t39*t6+0.5228853312*t60*t88-0.2109244407E1*t10*t39*t6
-0.3485902208*t69*t88;
    return(-0.5*t6*(2.0*t1*A*(4.0*t11+3.0*(2.0*alpha-2.0)*t10+2.0*(-3.0*alpha+
1.0)*t9)+t28*(0.8575855211E-1*t33*t31*t30-0.7110150146E-2*t3+0.464199863E-2)-
t28*(0.7673946604E-1+0.1406162938*t30*t39+0.4318652135E-1*t3)+0.1E1*t1*t3*t5)+
t30*t101);
  }
}

#include <math.h>
double rhosol(double x)
{
  double t3;
  {
    t3 = tanh(x/delta);
    return(0.247901727*t3+0.3544783819);
  }
}

#include <math.h>
double gradrho(double x)
{
  double t1;
  double t3;
  double t4;
  {
    t1 = 1/delta;
    t3 = tanh(x*t1);
    t4 = t3*t3;
    return(0.247901727*(1.0-t4)*t1);
  }
}

#include <math.h>
double gradphi(double x)
{
  {
    return(0.5*(1.0-pow(tanh(x/delta),2.0))/delta);
  }
}

#include <math.h>
double thetasol(double x)
{
  double t1;
  double t21;
  double t25;
  double t27;
  double t28;
  double t30;
  double t36;
  double t4;
  double t43;
  double t6;
  double t7;
  double t8;
  {
    t1 = 1/delta;
    t4 = tanh(x*t1);
    t6 = 0.5*t4+0.5;
    t7 = t6*t6;
    t8 = t7*t6;
    t21 = t7*t7;
    t25 = 30.0*t21-60.0*t8+30.0*t7;
    t27 = 0.247901727*t4+0.3544783819;
    t28 = t27*t27;
    t30 = t28*t28;
    t36 = log(t27);
    t43 = t4*t4;
    return(2.0*t1*A*(4.0*t8+3.0*(2.0*alpha-2.0)*t7+2.0*(-3.0*alpha+1.0)*t6)+t25
*(0.8575855211E-1*t30*t28*t27-0.7110150146E-2*t4+0.464199863E-2)-t25*(
0.7680326859E-1+0.1407332047*t27*t36+0.4322242738E-1*t4)+0.1E1*t1*t4*(1.0-t43)
);
  }
}

#include <math.h>
double phiSource(double x)
{
  double t1;
  double t11;
  double t12;
  double t13;
  double t26;
  double t3;
  double t30;
  double t32;
  double t33;
  double t35;
  double t41;
  double t48;
  double t5;
  double t9;
  {
    t1 = 1/delta;
    t3 = tanh(x*t1);
    t5 = 0.247901727*t3+0.3544783819;
    t9 = tanh(t5*t1);
    t11 = 0.5*t9+0.5;
    t12 = t11*t11;
    t13 = t12*t11;
    t26 = t12*t12;
    t30 = 30.0*t26-60.0*t13+30.0*t12;
    t32 = 0.247901727*t9+0.3544783819;
    t33 = t32*t32;
    t35 = t33*t33;
    t41 = log(t32);
    t48 = t9*t9;
    return(1/t5*(2.0*t1*A*(4.0*t13+3.0*(2.0*alpha-2.0)*t12+2.0*(-3.0*alpha+1.0)
*t11)+t30*(0.8575855211E-1*t35*t33*t32-0.7110150146E-2*t9+0.464199863E-2)-t30*(
0.7680326859E-1+0.1407332047*t32*t41+0.4322242738E-1*t9)+0.1E1*t1*t9*(1.0-t48))
);
  }
}

#include <math.h>
double musol(double x)
{
  double t10;
  double t11;
  double t12;
  double t13;
  double t20;
  double t24;
  double t3;
  double t5;
  double t6;
  double t7;
  double t8;
  {
    t3 = tanh(x/delta);
    t5 = 0.5*t3+0.5;
    t6 = t5*t5;
    t7 = t6*t6;
    t8 = t7*t5;
    t10 = 0.247901727*t3+0.3544783819;
    t11 = t10*t10;
    t12 = t11*t11;
    t13 = t12*t11;
    t20 = t6*t5;
    t24 = log(t10);
    return(0.3601859189E1*t8*t13-0.2062605607E1*t8-0.9004647972E1*t7*t13+
0.5156514018E1*t7+0.6003098648E1*t20*t13-0.3437676012E1*t20+0.3150862756+
0.1407332047*t24-0.8443992282*t8*t24+0.211099807E1*t7*t24-0.1407332047E1*t20*
t24);
  }
}

#include <math.h>
double veloSource(double x)
{
  double t1;
  double t10;
  double t102;
  double t11;
  double t24;
  double t28;
  double t3;
  double t30;
  double t31;
  double t33;
  double t39;
  double t4;
  double t5;
  double t51;
  double t55;
  double t56;
  double t6;
  double t60;
  double t81;
  double t88;
  double t9;
  {
    t1 = 1/delta;
    t3 = tanh(x*t1);
    t4 = t3*t3;
    t5 = 1.0-t4;
    t6 = t5*t1;
    t9 = 0.5*t3+0.5;
    t10 = t9*t9;
    t11 = t10*t9;
    t24 = t10*t10;
    t28 = 30.0*t24-60.0*t11+30.0*t10;
    t30 = 0.247901727*t3+0.3544783819;
    t31 = t30*t30;
    t33 = t31*t31;
    t39 = log(t30);
    t51 = t33*t31;
    t55 = t24*t9;
    t56 = t33*t30;
    t60 = t24*t5;
    t81 = 1/t30;
    t88 = t1*t81;
    t102 = 0.9004647972E1*t24*t51*t6+0.535744268E1*t55*t56*t6-0.5156514018E1*
t60*t1-0.1800929594E2*t11*t51*t6-0.133936067E2*t24*t56*t6+0.1031302804E2*t11*t5
*t1+0.9004647972E1*t10*t51*t6+0.8929071133E1*t11*t56*t6-0.5156514018E1*t10*t5*
t1+0.3488800449E-1*t6*t81-0.211099807E1*t24*t39*t6-0.2093280269*t55*t5*t88+
0.422199614E1*t11*t39*t6+0.5233200672*t60*t88-0.3488800449*t6*t81*t11
-0.211099807E1*t39*t10*t6;
    return(-0.5*t6*(2.0*t1*A*(4.0*t11+3.0*(2.0*alpha-2.0)*t10+2.0*(-3.0*alpha+
1.0)*t9)+t28*(0.8575855211E-1*t33*t31*t30-0.7110150146E-2*t3+0.464199863E-2)-
t28*(0.7680326859E-1+0.1407332047*t30*t39+0.4322242738E-1*t3)+0.1E1*t1*t3*t5)+
t30*t102);
  }
}

#include <math.h>
double rhosol(double x)
{
  double t3;
  {
    t3 = tanh(x/delta);
    return(0.247901727*t3+0.3544783819);
  }
}

#include <math.h>
double gradrho(double x)
{
  double t1;
  double t3;
  double t4;
  {
    t1 = 1/delta;
    t3 = tanh(x*t1);
    t4 = t3*t3;
    return(0.247901727*(1.0-t4)*t1);
  }
}

#include <math.h>
double gradphi(double x)
{
  {
    return(0.5*(1.0-pow(tanh(x/delta),2.0))/delta);
  }
}

#include <math.h>
double thetasol(double x)
{
  double t1;
  double t21;
  double t25;
  double t26;
  double t27;
  double t28;
  double t30;
  double t37;
  double t39;
  double t4;
  double t41;
  double t49;
  double t6;
  double t7;
  double t8;
  {
    t1 = 1/delta;
    t4 = tanh(x*t1);
    t6 = 0.5*t4+0.5;
    t7 = t6*t6;
    t8 = t7*t6;
    t21 = t7*t7;
    t25 = 30.0*t21-60.0*t8+30.0*t7;
    t26 = 0.247901727*t4;
    t27 = t26+0.3544783819;
    t28 = t27*t27;
    t30 = t28*t28;
    t37 = e(t26+0.247901727);
    t39 = pow(t37-5.0,2.0);
    t41 = log(t27);
    t49 = t4*t4;
    return(2.0*t1*A*(4.0*t8+3.0*(2.0*alpha-2.0)*t7+2.0*(-3.0*alpha+1.0)*t6)+t25
*(0.8575855211E-1*t30*t28*t27-0.7110150146E-2*t4+0.464199863E-2)-t25*(t39+av*
t27*t41+(bv-av)*t27+0.7179291685E-1-0.1135858337E-1*e)+0.1E1*t1*t4*(1.0-t49));
  }
}

#include <math.h>
double phiSource(double x)
{
  double t1;
  double t11;
  double t12;
  double t13;
  double t26;
  double t3;
  double t30;
  double t31;
  double t32;
  double t33;
  double t35;
  double t42;
  double t44;
  double t46;
  double t5;
  double t54;
  double t9;
  {
    t1 = 1/delta;
    t3 = tanh(x*t1);
    t5 = 0.247901727*t3+0.3544783819;
    t9 = tanh(t5*t1);
    t11 = 0.5*t9+0.5;
    t12 = t11*t11;
    t13 = t12*t11;
    t26 = t12*t12;
    t30 = 30.0*t26-60.0*t13+30.0*t12;
    t31 = 0.247901727*t9;
    t32 = t31+0.3544783819;
    t33 = t32*t32;
    t35 = t33*t33;
    t42 = e(t31+0.247901727);
    t44 = pow(t42-5.0,2.0);
    t46 = log(t32);
    t54 = t9*t9;
    return(1/t5*(2.0*t1*A*(4.0*t13+3.0*(2.0*alpha-2.0)*t12+2.0*(-3.0*alpha+1.0)
*t11)+t30*(0.8575855211E-1*t35*t33*t32-0.7110150146E-2*t9+0.464199863E-2)-t30*(
t44+av*t32*t46+(bv-av)*t32+0.7179291685E-1-0.1135858337E-1*e)+0.1E1*t1*t9*(1.0-
t54)));
  }
}



#include <math.h>
double rhosol(double x)
{
  double t3;
  {
    t3 = tanh(x/delta);
    return(0.247901727*t3+0.3544783819);
  }
}

#include <math.h>
double gradrho(double x)
{
  double t1;
  double t3;
  double t4;
  {
    t1 = 1/delta;
    t3 = tanh(x*t1);
    t4 = t3*t3;
    return(0.247901727*(1.0-t4)*t1);
  }
}

#include <math.h>
double gradphi(double x)
{
  {
    return(0.5*(1.0-pow(tanh(x/delta),2.0))/delta);
  }
}

#include <math.h>
double thetasol(double x)
{
  double t1;
  double t21;
  double t25;
  double t27;
  double t28;
  double t30;
  double t36;
  double t4;
  double t43;
  double t6;
  double t7;
  double t8;
  {
    t1 = 1/delta;
    t4 = tanh(x*t1);
    t6 = 0.5*t4+0.5;
    t7 = t6*t6;
    t8 = t7*t6;
    t21 = t7*t7;
    t25 = 30.0*t21-60.0*t8+30.0*t7;
    t27 = 0.247901727*t4+0.3544783819;
    t28 = t27*t27;
    t30 = t28*t28;
    t36 = log(t27);
    t43 = t4*t4;
    return(2.0*t1*A*(4.0*t8+3.0*(2.0*alpha-2.0)*t7+2.0*(-3.0*alpha+1.0)*t6)+t25
*(0.8575855211E-1*t30*t28*t27-0.7110150146E-2*t4+0.464199863E-2)-t25*(
0.7680903369E-1+0.1407437685*t27*t36+0.432256718E-1*t4)+0.1E1*t1*t4*(1.0-t43));
  }
}

#include <math.h>
double phiSource(double x)
{
  double t1;
  double t11;
  double t12;
  double t13;
  double t26;
  double t3;
  double t30;
  double t32;
  double t33;
  double t35;
  double t41;
  double t48;
  double t5;
  double t9;
  {
    t1 = 1/delta;
    t3 = tanh(x*t1);
    t5 = 0.247901727*t3+0.3544783819;
    t9 = tanh(t5*t1);
    t11 = 0.5*t9+0.5;
    t12 = t11*t11;
    t13 = t12*t11;
    t26 = t12*t12;
    t30 = 30.0*t26-60.0*t13+30.0*t12;
    t32 = 0.247901727*t9+0.3544783819;
    t33 = t32*t32;
    t35 = t33*t33;
    t41 = log(t32);
    t48 = t9*t9;
    return(1/t5*(2.0*t1*A*(4.0*t13+3.0*(2.0*alpha-2.0)*t12+2.0*(-3.0*alpha+1.0)
*t11)+t30*(0.8575855211E-1*t35*t33*t32-0.7110150146E-2*t9+0.464199863E-2)-t30*(
0.7680903369E-1+0.1407437685*t32*t41+0.432256718E-1*t9)+0.1E1*t1*t9*(1.0-t48))
);
  }
}

#include <math.h>
double musol(double x)
{
  double t10;
  double t11;
  double t12;
  double t13;
  double t20;
  double t24;
  double t3;
  double t5;
  double t6;
  double t7;
  double t8;
  {
    t3 = tanh(x/delta);
    t5 = 0.5*t3+0.5;
    t6 = t5*t5;
    t7 = t6*t6;
    t8 = t7*t5;
    t10 = 0.247901727*t3+0.3544783819;
    t11 = t10*t10;
    t12 = t11*t11;
    t13 = t12*t11;
    t20 = t6*t5;
    t24 = log(t10);
    return(0.3601859189E1*t8*t13-0.2062747515E1*t8-0.9004647972E1*t7*t13+
0.5156868787E1*t7+0.6003098648E1*t20*t13-0.3437912525E1*t20+0.3151099269+
0.1407437685*t24-0.844462611*t8*t24+0.2111156528E1*t7*t24-0.1407437685E1*t20*
t24);
  }
}

#include <math.h>
double veloSource(double x)
{
  double t1;
  double t10;
  double t101;
  double t11;
  double t24;
  double t28;
  double t3;
  double t30;
  double t31;
  double t33;
  double t39;
  double t4;
  double t5;
  double t51;
  double t55;
  double t56;
  double t6;
  double t60;
  double t69;
  double t81;
  double t88;
  double t9;
  {
    t1 = 1/delta;
    t3 = tanh(x*t1);
    t4 = t3*t3;
    t5 = 1.0-t4;
    t6 = t5*t1;
    t9 = 0.5*t3+0.5;
    t10 = t9*t9;
    t11 = t10*t9;
    t24 = t10*t10;
    t28 = 30.0*t24-60.0*t11+30.0*t10;
    t30 = 0.247901727*t3+0.3544783819;
    t31 = t30*t30;
    t33 = t31*t31;
    t39 = log(t30);
    t51 = t33*t31;
    t55 = t24*t9;
    t56 = t33*t30;
    t60 = t24*t5;
    t69 = t11*t5;
    t81 = 1/t30;
    t88 = t1*t81;
    t101 = 0.9004647972E1*t24*t51*t6+0.535744268E1*t55*t56*t6-0.5156868788E1*
t60*t1-0.1800929594E2*t11*t51*t6-0.133936067E2*t24*t56*t6+0.1031373757E2*t69*t1
+0.9004647972E1*t10*t51*t6+0.8929071133E1*t11*t56*t6-0.5156868788E1*t10*t5*t1+
0.3489062328E-1*t6*t81-0.2111156528E1*t24*t39*t6-0.2093437397*t55*t5*t88+
0.4222313056E1*t11*t39*t6+0.5233593493*t60*t88-0.2111156528E1*t10*t39*t6
-0.3489062328*t69*t88;
    return(-0.5*t6*(2.0*t1*A*(4.0*t11+3.0*(2.0*alpha-2.0)*t10+2.0*(-3.0*alpha+
1.0)*t9)+t28*(0.8575855211E-1*t33*t31*t30-0.7110150146E-2*t3+0.464199863E-2)-
t28*(0.7680903369E-1+0.1407437685*t30*t39+0.432256718E-1*t3)+0.1E1*t1*t3*t5)+
t30*t101);
  }
}

#include <math.h>
double rhosol(double x)
{
  double t3;
  {
    t3 = tanh(x/delta);
    return(0.247901727*t3+0.3544783819);
  }
}

#include <math.h>
double gradrho(double x)
{
  double t1;
  double t3;
  double t4;
  {
    t1 = 1/delta;
    t3 = tanh(x*t1);
    t4 = t3*t3;
    return(0.247901727*(1.0-t4)*t1);
  }
}

#include <math.h>
double gradphi(double x)
{
  {
    return(0.5*(1.0-pow(tanh(x/delta),2.0))/delta);
  }
}

#include <math.h>
double thetasol(double x)
{
  double t1;
  double t21;
  double t25;
  double t27;
  double t28;
  double t30;
  double t37;
  double t4;
  double t44;
  double t6;
  double t7;
  double t8;
  {
    t1 = 1/delta;
    t4 = tanh(x*t1);
    t6 = 0.5*t4+0.5;
    t7 = t6*t6;
    t8 = t7*t6;
    t21 = t7*t7;
    t25 = 30.0*t21-60.0*t8+30.0*t7;
    t27 = 0.247901727*t4+0.3544783819;
    t28 = t27*t27;
    t30 = t28*t28;
    t37 = log(t27);
    t44 = t4*t4;
    return(2.0*t1*A*(4.0*t8+3.0*(2.0*alpha-2.0)*t7+2.0*(-3.0*alpha+1.0)*t6)+t25
*(0.8575855211E-1*t30*t28*t27-0.7110150146E-2*t4+0.464199863E-2)-t25*(av*t27*
t37+(bv-av)*t27+cv)+0.1E1*t1*t4*(1.0-t44));
  }
}

#include <math.h>
double phiSource(double x)
{
  double t1;
  double t11;
  double t12;
  double t13;
  double t26;
  double t3;
  double t30;
  double t32;
  double t33;
  double t35;
  double t42;
  double t49;
  double t5;
  double t9;
  {
    t1 = 1/delta;
    t3 = tanh(x*t1);
    t5 = 0.247901727*t3+0.3544783819;
    t9 = tanh(t5*t1);
    t11 = 0.5*t9+0.5;
    t12 = t11*t11;
    t13 = t12*t11;
    t26 = t12*t12;
    t30 = 30.0*t26-60.0*t13+30.0*t12;
    t32 = 0.247901727*t9+0.3544783819;
    t33 = t32*t32;
    t35 = t33*t33;
    t42 = log(t32);
    t49 = t9*t9;
    return(1/t5*(2.0*t1*A*(4.0*t13+3.0*(2.0*alpha-2.0)*t12+2.0*(-3.0*alpha+1.0)
*t11)+t30*(0.8575855211E-1*t35*t33*t32-0.7110150146E-2*t9+0.464199863E-2)-t30*(
av*t32*t42+(bv-av)*t32+cv)+0.1E1*t1*t9*(1.0-t49)));
  }
}

#include <math.h>
double musol(double x)
{
  double t10;
  double t11;
  double t12;
  double t13;
  double t20;
  double t24;
  double t3;
  double t5;
  double t6;
  double t7;
  double t8;
  {
    t3 = tanh(x/delta);
    t5 = 0.5*t3+0.5;
    t6 = t5*t5;
    t7 = t6*t6;
    t8 = t7*t5;
    t10 = 0.247901727*t3+0.3544783819;
    t11 = t10*t10;
    t12 = t11*t11;
    t13 = t12*t11;
    t20 = t6*t5;
    t24 = log(t10);
    return(0.3601859189E1*t8*t13-0.1720879535*t8-0.9004647972E1*t7*t13+
0.4302198838*t7+0.6003098648E1*t20*t13-0.2868132559*t20+av*t24+bv-6.0*t8*av*t24
-6.0*t8*bv+15.0*t7*av*t24+15.0*t7*bv-10.0*t20*av*t24-10.0*t20*bv);
  }
}

#include <math.h>
double veloSource(double x)
{
  double t1;
  double t10;
  double t11;
  double t115;
  double t24;
  double t28;
  double t3;
  double t30;
  double t31;
  double t33;
  double t4;
  double t40;
  double t5;
  double t52;
  double t56;
  double t57;
  double t6;
  double t83;
  double t87;
  double t89;
  double t9;
  double t93;
  double t99;
  {
    t1 = 1/delta;
    t3 = tanh(x*t1);
    t4 = t3*t3;
    t5 = 1.0-t4;
    t6 = t5*t1;
    t9 = 0.5*t3+0.5;
    t10 = t9*t9;
    t11 = t10*t9;
    t24 = t10*t10;
    t28 = 30.0*t24-60.0*t11+30.0*t10;
    t30 = 0.247901727*t3+0.3544783819;
    t31 = t30*t30;
    t33 = t31*t31;
    t40 = log(t30);
    t52 = t33*t31;
    t56 = t24*t9;
    t57 = t33*t30;
    t83 = 1/t30;
    t87 = av*t24;
    t89 = t40*t5*t1;
    t93 = t6*t83;
    t99 = t11*av;
    t115 = 0.9004647972E1*t24*t52*t6+0.535744268E1*t56*t57*t6-0.4302198838*t24*
t5*t1-0.1800929594E2*t11*t52*t6-0.133936067E2*t24*t57*t6+0.8604397676*t11*t5*t1
+0.9004647972E1*t10*t52*t6+0.8929071133E1*t11*t57*t6-0.4302198838*t10*t5*t1+
0.247901727*av*t5*t1*t83-0.15E2*t87*t89-0.1487410362E1*t56*av*t93-0.15E2*t24*bv
*t6+0.3E2*t99*t89+0.3718525905E1*t87*t93+0.3E2*t11*bv*t6-0.15E2*t10*av*t89
-0.247901727E1*t99*t93-0.15E2*t10*bv*t6;
    return(-0.5*t6*(2.0*t1*A*(4.0*t11+3.0*(2.0*alpha-2.0)*t10+2.0*(-3.0*alpha+
1.0)*t9)+t28*(0.8575855211E-1*t33*t31*t30-0.7110150146E-2*t3+0.464199863E-2)-
t28*(av*t30*t40+(bv-av)*t30+cv)+0.1E1*t1*t3*t5)+t30*t115);
  }
}

#include <math.h>
double rhosol(double x)
{
  double t3;
  {
    t3 = tanh(x/delta);
    return(0.247901727*t3+0.3544783819);
  }
}

#include <math.h>
double gradrho(double x)
{
  double t1;
  double t3;
  double t4;
  {
    t1 = 1/delta;
    t3 = tanh(x*t1);
    t4 = t3*t3;
    return(0.247901727*(1.0-t4)*t1);
  }
}

#include <math.h>
double gradphi(double x)
{
  {
    return(0.5*(1.0-pow(tanh(x/delta),2.0))/delta);
  }
}

#include <math.h>
double thetasol(double x)
{
  double t1;
  double t21;
  double t25;
  double t27;
  double t28;
  double t30;
  double t36;
  double t4;
  double t43;
  double t6;
  double t7;
  double t8;
  {
    t1 = 1/delta;
    t4 = tanh(x*t1);
    t6 = 0.5*t4+0.5;
    t7 = t6*t6;
    t8 = t7*t6;
    t21 = t7*t7;
    t25 = 30.0*t21-60.0*t8+30.0*t7;
    t27 = 0.247901727*t4+0.3544783819;
    t28 = t27*t27;
    t30 = t28*t28;
    t36 = log(t27);
    t43 = t4*t4;
    return(2.0*t1*A*(4.0*t8+3.0*(2.0*alpha-2.0)*t7+2.0*(-3.0*alpha+1.0)*t6)+t25
*(0.8575855211E-1*t30*t28*t27-0.7110150146E-2*t4+0.464199863E-2)-t25*(
0.1638030767E-1+0.3001504024E-1*t27*t36+0.9218314185E-2*t4)+0.1E1*t1*t4*(1.0-
t43));
  }
}

#include <math.h>
double phiSource(double x)
{
  double t1;
  double t11;
  double t12;
  double t13;
  double t26;
  double t3;
  double t30;
  double t32;
  double t33;
  double t35;
  double t41;
  double t48;
  double t5;
  double t9;
  {
    t1 = 1/delta;
    t3 = tanh(x*t1);
    t5 = 0.247901727*t3+0.3544783819;
    t9 = tanh(t5*t1);
    t11 = 0.5*t9+0.5;
    t12 = t11*t11;
    t13 = t12*t11;
    t26 = t12*t12;
    t30 = 30.0*t26-60.0*t13+30.0*t12;
    t32 = 0.247901727*t9+0.3544783819;
    t33 = t32*t32;
    t35 = t33*t33;
    t41 = log(t32);
    t48 = t9*t9;
    return(1/t5*(2.0*t1*A*(4.0*t13+3.0*(2.0*alpha-2.0)*t12+2.0*(-3.0*alpha+1.0)
*t11)+t30*(0.8575855211E-1*t35*t33*t32-0.7110150146E-2*t9+0.464199863E-2)-t30*(
0.1638030767E-1+0.3001504024E-1*t32*t41+0.9218314185E-2*t9)+0.1E1*t1*t9*(1.0-
t48)));
  }
}

#include <math.h>
double musol(double x)
{
  double t10;
  double t11;
  double t12;
  double t13;
  double t20;
  double t24;
  double t3;
  double t5;
  double t6;
  double t7;
  double t8;
  {
    t3 = tanh(x/delta);
    t5 = 0.5*t3+0.5;
    t6 = t5*t5;
    t7 = t6*t6;
    t8 = t7*t5;
    t10 = 0.247901727*t3+0.3544783819;
    t11 = t10*t10;
    t12 = t11*t11;
    t13 = t12*t11;
    t20 = t6*t5;
    t24 = log(t10);
    return(0.3601859189E1*t8*t13-0.5752903361*t8-0.9004647972E1*t7*t13+
0.143822584E1*t7+0.6003098648E1*t20*t13-0.9588172269*t20+0.672003971E-1+
0.3001504024E-1*t24-0.1800902414*t8*t24+0.4502256036*t7*t24-0.3001504024*t20*
t24);
  }
}

#include <math.h>
double veloSource(double x)
{
  double t1;
  double t10;
  double t102;
  double t11;
  double t24;
  double t28;
  double t3;
  double t30;
  double t31;
  double t33;
  double t39;
  double t4;
  double t5;
  double t51;
  double t55;
  double t56;
  double t6;
  double t60;
  double t81;
  double t88;
  double t9;
  {
    t1 = 1/delta;
    t3 = tanh(x*t1);
    t4 = t3*t3;
    t5 = 1.0-t4;
    t6 = t5*t1;
    t9 = 0.5*t3+0.5;
    t10 = t9*t9;
    t11 = t10*t9;
    t24 = t10*t10;
    t28 = 30.0*t24-60.0*t11+30.0*t10;
    t30 = 0.247901727*t3+0.3544783819;
    t31 = t30*t30;
    t33 = t31*t31;
    t39 = log(t30);
    t51 = t33*t31;
    t55 = t24*t9;
    t56 = t33*t30;
    t60 = t24*t5;
    t81 = 1/t30;
    t88 = t1*t81;
    t102 = 0.9004647972E1*t24*t51*t6+0.535744268E1*t55*t56*t6-0.143822584E1*t60
*t1-0.1800929594E2*t11*t51*t6-0.133936067E2*t24*t56*t6+0.287645168E1*t11*t5*t1+
0.9004647972E1*t10*t51*t6+0.8929071133E1*t11*t56*t6-0.143822584E1*t10*t5*t1+
0.7440780311E-2*t6*t81-0.4502256035*t24*t39*t6-0.4464468186E-1*t55*t5*t88+
0.9004512072*t11*t39*t6+0.1116117047*t60*t88-0.7440780311E-1*t6*t81*t11
-0.4502256036*t39*t10*t6;
    return(-0.5*t6*(2.0*t1*A*(4.0*t11+3.0*(2.0*alpha-2.0)*t10+2.0*(-3.0*alpha+
1.0)*t9)+t28*(0.8575855211E-1*t33*t31*t30-0.7110150146E-2*t3+0.464199863E-2)-
t28*(0.1638030767E-1+0.3001504024E-1*t30*t39+0.9218314185E-2*t3)+0.1E1*t1*t3*t5
)+t30*t102);
  }
}

#include <math.h>
double rhosol(double x)
{
  double t3;
  {
    t3 = tanh(x/delta);
    return(0.247901727*t3+0.3544783819);
  }
}

#include <math.h>
double gradrho(double x)
{
  double t1;
  double t3;
  double t4;
  {
    t1 = 1/delta;
    t3 = tanh(x*t1);
    t4 = t3*t3;
    return(0.247901727*(1.0-t4)*t1);
  }
}

#include <math.h>
double gradphi(double x)
{
  {
    return(0.5*(1.0-pow(tanh(x/delta),2.0))/delta);
  }
}

#include <math.h>
double thetasol(double x)
{
  double t1;
  double t21;
  double t25;
  double t27;
  double t28;
  double t30;
  double t36;
  double t4;
  double t43;
  double t6;
  double t7;
  double t8;
  {
    t1 = 1/delta;
    t4 = tanh(x*t1);
    t6 = 0.5*t4+0.5;
    t7 = t6*t6;
    t8 = t7*t6;
    t21 = t7*t7;
    t25 = 30.0*t21-60.0*t8+30.0*t7;
    t27 = 0.247901727*t4+0.3544783819;
    t28 = t27*t27;
    t30 = t28*t28;
    t36 = log(t27);
    t43 = t4*t4;
    return(2.0*t1*A*(4.0*t8+3.0*(2.0*alpha-2.0)*t7+2.0*(-3.0*alpha+1.0)*t6)+t25
*(0.8575855211E-1*t30*t28*t27-0.7110150146E-2*t4+0.464199863E-2)-t25*(
0.3001504024E-1*t27*t36+0.9218314185E-2*t4+0.1638030767E-1)+0.1E1*t1*t4*(1.0-
t43));
  }
}

#include <math.h>
double phiSource(double x)
{
  double t1;
  double t11;
  double t12;
  double t13;
  double t26;
  double t3;
  double t30;
  double t32;
  double t33;
  double t35;
  double t41;
  double t48;
  double t5;
  double t9;
  {
    t1 = 1/delta;
    t3 = tanh(x*t1);
    t5 = 0.247901727*t3+0.3544783819;
    t9 = tanh(t5*t1);
    t11 = 0.5*t9+0.5;
    t12 = t11*t11;
    t13 = t12*t11;
    t26 = t12*t12;
    t30 = 30.0*t26-60.0*t13+30.0*t12;
    t32 = 0.247901727*t9+0.3544783819;
    t33 = t32*t32;
    t35 = t33*t33;
    t41 = log(t32);
    t48 = t9*t9;
    return(1/t5*(2.0*t1*A*(4.0*t13+3.0*(2.0*alpha-2.0)*t12+2.0*(-3.0*alpha+1.0)
*t11)+t30*(0.8575855211E-1*t35*t33*t32-0.7110150146E-2*t9+0.464199863E-2)-t30*(
0.3001504024E-1*t32*t41+0.9218314185E-2*t9+0.1638030767E-1)+0.1E1*t1*t9*(1.0-
t48)));
  }
}

#include <math.h>
double musol(double x)
{
  double t10;
  double t11;
  double t12;
  double t13;
  double t20;
  double t24;
  double t3;
  double t5;
  double t6;
  double t7;
  double t8;
  {
    t3 = tanh(x/delta);
    t5 = 0.5*t3+0.5;
    t6 = t5*t5;
    t7 = t6*t6;
    t8 = t7*t5;
    t10 = 0.247901727*t3+0.3544783819;
    t11 = t10*t10;
    t12 = t11*t11;
    t13 = t12*t11;
    t20 = t6*t5;
    t24 = log(t10);
    return(0.3601859189E1*t8*t13-0.5752903361*t8-0.9004647972E1*t7*t13+
0.143822584E1*t7+0.6003098648E1*t20*t13-0.9588172269*t20+0.3001504024E-1*t24+
0.672003971E-1-0.1800902414*t8*t24+0.4502256036*t7*t24-0.3001504024*t20*t24);
  }
}

#include <math.h>
double veloSource(double x)
{
  double t1;
  double t10;
  double t102;
  double t11;
  double t24;
  double t28;
  double t3;
  double t30;
  double t31;
  double t33;
  double t39;
  double t4;
  double t5;
  double t51;
  double t55;
  double t56;
  double t6;
  double t60;
  double t81;
  double t88;
  double t9;
  {
    t1 = 1/delta;
    t3 = tanh(x*t1);
    t4 = t3*t3;
    t5 = 1.0-t4;
    t6 = t5*t1;
    t9 = 0.5*t3+0.5;
    t10 = t9*t9;
    t11 = t10*t9;
    t24 = t10*t10;
    t28 = 30.0*t24-60.0*t11+30.0*t10;
    t30 = 0.247901727*t3+0.3544783819;
    t31 = t30*t30;
    t33 = t31*t31;
    t39 = log(t30);
    t51 = t33*t31;
    t55 = t24*t9;
    t56 = t33*t30;
    t60 = t24*t5;
    t81 = 1/t30;
    t88 = t1*t81;
    t102 = 0.9004647972E1*t24*t51*t6+0.535744268E1*t55*t56*t6-0.143822584E1*t60
*t1-0.1800929594E2*t11*t51*t6-0.133936067E2*t24*t56*t6+0.287645168E1*t11*t5*t1+
0.9004647972E1*t10*t51*t6+0.8929071133E1*t11*t56*t6-0.143822584E1*t10*t5*t1+
0.7440780311E-2*t6*t81-0.4502256035*t24*t39*t6-0.4464468186E-1*t55*t5*t88+
0.9004512072*t11*t39*t6+0.1116117047*t60*t88-0.7440780311E-1*t6*t81*t11
-0.4502256036*t39*t10*t6;
    return(-0.5*t6*(2.0*t1*A*(4.0*t11+3.0*(2.0*alpha-2.0)*t10+2.0*(-3.0*alpha+
1.0)*t9)+t28*(0.8575855211E-1*t33*t31*t30-0.7110150146E-2*t3+0.464199863E-2)-
t28*(0.3001504024E-1*t30*t39+0.9218314185E-2*t3+0.1638030767E-1)+0.1E1*t1*t3*t5
)+t30*t102);
  }
}


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
  double t45;
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
    t45 = t4*t4;
    return(2.0*t1*A*(4.0*t8+3.0*(2.0*alpha-2.0)*t7+2.0*(-3.0*alpha+1.0)*t6)+
beta*(t25*(0.35*t30*t28*t27-0.2901812679E-1*t4+0.189450439E-1)-t25*(
0.1753186564*t27*t36-0.6186232473E-1*t4-0.8845786205E-1))+0.1E1*t1*t4*(1.0-t45)
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
  double t5;
  double t50;
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
    t50 = t9*t9;
    return(1/t5*(2.0*t1*A*(4.0*t13+3.0*(2.0*alpha-2.0)*t12+2.0*(-3.0*alpha+1.0)
*t11)+beta*(t30*(0.35*t35*t33*t32-0.2901812679E-1*t9+0.189450439E-1)-t30*(
0.1753186564*t32*t41-0.6186232473E-1*t9-0.8845786205E-1))+0.1E1*t1*t9*(1.0-t50)
));
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
  double t32;
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
    t32 = -0.735E12*t8*t13+0.1284896222E11*t8+0.18375E13*t7*t13-0.3212240554E11
*t7-0.1225E13*t20*t13+0.2141493703E11*t20-8765932820.0*t24+3711254307.0+
0.5259559692E11*t8*t24-0.1314889923E12*t7*t24+0.876593282E11*t20*t24;
    return(-0.2E-10*beta*t32);
  }
}

#include <math.h>
double veloSource(double x)
{
  double t1;
  double t10;
  double t104;
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
  double t54;
  double t58;
  double t59;
  double t6;
  double t63;
  double t72;
  double t84;
  double t9;
  double t91;
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
    t54 = t33*t31;
    t58 = t24*t9;
    t59 = t33*t30;
    t63 = t24*t5;
    t72 = t11*t5;
    t84 = 1/t30;
    t91 = t1*t84;
    t104 = -0.18375E13*t24*t54*t6-0.1093246616E13*t58*t59*t6+0.3212240555E11*
t63*t1+0.3675E13*t11*t54*t6+0.273311654E13*t24*t59*t6-0.6424481108E11*t72*t1
-0.18375E13*t10*t54*t6-0.1822077693E13*t11*t59*t6+0.3212240554E11*t10*t5*t1
-2173089885.0*t6*t84+0.1314889923E12*t24*t39*t6+0.1303853931E11*t58*t5*t91
-0.2629779846E12*t11*t39*t6-0.3259634827E11*t63*t91+0.1314889923E12*t10*t39*t6+
0.2173089885E11*t72*t91;
    return(-0.5*t6*(2.0*t1*A*(4.0*t11+3.0*(2.0*alpha-2.0)*t10+2.0*(-3.0*alpha+
1.0)*t9)+beta*(t28*(0.35*t33*t31*t30-0.2901812679E-1*t3+0.189450439E-1)-t28*(
0.1753186564*t30*t39-0.6186232473E-1*t3-0.8845786205E-1))+0.1E1*t1*t3*t5)
-0.2E-10*t30*beta*t104);
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
  double t45;
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
    t45 = t4*t4;
    return(2.0*t1*A*(4.0*t8+3.0*(2.0*alpha-2.0)*t7+2.0*(-3.0*alpha+1.0)*t6)+
beta*(t25*(0.35*t30*t28*t27-0.2901812679E-1*t4+0.189450439E-1)-t25*(
0.1753186564*t27*t36-0.6186232473E-1*t4-0.8845786205E-1))+0.1E1*t1*t4*(1.0-t45)
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
  double t5;
  double t50;
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
    t50 = t9*t9;
    return(1/t5*(2.0*t1*A*(4.0*t13+3.0*(2.0*alpha-2.0)*t12+2.0*(-3.0*alpha+1.0)
*t11)+beta*(t30*(0.35*t35*t33*t32-0.2901812679E-1*t9+0.189450439E-1)-t30*(
0.1753186564*t32*t41-0.6186232473E-1*t9-0.8845786205E-1))+0.1E1*t1*t9*(1.0-t50)
));
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
  double t32;
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
    t32 = -0.735E12*t8*t13+0.1284896222E11*t8+0.18375E13*t7*t13-0.3212240554E11
*t7-0.1225E13*t20*t13+0.2141493703E11*t20-8765932820.0*t24+3711254307.0+
0.5259559692E11*t8*t24-0.1314889923E12*t7*t24+0.876593282E11*t20*t24;
    return(-0.2E-10*beta*t32);
  }
}

#include <math.h>
double veloSource(double x)
{
  double t1;
  double t10;
  double t104;
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
  double t54;
  double t58;
  double t59;
  double t6;
  double t63;
  double t72;
  double t84;
  double t9;
  double t91;
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
    t54 = t33*t31;
    t58 = t24*t9;
    t59 = t33*t30;
    t63 = t24*t5;
    t72 = t11*t5;
    t84 = 1/t30;
    t91 = t1*t84;
    t104 = -0.18375E13*t24*t54*t6-0.1093246616E13*t58*t59*t6+0.3212240555E11*
t63*t1+0.3675E13*t11*t54*t6+0.273311654E13*t24*t59*t6-0.6424481108E11*t72*t1
-0.18375E13*t10*t54*t6-0.1822077693E13*t11*t59*t6+0.3212240554E11*t10*t5*t1
-2173089885.0*t6*t84+0.1314889923E12*t24*t39*t6+0.1303853931E11*t58*t5*t91
-0.2629779846E12*t11*t39*t6-0.3259634827E11*t63*t91+0.1314889923E12*t10*t39*t6+
0.2173089885E11*t72*t91;
    return(-0.5*t6*(2.0*t1*A*(4.0*t11+3.0*(2.0*alpha-2.0)*t10+2.0*(-3.0*alpha+
1.0)*t9)+beta*(t28*(0.35*t33*t31*t30-0.2901812679E-1*t3+0.189450439E-1)-t28*(
0.1753186564*t30*t39-0.6186232473E-1*t3-0.8845786205E-1))+0.1E1*t1*t3*t5)
-0.2E-10*t30*beta*t104);
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
  double t45;
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
    t45 = t4*t4;
    return(2.0*t1*A*(4.0*t8+3.0*(2.0*alpha-2.0)*t7+2.0*(-3.0*alpha+1.0)*t6)+
beta*(t25*(0.35*t30*t28*t27-0.2901812679E-1*t4+0.189450439E-1)-t25*(
0.1753186564*t27*t36-0.6186232473E-1*t4-0.8845786205E-1))+0.1E1*t1*t4*(1.0-t45)
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
  double t5;
  double t50;
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
    t50 = t9*t9;
    return(1/t5*(2.0*t1*A*(4.0*t13+3.0*(2.0*alpha-2.0)*t12+2.0*(-3.0*alpha+1.0)
*t11)+beta*(t30*(0.35*t35*t33*t32-0.2901812679E-1*t9+0.189450439E-1)-t30*(
0.1753186564*t32*t41-0.6186232473E-1*t9-0.8845786205E-1))+0.1E1*t1*t9*(1.0-t50)
));
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
  double t32;
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
    t32 = -0.735E12*t8*t13+0.1284896222E11*t8+0.18375E13*t7*t13-0.3212240554E11
*t7-0.1225E13*t20*t13+0.2141493703E11*t20-8765932820.0*t24+3711254307.0+
0.5259559692E11*t8*t24-0.1314889923E12*t7*t24+0.876593282E11*t20*t24;
    return(-0.2E-10*beta*t32);
  }
}

#include <math.h>
double veloSource(double x)
{
  double t1;
  double t10;
  double t104;
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
  double t54;
  double t58;
  double t59;
  double t6;
  double t63;
  double t72;
  double t84;
  double t9;
  double t91;
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
    t54 = t33*t31;
    t58 = t24*t9;
    t59 = t33*t30;
    t63 = t24*t5;
    t72 = t11*t5;
    t84 = 1/t30;
    t91 = t1*t84;
    t104 = -0.18375E13*t24*t54*t6-0.1093246616E13*t58*t59*t6+0.3212240555E11*
t63*t1+0.3675E13*t11*t54*t6+0.273311654E13*t24*t59*t6-0.6424481108E11*t72*t1
-0.18375E13*t10*t54*t6-0.1822077693E13*t11*t59*t6+0.3212240554E11*t10*t5*t1
-2173089885.0*t6*t84+0.1314889923E12*t24*t39*t6+0.1303853931E11*t58*t5*t91
-0.2629779846E12*t11*t39*t6-0.3259634827E11*t63*t91+0.1314889923E12*t10*t39*t6+
0.2173089885E11*t72*t91;
    return(-0.5*t6*(2.0*t1*A*(4.0*t11+3.0*(2.0*alpha-2.0)*t10+2.0*(-3.0*alpha+
1.0)*t9)+beta*(t28*(0.35*t33*t31*t30-0.2901812679E-1*t3+0.189450439E-1)-t28*(
0.1753186564*t30*t39-0.6186232473E-1*t3-0.8845786205E-1))+0.1E1*t1*t3*t5)
-0.2E-10*t30*beta*t104);
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
  double t45;
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
    t45 = t4*t4;
    return(2.0*t1*A*(4.0*t8+3.0*(2.0*alpha-2.0)*t7+2.0*(-3.0*alpha+1.0)*t6)+
beta*(t25*(0.35*t30*t28*t27-0.2901812679E-1*t4+0.189450439E-1)-t25*(
0.4691458936E-9*t27*t36+0.144085572E-9*t4+0.2060301114E-9))+0.1E1*t1*t4*(1.0-
t45));
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
  double t5;
  double t50;
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
    t50 = t9*t9;
    return(1/t5*(2.0*t1*A*(4.0*t13+3.0*(2.0*alpha-2.0)*t12+2.0*(-3.0*alpha+1.0)
*t11)+beta*(t30*(0.35*t35*t33*t32-0.2901812679E-1*t9+0.189450439E-1)-t30*(
0.4691458936E-9*t32*t41+0.144085572E-9*t9+0.2060301114E-9))+0.1E1*t1*t9*(1.0-
t50)));
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
  double t32;
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
    t32 = -0.735E20*t8*t13+0.3511648838E19*t8+0.18375E21*t7*t13-0.8779122094E19
*t7-0.1225E21*t20*t13+0.5852748063E19*t20-2345729468.0*t24-5251832095.0+
0.1407437681E11*t8*t24-0.3518594202E11*t7*t24+0.2345729468E11*t20*t24;
    return(-0.2E-18*beta*t32);
  }
}

#include <math.h>
double veloSource(double x)
{
  double t1;
  double t10;
  double t104;
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
  double t54;
  double t58;
  double t59;
  double t6;
  double t63;
  double t72;
  double t84;
  double t9;
  double t91;
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
    t54 = t33*t31;
    t58 = t24*t9;
    t59 = t33*t30;
    t63 = t24*t5;
    t72 = t11*t5;
    t84 = 1/t30;
    t91 = t1*t84;
    t104 = -0.18375E21*t24*t54*t6-0.1093246616E21*t58*t59*t6+0.8779122095E19*
t63*t1+0.3675E21*t11*t54*t6+0.273311654E21*t24*t59*t6-0.1755824419E20*t72*t1
-0.18375E21*t10*t54*t6-0.1822077693E21*t11*t59*t6+0.8779122094E19*t10*t5*t1
-0.5815103862E9*t6*t84+0.3518594202E11*t24*t39*t6+3489062318.0*t58*t5*t91
-0.7037188404E11*t11*t39*t6-8722655793.0*t63*t91+0.3518594202E11*t10*t39*t6+
5815103862.0*t72*t91;
    return(-0.5*t6*(2.0*t1*A*(4.0*t11+3.0*(2.0*alpha-2.0)*t10+2.0*(-3.0*alpha+
1.0)*t9)+beta*(t28*(0.35*t33*t31*t30-0.2901812679E-1*t3+0.189450439E-1)-t28*(
0.4691458936E-9*t30*t39+0.144085572E-9*t3+0.2060301114E-9))+0.1E1*t1*t3*t5)
-0.2E-18*t30*beta*t104);
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
  double t45;
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
    t45 = t4*t4;
    return(2.0*t1*A*(4.0*t8+3.0*(2.0*alpha-2.0)*t7+2.0*(-3.0*alpha+1.0)*t6)+
beta*(t25*(0.35*t30*t28*t27-0.2901812679E-1*t4+0.189450439E-1)-t25*(
0.4691458936E-9*t27*t36+0.144085572E-9*t4+0.2060301114E-9))+0.1E1*t1*t4*(1.0-
t45));
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
  double t5;
  double t50;
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
    t50 = t9*t9;
    return(1/t5*(2.0*t1*A*(4.0*t13+3.0*(2.0*alpha-2.0)*t12+2.0*(-3.0*alpha+1.0)
*t11)+beta*(t30*(0.35*t35*t33*t32-0.2901812679E-1*t9+0.189450439E-1)-t30*(
0.4691458936E-9*t32*t41+0.144085572E-9*t9+0.2060301114E-9))+0.1E1*t1*t9*(1.0-
t50)));
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
  double t32;
  double t5;
  double t6;
  double t7;
  double t8;
  {
    t3 = tanh(x/delta);
    t5 = 0.5+0.5*t3;
    t6 = t5*t5;
    t7 = t6*t6;
    t8 = t7*t5;
    t10 = 0.247901727*t3+0.3544783819;
    t11 = t10*t10;
    t12 = t11*t11;
    t13 = t12*t11;
    t20 = t6*t5;
    t24 = log(t10);
    t32 = -5251832095.0-0.735E20*t8*t13+0.3511648838E19*t8+0.18375E21*t7*t13
-0.8779122094E19*t7-0.1225E21*t20*t13+0.5852748063E19*t20-2345729468.0*t24+
0.1407437681E11*t8*t24-0.3518594202E11*t7*t24+0.2345729468E11*t24*t20;
    return(-0.2E-18*beta*t32);
  }
}

#include <math.h>
double veloSource(double x)
{
  double t1;
  double t10;
  double t104;
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
  double t54;
  double t58;
  double t59;
  double t6;
  double t63;
  double t72;
  double t84;
  double t9;
  double t91;
  {
    t1 = 1/delta;
    t3 = tanh(x*t1);
    t4 = t3*t3;
    t5 = 1.0-t4;
    t6 = t5*t1;
    t9 = 0.5+0.5*t3;
    t10 = t9*t9;
    t11 = t10*t9;
    t24 = t10*t10;
    t28 = 30.0*t24-60.0*t11+30.0*t10;
    t30 = 0.247901727*t3+0.3544783819;
    t31 = t30*t30;
    t33 = t31*t31;
    t39 = log(t30);
    t54 = t33*t31;
    t58 = t24*t9;
    t59 = t33*t30;
    t63 = t24*t5;
    t72 = t11*t5;
    t84 = 1/t30;
    t91 = t1*t84;
    t104 = -0.18375E21*t24*t54*t6-0.1093246616E21*t58*t59*t6+0.8779122095E19*
t63*t1+0.3675E21*t11*t54*t6+0.273311654E21*t24*t59*t6-0.1755824419E20*t72*t1
-0.18375E21*t10*t54*t6-0.1822077693E21*t11*t59*t6+0.8779122094E19*t10*t5*t1
-0.5815103862E9*t6*t84+0.3518594202E11*t24*t39*t6+3489062318.0*t58*t5*t91
-0.7037188404E11*t11*t39*t6-8722655793.0*t63*t91+0.3518594202E11*t10*t39*t6+
5815103862.0*t72*t91;
    return(-0.5*t6*(2.0*t1*A*(4.0*t11+3.0*(2.0*alpha-2.0)*t10+2.0*(-3.0*alpha+
1.0)*t9)+beta*(t28*(0.35*t33*t31*t30-0.2901812679E-1*t3+0.189450439E-1)-t28*(
0.4691458936E-9*t30*t39+0.144085572E-9*t3+0.2060301114E-9))+0.1E1*t1*t3*t5)
-0.2E-18*t30*beta*t104);
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
  double t45;
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
    t45 = t4*t4;
    return(2.0*t1*A*(4.0*t8+3.0*(2.0*alpha-2.0)*t7+2.0*(-3.0*alpha+1.0)*t6)+
beta*(t25*(0.35*t30*t28*t27-0.2901812679E-1*t4+0.189450439E-1)-t25*(
0.4691458936E-9*t27*t36+0.144085572E-9*t4+0.2060301114E-9))+0.1E1*t1*t4*(1.0-
t45));
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
  double t5;
  double t50;
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
    t50 = t9*t9;
    return(1/t5*(2.0*t1*A*(4.0*t13+3.0*(2.0*alpha-2.0)*t12+2.0*(-3.0*alpha+1.0)
*t11)+beta*(t30*(0.35*t35*t33*t32-0.2901812679E-1*t9+0.189450439E-1)-t30*(
0.4691458936E-9*t32*t41+0.144085572E-9*t9+0.2060301114E-9))+0.1E1*t1*t9*(1.0-
t50)));
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
  double t32;
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
    t32 = -0.735E20*t8*t13+0.3511648838E19*t8+0.18375E21*t7*t13-0.8779122094E19
*t7-0.1225E21*t20*t13+0.5852748063E19*t20-2345729468.0*t24-5251832095.0+
0.1407437681E11*t8*t24-0.3518594202E11*t7*t24+0.2345729468E11*t20*t24;
    return(-0.2E-18*beta*t32);
  }
}

#include <math.h>
double veloSource(double x)
{
  double t1;
  double t10;
  double t104;
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
  double t54;
  double t58;
  double t59;
  double t6;
  double t63;
  double t72;
  double t84;
  double t9;
  double t91;
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
    t54 = t33*t31;
    t58 = t24*t9;
    t59 = t33*t30;
    t63 = t24*t5;
    t72 = t11*t5;
    t84 = 1/t30;
    t91 = t1*t84;
    t104 = -0.18375E21*t24*t54*t6-0.1093246616E21*t58*t59*t6+0.8779122095E19*
t63*t1+0.3675E21*t11*t54*t6+0.273311654E21*t24*t59*t6-0.1755824419E20*t72*t1
-0.18375E21*t10*t54*t6-0.1822077693E21*t11*t59*t6+0.8779122094E19*t10*t5*t1
-0.5815103862E9*t6*t84+0.3518594202E11*t24*t39*t6+3489062318.0*t58*t5*t91
-0.7037188404E11*t11*t39*t6-8722655793.0*t63*t91+0.3518594202E11*t10*t39*t6+
5815103862.0*t72*t91;
    return(-0.5*t6*(2.0*t1*A*(4.0*t11+3.0*(2.0*alpha-2.0)*t10+2.0*(-3.0*alpha+
1.0)*t9)+beta*(t28*(0.35*t33*t31*t30-0.2901812679E-1*t3+0.189450439E-1)-t28*(
0.4691458936E-9*t30*t39+0.144085572E-9*t3+0.2060301114E-9))+0.1E1*t1*t3*t5)
-0.2E-18*t30*beta*t104);
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
  double t45;
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
    t45 = t4*t4;
    return(2.0*t1*A*(4.0*t8+3.0*(2.0*alpha-2.0)*t7+2.0*(-3.0*alpha+1.0)*t6)+
beta*(t25*(0.35*t30*t28*t27-0.2901812679E-1*t4+0.189450439E-1)-t25*(
0.4691458936E-9*t27*t36+0.144085572E-9*t4+0.2060301114E-9))+0.1E1*t1*t4*(1.0-
t45));
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
  double t5;
  double t50;
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
    t50 = t9*t9;
    return(1/t5*(2.0*t1*A*(4.0*t13+3.0*(2.0*alpha-2.0)*t12+2.0*(-3.0*alpha+1.0)
*t11)+beta*(t30*(0.35*t35*t33*t32-0.2901812679E-1*t9+0.189450439E-1)-t30*(
0.4691458936E-9*t32*t41+0.144085572E-9*t9+0.2060301114E-9))+0.1E1*t1*t9*(1.0-
t50)));
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
  double t32;
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
    t32 = -0.735E20*t8*t13+0.3511648838E19*t8+0.18375E21*t7*t13-0.8779122094E19
*t7-0.1225E21*t20*t13+0.5852748063E19*t20-2345729468.0*t24-5251832095.0+
0.1407437681E11*t8*t24-0.3518594202E11*t7*t24+0.2345729468E11*t20*t24;
    return(-0.2E-18*beta*t32);
  }
}

#include <math.h>
double veloSource(double x)
{
  double t1;
  double t10;
  double t104;
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
  double t54;
  double t58;
  double t59;
  double t6;
  double t63;
  double t72;
  double t84;
  double t9;
  double t91;
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
    t54 = t33*t31;
    t58 = t24*t9;
    t59 = t33*t30;
    t63 = t24*t5;
    t72 = t11*t5;
    t84 = 1/t30;
    t91 = t1*t84;
    t104 = -0.18375E21*t24*t54*t6-0.1093246616E21*t58*t59*t6+0.8779122095E19*
t63*t1+0.3675E21*t11*t54*t6+0.273311654E21*t24*t59*t6-0.1755824419E20*t72*t1
-0.18375E21*t10*t54*t6-0.1822077693E21*t11*t59*t6+0.8779122094E19*t10*t5*t1
-0.5815103862E9*t6*t84+0.3518594202E11*t24*t39*t6+3489062318.0*t58*t5*t91
-0.7037188404E11*t11*t39*t6-8722655793.0*t63*t91+0.3518594202E11*t10*t39*t6+
5815103862.0*t72*t91;
    return(-0.5*t6*(2.0*t1*A*(4.0*t11+3.0*(2.0*alpha-2.0)*t10+2.0*(-3.0*alpha+
1.0)*t9)+beta*(t28*(0.35*t33*t31*t30-0.2901812679E-1*t3+0.189450439E-1)-t28*(
0.4691458936E-9*t30*t39+0.144085572E-9*t3+0.2060301114E-9))+0.1E1*t1*t3*t5)
-0.2E-18*t30*beta*t104);
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
  double t45;
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
    t45 = t4*t4;
    return(2.0*t1*A*(4.0*t8+3.0*(2.0*alpha-2.0)*t7+2.0*(-3.0*alpha+1.0)*t6)+
beta*(t25*(0.35*t30*t28*t27-0.2901812679E-1*t4+0.189450439E-1)-t25*(
0.4691458936E-9*t27*t36+0.144085572E-9*t4+0.2060301114E-9))+0.1E1*t1*t4*(1.0-
t45));
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
  double t5;
  double t50;
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
    t50 = t9*t9;
    return(1/t5*(2.0*t1*A*(4.0*t13+3.0*(2.0*alpha-2.0)*t12+2.0*(-3.0*alpha+1.0)
*t11)+beta*(t30*(0.35*t35*t33*t32-0.2901812679E-1*t9+0.189450439E-1)-t30*(
0.4691458936E-9*t32*t41+0.144085572E-9*t9+0.2060301114E-9))+0.1E1*t1*t9*(1.0-
t50)));
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
  double t32;
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
    t32 = -0.735E20*t8*t13+0.3511648838E19*t8+0.18375E21*t7*t13-0.8779122094E19
*t7-0.1225E21*t20*t13+0.5852748063E19*t20-2345729468.0*t24-5251832095.0+
0.1407437681E11*t8*t24-0.3518594202E11*t7*t24+0.2345729468E11*t20*t24;
    return(-0.2E-18*beta*t32);
  }
}

#include <math.h>
double veloSource(double x)
{
  double t1;
  double t10;
  double t104;
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
  double t54;
  double t58;
  double t59;
  double t6;
  double t63;
  double t72;
  double t84;
  double t9;
  double t91;
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
    t54 = t33*t31;
    t58 = t24*t9;
    t59 = t33*t30;
    t63 = t24*t5;
    t72 = t11*t5;
    t84 = 1/t30;
    t91 = t1*t84;
    t104 = -0.18375E21*t24*t54*t6-0.1093246616E21*t58*t59*t6+0.8779122095E19*
t63*t1+0.3675E21*t11*t54*t6+0.273311654E21*t24*t59*t6-0.1755824419E20*t72*t1
-0.18375E21*t10*t54*t6-0.1822077693E21*t11*t59*t6+0.8779122094E19*t10*t5*t1
-0.5815103862E9*t6*t84+0.3518594202E11*t24*t39*t6+3489062318.0*t58*t5*t91
-0.7037188404E11*t11*t39*t6-8722655793.0*t63*t91+0.3518594202E11*t10*t39*t6+
5815103862.0*t72*t91;
    return(-0.5*t6*(2.0*t1*A*(4.0*t11+3.0*(2.0*alpha-2.0)*t10+2.0*(-3.0*alpha+
1.0)*t9)+beta*(t28*(0.35*t33*t31*t30-0.2901812679E-1*t3+0.189450439E-1)-t28*(
0.4691458936E-9*t30*t39+0.144085572E-9*t3+0.2060301114E-9))+0.1E1*t1*t3*t5)
-0.2E-18*t30*beta*t104);
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
  double t45;
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
    t45 = t4*t4;
    return(2.0*t1*A*(4.0*t8+3.0*(2.0*alpha-2.0)*t7+2.0*(-3.0*alpha+1.0)*t6)+
beta*(t25*(0.35*t30*t28*t27-0.2901812679E-1*t4+0.189450439E-1)-t25*(
0.4691458936E-9*t27*t36+0.144085572E-9*t4+0.2060301114E-9))+0.1E1*t1*t4*(1.0-
t45));
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
  double t5;
  double t50;
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
    t50 = t9*t9;
    return(1/t5*(2.0*t1*A*(4.0*t13+3.0*(2.0*alpha-2.0)*t12+2.0*(-3.0*alpha+1.0)
*t11)+beta*(t30*(0.35*t35*t33*t32-0.2901812679E-1*t9+0.189450439E-1)-t30*(
0.4691458936E-9*t32*t41+0.144085572E-9*t9+0.2060301114E-9))+0.1E1*t1*t9*(1.0-
t50)));
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
  double t32;
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
    t32 = -0.735E20*t8*t13+0.3511648838E19*t8+0.18375E21*t7*t13-0.8779122094E19
*t7-0.1225E21*t20*t13+0.5852748063E19*t20-2345729468.0*t24-5251832095.0+
0.1407437681E11*t8*t24-0.3518594202E11*t7*t24+0.2345729468E11*t20*t24;
    return(-0.2E-18*beta*t32);
  }
}

#include <math.h>
double veloSource(double x)
{
  double t1;
  double t10;
  double t104;
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
  double t54;
  double t58;
  double t59;
  double t6;
  double t63;
  double t72;
  double t84;
  double t9;
  double t91;
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
    t54 = t33*t31;
    t58 = t24*t9;
    t59 = t33*t30;
    t63 = t24*t5;
    t72 = t11*t5;
    t84 = 1/t30;
    t91 = t1*t84;
    t104 = -0.18375E21*t24*t54*t6-0.1093246616E21*t58*t59*t6+0.8779122095E19*
t63*t1+0.3675E21*t11*t54*t6+0.273311654E21*t24*t59*t6-0.1755824419E20*t72*t1
-0.18375E21*t10*t54*t6-0.1822077693E21*t11*t59*t6+0.8779122094E19*t10*t5*t1
-0.5815103862E9*t6*t84+0.3518594202E11*t24*t39*t6+3489062318.0*t58*t5*t91
-0.7037188404E11*t11*t39*t6-8722655793.0*t63*t91+0.3518594202E11*t10*t39*t6+
5815103862.0*t72*t91;
    return(-0.5*t6*(2.0*t1*A*(4.0*t11+3.0*(2.0*alpha-2.0)*t10+2.0*(-3.0*alpha+
1.0)*t9)+beta*(t28*(0.35*t33*t31*t30-0.2901812679E-1*t3+0.189450439E-1)-t28*(
0.4691458936E-9*t30*t39+0.144085572E-9*t3+0.2060301114E-9))+0.1E1*t1*t3*t5)
-0.2E-18*t30*beta*t104);
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
  double t46;
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
    t46 = t4*t4;
    return(2.0*t1*A*(4.0*t8+3.0*(2.0*alpha-2.0)*t7+2.0*(-3.0*alpha+1.0)*t6)+
beta*(t25*(0.35*t30*t28*t27-0.2901812679E-1*t4+0.189450439E-1)-t25*(av*t27*t37+
(bv-av)*t27+cv))+0.1E1*t1*t4*(1.0-t46));
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
  double t5;
  double t51;
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
    t51 = t9*t9;
    return(1/t5*(2.0*t1*A*(4.0*t13+3.0*(2.0*alpha-2.0)*t12+2.0*(-3.0*alpha+1.0)
*t11)+beta*(t30*(0.35*t35*t33*t32-0.2901812679E-1*t9+0.189450439E-1)-t30*(av*
t32*t42+(bv-av)*t32+cv))+0.1E1*t1*t9*(1.0-t51)));
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
  double t43;
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
    t43 = -0.735E11*t8*t13+3511648806.0*t8+0.18375E12*t7*t13-8779122015.0*t7
-0.1225E12*t20*t13+5852748010.0*t20-5000000000.0*av*t24-5000000000.0*bv+0.3E11*
t8*av*t24+0.3E11*t8*bv-0.75E11*t7*av*t24-0.75E11*t7*bv+0.5E11*t20*av*t24+0.5E11
*t20*bv;
    return(-0.2E-9*beta*t43);
  }
}

#include <math.h>
double veloSource(double x)
{
  double t1;
  double t10;
  double t102;
  double t11;
  double t118;
  double t24;
  double t28;
  double t3;
  double t30;
  double t31;
  double t33;
  double t4;
  double t40;
  double t5;
  double t55;
  double t59;
  double t6;
  double t60;
  double t86;
  double t9;
  double t90;
  double t92;
  double t96;
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
    t55 = t33*t31;
    t59 = t24*t9;
    t60 = t33*t30;
    t86 = 1/t30;
    t90 = av*t24;
    t92 = t40*t5*t1;
    t96 = t6*t86;
    t102 = t11*av;
    t118 = -0.18375E12*t24*t55*t6-0.1093246616E12*t59*t60*t6+8779122015.0*t24*
t5*t1+0.3675E12*t11*t55*t6+0.273311654E12*t24*t60*t6-0.1755824403E11*t11*t5*t1
-0.18375E12*t10*t55*t6-0.1822077693E12*t11*t60*t6+8779122015.0*t10*t5*t1
-1239508635.0*av*t5*t1*t86+0.75E11*t90*t92+7437051810.0*t59*av*t96+0.75E11*t24*
bv*t6-0.15E12*t102*t92-0.1859262952E11*t90*t96-0.15E12*t11*bv*t6+0.75E11*t10*av
*t92+0.1239508635E11*t102*t96+0.75E11*t10*bv*t6;
    return(-0.5*t6*(2.0*t1*A*(4.0*t11+3.0*(2.0*alpha-2.0)*t10+2.0*(-3.0*alpha+
1.0)*t9)+beta*(t28*(0.35*t33*t31*t30-0.2901812679E-1*t3+0.189450439E-1)-t28*(av
*t30*t40+(bv-av)*t30+cv))+0.1E1*t1*t3*t5)-0.2E-9*t30*beta*t118);
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
  double t45;
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
    t45 = t4*t4;
    return(2.0*t1*A*(4.0*t8+3.0*(2.0*alpha-2.0)*t7+2.0*(-3.0*alpha+1.0)*t6)+
beta*(t25*(0.35*t30*t28*t27-0.2901812679E-1*t4+0.189450439E-1)-t25*(
0.1688925222*t27*t36+0.5187080609E-1*t4+0.9217084033E-1))+0.1E1*t1*t4*(1.0-t45)
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
  double t5;
  double t50;
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
    t50 = t9*t9;
    return(1/t5*(2.0*t1*A*(4.0*t13+3.0*(2.0*alpha-2.0)*t12+2.0*(-3.0*alpha+1.0)
*t11)+beta*(t30*(0.35*t35*t33*t32-0.2901812679E-1*t9+0.189450439E-1)-t30*(
0.1688925222*t32*t41+0.5187080609E-1*t9+0.9217084033E-1))+0.1E1*t1*t9*(1.0-t50)
));
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
  double t32;
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
    t32 = -0.735E11*t8*t13+0.1485560617E11*t8+0.18375E12*t7*t13-0.3713901542E11
*t7-0.1225E12*t20*t13+0.2475934361E11*t20-844462611.0*t24-1890659560.0+
5066775666.0*t8*t24-0.1266693916E11*t7*t24+8444626110.0*t20*t24;
    return(-0.2E-9*beta*t32);
  }
}

#include <math.h>
double veloSource(double x)
{
  double t1;
  double t10;
  double t104;
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
  double t54;
  double t58;
  double t59;
  double t6;
  double t63;
  double t72;
  double t84;
  double t9;
  double t91;
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
    t54 = t33*t31;
    t58 = t24*t9;
    t59 = t33*t30;
    t63 = t24*t5;
    t72 = t11*t5;
    t84 = 1/t30;
    t91 = t1*t84;
    t104 = -0.18375E12*t24*t54*t6-0.1093246616E12*t58*t59*t6+0.3713901542E11*
t63*t1+0.3675E12*t11*t54*t6+0.273311654E12*t24*t59*t6-0.7427803084E11*t72*t1
-0.18375E12*t10*t54*t6-0.1822077693E12*t11*t59*t6+0.3713901542E11*t10*t5*t1
-0.2093437397E9*t6*t84+0.1266693916E11*t24*t39*t6+1256062438.0*t58*t5*t91
-0.2533387832E11*t11*t39*t6-3140156094.0*t63*t91+0.1266693916E11*t10*t39*t6+
2093437397.0*t72*t91;
    return(-0.5*t6*(2.0*t1*A*(4.0*t11+3.0*(2.0*alpha-2.0)*t10+2.0*(-3.0*alpha+
1.0)*t9)+beta*(t28*(0.35*t33*t31*t30-0.2901812679E-1*t3+0.189450439E-1)-t28*(
0.1688925222*t30*t39+0.5187080609E-1*t3+0.9217084033E-1))+0.1E1*t1*t3*t5)
-0.2E-9*t30*beta*t104);
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
  double t45;
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
    t45 = t4*t4;
    return(2.0*t1*A*(4.0*t8+3.0*(2.0*alpha-2.0)*t7+2.0*(-3.0*alpha+1.0)*t6)+
beta*(t25*(0.35*t30*t28*t27-0.2901812679E-1*t4+0.189450439E-1)-t25*(
0.1688925222*t27*t36+0.5187080609E-1*t4+0.9217084033E-1))+0.1E1*t1*t4*(1.0-t45)
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
  double t5;
  double t50;
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
    t50 = t9*t9;
    return(1/t5*(2.0*t1*A*(4.0*t13+3.0*(2.0*alpha-2.0)*t12+2.0*(-3.0*alpha+1.0)
*t11)+beta*(t30*(0.35*t35*t33*t32-0.2901812679E-1*t9+0.189450439E-1)-t30*(
0.1688925222*t32*t41+0.5187080609E-1*t9+0.9217084033E-1))+0.1E1*t1*t9*(1.0-t50)
));
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
  double t32;
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
    t32 = -0.735E11*t8*t13+0.1485560617E11*t8+0.18375E12*t7*t13-0.3713901542E11
*t7-0.1225E12*t20*t13+0.2475934361E11*t20-844462611.0*t24-1890659560.0+
5066775666.0*t8*t24-0.1266693916E11*t7*t24+8444626110.0*t20*t24;
    return(-0.2E-9*beta*t32);
  }
}

#include <math.h>
double veloSource(double x)
{
  double t1;
  double t10;
  double t104;
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
  double t54;
  double t58;
  double t59;
  double t6;
  double t63;
  double t72;
  double t84;
  double t9;
  double t91;
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
    t54 = t33*t31;
    t58 = t24*t9;
    t59 = t33*t30;
    t63 = t24*t5;
    t72 = t11*t5;
    t84 = 1/t30;
    t91 = t1*t84;
    t104 = -0.18375E12*t24*t54*t6-0.1093246616E12*t58*t59*t6+0.3713901542E11*
t63*t1+0.3675E12*t11*t54*t6+0.273311654E12*t24*t59*t6-0.7427803084E11*t72*t1
-0.18375E12*t10*t54*t6-0.1822077693E12*t11*t59*t6+0.3713901542E11*t10*t5*t1
-0.2093437397E9*t6*t84+0.1266693916E11*t24*t39*t6+1256062438.0*t58*t5*t91
-0.2533387832E11*t11*t39*t6-3140156094.0*t63*t91+0.1266693916E11*t10*t39*t6+
2093437397.0*t72*t91;
    return(-0.5*t6*(2.0*t1*A*(4.0*t11+3.0*(2.0*alpha-2.0)*t10+2.0*(-3.0*alpha+
1.0)*t9)+beta*(t28*(0.35*t33*t31*t30-0.2901812679E-1*t3+0.189450439E-1)-t28*(
0.1688925222*t30*t39+0.5187080609E-1*t3+0.9217084033E-1))+0.1E1*t1*t3*t5)
-0.2E-9*t30*beta*t104);
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
  double t45;
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
    t45 = t4*t4;
    return(2.0*t1*A*(4.0*t8+3.0*(2.0*alpha-2.0)*t7+2.0*(-3.0*alpha+1.0)*t6)+
beta*(t25*(0.35*t30*t28*t27-0.2901812679E-1*t4+0.189450439E-1)-t25*(
0.9382917919E-1*t27*t36+0.2881711456E-1*t4+0.5120602249E-1))+0.1E1*t1*t4*(1.0-
t45));
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
  double t5;
  double t50;
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
    t50 = t9*t9;
    return(1/t5*(2.0*t1*A*(4.0*t13+3.0*(2.0*alpha-2.0)*t12+2.0*(-3.0*alpha+1.0)
*t11)+beta*(t30*(0.35*t35*t33*t32-0.2901812679E-1*t9+0.189450439E-1)-t30*(
0.9382917919E-1*t32*t41+0.2881711456E-1*t9+0.5120602249E-1))+0.1E1*t1*t9*(1.0-
t50)));
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
  double t32;
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
    t32 = -0.147E13*t8*t13+0.1962769471E12*t8+0.3675E13*t7*t13-0.4906923676E12*
t7-0.245E13*t20*t13+0.3271282451E12*t20-9382917919.0*t24-0.2100732849E11+
0.5629750751E11*t8*t24-0.1407437688E12*t7*t24+0.9382917919E11*t20*t24;
    return(-0.1E-10*beta*t32);
  }
}

#include <math.h>
double veloSource(double x)
{
  double t1;
  double t10;
  double t104;
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
  double t54;
  double t58;
  double t59;
  double t6;
  double t63;
  double t72;
  double t84;
  double t9;
  double t91;
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
    t54 = t33*t31;
    t58 = t24*t9;
    t59 = t33*t30;
    t63 = t24*t5;
    t72 = t11*t5;
    t84 = 1/t30;
    t91 = t1*t84;
    t104 = -0.3675E13*t24*t54*t6-0.2186493232E13*t58*t59*t6+0.4906923678E12*t63
*t1+0.735E13*t11*t54*t6+0.546623308E13*t24*t59*t6-0.9813847352E12*t72*t1
-0.3675E13*t10*t54*t6-0.3644155387E13*t11*t59*t6+0.4906923676E12*t10*t5*t1
-2326041556.0*t6*t84+0.1407437688E12*t24*t39*t6+0.1395624934E11*t58*t5*t91
-0.2814875376E12*t11*t39*t6-0.3489062335E11*t63*t91+0.1407437688E12*t10*t39*t6+
0.2326041556E11*t72*t91;
    return(-0.5*t6*(2.0*t1*A*(4.0*t11+3.0*(2.0*alpha-2.0)*t10+2.0*(-3.0*alpha+
1.0)*t9)+beta*(t28*(0.35*t33*t31*t30-0.2901812679E-1*t3+0.189450439E-1)-t28*(
0.9382917919E-1*t30*t39+0.2881711456E-1*t3+0.5120602249E-1))+0.1E1*t1*t3*t5)
-0.1E-10*t30*beta*t104);
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
  double t45;
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
    t45 = t4*t4;
    return(2.0*t1*A*(4.0*t8+3.0*(2.0*alpha-2.0)*t7+2.0*(-3.0*alpha+1.0)*t6)+
beta*(t25*(0.35*t30*t28*t27-0.2901812679E-1*t4+0.189450439E-1)-t25*(
0.5670894522*t27*t36+0.1741663078*t4+0.3094815015))+0.1E1*t1*t4*(1.0-t45));
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
  double t5;
  double t50;
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
    t50 = t9*t9;
    return(1/t5*(2.0*t1*A*(4.0*t13+3.0*(2.0*alpha-2.0)*t12+2.0*(-3.0*alpha+1.0)
*t11)+beta*(t30*(0.35*t35*t33*t32-0.2901812679E-1*t9+0.189450439E-1)-t30*(
0.5670894522*t32*t41+0.1741663078*t9+0.3094815015))+0.1E1*t1*t9*(1.0-t50)));
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
  double t32;
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
    t32 = -0.735E11*t8*t13+0.4160118931E11*t8+0.18375E12*t7*t13-0.1040029733E12
*t7-0.1225E12*t20*t13+0.6933531551E11*t20-2835447261.0*t24-6348256750.0+
0.1701268357E11*t8*t24-0.4253170892E11*t7*t24+0.2835447261E11*t20*t24;
    return(-0.2E-9*beta*t32);
  }
}

#include <math.h>
double veloSource(double x)
{
  double t1;
  double t10;
  double t104;
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
  double t54;
  double t58;
  double t59;
  double t6;
  double t63;
  double t72;
  double t84;
  double t9;
  double t91;
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
    t54 = t33*t31;
    t58 = t24*t9;
    t59 = t33*t30;
    t63 = t24*t5;
    t72 = t11*t5;
    t84 = 1/t30;
    t91 = t1*t84;
    t104 = -0.18375E12*t24*t54*t6-0.1093246616E12*t58*t59*t6+0.1040029733E12*
t63*t1+0.3675E12*t11*t54*t6+0.273311654E12*t24*t59*t6-0.2080059466E12*t72*t1
-0.18375E12*t10*t54*t6-0.1822077693E12*t11*t59*t6+0.1040029733E12*t10*t5*t1
-0.7029122728E9*t6*t84+0.4253170892E11*t24*t39*t6+4217473638.0*t58*t5*t91
-0.8506341784E11*t11*t39*t6-0.1054368409E11*t63*t91+0.4253170892E11*t10*t39*t6+
7029122728.0*t72*t91;
    return(-0.5*t6*(2.0*t1*A*(4.0*t11+3.0*(2.0*alpha-2.0)*t10+2.0*(-3.0*alpha+
1.0)*t9)+beta*(t28*(0.35*t33*t31*t30-0.2901812679E-1*t3+0.189450439E-1)-t28*(
0.5670894522*t30*t39+0.1741663078*t3+0.3094815015))+0.1E1*t1*t3*t5)-0.2E-9*t30*
beta*t104);
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
  double t45;
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
    t45 = t4*t4;
    return(2.0*t1*A*(4.0*t8+3.0*(2.0*alpha-2.0)*t7+2.0*(-3.0*alpha+1.0)*t6)+
beta*(t25*(0.35*t30*t28*t27-0.2901812679E-1*t4+0.189450439E-1)-t25*(
0.1753186565*t27*t36+0.5384442083E-1*t4+0.9567781742E-1))+0.1E1*t1*t4*(1.0-t45)
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
  double t5;
  double t50;
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
    t50 = t9*t9;
    return(1/t5*(2.0*t1*A*(4.0*t13+3.0*(2.0*alpha-2.0)*t12+2.0*(-3.0*alpha+1.0)
*t11)+beta*(t30*(0.35*t35*t33*t32-0.2901812679E-1*t9+0.189450439E-1)-t30*(
0.1753186565*t32*t41+0.5384442083E-1*t9+0.9567781742E-1))+0.1E1*t1*t9*(1.0-t50)
));
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
  double t32;
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
    t32 = -0.294E11*t8*t13+6114891422.0*t8+0.735E11*t7*t13-0.1528722856E11*t7
-0.49E11*t20*t13+0.101914857E11*t20-350637313.0*t24-785038650.0+2103823878.0*t8
*t24-5259559696.0*t7*t24+3506373130.0*t20*t24;
    return(-0.5E-9*beta*t32);
  }
}

#include <math.h>
double veloSource(double x)
{
  double t1;
  double t10;
  double t104;
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
  double t54;
  double t58;
  double t59;
  double t6;
  double t63;
  double t72;
  double t84;
  double t9;
  double t91;
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
    t54 = t33*t31;
    t58 = t24*t9;
    t59 = t33*t30;
    t63 = t24*t5;
    t72 = t11*t5;
    t84 = 1/t30;
    t91 = t1*t84;
    t104 = -0.735E11*t24*t54*t6-0.4372986464E11*t58*t59*t6+0.1528722856E11*t63*
t1+0.147E12*t11*t54*t6+0.1093246616E12*t24*t59*t6-0.3057445712E11*t72*t1
-0.735E11*t10*t54*t6-0.7288310774E11*t11*t59*t6+0.1528722855E11*t10*t5*t1
-0.8692359544E8*t6*t84+5259559695.0*t24*t39*t6+0.5215415727E9*t58*t5*t91
-0.1051911939E11*t11*t39*t6-1303853932.0*t63*t91+5259559695.0*t10*t39*t6+
0.8692359544E9*t72*t91;
    return(-0.5*t6*(2.0*t1*A*(4.0*t11+3.0*(2.0*alpha-2.0)*t10+2.0*(-3.0*alpha+
1.0)*t9)+beta*(t28*(0.35*t33*t31*t30-0.2901812679E-1*t3+0.189450439E-1)-t28*(
0.1753186565*t30*t39+0.5384442083E-1*t3+0.9567781742E-1))+0.1E1*t1*t3*t5)
-0.5E-9*t30*beta*t104);
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
  double t45;
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
    t45 = t4*t4;
    return(2.0*t1*A*(4.0*t8+3.0*(2.0*alpha-2.0)*t7+2.0*(-3.0*alpha+1.0)*t6)+
beta*(t25*(0.3*t30*t28*t27-0.248726801E-1*t4+0.1623860903E-1)-t25*(0.1753186565
*t27*t36+0.5384442083E-1*t4+0.9567781742E-1))+0.1E1*t1*t4*(1.0-t45));
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
  double t5;
  double t50;
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
    t50 = t9*t9;
    return(1/t5*(2.0*t1*A*(4.0*t13+3.0*(2.0*alpha-2.0)*t12+2.0*(-3.0*alpha+1.0)
*t11)+beta*(t30*(0.3*t35*t33*t32-0.248726801E-1*t9+0.1623860903E-1)-t30*(
0.1753186565*t32*t41+0.5384442083E-1*t9+0.9567781742E-1))+0.1E1*t1*t9*(1.0-t50)
));
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
  double t32;
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
    t32 = -0.252E11*t8*t13+5914225776.0*t8+0.63E11*t7*t13-0.1478556444E11*t7
-0.42E11*t20*t13+9857042960.0*t20-350637313.0*t24-785038650.0+2103823878.0*t8*
t24-5259559695.0*t7*t24+3506373130.0*t20*t24;
    return(-0.5E-9*beta*t32);
  }
}

#include <math.h>
double veloSource(double x)
{
  double t1;
  double t10;
  double t104;
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
  double t54;
  double t58;
  double t59;
  double t6;
  double t63;
  double t72;
  double t84;
  double t9;
  double t91;
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
    t54 = t33*t31;
    t58 = t24*t9;
    t59 = t33*t30;
    t63 = t24*t5;
    t72 = t11*t5;
    t84 = 1/t30;
    t91 = t1*t84;
    t104 = -0.63E11*t24*t54*t6-0.3748274112E11*t58*t59*t6+0.1478556444E11*t63*
t1+0.126E12*t11*t54*t6+0.9370685281E11*t24*t59*t6-0.2957112888E11*t72*t1
-0.63E11*t10*t54*t6-0.624712352E11*t11*t59*t6+0.1478556444E11*t10*t5*t1
-0.8692359544E8*t6*t84+5259559695.0*t24*t39*t6+0.5215415727E9*t58*t5*t91
-0.1051911939E11*t11*t39*t6-1303853932.0*t63*t91+5259559695.0*t10*t39*t6+
0.8692359544E9*t72*t91;
    return(-0.5*t6*(2.0*t1*A*(4.0*t11+3.0*(2.0*alpha-2.0)*t10+2.0*(-3.0*alpha+
1.0)*t9)+beta*(t28*(0.3*t33*t31*t30-0.248726801E-1*t3+0.1623860903E-1)-t28*(
0.1753186565*t30*t39+0.5384442083E-1*t3+0.9567781742E-1))+0.1E1*t1*t3*t5)
-0.5E-9*t30*beta*t104);
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
  double t45;
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
    t45 = t4*t4;
    return(2.0*t1*A*(4.0*t8+3.0*(2.0*alpha-2.0)*t7+2.0*(-3.0*alpha+1.0)*t6)+
beta*(t25*(0.35*t30*t28*t27-0.2901812679E-1*t4+0.189450439E-1)-t25*(
0.1753186565*t27*t36+0.5384442083E-1*t4+0.9567781742E-1))+0.1E1*t1*t4*(1.0-t45)
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
  double t5;
  double t50;
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
    t50 = t9*t9;
    return(1/t5*(2.0*t1*A*(4.0*t13+3.0*(2.0*alpha-2.0)*t12+2.0*(-3.0*alpha+1.0)
*t11)+beta*(t30*(0.35*t35*t33*t32-0.2901812679E-1*t9+0.189450439E-1)-t30*(
0.1753186565*t32*t41+0.5384442083E-1*t9+0.9567781742E-1))+0.1E1*t1*t9*(1.0-t50)
));
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
  double t32;
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
    t32 = -0.294E11*t8*t13+6114891422.0*t8+0.735E11*t7*t13-0.1528722856E11*t7
-0.49E11*t20*t13+0.101914857E11*t20-350637313.0*t24-785038650.0+2103823878.0*t8
*t24-5259559696.0*t7*t24+3506373130.0*t20*t24;
    return(-0.5E-9*beta*t32);
  }
}

#include <math.h>
double veloSource(double x)
{
  double t1;
  double t10;
  double t105;
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
  double t54;
  double t58;
  double t59;
  double t6;
  double t63;
  double t84;
  double t9;
  double t91;
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
    t54 = t33*t31;
    t58 = t24*t9;
    t59 = t33*t30;
    t63 = t24*t5;
    t84 = 1/t30;
    t91 = t1*t84;
    t105 = -0.735E11*t24*t54*t6-0.4372986464E11*t58*t59*t6+0.1528722856E11*t63*
t1+0.147E12*t11*t54*t6+0.1093246616E12*t24*t59*t6-0.3057445712E11*t11*t5*t1
-0.735E11*t10*t54*t6-0.7288310774E11*t11*t59*t6+0.1528722855E11*t10*t5*t1
-0.8692359544E8*t6*t84+5259559695.0*t24*t39*t6+0.5215415727E9*t58*t5*t91
-0.1051911939E11*t11*t39*t6-1303853932.0*t63*t91+0.8692359544E9*t6*t84*t11+
5259559695.0*t39*t10*t6;
    return(-0.5*t6*(2.0*t1*A*(4.0*t11+3.0*(2.0*alpha-2.0)*t10+2.0*(-3.0*alpha+
1.0)*t9)+beta*(t28*(0.35*t33*t31*t30-0.2901812679E-1*t3+0.189450439E-1)-t28*(
0.1753186565*t30*t39+0.5384442083E-1*t3+0.9567781742E-1))+0.1E1*t1*t3*t5)
-0.5E-9*t30*beta*t105);
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
  double t45;
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
    t45 = t4*t4;
    return(2.0*t1*A*(4.0*t8+3.0*(2.0*alpha-2.0)*t7+2.0*(-3.0*alpha+1.0)*t6)+
beta*(t25*(0.35*t30*t28*t27-0.2901812679E-1*t4+0.189450439E-1)-t25*(
0.1753186565*t27*t36+0.5384442083E-1*t4+0.9567781742E-1))+0.1E1*t1*t4*(1.0-t45)
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
  double t5;
  double t50;
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
    t50 = t9*t9;
    return(1/t5*(2.0*t1*A*(4.0*t13+3.0*(2.0*alpha-2.0)*t12+2.0*(-3.0*alpha+1.0)
*t11)+beta*(t30*(0.35*t35*t33*t32-0.2901812679E-1*t9+0.189450439E-1)-t30*(
0.1753186565*t32*t41+0.5384442083E-1*t9+0.9567781742E-1))+0.1E1*t1*t9*(1.0-t50)
));
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
  double t32;
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
    t32 = -0.294E11*t8*t13+6114891422.0*t8+0.735E11*t7*t13-0.1528722856E11*t7
-0.49E11*t20*t13+0.101914857E11*t20-350637313.0*t24-785038650.0+2103823878.0*t8
*t24-5259559696.0*t7*t24+3506373130.0*t20*t24;
    return(-0.5E-9*beta*t32);
  }
}

#include <math.h>
double veloSource(double x)
{
  double t1;
  double t10;
  double t104;
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
  double t54;
  double t58;
  double t59;
  double t6;
  double t63;
  double t72;
  double t84;
  double t9;
  double t91;
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
    t54 = t33*t31;
    t58 = t24*t9;
    t59 = t33*t30;
    t63 = t24*t5;
    t72 = t11*t5;
    t84 = 1/t30;
    t91 = t1*t84;
    t104 = -0.735E11*t24*t54*t6-0.4372986464E11*t58*t59*t6+0.1528722856E11*t63*
t1+0.147E12*t11*t54*t6+0.1093246616E12*t24*t59*t6-0.3057445712E11*t72*t1
-0.735E11*t10*t54*t6-0.7288310774E11*t11*t59*t6+0.1528722855E11*t10*t5*t1
-0.8692359544E8*t6*t84+5259559695.0*t24*t39*t6+0.5215415727E9*t58*t5*t91
-0.1051911939E11*t11*t39*t6-1303853932.0*t63*t91+5259559695.0*t10*t39*t6+
0.8692359544E9*t72*t91;
    return(-0.5*t6*(2.0*t1*A*(4.0*t11+3.0*(2.0*alpha-2.0)*t10+2.0*(-3.0*alpha+
1.0)*t9)+beta*(t28*(0.35*t33*t31*t30-0.2901812679E-1*t3+0.189450439E-1)-t28*(
0.1753186565*t30*t39+0.5384442083E-1*t3+0.9567781742E-1))+0.1E1*t1*t3*t5)
-0.5E-9*t30*beta*t104);
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
  double t45;
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
    t45 = t4*t4;
    return(2.0*t1*A*(4.0*t8+3.0*(2.0*alpha-2.0)*t7+2.0*(-3.0*alpha+1.0)*t6)+
beta*(t25*(0.35*t30*t28*t27-0.2901812679E-1*t4+0.189450439E-1)-t25*(
0.1753186565*t27*t36+0.5384442083E-1*t4+0.9567781742E-1))+0.1E1*t1*t4*(1.0-t45)
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
  double t5;
  double t50;
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
    t50 = t9*t9;
    return(1/t5*(2.0*t1*A*(4.0*t13+3.0*(2.0*alpha-2.0)*t12+2.0*(-3.0*alpha+1.0)
*t11)+beta*(t30*(0.35*t35*t33*t32-0.2901812679E-1*t9+0.189450439E-1)-t30*(
0.1753186565*t32*t41+0.5384442083E-1*t9+0.9567781742E-1))+0.1E1*t1*t9*(1.0-t50)
));
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
  double t32;
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
    t32 = -0.294E11*t8*t13+6114891422.0*t8+0.735E11*t7*t13-0.1528722856E11*t7
-0.49E11*t20*t13+0.101914857E11*t20-350637313.0*t24-785038650.0+2103823878.0*t8
*t24-5259559696.0*t7*t24+3506373130.0*t20*t24;
    return(-0.5E-9*beta*t32);
  }
}

#include <math.h>
double veloSource(double x)
{
  double t1;
  double t10;
  double t104;
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
  double t54;
  double t58;
  double t59;
  double t6;
  double t63;
  double t72;
  double t84;
  double t9;
  double t91;
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
    t54 = t33*t31;
    t58 = t24*t9;
    t59 = t33*t30;
    t63 = t24*t5;
    t72 = t11*t5;
    t84 = 1/t30;
    t91 = t1*t84;
    t104 = -0.735E11*t24*t54*t6-0.4372986464E11*t58*t59*t6+0.1528722856E11*t63*
t1+0.147E12*t11*t54*t6+0.1093246616E12*t24*t59*t6-0.3057445712E11*t72*t1
-0.735E11*t10*t54*t6-0.7288310774E11*t11*t59*t6+0.1528722855E11*t10*t5*t1
-0.8692359544E8*t6*t84+5259559695.0*t24*t39*t6+0.5215415727E9*t58*t5*t91
-0.1051911939E11*t11*t39*t6-1303853932.0*t63*t91+5259559695.0*t10*t39*t6+
0.8692359544E9*t72*t91;
    return(-0.5*t6*(2.0*t1*A*(4.0*t11+3.0*(2.0*alpha-2.0)*t10+2.0*(-3.0*alpha+
1.0)*t9)+beta*(t28*(0.35*t33*t31*t30-0.2901812679E-1*t3+0.189450439E-1)-t28*(
0.1753186565*t30*t39+0.5384442083E-1*t3+0.9567781742E-1))+0.1E1*t1*t3*t5)
-0.5E-9*t30*beta*t104);
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
    t102 = 0.3675E2*t24*t51*t6+0.2186493232E2*t55*t56*t6-0.7643614278E1*t60*t1
-0.735E2*t11*t51*t6-0.546623308E2*t24*t56*t6+0.1528722856E2*t11*t5*t1+0.3675E2*
t10*t51*t6+0.3644155387E2*t11*t56*t6-0.7643614278E1*t10*t5*t1+0.4346179772E-1*
t6*t81-0.2629779848E1*t24*t39*t6-0.2607707863*t55*t5*t88+0.5259559696E1*t11*t39
*t6+0.6519269659*t60*t88-0.4346179772*t6*t81*t11-0.2629779848E1*t39*t10*t6;
    return(-0.5*t6*(2.0*t1*A*(4.0*t11+3.0*(2.0*alpha-2.0)*t10+2.0*(-3.0*alpha+
1.0)*t9)+t28*(0.35*t33*t31*t30-0.2901812679E-1*t3+0.189450439E-1)-t28*(
0.1753186565*t30*t39+0.5384442083E-1*t3+0.9567781742E-1)+0.1E1*t1*t3*t5)+t30*
t102);
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


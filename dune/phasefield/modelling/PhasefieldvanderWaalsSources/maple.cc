double evalRho(double x)
{
  {
    return(0.4958034541*x+0.1065766549);
  }
}


inline double helmholtz ( double rho ,double phi ) const
{
  double t1;
  double t16;
  double t17;
  double t18;
  double t2;
  double t20;
  double t21;
  double t5;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t5 = t1*phi;
    t16 = 6.0*t2*phi;
    t17 = 15.0*t2;
    t18 = 10.0*t5;
    t20 = log(rho);
    t21 = rho*t20;
    return(2.0*0*(t2+(2.0*alpha_-2.0)*t5+(-3.0*alpha_+1.0)*t1+alpha_)/delta_+(t16-
t17+t18)*(2.0*t21-0.9862667552*rho+0.1204760218E1)+(1.0-t16+t17-t18)*(
0.1753186565*t21+0.2172006685*rho+0.186848759E-1));
  }
}


inline double reactionSource ( double rho ,double phi ) const
{
  double t1;
  double t17;
  double t2;
  double t21;
  double t22;
  double t23;
  {
    t1 = phi*phi;
    t2 = t1*phi;
    t17 = t1*t1;
    t21 = 30.0*t17-60.0*t2+30.0*t1;
    t22 = log(rho);
    t23 = rho*t22;
    return(2.0*A_*(4.0*t2+3.0*(2.0*alpha_-2.0)*t1+2.0*(-3.0*alpha_+1.0)*phi)/delta_
+t21*(2.0*t23-0.9862667552*rho+0.1204760218E1)-t21*(0.1753186565*t23+
0.2172006685*rho+0.186848759E-1));
  }
}


inline double dphireactionSource ( double rho ,double phi ) const
{
  double t1;
  double t17;
  double t18;
  double t19;
  {
    t1 = phi*phi;
    t17 = 120.0*t1*phi-180.0*t1+60.0*phi;
    t18 = log(rho);
    t19 = rho*t18;
    return(2.0*A_*(2.0+12.0*t1+6.0*(2.0*alpha_-2.0)*phi-6.0*alpha_)/delta_+t17*(
0.1204760218E1+2.0*t19-0.9862667552*rho)-t17*(0.186848759E-1+0.1753186565*t19+
0.2172006685*rho));
  }
}


inline double chemicalPotential ( double rho ,double phi ) const
{
  double t1;
  double t11;
  double t2;
  double t3;
  double t4;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = t2*phi;
    t4 = log(rho);
    t11 = t1*phi;
    return(0.1094808806E2*t3*t4+0.372728352E1*t3-0.2737022015E2*t2*t4
-0.93182088E1*t2+0.1824681344E2*t11*t4+0.62121392E1*t11+0.1753186565*t4+
0.392519325);
  }
}

inline double drhochemicalPotential ( double rho ,double phi ) const
{
  double t1;
  double t2;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    return(0.5E-9*(0.2189617612E11*t2*phi-0.547404403E11*t2+0.3649362687E11*t1*
phi+350637313.0)/rho);
  }
}


inline double dphichemicalPotential ( double rho ,double phi ) const
{
  double t1;
  double t2;
  double t3;
  double t7;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = log(rho);
    t7 = t1*phi;
    return(0.547404403E2*t2*t3+0.186364176E2*t2-0.1094808806E3*t7*t3
-0.372728352E2*t7+0.547404403E2*t1*t3+0.186364176E2*t1);
  }
}


inline double pressure ( double rho ,double phi ) const
{
  double t1;
  double t10;
  double t2;
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
    t9 = log(rho);
    t10 = rho*t9;
    return((t4-t5+t7)*(-2.0*t10+0.9862667552*rho-0.1204760218E1+rho*(2.0*t9+
0.1013733245E1))+(1.0-t4+t5-t7)*(-0.1753186565*t10-0.2172006685*rho
-0.186848759E-1+rho*(0.1753186565*t9+0.392519325))-2.0*A_*(t2+(2.0*alpha_-2.0)*t6
+(-3.0*alpha_+1.0)*t1+alpha_)/delta_);
  }
}


inline double a ( double rho ,double phi ) const
{
  double t1;
  double t2;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t9 = sqrt(0.1094808806E12*t2*phi-0.2737022015E12*t2+0.1824681344E12*t1*phi+
1753186565.0);
    return(0.1E-4*t9);
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
  double t21;
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
    t20 = log(rho);
    t21 = rho*t20;
    return(2.0/delta*A*(t4+(2.0*alpha-2.0)*t7+(-3.0*alpha+1.0)*t3+alpha)+(t16-
t17+t18)*(2.0*t21-0.9862667552*rho+0.1204760218E1)+(1.0-t16+t17-t18)*(
0.1753186565*t21+0.2172006685*rho+0.186848759E-1));
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t17;
  double t21;
  double t22;
  double t23;
  double t3;
  double t4;
  {
    t3 = phi*phi;
    t4 = t3*phi;
    t17 = t3*t3;
    t21 = 30.0*t17-60.0*t4+30.0*t3;
    t22 = log(rho);
    t23 = rho*t22;
    return(2.0/delta*A*(4.0*t4+3.0*(2.0*alpha-2.0)*t3+2.0*(-3.0*alpha+1.0)*phi)
+t21*(2.0*t23-0.9862667552*rho+0.1204760218E1)-t21*(0.1753186565*t23+
0.2172006685*rho+0.186848759E-1));
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  double t17;
  double t18;
  double t19;
  double t3;
  {
    t3 = phi*phi;
    t17 = 120.0*t3*phi-180.0*t3+60.0*phi;
    t18 = log(rho);
    t19 = rho*t18;
    return(2.0/delta*A*(12.0*t3+6.0*(2.0*alpha-2.0)*phi-6.0*alpha+2.0)+t17*(2.0
*t19-0.9862667552*rho+0.1204760218E1)-t17*(0.1753186565*t19+0.2172006685*rho+
0.186848759E-1));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t11;
  double t2;
  double t3;
  double t4;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = t2*phi;
    t4 = log(rho);
    t11 = t1*phi;
    return(0.1094808806E2*t3*t4+0.372728352E1*t3-0.2737022015E2*t2*t4
-0.93182088E1*t2+0.1824681344E2*t11*t4+0.62121392E1*t11+0.1753186565*t4+
0.392519325);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  double t1;
  double t2;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    return(0.5E-9*(0.2189617612E11*t2*phi-0.547404403E11*t2+0.3649362687E11*t1*
phi+350637313.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t2;
  double t3;
  double t7;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = log(rho);
    t7 = t1*phi;
    return(0.547404403E2*t2*t3+0.186364176E2*t2-0.1094808806E3*t7*t3
-0.372728352E2*t7+0.547404403E2*t1*t3+0.186364176E2*t1);
  }
}

#include <math.h>
double pressure(double rho,double phi)
{
  double t1;
  double t10;
  double t2;
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
    t9 = log(rho);
    t10 = rho*t9;
    return((t4-t5+t7)*(-2.0*t10+0.9862667552*rho-0.1204760218E1+rho*(2.0*t9+
0.1013733245E1))+(1.0-t4+t5-t7)*(-0.1753186565*t10-0.2172006685*rho
-0.186848759E-1+rho*(0.1753186565*t9+0.392519325))-2.0/delta*A*(t2+(2.0*alpha
-2.0)*t6+(-3.0*alpha+1.0)*t1+alpha));
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t1;
  double t2;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t9 = sqrt(0.1094808806E12*t2*phi-0.2737022015E12*t2+0.1824681344E12*t1*phi+
1753186565.0);
    return(0.1E-4*t9);
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t16;
  double t17;
  double t18;
  double t20;
  double t21;
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
    t20 = log(rho);
    t21 = rho*t20;
    return(2.0/delta*A*(t4+(2.0*alpha-2.0)*t7+(-3.0*alpha+1.0)*t3+alpha)+(t16-
t17+t18)*(2.0*t21-0.986266755*rho+0.1204760218E1)+(1.0-t16+t17-t18)*(
0.1753186565*t21+0.2172006685*rho+0.186848759E-1));
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t17;
  double t21;
  double t22;
  double t23;
  double t3;
  double t4;
  {
    t3 = phi*phi;
    t4 = t3*phi;
    t17 = t3*t3;
    t21 = 30.0*t17-60.0*t4+30.0*t3;
    t22 = log(rho);
    t23 = rho*t22;
    return(2.0/delta*A*(4.0*t4+3.0*(2.0*alpha-2.0)*t3+2.0*(-3.0*alpha+1.0)*phi)
+t21*(2.0*t23-0.986266755*rho+0.1204760218E1)-t21*(0.1753186565*t23+
0.2172006685*rho+0.186848759E-1));
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  double t17;
  double t18;
  double t19;
  double t3;
  {
    t3 = phi*phi;
    t17 = 120.0*t3*phi-180.0*t3+60.0*phi;
    t18 = log(rho);
    t19 = rho*t18;
    return(2.0/delta*A*(12.0*t3+6.0*(2.0*alpha-2.0)*phi-6.0*alpha+2.0)+t17*(2.0
*t19-0.986266755*rho+0.1204760218E1)-t17*(0.1753186565*t19+0.2172006685*rho+
0.186848759E-1));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t11;
  double t2;
  double t3;
  double t4;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = t2*phi;
    t4 = log(rho);
    t11 = t1*phi;
    return(0.1094808806E2*t3*t4+0.372728352E1*t3-0.2737022015E2*t2*t4
-0.93182088E1*t2+0.1824681344E2*t11*t4+0.62121392E1*t11+0.1753186565*t4+
0.392519325);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  double t1;
  double t2;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    return(0.5E-9*(0.2189617612E11*t2*phi-0.547404403E11*t2+0.3649362687E11*t1*
phi+350637313.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t2;
  double t3;
  double t7;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = log(rho);
    t7 = t1*phi;
    return(0.547404403E2*t2*t3+0.186364176E2*t2-0.1094808806E3*t7*t3
-0.372728352E2*t7+0.547404403E2*t1*t3+0.186364176E2*t1);
  }
}

#include <math.h>
double pressure(double rho,double phi)
{
  double t1;
  double t10;
  double t2;
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
    t9 = log(rho);
    t10 = rho*t9;
    return((t4-t5+t7)*(-2.0*t10+0.986266755*rho-0.1204760218E1+rho*(2.0*t9+
0.1013733245E1))+(1.0-t4+t5-t7)*(-0.1753186565*t10-0.2172006685*rho
-0.186848759E-1+rho*(0.1753186565*t9+0.392519325))-2.0/delta*A*(t2+(2.0*alpha
-2.0)*t6+(-3.0*alpha+1.0)*t1+alpha));
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t1;
  double t2;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t9 = sqrt(0.1094808806E12*t2*phi-0.2737022015E12*t2+0.1824681344E12*t1*phi+
1753186565.0);
    return(0.1E-4*t9);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(solrho(x));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t16;
  double t17;
  double t18;
  double t20;
  double t21;
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
    t20 = log(rho);
    t21 = rho*t20;
    return(2.0/delta*A*(t4+(2.0*alpha-2.0)*t7+(-3.0*alpha+1.0)*t3+alpha)+(t16-
t17+t18)*(2.0*t21-0.986266755*rho+0.1204760218E1)+(1.0-t16+t17-t18)*(
0.1753186565*t21+0.2172006685*rho+0.186848759E-1));
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t17;
  double t21;
  double t22;
  double t23;
  double t3;
  double t4;
  {
    t3 = phi*phi;
    t4 = t3*phi;
    t17 = t3*t3;
    t21 = 30.0*t17-60.0*t4+30.0*t3;
    t22 = log(rho);
    t23 = rho*t22;
    return(2.0/delta*A*(4.0*t4+3.0*(2.0*alpha-2.0)*t3+2.0*(-3.0*alpha+1.0)*phi)
+t21*(2.0*t23-0.986266755*rho+0.1204760218E1)-t21*(0.1753186565*t23+
0.2172006685*rho+0.186848759E-1));
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  double t17;
  double t18;
  double t19;
  double t3;
  {
    t3 = phi*phi;
    t17 = 120.0*t3*phi-180.0*t3+60.0*phi;
    t18 = log(rho);
    t19 = rho*t18;
    return(2.0/delta*A*(12.0*t3+6.0*(2.0*alpha-2.0)*phi-6.0*alpha+2.0)+t17*(2.0
*t19-0.986266755*rho+0.1204760218E1)-t17*(0.1753186565*t19+0.2172006685*rho+
0.186848759E-1));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t11;
  double t2;
  double t3;
  double t4;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = t2*phi;
    t4 = log(rho);
    t11 = t1*phi;
    return(0.1094808806E2*t3*t4+0.372728352E1*t3-0.2737022015E2*t2*t4
-0.93182088E1*t2+0.1824681344E2*t11*t4+0.62121392E1*t11+0.1753186565*t4+
0.392519325);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  double t1;
  double t2;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    return(0.5E-9*(0.2189617612E11*t2*phi-0.547404403E11*t2+0.3649362687E11*t1*
phi+350637313.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t2;
  double t3;
  double t7;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = log(rho);
    t7 = t1*phi;
    return(0.547404403E2*t2*t3+0.186364176E2*t2-0.1094808806E3*t7*t3
-0.372728352E2*t7+0.547404403E2*t1*t3+0.186364176E2*t1);
  }
}

#include <math.h>
double pressure(double rho,double phi)
{
  double t1;
  double t10;
  double t2;
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
    t9 = log(rho);
    t10 = rho*t9;
    return((t4-t5+t7)*(-2.0*t10+0.986266755*rho-0.1204760218E1+rho*(2.0*t9+
0.1013733245E1))+(1.0-t4+t5-t7)*(-0.1753186565*t10-0.2172006685*rho
-0.186848759E-1+rho*(0.1753186565*t9+0.392519325))-2.0/delta*A*(t2+(2.0*alpha
-2.0)*t6+(-3.0*alpha+1.0)*t1+alpha));
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t1;
  double t2;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t9 = sqrt(0.1094808806E12*t2*phi-0.2737022015E12*t2+0.1824681344E12*t1*phi+
1753186565.0);
    return(0.1E-4*t9);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(solrho(x));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t16;
  double t17;
  double t18;
  double t20;
  double t21;
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
    t20 = log(rho);
    t21 = rho*t20;
    return(2.0/delta*A*(t4+(2.0*alpha-2.0)*t7+(-3.0*alpha+1.0)*t3+alpha)+(t16-
t17+t18)*(2.0*t21-0.986266755*rho+0.1204760218E1)+(1.0-t16+t17-t18)*(
0.1753186565*t21+0.2172006685*rho+0.186848759E-1));
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t17;
  double t21;
  double t22;
  double t23;
  double t3;
  double t4;
  {
    t3 = phi*phi;
    t4 = t3*phi;
    t17 = t3*t3;
    t21 = 30.0*t17-60.0*t4+30.0*t3;
    t22 = log(rho);
    t23 = rho*t22;
    return(2.0/delta*A*(4.0*t4+3.0*(2.0*alpha-2.0)*t3+2.0*(-3.0*alpha+1.0)*phi)
+t21*(2.0*t23-0.986266755*rho+0.1204760218E1)-t21*(0.1753186565*t23+
0.2172006685*rho+0.186848759E-1));
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  double t17;
  double t18;
  double t19;
  double t3;
  {
    t3 = phi*phi;
    t17 = 120.0*t3*phi-180.0*t3+60.0*phi;
    t18 = log(rho);
    t19 = rho*t18;
    return(2.0/delta*A*(12.0*t3+6.0*(2.0*alpha-2.0)*phi-6.0*alpha+2.0)+t17*(2.0
*t19-0.986266755*rho+0.1204760218E1)-t17*(0.1753186565*t19+0.2172006685*rho+
0.186848759E-1));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t11;
  double t2;
  double t3;
  double t4;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = t2*phi;
    t4 = log(rho);
    t11 = t1*phi;
    return(0.1094808806E2*t3*t4+0.372728352E1*t3-0.2737022015E2*t2*t4
-0.93182088E1*t2+0.1824681344E2*t11*t4+0.62121392E1*t11+0.1753186565*t4+
0.392519325);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  double t1;
  double t2;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    return(0.5E-9*(0.2189617612E11*t2*phi-0.547404403E11*t2+0.3649362687E11*t1*
phi+350637313.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t2;
  double t3;
  double t7;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = log(rho);
    t7 = t1*phi;
    return(0.547404403E2*t2*t3+0.186364176E2*t2-0.1094808806E3*t7*t3
-0.372728352E2*t7+0.547404403E2*t1*t3+0.186364176E2*t1);
  }
}

#include <math.h>
double pressure(double rho,double phi)
{
  double t1;
  double t10;
  double t2;
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
    t9 = log(rho);
    t10 = rho*t9;
    return((t4-t5+t7)*(-2.0*t10+0.986266755*rho-0.1204760218E1+rho*(2.0*t9+
0.1013733245E1))+(1.0-t4+t5-t7)*(-0.1753186565*t10-0.2172006685*rho
-0.186848759E-1+rho*(0.1753186565*t9+0.392519325))-2.0/delta*A*(t2+(2.0*alpha
-2.0)*t6+(-3.0*alpha+1.0)*t1+alpha));
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t1;
  double t2;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t9 = sqrt(0.1094808806E12*t2*phi-0.2737022015E12*t2+0.1824681344E12*t1*phi+
1753186565.0);
    return(0.1E-4*t9);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(solrho(x));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t16;
  double t17;
  double t18;
  double t20;
  double t21;
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
    t20 = log(rho);
    t21 = rho*t20;
    return(2.0/delta*A*(t4+(2.0*alpha-2.0)*t7+(-3.0*alpha+1.0)*t3+alpha)+(t16-
t17+t18)*(2.0*t21-0.986266755*rho+0.1204760218E1)+(1.0-t16+t17-t18)*(
0.1753186565*t21+0.2172006685*rho+0.186848759E-1));
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t17;
  double t21;
  double t22;
  double t23;
  double t3;
  double t4;
  {
    t3 = phi*phi;
    t4 = t3*phi;
    t17 = t3*t3;
    t21 = 30.0*t17-60.0*t4+30.0*t3;
    t22 = log(rho);
    t23 = rho*t22;
    return(2.0/delta*A*(4.0*t4+3.0*(2.0*alpha-2.0)*t3+2.0*(-3.0*alpha+1.0)*phi)
+t21*(2.0*t23-0.986266755*rho+0.1204760218E1)-t21*(0.1753186565*t23+
0.2172006685*rho+0.186848759E-1));
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  double t17;
  double t18;
  double t19;
  double t3;
  {
    t3 = phi*phi;
    t17 = 120.0*t3*phi-180.0*t3+60.0*phi;
    t18 = log(rho);
    t19 = rho*t18;
    return(2.0/delta*A*(12.0*t3+6.0*(2.0*alpha-2.0)*phi-6.0*alpha+2.0)+t17*(2.0
*t19-0.986266755*rho+0.1204760218E1)-t17*(0.1753186565*t19+0.2172006685*rho+
0.186848759E-1));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t11;
  double t2;
  double t3;
  double t4;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = t2*phi;
    t4 = log(rho);
    t11 = t1*phi;
    return(0.1094808806E2*t3*t4+0.372728352E1*t3-0.2737022015E2*t2*t4
-0.93182088E1*t2+0.1824681344E2*t11*t4+0.62121392E1*t11+0.1753186565*t4+
0.392519325);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  double t1;
  double t2;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    return(0.5E-9*(0.2189617612E11*t2*phi-0.547404403E11*t2+0.3649362687E11*t1*
phi+350637313.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t2;
  double t3;
  double t7;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = log(rho);
    t7 = t1*phi;
    return(0.547404403E2*t2*t3+0.186364176E2*t2-0.1094808806E3*t7*t3
-0.372728352E2*t7+0.547404403E2*t1*t3+0.186364176E2*t1);
  }
}

#include <math.h>
double pressure(double rho,double phi)
{
  double t1;
  double t10;
  double t2;
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
    t9 = log(rho);
    t10 = rho*t9;
    return((t4-t5+t7)*(-2.0*t10+0.986266755*rho-0.1204760218E1+rho*(2.0*t9+
0.1013733245E1))+(1.0-t4+t5-t7)*(-0.1753186565*t10-0.2172006685*rho
-0.186848759E-1+rho*(0.1753186565*t9+0.392519325))-2.0/delta*A*(t2+(2.0*alpha
-2.0)*t6+(-3.0*alpha+1.0)*t1+alpha));
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t1;
  double t2;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t9 = sqrt(0.1094808806E12*t2*phi-0.2737022015E12*t2+0.1824681344E12*t1*phi+
1753186565.0);
    return(0.1E-4*t9);
  }
}

#include <math.h>
double evalRho(double x)
{
  double t1;
  double t3;
  {
    t1 = log(sol1);
    t3 = nu(x);
    return(exp((0.62121392+a*t1-0.62121392*t3)/(0.1753186565+t3*(a-0.1753186565
))));
  }
}

#include <math.h>
double evalRho(double x)
{
  double t1;
  double t3;
  double t4;
  double t5;
  double t8;
  {
    t1 = log(sol1);
    t3 = x*x;
    t4 = t3*t3;
    t5 = t4*x;
    t8 = t3*x;
    return(exp((2.0*t1+0.62121392-0.372728352E1*t5+0.93182088E1*t4-0.62121392E1
*t8)/(0.1094808806E2*t5-0.2737022016E2*t4+0.1824681344E2*t8+0.1753186565)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t16;
  double t17;
  double t18;
  double t20;
  double t21;
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
    t20 = log(rho);
    t21 = rho*t20;
    return(2.0/delta*A*(t4+(2.0*alpha-2.0)*t7+(-3.0*alpha+1.0)*t3+alpha)+(t16-
t17+t18)*(2.0*t21-0.986266755*rho+0.1204760218E1)+(1.0-t16+t17-t18)*(
0.1753186565*t21+0.2172006685*rho+0.186848759E-1));
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t17;
  double t21;
  double t22;
  double t23;
  double t3;
  double t4;
  {
    t3 = phi*phi;
    t4 = t3*phi;
    t17 = t3*t3;
    t21 = 30.0*t17-60.0*t4+30.0*t3;
    t22 = log(rho);
    t23 = rho*t22;
    return(2.0/delta*A*(4.0*t4+3.0*(2.0*alpha-2.0)*t3+2.0*(-3.0*alpha+1.0)*phi)
+t21*(2.0*t23-0.986266755*rho+0.1204760218E1)-t21*(0.1753186565*t23+
0.2172006685*rho+0.186848759E-1));
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  double t17;
  double t18;
  double t19;
  double t3;
  {
    t3 = phi*phi;
    t17 = 120.0*t3*phi-180.0*t3+60.0*phi;
    t18 = log(rho);
    t19 = rho*t18;
    return(2.0/delta*A*(12.0*t3+6.0*(2.0*alpha-2.0)*phi-6.0*alpha+2.0)+t17*(2.0
*t19-0.986266755*rho+0.1204760218E1)-t17*(0.1753186565*t19+0.2172006685*rho+
0.186848759E-1));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t11;
  double t2;
  double t3;
  double t4;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = t2*phi;
    t4 = log(rho);
    t11 = t1*phi;
    return(0.1094808806E2*t3*t4+0.372728352E1*t3-0.2737022015E2*t2*t4
-0.93182088E1*t2+0.1824681344E2*t11*t4+0.62121392E1*t11+0.1753186565*t4+
0.392519325);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  double t1;
  double t2;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    return(0.5E-9*(0.2189617612E11*t2*phi-0.547404403E11*t2+0.3649362687E11*t1*
phi+350637313.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t2;
  double t3;
  double t7;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = log(rho);
    t7 = t1*phi;
    return(0.547404403E2*t2*t3+0.186364176E2*t2-0.1094808806E3*t7*t3
-0.372728352E2*t7+0.547404403E2*t1*t3+0.186364176E2*t1);
  }
}

#include <math.h>
double pressure(double rho,double phi)
{
  double t1;
  double t10;
  double t2;
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
    t9 = log(rho);
    t10 = rho*t9;
    return((t4-t5+t7)*(-2.0*t10+0.986266755*rho-0.1204760218E1+rho*(2.0*t9+
0.1013733245E1))+(1.0-t4+t5-t7)*(-0.1753186565*t10-0.2172006685*rho
-0.186848759E-1+rho*(0.1753186565*t9+0.392519325))-2.0/delta*A*(t2+(2.0*alpha
-2.0)*t6+(-3.0*alpha+1.0)*t1+alpha));
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t1;
  double t2;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t9 = sqrt(0.1094808806E12*t2*phi-0.2737022015E12*t2+0.1824681344E12*t1*phi+
1753186565.0);
    return(0.1E-4*t9);
  }
}

#include <math.h>
double evalRho(double x)
{
  double t1;
  double t3;
  double t4;
  double t5;
  double t8;
  {
    t1 = log(sol1);
    t3 = x*x;
    t4 = t3*t3;
    t5 = t4*x;
    t8 = t3*x;
    return(exp((a*t1+0.62121392-0.372728352E1*t5+0.93182088E1*t4-0.62121392E1*
t8)/((6.0*t5-15.0*t4+10.0*t8)*(a-0.1753186565)+0.1753186565)));
  }
}

#include <math.h>
double evalRho(double x)
{
  double t1;
  double t3;
  double t4;
  double t5;
  double t8;
  {
    t1 = log(sol1);
    t3 = x*x;
    t4 = t3*t3;
    t5 = t4*x;
    t8 = t3*x;
    return(exp((2.0*t1+0.62121392-0.372728352E1*t5+0.93182088E1*t4-0.62121392E1
*t8)/(0.1094808806E2*t5-0.2737022016E2*t4+0.1824681344E2*t8+0.1753186565)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t16;
  double t17;
  double t18;
  double t20;
  double t21;
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
    t20 = log(rho);
    t21 = rho*t20;
    return(2.0/delta*A*(t4+(2.0*alpha-2.0)*t7+(-3.0*alpha+1.0)*t3+alpha)+(t16-
t17+t18)*(2.0*t21-0.986266755*rho+0.1204760218E1)+(1.0-t16+t17-t18)*(
0.1753186565*t21+0.2172006685*rho+0.186848759E-1));
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t17;
  double t21;
  double t22;
  double t23;
  double t3;
  double t4;
  {
    t3 = phi*phi;
    t4 = t3*phi;
    t17 = t3*t3;
    t21 = 30.0*t17-60.0*t4+30.0*t3;
    t22 = log(rho);
    t23 = rho*t22;
    return(2.0/delta*A*(4.0*t4+3.0*(2.0*alpha-2.0)*t3+2.0*(-3.0*alpha+1.0)*phi)
+t21*(2.0*t23-0.986266755*rho+0.1204760218E1)-t21*(0.1753186565*t23+
0.2172006685*rho+0.186848759E-1));
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  double t17;
  double t18;
  double t19;
  double t3;
  {
    t3 = phi*phi;
    t17 = 120.0*t3*phi-180.0*t3+60.0*phi;
    t18 = log(rho);
    t19 = rho*t18;
    return(2.0/delta*A*(12.0*t3+6.0*(2.0*alpha-2.0)*phi-6.0*alpha+2.0)+t17*(2.0
*t19-0.986266755*rho+0.1204760218E1)-t17*(0.1753186565*t19+0.2172006685*rho+
0.186848759E-1));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t11;
  double t2;
  double t3;
  double t4;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = t2*phi;
    t4 = log(rho);
    t11 = t1*phi;
    return(0.1094808806E2*t3*t4+0.372728352E1*t3-0.2737022015E2*t2*t4
-0.93182088E1*t2+0.1824681344E2*t11*t4+0.62121392E1*t11+0.1753186565*t4+
0.392519325);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  double t1;
  double t2;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    return(0.5E-9*(0.2189617612E11*t2*phi-0.547404403E11*t2+0.3649362687E11*t1*
phi+350637313.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t2;
  double t3;
  double t7;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = log(rho);
    t7 = t1*phi;
    return(0.547404403E2*t2*t3+0.186364176E2*t2-0.1094808806E3*t7*t3
-0.372728352E2*t7+0.547404403E2*t1*t3+0.186364176E2*t1);
  }
}

#include <math.h>
double pressure(double rho,double phi)
{
  double t1;
  double t10;
  double t2;
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
    t9 = log(rho);
    t10 = rho*t9;
    return((t4-t5+t7)*(-2.0*t10+0.986266755*rho-0.1204760218E1+rho*(2.0*t9+
0.1013733245E1))+(1.0-t4+t5-t7)*(-0.1753186565*t10-0.2172006685*rho
-0.186848759E-1+rho*(0.1753186565*t9+0.392519325))-2.0/delta*A*(t2+(2.0*alpha
-2.0)*t6+(-3.0*alpha+1.0)*t1+alpha));
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t1;
  double t2;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t9 = sqrt(0.1094808806E12*t2*phi-0.2737022015E12*t2+0.1824681344E12*t1*phi+
1753186565.0);
    return(0.1E-4*t9);
  }
}

#include <math.h>
double evalRho(double x)
{
  double t1;
  double t3;
  double t4;
  double t5;
  double t8;
  {
    t1 = log(sol1);
    t3 = x*x;
    t4 = t3*t3;
    t5 = t4*x;
    t8 = t3*x;
    return(exp((0.62121392+a*t1-0.372728352E1*t5+0.93182088E1*t4-0.62121392E1*
t8)/(0.1753186565+(6.0*t5-15.0*t4+10.0*t8)*(a-0.1753186565))));
  }
}

#include <math.h>
double evalRho(double x)
{
  double t1;
  double t3;
  double t4;
  double t5;
  double t8;
  {
    t1 = log(sol1);
    t3 = x*x;
    t4 = t3*t3;
    t5 = t4*x;
    t8 = t3*x;
    return(exp((0.62121392+a*t1-0.372728352E1*t5+0.93182088E1*t4-0.62121392E1*
t8)/(0.1753186565+(6.0*t5-15.0*t4+10.0*t8)*(a-0.1753186565))));
  }
}

#include <math.h>
double evalRho(double x)
{
  double t1;
  double t3;
  double t4;
  double t5;
  double t8;
  {
    t1 = log(sol1);
    t3 = x*x;
    t4 = t3*t3;
    t5 = t4*x;
    t8 = t3*x;
    return(exp((2.0*t1+0.62121392-0.372728352E1*t5+0.93182088E1*t4-0.62121392E1
*t8)/(0.1094808806E2*t5-0.2737022016E2*t4+0.1824681344E2*t8+0.1753186565)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t16;
  double t17;
  double t18;
  double t20;
  double t21;
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
    t20 = log(rho);
    t21 = rho*t20;
    return(2.0/delta*A*(t4+(2.0*alpha-2.0)*t7+(-3.0*alpha+1.0)*t3+alpha)+(t16-
t17+t18)*(2.0*t21-0.986266755*rho+0.1204760218E1)+(1.0-t16+t17-t18)*(
0.1753186565*t21+0.2172006685*rho+0.186848759E-1));
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t17;
  double t21;
  double t22;
  double t23;
  double t3;
  double t4;
  {
    t3 = phi*phi;
    t4 = t3*phi;
    t17 = t3*t3;
    t21 = 30.0*t17-60.0*t4+30.0*t3;
    t22 = log(rho);
    t23 = rho*t22;
    return(2.0/delta*A*(4.0*t4+3.0*(2.0*alpha-2.0)*t3+2.0*(-3.0*alpha+1.0)*phi)
+t21*(2.0*t23-0.986266755*rho+0.1204760218E1)-t21*(0.1753186565*t23+
0.2172006685*rho+0.186848759E-1));
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  double t17;
  double t18;
  double t19;
  double t3;
  {
    t3 = phi*phi;
    t17 = 120.0*t3*phi-180.0*t3+60.0*phi;
    t18 = log(rho);
    t19 = rho*t18;
    return(2.0/delta*A*(12.0*t3+6.0*(2.0*alpha-2.0)*phi-6.0*alpha+2.0)+t17*(2.0
*t19-0.986266755*rho+0.1204760218E1)-t17*(0.1753186565*t19+0.2172006685*rho+
0.186848759E-1));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t11;
  double t2;
  double t3;
  double t4;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = t2*phi;
    t4 = log(rho);
    t11 = t1*phi;
    return(0.1094808806E2*t3*t4+0.372728352E1*t3-0.2737022015E2*t2*t4
-0.93182088E1*t2+0.1824681344E2*t11*t4+0.62121392E1*t11+0.1753186565*t4+
0.392519325);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  double t1;
  double t2;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    return(0.5E-9*(0.2189617612E11*t2*phi-0.547404403E11*t2+0.3649362687E11*t1*
phi+350637313.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t2;
  double t3;
  double t7;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = log(rho);
    t7 = t1*phi;
    return(0.547404403E2*t2*t3+0.186364176E2*t2-0.1094808806E3*t7*t3
-0.372728352E2*t7+0.547404403E2*t1*t3+0.186364176E2*t1);
  }
}

#include <math.h>
double pressure(double rho,double phi)
{
  double t1;
  double t10;
  double t2;
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
    t9 = log(rho);
    t10 = rho*t9;
    return((t4-t5+t7)*(-2.0*t10+0.986266755*rho-0.1204760218E1+rho*(2.0*t9+
0.1013733245E1))+(1.0-t4+t5-t7)*(-0.1753186565*t10-0.2172006685*rho
-0.186848759E-1+rho*(0.1753186565*t9+0.392519325))-2.0/delta*A*(t2+(2.0*alpha
-2.0)*t6+(-3.0*alpha+1.0)*t1+alpha));
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t1;
  double t2;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t9 = sqrt(0.1094808806E12*t2*phi-0.2737022015E12*t2+0.1824681344E12*t1*phi+
1753186565.0);
    return(0.1E-4*t9);
  }
}

#include <math.h>
double evalRho(double x)
{
  double t1;
  double t3;
  double t4;
  double t5;
  double t8;
  {
    t1 = log(sol1);
    t3 = x*x;
    t4 = t3*t3;
    t5 = t4*x;
    t8 = t3*x;
    return(exp((2.0*t1+0.62121392-0.372728352E1*t5+0.93182088E1*t4-0.62121392E1
*t8)/(0.1094808806E2*t5-0.2737022016E2*t4+0.1824681344E2*t8+0.1753186565)));
  }
}

#include <math.h>
double evalRho(double x)
{
  double t1;
  double t3;
  double t4;
  double t5;
  double t8;
  {
    t1 = log(sol1);
    t3 = x*x;
    t4 = t3*t3;
    t5 = t4*x;
    t8 = t3*x;
    return(exp((2.0*t1+0.62121392-0.372728352E1*t5+0.93182088E1*t4-0.62121392E1
*t8)/(0.1094808806E2*t5-0.2737022016E2*t4+0.1824681344E2*t8+0.1753186565)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t16;
  double t17;
  double t18;
  double t20;
  double t21;
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
    t20 = log(rho);
    t21 = rho*t20;
    return(2.0/delta*A*(t4+(2.0*alpha-2.0)*t7+(-3.0*alpha+1.0)*t3+alpha)+(t16-
t17+t18)*(2.0*t21-0.986266755*rho+0.1204760218E1)+(1.0-t16+t17-t18)*(
0.1753186565*t21+0.2172006685*rho+0.186848759E-1));
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t17;
  double t21;
  double t22;
  double t23;
  double t3;
  double t4;
  {
    t3 = phi*phi;
    t4 = t3*phi;
    t17 = t3*t3;
    t21 = 30.0*t17-60.0*t4+30.0*t3;
    t22 = log(rho);
    t23 = rho*t22;
    return(2.0/delta*A*(4.0*t4+3.0*(2.0*alpha-2.0)*t3+2.0*(-3.0*alpha+1.0)*phi)
+t21*(2.0*t23-0.986266755*rho+0.1204760218E1)-t21*(0.1753186565*t23+
0.2172006685*rho+0.186848759E-1));
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  double t17;
  double t18;
  double t19;
  double t3;
  {
    t3 = phi*phi;
    t17 = 120.0*t3*phi-180.0*t3+60.0*phi;
    t18 = log(rho);
    t19 = rho*t18;
    return(2.0/delta*A*(12.0*t3+6.0*(2.0*alpha-2.0)*phi-6.0*alpha+2.0)+t17*(2.0
*t19-0.986266755*rho+0.1204760218E1)-t17*(0.1753186565*t19+0.2172006685*rho+
0.186848759E-1));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t11;
  double t2;
  double t3;
  double t4;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = t2*phi;
    t4 = log(rho);
    t11 = t1*phi;
    return(0.1094808806E2*t3*t4+0.372728352E1*t3-0.2737022015E2*t2*t4
-0.93182088E1*t2+0.1824681344E2*t11*t4+0.62121392E1*t11+0.1753186565*t4+
0.392519325);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  double t1;
  double t2;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    return(0.5E-9*(0.2189617612E11*t2*phi-0.547404403E11*t2+0.3649362687E11*t1*
phi+350637313.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t2;
  double t3;
  double t7;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = log(rho);
    t7 = t1*phi;
    return(0.547404403E2*t2*t3+0.186364176E2*t2-0.1094808806E3*t7*t3
-0.372728352E2*t7+0.547404403E2*t1*t3+0.186364176E2*t1);
  }
}

#include <math.h>
double pressure(double rho,double phi)
{
  double t1;
  double t10;
  double t2;
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
    t9 = log(rho);
    t10 = rho*t9;
    return((t4-t5+t7)*(-2.0*t10+0.986266755*rho-0.1204760218E1+rho*(2.0*t9+
0.1013733245E1))+(1.0-t4+t5-t7)*(-0.1753186565*t10-0.2172006685*rho
-0.186848759E-1+rho*(0.1753186565*t9+0.392519325))-2.0/delta*A*(t2+(2.0*alpha
-2.0)*t6+(-3.0*alpha+1.0)*t1+alpha));
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t1;
  double t2;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t9 = sqrt(0.1094808806E12*t2*phi-0.2737022015E12*t2+0.1824681344E12*t1*phi+
1753186565.0);
    return(0.1E-4*t9);
  }
}

#include <math.h>
double evalRho(double x)
{
  double t1;
  double t3;
  double t4;
  double t5;
  double t8;
  {
    t1 = log(sol1);
    t3 = x*x;
    t4 = t3*t3;
    t5 = t4*x;
    t8 = t3*x;
    return(exp((a*t1+0.62121392-0.372728352E1*t5+0.93182088E1*t4-0.62121392E1*
t8)/((6.0*t5-15.0*t4+10.0*t8)*(a-0.1753186565)+0.1753186565)));
  }
}

#include <math.h>
double evalRho(double x)
{
  double t2;
  double t3;
  double t4;
  double t7;
  {
    t2 = x*x;
    t3 = t2*t2;
    t4 = t3*x;
    t7 = t2*x;
    return(exp((-0.2238890788E1*a+0.62121392-0.372728352E1*t4+0.93182088E1*t3
-0.62121392E1*t7)/((6.0*t4-15.0*t3+10.0*t7)*(a-0.1753186565)+0.1753186565)));
  }
}

#include <math.h>
double evalRho(double x)
{
  double t1;
  double t2;
  double t3;
  double t6;
  {
    t1 = x*x;
    t2 = t1*t1;
    t3 = t2*x;
    t6 = t1*x;
    return(exp((-0.3856567656E1-0.372728352E1*t3+0.93182088E1*t2-0.62121392E1*
t6)/(0.1094808806E2*t3-0.2737022016E2*t2+0.1824681344E2*t6+0.1753186565)));
  }
}

#include <math.h>
double evalRho(double x)
{
  double t1;
  double t2;
  double t3;
  double t6;
  {
    t1 = x*x;
    t2 = t1*t1;
    t3 = t2*x;
    t6 = t1*x;
    return(exp((-0.3856567656E1-0.372728352E1*t3+0.93182088E1*t2-0.62121392E1*
t6)/(0.1094808806E2*t3-0.2737022016E2*t2+0.1824681344E2*t6+0.1753186565)));
  }
}

#include <math.h>
double evalRho(double x)
{
  double t1;
  double t2;
  double t3;
  double t6;
  {
    t1 = x*x;
    t2 = t1*t1;
    t3 = t2*x;
    t6 = t1*x;
    return(exp((-0.392519325-0.372728352E1*t3+0.93182088E1*t2-0.62121392E1*t6)/
(0.1094808806E2*t3-0.2737022016E2*t2+0.1824681344E2*t6+0.1753186565)));
  }
}

#include <math.h>
double evalRho(double x)
{
  double t1;
  double t2;
  double t3;
  double t6;
  {
    t1 = x*x;
    t2 = t1*t1;
    t3 = t2*x;
    t6 = t1*x;
    return(exp((-0.3856567656E1-0.372728352E1*t3+0.93182088E1*t2-0.62121392E1*
t6)/(0.1094808806E2*t3-0.2737022016E2*t2+0.1824681344E2*t6+0.1753186565)));
  }
}

#include <math.h>
double evalRho(double x)
{
  double t1;
  double t2;
  double t3;
  double t6;
  {
    t1 = x*x;
    t2 = t1*t1;
    t3 = t2*x;
    t6 = t1*x;
    return(exp((-0.392519325-0.372728352E1*t3+0.93182088E1*t2-0.62121392E1*t6)/
(0.1094808806E2*t3-0.2737022016E2*t2+0.1824681344E2*t6+0.1753186565)));
  }
}


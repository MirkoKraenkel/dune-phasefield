#include <math.h>
double evalRho(double x)
{
  double t1;
  {
    t1 = nu(x);
    return(exp((0.2079441542E1-t1*(b-f)-f)/(-0.15E1*t1+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  {
    return(F(rho,phi));
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  {
    return(diff(F(rho,phi),phi));
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(diff(diff(F(rho,phi),phi),phi));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  {
    return(Potential(rho,phi));
  }
}

#include <math.h>
double drhochemicalPotential(double rho,double phi)
{
  {
    return(diff(Potential(rho,phi),rho));
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  {
    return(diff(Potential(rho,phi),phi));
  }
}

#include <math.h>
double pressure(double rho,double phi)
{
  {
    return(Pressure(rho,phi));
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t1;
  {
    t1 = diff(Pressure(rho,phi),rho);
    return(sqrt(t1));
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441542E1-x*(b-f)-f)/(-0.15E1*x+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t1;
  double t11;
  double t3;
  double t4;
  double t6;
  {
    t1 = log(phi);
    t3 = 1.0-phi;
    t4 = log(t3);
    t6 = phi*phi;
    t11 = log(rho);
    return(rho*(phi*t1+t3*t4-t6/2.0)+0.15E1*phi*rho*t11+3.0*t3*rho*t11);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t1;
  double t3;
  double t6;
  {
    t1 = log(phi);
    t3 = log(1.0-phi);
    t6 = log(rho);
    return(rho*(t1-t3-phi)-0.15E1*rho*t6);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho*(1/phi+1/(1.0-phi)-1.0));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = phi*phi;
    t10 = log(rho);
    return(phi*t1+t5-1.0*t5*phi-0.5*t8-0.15E1*phi*t10-0.15E1*phi+3.0*t10+3.0);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t2;
  double t4;
  double t6;
  {
    t1 = log(phi);
    t2 = 1.0*phi;
    t4 = log(1.0-t2);
    t6 = log(rho);
    return(t1-0.15E1-1.0*t4-t2-0.15E1*t6);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441542E1-x*(b-f)-f)/(-0.15E1*x+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t1;
  double t11;
  double t3;
  double t4;
  double t6;
  {
    t1 = log(phi);
    t3 = 1.0-phi;
    t4 = log(t3);
    t6 = phi*phi;
    t11 = log(rho);
    return(rho*(phi*t1+t3*t4-t6/2.0)+0.15E1*phi*rho*t11+3.0*t3*rho*t11);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t1;
  double t3;
  double t6;
  {
    t1 = log(phi);
    t3 = log(1.0-phi);
    t6 = log(rho);
    return(rho*(t1-t3-phi)-0.15E1*rho*t6);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho*(1/phi+1/(1.0-phi)-1.0));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = phi*phi;
    t10 = log(rho);
    return(phi*t1+t5-1.0*t5*phi-0.5*t8-0.15E1*phi*t10-0.15E1*phi+3.0*t10+3.0);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t2;
  double t4;
  double t6;
  {
    t1 = log(phi);
    t2 = 1.0*phi;
    t4 = log(1.0-t2);
    t6 = log(rho);
    return(t1-0.15E1-1.0*t4-t2-0.15E1*t6);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441542E1-x*(b-f)-f)/(-0.15E1*x+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t1;
  double t12;
  double t3;
  double t4;
  double t6;
  {
    t1 = log(phi);
    t3 = 1.0-phi;
    t4 = log(t3);
    t6 = phi*phi;
    t12 = log(rho);
    return(rho*(phi*t1+t3*t4-t6/2.0)/2.0+0.15E1*phi*rho*t12+3.0*t3*rho*t12);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t1;
  double t3;
  double t7;
  {
    t1 = log(phi);
    t3 = log(1.0-phi);
    t7 = log(rho);
    return(rho*(t1-t3-phi)/2.0-0.15E1*rho*t7);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho*(1/phi+1/(1.0-phi)-1.0)/2.0);
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t12;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t12 = log(rho);
    return(0.5*phi*t1+0.5*t6-0.5*t6*phi-0.25*t10-0.15E1*phi*t12-0.15E1*phi+3.0*
t12+3.0);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(0.5*t1-0.15E1-0.5*t5-0.5*phi-0.15E1*t8);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441542E1-x*(b-f)-f)/(-0.15E1*x+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t1;
  double t12;
  double t3;
  double t4;
  double t6;
  {
    t1 = log(phi);
    t3 = 1.0-phi;
    t4 = log(t3);
    t6 = phi*phi;
    t12 = log(rho);
    return(rho*(phi*t1+t3*t4-t6/2.0)/2.0+0.15E1*phi*rho*t12+3.0*t3*rho*t12);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t1;
  double t3;
  double t7;
  {
    t1 = log(phi);
    t3 = log(1.0-phi);
    t7 = log(rho);
    return(rho*(t1-t3-phi)/2.0-0.15E1*rho*t7);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho*(1/phi+1/(1.0-phi)-1.0)/2.0);
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t12;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t12 = log(rho);
    return(0.5*phi*t1+0.5*t6-0.5*t6*phi-0.25*t10-0.15E1*phi*t12-0.15E1*phi+3.0*
t12+3.0);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(0.5*t1-0.15E1-0.5*t5-0.5*phi-0.15E1*t8);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441542E1-x*(b-f)-f)/(-0.15E1*x+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t1;
  double t12;
  double t3;
  double t4;
  double t6;
  {
    t1 = log(phi);
    t3 = 1.0-phi;
    t4 = log(t3);
    t6 = phi*phi;
    t12 = log(rho);
    return(rho*(phi*t1+t3*t4-t6/2.0)/2.0+0.15E1*phi*rho*t12+3.0*t3*rho*t12);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t1;
  double t3;
  double t7;
  {
    t1 = log(phi);
    t3 = log(1.0-phi);
    t7 = log(rho);
    return(rho*(t1-t3-phi)/2.0-0.15E1*rho*t7);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho*(1/phi+1/(1.0-phi)-1.0)/2.0);
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t12;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t12 = log(rho);
    return(0.5*phi*t1+0.5*t6-0.5*t6*phi-0.25*t10-0.15E1*phi*t12-0.15E1*phi+3.0*
t12+3.0);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(0.5*t1-0.15E1-0.5*t5-0.5*phi-0.15E1*t8);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441542E1-x*(b-f)-f)/(-0.15E1*x+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t1;
  double t12;
  double t3;
  double t4;
  double t6;
  {
    t1 = log(phi);
    t3 = 1.0-phi;
    t4 = log(t3);
    t6 = phi*phi;
    t12 = log(rho);
    return(rho*(phi*t1+t3*t4-t6/2.0)/2.0+0.15E1*phi*rho*t12+3.0*t3*rho*t12);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t1;
  double t3;
  double t7;
  {
    t1 = log(phi);
    t3 = log(1.0-phi);
    t7 = log(rho);
    return(rho*(t1-t3-phi)/2.0-0.15E1*rho*t7);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho*(1/phi+1/(1.0-phi)-1.0)/2.0);
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t12;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t12 = log(rho);
    return(0.5*phi*t1+0.5*t6-0.5*t6*phi-0.25*t10-0.15E1*phi*t12-0.15E1*phi+3.0*
t12+3.0);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(0.5*t1-0.15E1-0.5*t5-0.5*phi-0.15E1*t8);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441542E1-x*(b-f)-f)/(-0.15E1*x+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t1;
  double t12;
  double t3;
  double t4;
  double t6;
  {
    t1 = log(phi);
    t3 = 1.0-phi;
    t4 = log(t3);
    t6 = phi*phi;
    t12 = log(rho);
    return(rho*(phi*t1+t3*t4-t6/2.0)/2.0+0.15E1*phi*rho*t12+3.0*t3*rho*t12);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t1;
  double t3;
  double t7;
  {
    t1 = log(phi);
    t3 = log(1.0-phi);
    t7 = log(rho);
    return(rho*(t1-t3-phi)/2.0-0.15E1*rho*t7);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho*(1/phi+1/(1.0-phi)-1.0)/2.0);
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t12;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t12 = log(rho);
    return(0.5*phi*t1+0.5*t6-0.5*t6*phi-0.25*t10-0.15E1*phi*t12-0.15E1*phi+3.0*
t12+3.0);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(0.5*t1-0.15E1-0.5*t5-0.5*phi-0.15E1*t8);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441542E1-x*(b-f)-f)/(-0.15E1*x+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t1;
  double t12;
  double t3;
  double t4;
  double t6;
  {
    t1 = log(phi);
    t3 = 1.0-phi;
    t4 = log(t3);
    t6 = phi*phi;
    t12 = log(rho);
    return(rho*(phi*t1+t3*t4-t6/2.0)/10.0+0.15E1*phi*rho*t12+3.0*t3*rho*t12);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t1;
  double t3;
  double t7;
  {
    t1 = log(phi);
    t3 = log(1.0-phi);
    t7 = log(rho);
    return(rho*(t1-t3-phi)/10.0-0.15E1*rho*t7);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho*(1/phi+1/(1.0-phi)-1.0)/10.0);
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t12;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t12 = log(rho);
    return(0.1*phi*t1+0.1*t6-0.1*t6*phi-0.5E-1*t10-0.15E1*phi*t12-0.15E1*phi+
3.0*t12+3.0);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(0.1*t1-0.15E1-0.1*t5-0.1*phi-0.15E1*t8);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441542E1-x*(b-f)-f)/(-0.15E1*x+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t1;
  double t12;
  double t3;
  double t4;
  double t6;
  {
    t1 = log(phi);
    t3 = 1.0-phi;
    t4 = log(t3);
    t6 = phi*phi;
    t12 = log(rho);
    return(0.1E2*rho*(phi*t1+t3*t4-t6/2.0)+0.15E1*phi*rho*t12+3.0*t3*rho*t12);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t1;
  double t3;
  double t7;
  {
    t1 = log(phi);
    t3 = log(1.0-phi);
    t7 = log(rho);
    return(0.1E2*rho*(t1-t3-phi)-0.15E1*rho*t7);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(0.1E2*rho*(1/phi+1/(1.0-phi)-1.0));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t12;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t12 = log(rho);
    return(10.0*phi*t1+10.0*t6-10.0*t6*phi-5.0*t10-0.15E1*phi*t12-0.15E1*phi+
3.0*t12+3.0);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(10.0*t1-0.15E1-10.0*t5-10.0*phi-0.15E1*t8);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441542E1-x*(b-f)-f)/(-0.15E1*x+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t2;
  {
    t2 = log(rho);
    return(0.15E1*phi*rho*t2+3.0*(1.0-phi)*rho*t2);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t1;
  {
    t1 = log(rho);
    return(-0.15E1*rho*t1);
  }
}

int dphireactionSource(double rho,double phi)
{
  {
    return(0);
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t12;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t12 = log(rho);
    return(10.0*phi*t1+10.0*t6-10.0*t6*phi-5.0*t10-0.15E1*phi*t12-0.15E1*phi+
3.0*t12+3.0);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(10.0*t1-0.15E1-10.0*t5-10.0*phi-0.15E1*t8);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.1386294361E1*a-x*(b-f)-f)/(x*(a-3.0)+3.0)));
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441542E1-x*(b-f)-f)/(-0.15E1*x+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t1;
  double t12;
  double t3;
  double t4;
  double t6;
  {
    t1 = log(phi);
    t3 = 1.0-phi;
    t4 = log(t3);
    t6 = phi*phi;
    t12 = log(rho);
    return(0.1E2*rho*(phi*t1+t3*t4-t6/2.0)+0.15E1*phi*rho*t12+3.0*t3*rho*t12);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t1;
  double t3;
  double t7;
  {
    t1 = log(phi);
    t3 = log(1.0-phi);
    t7 = log(rho);
    return(0.1E2*rho*(t1-t3-phi)-0.15E1*rho*t7);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(0.1E2*rho*(1/phi+1/(1.0-phi)-1.0));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t12;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t12 = log(rho);
    return(10.0*phi*t1+10.0*t6-10.0*t6*phi-5.0*t10-0.15E1*phi*t12-0.15E1*phi+
3.0*t12+3.0);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(10.0*t1-0.15E1-10.0*t5-10.0*phi-0.15E1*t8);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t12;
  double t13;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t12 = phi*a;
    t13 = log(rho);
    return(10.0*phi*t1+10.0*t6-10.0*t6*phi-5.0*t10+t12*t13+t12+3.0*t13+3.0-3.0*
phi*t13-3.0*phi);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441542E1-x*(b-f)-f)/(-0.15E1*x+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t1;
  double t12;
  double t3;
  double t4;
  double t6;
  {
    t1 = log(phi);
    t3 = 1.0-phi;
    t4 = log(t3);
    t6 = phi*phi;
    t12 = log(rho);
    return(rho*(phi*t1+t3*t4-t6/2.0)/10.0+0.15E1*phi*rho*t12+3.0*t3*rho*t12);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t1;
  double t3;
  double t7;
  {
    t1 = log(phi);
    t3 = log(1.0-phi);
    t7 = log(rho);
    return(rho*(t1-t3-phi)/10.0-0.15E1*rho*t7);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho*(1/phi+1/(1.0-phi)-1.0)/10.0);
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t12;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t12 = log(rho);
    return(0.1*phi*t1+0.1*t6-0.1*t6*phi-0.5E-1*t10-0.15E1*phi*t12-0.15E1*phi+
3.0*t12+3.0);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(0.1*t1-0.15E1-0.1*t5-0.1*phi-0.15E1*t8);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t1;
  double t3;
  double t8;
  {
    t1 = log(phi);
    t3 = log(1.0-phi);
    t8 = log(rho);
    return(rho*(t1-t3-phi)/10.0+a*rho*t8-3.0*rho*t8);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t1;
  double t3;
  double t8;
  {
    t1 = log(phi);
    t3 = log(1.0-phi);
    t8 = log(rho);
    return(rho*(t1-t3-phi)/10.0+a*rho*t8-3.0*rho*t8);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t1;
  double t3;
  double t8;
  {
    t1 = log(phi);
    t3 = log(1.0-phi);
    t8 = log(rho);
    return(rho*(t1-t3-phi)/10.0+a*rho*t8-3.0*rho*t8);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t1;
  double t3;
  double t8;
  {
    t1 = log(phi);
    t3 = log(1.0-phi);
    t8 = log(rho);
    return(rho*(t1-t3-phi)/10.0+a*rho*t8-3.0*rho*t8);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t1;
  double t3;
  double t8;
  {
    t1 = log(phi);
    t3 = log(1.0-phi);
    t8 = log(rho);
    return(rho*(t1-t3-phi)/10.0+a*rho*t8-3.0*rho*t8);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441542E1-x*(b-f)-f)/(-0.15E1*x+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t1;
  double t12;
  double t3;
  double t4;
  double t6;
  {
    t1 = log(phi);
    t3 = 1.0-phi;
    t4 = log(t3);
    t6 = phi*phi;
    t12 = log(rho);
    return(rho*(phi*t1+t3*t4-t6/2.0)/10.0+0.15E1*phi*rho*t12+3.0*t3*rho*t12);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t1;
  double t3;
  double t7;
  {
    t1 = log(phi);
    t3 = log(1.0-phi);
    t7 = log(rho);
    return(rho*(t1-t3-phi)/10.0-0.15E1*rho*t7);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho*(1/phi+1/(1.0-phi)-1.0)/10.0);
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t12;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t12 = log(rho);
    return(0.1*phi*t1+0.1*t6-0.1*t6*phi-0.5E-1*t10-0.15E1*phi*t12-0.15E1*phi+
3.0*t12+3.0);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(0.1*t1-0.15E1-0.1*t5-0.1*phi-0.15E1*t8);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((Const-x*(b-f)-f)/(x*(a-e)+e)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t1;
  double t12;
  double t13;
  double t3;
  double t4;
  double t6;
  {
    t1 = log(phi);
    t3 = 1.0-phi;
    t4 = log(t3);
    t6 = phi*phi;
    t12 = log(rho);
    t13 = rho*t12;
    return(rho*(phi*t1+t3*t4-t6/2.0)/10.0+phi*a*t13+t3*e*t13);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t1;
  double t3;
  double t8;
  {
    t1 = log(phi);
    t3 = log(1.0-phi);
    t8 = log(rho);
    return(rho*(t1-t3-phi)/10.0+a*rho*t8-e*rho*t8);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho*(1/phi+1/(1.0-phi)-1.0)/10.0);
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t12;
  double t13;
  double t16;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t12 = phi*a;
    t13 = log(rho);
    t16 = phi*e;
    return(0.1*phi*t1+0.1*t6-0.1*t6*phi-0.5E-1*t10+t12*t13+t12+e*t13+e-1.0*t16*
t13-1.0*t16);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return((phi*a+e-phi*e)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(0.1*t1-0.1*t5-0.1*phi+a*t8+a-1.0*e*t8-1.0*e);
  }
}

double pressure(double rho,double phi)
{
  {
    return(rho*(phi*a+e-phi*e));
  }
}

#include <math.h>
double a(double rho,double phi)
{
  {
    return(sqrt(phi*a+e-phi*e));
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho*(1/phi+1/(1.0-phi)-1.0)/10.0);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho*(1/phi+1/(1.0-phi)-1.0)/10.0);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho*(1/phi+1/(1.0-phi)-1.0)/10.0);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho*(1/phi+1/(1.0-phi)-1.0)/10.0);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho*(1/phi+1/(1.0-phi)-1.0)/10.0);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho*(1/phi+1/(1.0-phi)-1.0)/10.0);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho*(1/phi+1/(1.0-phi)-1.0)/10.0);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho*(1/phi+1/(1.0-phi)-1.0)/10.0);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho*(1/phi+1/(1.0-phi)-1.0)/10.0);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho*(1/phi+1/(1.0-phi)-1.0)/10.0);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho*(1/phi+1/(1.0-phi)-1.0)/10.0);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho*(1/phi+1/(1.0-phi)-1.0)/10.0);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho*(1/phi+1/(1.0-phi)-1.0)/10.0);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho*(1/phi+1/(1.0-phi)-1.0)/10.0);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441542E1-x*(b-f)-f)/(-0.15E1*x+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t12;
  double t2;
  double t4;
  double t5;
  double t7;
  {
    t2 = log(phi);
    t4 = 1.0-phi;
    t5 = log(t4);
    t7 = phi*phi;
    t12 = log(rho);
    return(rho*theta*(phi*t2+t4*t5-t7/2.0)+0.15E1*phi*rho*t12+3.0*t4*rho*t12);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t2;
  double t4;
  double t7;
  {
    t2 = log(phi);
    t4 = log(1.0-phi);
    t7 = log(rho);
    return(rho*theta*(t2-t4-phi)-0.15E1*rho*t7);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho*theta*(1/phi+1/(1.0-phi)-1.0));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t10;
  double t13;
  double t2;
  double t6;
  double t7;
  {
    t2 = log(phi);
    t6 = log(1.0-1.0*phi);
    t7 = theta*t6;
    t10 = phi*phi;
    t13 = log(rho);
    return(theta*phi*t2+t7-1.0*t7*phi-0.5*theta*t10-0.15E1*phi*t13-0.15E1*phi+
3.0*t13+3.0);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t5;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t10 = log(rho);
    return(theta*t1-1.0*theta*t5-1.0*theta*phi-0.15E1*t10-0.15E1);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441542E1-x*(b-f)-f)/(-0.15E1*x+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t13;
  double t3;
  double t5;
  double t6;
  double t8;
  {
    t3 = log(phi);
    t5 = 1.0-phi;
    t6 = log(t5);
    t8 = phi*phi;
    t13 = log(rho);
    return(rho/delta*(phi*t3+t5*t6-t8/2.0)+0.15E1*phi*rho*t13+3.0*t5*rho*t13);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t3;
  double t5;
  double t8;
  {
    t3 = log(phi);
    t5 = log(1.0-phi);
    t8 = log(rho);
    return(rho/delta*(t3-t5-phi)-0.15E1*rho*t8);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho/delta*(1/phi+1/(1.0-phi)-1.0));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t12;
  double t13;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t12 = phi*delta;
    t13 = log(rho);
    return(0.5*(2.0*phi*t1+2.0*t6-2.0*t6*phi-1.0*t10-3.0*t12*t13-3.0*t12+6.0*
delta*t13+6.0*delta)/delta);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(0.5*(2.0*t1-2.0*t5-2.0*phi-3.0*delta*t8-3.0*delta)/delta);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441542E1-x*(b-f)-f)/(-0.15E1*x+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t13;
  double t3;
  double t5;
  double t6;
  double t8;
  {
    t3 = log(phi);
    t5 = 1.0-phi;
    t6 = log(t5);
    t8 = phi*phi;
    t13 = log(rho);
    return(rho/delta*(phi*t3+t5*t6-t8/2.0)+0.15E1*phi*rho*t13+3.0*t5*rho*t13);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t3;
  double t5;
  double t8;
  {
    t3 = log(phi);
    t5 = log(1.0-phi);
    t8 = log(rho);
    return(rho/delta*(t3-t5-phi)-0.15E1*rho*t8);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho/delta*(1/phi+1/(1.0-phi)-1.0));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t12;
  double t13;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t12 = phi*delta;
    t13 = log(rho);
    return(0.5*(2.0*phi*t1+2.0*t6-2.0*t6*phi-1.0*t10-3.0*t12*t13-3.0*t12+6.0*
delta*t13+6.0*delta)/delta);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(0.5*(2.0*t1-2.0*t5-2.0*phi-3.0*t8*delta-3.0*delta)/delta);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441542E1-x*(b-f)-f)/(-0.15E1*x+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t13;
  double t3;
  double t5;
  double t6;
  double t8;
  {
    t3 = log(phi);
    t5 = 1.0-phi;
    t6 = log(t5);
    t8 = phi*phi;
    t13 = log(rho);
    return(rho/delta*(phi*t3+t5*t6-t8/2.0)+0.15E1*phi*rho*t13+3.0*t5*rho*t13);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t3;
  double t5;
  double t8;
  {
    t3 = log(phi);
    t5 = log(1.0-phi);
    t8 = log(rho);
    return(rho/delta*(t3-t5-phi)-0.15E1*rho*t8);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho/delta*(1/phi+1/(1.0-phi)-1.0));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t12;
  double t13;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t12 = phi*delta;
    t13 = log(rho);
    return(0.5*(2.0*phi*t1+2.0*t6-2.0*t6*phi-1.0*t10-3.0*t12*t13-3.0*t12+6.0*
delta*t13+6.0*delta)/delta);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(0.5*(2.0*t1-2.0*t5-2.0*phi-3.0*delta*t8-3.0*delta)/delta);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441542E1-x*(b-f)-f)/(-0.15E1*x+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t13;
  double t3;
  double t5;
  double t6;
  double t8;
  {
    t3 = log(phi);
    t5 = 1.0-phi;
    t6 = log(t5);
    t8 = phi*phi;
    t13 = log(rho);
    return(rho/delta*(phi*t3+t5*t6-t8/2.0)+0.15E1*phi*rho*t13+3.0*t5*rho*t13);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t3;
  double t5;
  double t8;
  {
    t3 = log(phi);
    t5 = log(1.0-phi);
    t8 = log(rho);
    return(rho/delta*(t3-t5-phi)-0.15E1*rho*t8);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho/delta*(1/phi+1/(1.0-phi)-1.0));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t12;
  double t13;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t12 = phi*delta;
    t13 = log(rho);
    return(0.5*(2.0*phi*t1+2.0*t6-2.0*t6*phi-1.0*t10-3.0*t12*t13-3.0*t12+6.0*
delta*t13+6.0*delta)/delta);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(0.5*(2.0*t1-2.0*t5-2.0*phi-3.0*delta*t8-3.0*delta)/delta);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441542E1-x*(b-f)-f)/(-0.15E1*x+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t13;
  double t3;
  double t5;
  double t6;
  double t8;
  {
    t3 = log(phi);
    t5 = 1.0-phi;
    t6 = log(t5);
    t8 = phi*phi;
    t13 = log(rho);
    return(rho/delta*(phi*t3+t5*t6-t8/2.0)+0.15E1*phi*rho*t13+3.0*t5*rho*t13);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t3;
  double t5;
  double t8;
  {
    t3 = log(phi);
    t5 = log(1.0-phi);
    t8 = log(rho);
    return(rho/delta*(t3-t5-phi)-0.15E1*rho*t8);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho/delta*(1/phi+1/(1.0-phi)-1.0));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t12;
  double t13;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t12 = phi*delta;
    t13 = log(rho);
    return(0.5*(2.0*phi*t1+2.0*t6-2.0*t6*phi-1.0*t10-3.0*t12*t13-3.0*t12+6.0*
delta*t13+6.0*delta)/delta);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(0.5*(2.0*t1-2.0*t5-2.0*phi-3.0*delta*t8-3.0*delta)/delta);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441542E1-x*(b-f)-f)/(-0.15E1*x+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t13;
  double t3;
  double t5;
  double t6;
  double t8;
  {
    t3 = log(phi);
    t5 = 1.0-phi;
    t6 = log(t5);
    t8 = phi*phi;
    t13 = log(rho);
    return(rho/delta*(phi*t3+t5*t6-t8/2.0)+0.15E1*phi*rho*t13+3.0*t5*rho*t13);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t3;
  double t5;
  double t8;
  {
    t3 = log(phi);
    t5 = log(1.0-phi);
    t8 = log(rho);
    return(rho/delta*(t3-t5-phi)-0.15E1*rho*t8);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho/delta*(1/phi+1/(1.0-phi)-1.0));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t12;
  double t13;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t12 = phi*delta;
    t13 = log(rho);
    return(0.5*(2.0*phi*t1+2.0*t6-2.0*t6*phi-1.0*t10-3.0*t12*t13-3.0*t12+6.0*
delta*t13+6.0*delta)/delta);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(0.5*(2.0*t1-2.0*t5-2.0*phi-3.0*delta*t8-3.0*delta)/delta);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441542E1-x*(b-f)-f)/(-0.15E1*x+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t13;
  double t3;
  double t5;
  double t6;
  double t8;
  {
    t3 = log(phi);
    t5 = 1.0-phi;
    t6 = log(t5);
    t8 = phi*phi;
    t13 = log(rho);
    return(rho/delta*(phi*t3+t5*t6-t8/2.0)+0.15E1*phi*rho*t13+3.0*t5*rho*t13);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t3;
  double t5;
  double t8;
  {
    t3 = log(phi);
    t5 = log(1.0-phi);
    t8 = log(rho);
    return(rho/delta*(t3-t5-phi)-0.15E1*rho*t8);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho/delta*(1/phi+1/(1.0-phi)-1.0));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t11;
  double t12;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t11 = phi*delta;
    t12 = log(rho);
    return(-0.5*(-2.0*phi*t1-2.0*t6+2.0*t6*phi+t10+3.0*t11*t12+3.0*t11-6.0*
delta*t12-6.0*delta)/delta);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(-0.5*(-2.0*t1+2.0*t5+2.0*phi+3.0*delta*t8+3.0*delta)/delta);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441542E1-x*(0.15E1-f)-f)/(x*(a-3.0)+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t13;
  double t3;
  double t5;
  double t6;
  double t8;
  {
    t3 = log(phi);
    t5 = 1.0-phi;
    t6 = log(t5);
    t8 = phi*phi;
    t13 = log(rho);
    return(rho/delta*(phi*t3+t5*t6-t8/2.0)+0.15E1*phi*rho*t13+3.0*t5*rho*t13);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t3;
  double t5;
  double t8;
  {
    t3 = log(phi);
    t5 = log(1.0-phi);
    t8 = log(rho);
    return(rho/delta*(t3-t5-phi)-0.15E1*rho*t8);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho/delta*(1/phi+1/(1.0-phi)-1.0));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t11;
  double t12;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t11 = phi*delta;
    t12 = log(rho);
    return(-0.5*(-2.0*phi*t1-2.0*t6+2.0*t6*phi+t10+3.0*t11*t12+3.0*t11-6.0*
delta*t12-6.0*delta)/delta);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(-0.5*(-2.0*t1+2.0*t5+2.0*phi+3.0*delta*t8+3.0*delta)/delta);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441542E1-x*(0.15E1-f)-f)/(x*(a-3.0)+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t13;
  double t3;
  double t5;
  double t6;
  double t8;
  {
    t3 = log(phi);
    t5 = 1.0-phi;
    t6 = log(t5);
    t8 = phi*phi;
    t13 = log(rho);
    return(rho/delta*(phi*t3+t5*t6-t8/2.0)+0.15E1*phi*rho*t13+3.0*t5*rho*t13);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t3;
  double t5;
  double t8;
  {
    t3 = log(phi);
    t5 = log(1.0-phi);
    t8 = log(rho);
    return(rho/delta*(t3-t5-phi)-0.15E1*rho*t8);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho/delta*(1/phi+1/(1.0-phi)-1.0));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t12;
  double t13;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t12 = phi*delta;
    t13 = log(rho);
    return(0.5*(2.0*phi*t1+2.0*t6-2.0*t6*phi-1.0*t10-3.0*t12*t13-3.0*t12+6.0*
delta*t13+6.0*delta)/delta);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(0.5*(2.0*t1-2.0*t5-2.0*phi-3.0*t8*delta-3.0*delta)/delta);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441542E1-x*(0.15E1-f)-f)/(x*(a-3.0)+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t13;
  double t3;
  double t5;
  double t6;
  double t8;
  {
    t3 = log(phi);
    t5 = 1.0-phi;
    t6 = log(t5);
    t8 = phi*phi;
    t13 = log(rho);
    return(rho/delta*(phi*t3+t5*t6-t8/2.0)+0.15E1*phi*rho*t13+3.0*t5*rho*t13);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t3;
  double t5;
  double t8;
  {
    t3 = log(phi);
    t5 = log(1.0-phi);
    t8 = log(rho);
    return(rho/delta*(t3-t5-phi)-0.15E1*rho*t8);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho/delta*(1/phi+1/(1.0-phi)-1.0));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t12;
  double t13;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t12 = phi*delta;
    t13 = log(rho);
    return(0.5*(2.0*phi*t1+2.0*t6-2.0*t6*phi-1.0*t10-3.0*t12*t13-3.0*t12+6.0*
delta*t13+6.0*delta)/delta);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(0.5*(2.0*t1-2.0*t5-2.0*phi-3.0*t8*delta-3.0*delta)/delta);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441542E1-x*(0.15E1-f)-f)/(x*(a-3.0)+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t13;
  double t3;
  double t5;
  double t6;
  double t8;
  {
    t3 = log(phi);
    t5 = 1.0-phi;
    t6 = log(t5);
    t8 = phi*phi;
    t13 = log(rho);
    return(rho/delta*(phi*t3+t5*t6-t8/2.0)+0.15E1*phi*rho*t13+3.0*t5*rho*t13);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t3;
  double t5;
  double t8;
  {
    t3 = log(phi);
    t5 = log(1.0-phi);
    t8 = log(rho);
    return(rho/delta*(t3-t5-phi)-0.15E1*rho*t8);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho/delta*(1/phi+1/(1.0-phi)-1.0));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t12;
  double t13;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t12 = phi*delta;
    t13 = log(rho);
    return(0.5*(2.0*phi*t1+2.0*t6-2.0*t6*phi-1.0*t10-3.0*t12*t13-3.0*t12+6.0*
delta*t13+6.0*delta)/delta);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(0.5*(2.0*t1-2.0*t5-2.0*phi-3.0*t8*delta-3.0*delta)/delta);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441542E1-x*(0.15E1-f)-f)/(x*(a-3.0)+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t13;
  double t3;
  double t5;
  double t6;
  double t8;
  {
    t3 = log(phi);
    t5 = 1.0-phi;
    t6 = log(t5);
    t8 = phi*phi;
    t13 = log(rho);
    return(rho/delta*(phi*t3+t5*t6-t8/2.0)+0.15E1*phi*rho*t13+3.0*t5*rho*t13);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t3;
  double t5;
  double t8;
  {
    t3 = log(phi);
    t5 = log(1.0-phi);
    t8 = log(rho);
    return(rho/delta*(t3-t5-phi)-0.15E1*rho*t8);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho/delta*(1/phi+1/(1.0-phi)-1.0));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t12;
  double t13;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t12 = phi*delta;
    t13 = log(rho);
    return(0.5*(2.0*phi*t1+2.0*t6-2.0*t6*phi-1.0*t10-3.0*t12*t13-3.0*t12+6.0*
delta*t13+6.0*delta)/delta);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(0.5*(2.0*t1-2.0*t5-2.0*phi-3.0*delta*t8-3.0*delta)/delta);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441542E1-x*(0.15E1-f)-f)/(x*(a-3.0)+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t13;
  double t3;
  double t5;
  double t6;
  double t8;
  {
    t3 = log(phi);
    t5 = 1.0-phi;
    t6 = log(t5);
    t8 = phi*phi;
    t13 = log(rho);
    return(rho/delta*(phi*t3+t5*t6-t8/2.0)+0.15E1*phi*rho*t13+3.0*t5*rho*t13);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t3;
  double t5;
  double t8;
  {
    t3 = log(phi);
    t5 = log(1.0-phi);
    t8 = log(rho);
    return(rho/delta*(t3-t5-phi)-0.15E1*rho*t8);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho/delta*(1/phi+1/(1.0-phi)-1.0));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t12;
  double t13;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t12 = phi*delta;
    t13 = log(rho);
    return(0.5*(2.0*phi*t1+2.0*t6-2.0*t6*phi-1.0*t10-3.0*t12*t13-3.0*t12+6.0*
delta*t13+6.0*delta)/delta);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(0.5*(2.0*t1-2.0*t5-2.0*phi-3.0*delta*t8-3.0*delta)/delta);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441542E1-x*(0.15E1-f)-f)/(x*(a-3.0)+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t13;
  double t3;
  double t5;
  double t6;
  double t8;
  {
    t3 = log(phi);
    t5 = 1.0-phi;
    t6 = log(t5);
    t8 = phi*phi;
    t13 = log(rho);
    return(rho/delta*(phi*t3+t5*t6-t8/2.0)+0.15E1*phi*rho*t13+3.0*t5*rho*t13);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t3;
  double t5;
  double t8;
  {
    t3 = log(phi);
    t5 = log(1.0-phi);
    t8 = log(rho);
    return(rho/delta*(t3-t5-phi)-0.15E1*rho*t8);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho/delta*(1/phi+1/(1.0-phi)-1.0));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t12;
  double t13;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t12 = phi*delta;
    t13 = log(rho);
    return(0.5*(2.0*phi*t1+2.0*t6-2.0*t6*phi-1.0*t10-3.0*t12*t13-3.0*t12+6.0*
delta*t13+6.0*delta)/delta);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(0.5*(2.0*t1-2.0*t5-2.0*phi-3.0*t8*delta-3.0*delta)/delta);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441542E1-x*(0.15E1-f)-f)/(x*(a-3.0)+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t13;
  double t3;
  double t5;
  double t6;
  double t8;
  {
    t3 = log(phi);
    t5 = 1.0-phi;
    t6 = log(t5);
    t8 = phi*phi;
    t13 = log(rho);
    return(rho/delta*(phi*t3+t5*t6-t8/2.0)+0.15E1*phi*rho*t13+3.0*t5*rho*t13);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t3;
  double t5;
  double t8;
  {
    t3 = log(phi);
    t5 = log(1.0-phi);
    t8 = log(rho);
    return(rho/delta*(t3-t5-phi)-0.15E1*rho*t8);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho/delta*(1/phi+1/(1.0-phi)-1.0));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t12;
  double t13;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t12 = phi*delta;
    t13 = log(rho);
    return(0.5*(2.0*phi*t1+2.0*t6-2.0*t6*phi-1.0*t10-3.0*t12*t13-3.0*t12+6.0*
delta*t13+6.0*delta)/delta);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(0.5*(2.0*t1-2.0*t5-2.0*phi-3.0*delta*t8-3.0*delta)/delta);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441542E1-x*(0.15E1-f)-f)/(x*(a-3.0)+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t13;
  double t3;
  double t5;
  double t6;
  double t8;
  {
    t3 = log(phi);
    t5 = 1.0-phi;
    t6 = log(t5);
    t8 = phi*phi;
    t13 = log(rho);
    return(rho/delta*(phi*t3+t5*t6-t8/2.0)+0.15E1*phi*rho*t13+3.0*t5*rho*t13);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t3;
  double t5;
  double t8;
  {
    t3 = log(phi);
    t5 = log(1.0-phi);
    t8 = log(rho);
    return(rho/delta*(t3-t5-phi)-0.15E1*rho*t8);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho/delta*(1/phi+1/(1.0-phi)-1.0));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t12;
  double t13;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t12 = phi*delta;
    t13 = log(rho);
    return(0.5*(2.0*phi*t1+2.0*t6-2.0*t6*phi-1.0*t10-3.0*t12*t13-3.0*t12+6.0*
delta*t13+6.0*delta)/delta);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(0.5*(2.0*t1-2.0*t5-2.0*phi-3.0*delta*t8-3.0*delta)/delta);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2579941723E1-x*(0.15E1-f)-f)/(x*(a-5.0)+5.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t13;
  double t3;
  double t5;
  double t6;
  double t8;
  {
    t3 = log(phi);
    t5 = 1.0-phi;
    t6 = log(t5);
    t8 = phi*phi;
    t13 = log(rho);
    return(rho/delta*(phi*t3+t5*t6-t8/2.0)+0.15E1*phi*rho*t13+5.0*t5*rho*t13);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t3;
  double t5;
  double t8;
  {
    t3 = log(phi);
    t5 = log(1.0-phi);
    t8 = log(rho);
    return(rho/delta*(t3-t5-phi)-0.35E1*rho*t8);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho/delta*(1/phi+1/(1.0-phi)-1.0));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t12;
  double t13;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t12 = phi*delta;
    t13 = log(rho);
    return(0.5*(2.0*phi*t1+2.0*t6-2.0*t6*phi-1.0*t10-7.0*t12*t13-7.0*t12+10.0*
delta*t13+10.0*delta)/delta);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.5*(7.0*phi-10.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(0.5*(2.0*t1-2.0*t5-2.0*phi-7.0*t8*delta-7.0*delta)/delta);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.35E1*phi*rho+5.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-14.0*phi+20.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2579941723E1-x*(0.15E1-f)-f)/(x*(a-5.0)+5.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t13;
  double t3;
  double t5;
  double t6;
  double t8;
  {
    t3 = log(phi);
    t5 = 1.0-phi;
    t6 = log(t5);
    t8 = phi*phi;
    t13 = log(rho);
    return(rho/delta*(phi*t3+t5*t6-t8/2.0)+0.15E1*phi*rho*t13+5.0*t5*rho*t13);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t3;
  double t5;
  double t8;
  {
    t3 = log(phi);
    t5 = log(1.0-phi);
    t8 = log(rho);
    return(rho/delta*(t3-t5-phi)-0.35E1*rho*t8);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho/delta*(1/phi+1/(1.0-phi)-1.0));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t12;
  double t13;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t12 = phi*delta;
    t13 = log(rho);
    return(0.5*(2.0*phi*t1+2.0*t6-2.0*t6*phi-1.0*t10-7.0*t12*t13-7.0*t12+10.0*
delta*t13+10.0*delta)/delta);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.5*(7.0*phi-10.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(0.5*(2.0*t1-2.0*t5-2.0*phi-7.0*delta*t8-7.0*delta)/delta);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.35E1*phi*rho+5.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-14.0*phi+20.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2579941723E1-x*(0.15E1-f)-f)/(x*(a-5.0)+5.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t13;
  double t3;
  double t5;
  double t6;
  double t8;
  {
    t3 = log(phi);
    t5 = 1.0-phi;
    t6 = log(t5);
    t8 = phi*phi;
    t13 = log(rho);
    return(rho/delta*(phi*t3+t5*t6-t8/2.0)+0.15E1*phi*rho*t13+5.0*t5*rho*t13);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t3;
  double t5;
  double t8;
  {
    t3 = log(phi);
    t5 = log(1.0-phi);
    t8 = log(rho);
    return(rho/delta*(t3-t5-phi)-0.35E1*rho*t8);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho/delta*(1/phi+1/(1.0-phi)-1.0));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t11;
  double t12;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t11 = phi*delta;
    t12 = log(rho);
    return(-0.5*(-2.0*phi*t1-2.0*t6+2.0*t6*phi+t10+7.0*t11*t12+7.0*t11-10.0*
delta*t12-10.0*delta)/delta);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.5*(7.0*phi-10.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(-0.5*(-2.0*t1+2.0*t5+2.0*phi+7.0*t8*delta+7.0*delta)/delta);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.35E1*phi*rho+5.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-14.0*phi+20.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2579941723E1-x*(0.15E1-f)-f)/(x*(a-5.0)+5.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t13;
  double t3;
  double t5;
  double t6;
  double t8;
  {
    t3 = log(phi);
    t5 = 1.0-phi;
    t6 = log(t5);
    t8 = phi*phi;
    t13 = log(rho);
    return(rho/delta*(phi*t3+t5*t6-t8/2.0)+0.15E1*phi*rho*t13+5.0*t5*rho*t13);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t3;
  double t5;
  double t8;
  {
    t3 = log(phi);
    t5 = log(1.0-phi);
    t8 = log(rho);
    return(rho/delta*(t3-t5-phi)-0.35E1*rho*t8);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho/delta*(1/phi+1/(1.0-phi)-1.0));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t12;
  double t13;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t12 = phi*delta;
    t13 = log(rho);
    return(0.5*(2.0*phi*t1+2.0*t6-2.0*t6*phi-1.0*t10-7.0*t12*t13-7.0*t12+10.0*
delta*t13+10.0*delta)/delta);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.5*(7.0*phi-10.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(0.5*(2.0*t1-2.0*t5-2.0*phi-7.0*delta*t8-7.0*delta)/delta);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.35E1*phi*rho+5.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-14.0*phi+20.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((Const-x*(0.15E1-f)-f)/(x*(a-15.0)+15.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t13;
  double t3;
  double t5;
  double t6;
  double t8;
  {
    t3 = log(phi);
    t5 = 1.0-phi;
    t6 = log(t5);
    t8 = phi*phi;
    t13 = log(rho);
    return(rho/delta*(phi*t3+t5*t6-t8/2.0)+0.15E1*phi*rho*t13+15.0*t5*rho*t13);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t3;
  double t5;
  double t8;
  {
    t3 = log(phi);
    t5 = log(1.0-phi);
    t8 = log(rho);
    return(rho/delta*(t3-t5-phi)-0.135E2*rho*t8);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho/delta*(1/phi+1/(1.0-phi)-1.0));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t11;
  double t12;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t11 = phi*delta;
    t12 = log(rho);
    return(-0.5*(-2.0*phi*t1-2.0*t6+2.0*t6*phi+t10+27.0*t11*t12+27.0*t11-30.0*
delta*t12-30.0*delta)/delta);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(9.0*phi-10.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(-0.5*(-2.0*t1+2.0*t5+2.0*phi+27.0*delta*t8+27.0*delta)/delta);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.135E2*phi*rho+15.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-54.0*phi+60.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441542E1-x*(0.15E1-f)-f)/(x*(a-3.0)+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t13;
  double t3;
  double t5;
  double t6;
  double t8;
  {
    t3 = log(phi);
    t5 = 1.0-phi;
    t6 = log(t5);
    t8 = phi*phi;
    t13 = log(rho);
    return(rho/delta*(phi*t3+t5*t6-t8/2.0)+0.15E1*phi*rho*t13+3.0*t5*rho*t13);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t3;
  double t5;
  double t8;
  {
    t3 = log(phi);
    t5 = log(1.0-phi);
    t8 = log(rho);
    return(rho/delta*(t3-t5-phi)-0.15E1*rho*t8);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho/delta*(1/phi+1/(1.0-phi)-1.0));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t12;
  double t13;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t12 = phi*delta;
    t13 = log(rho);
    return(0.5*(2.0*phi*t1+2.0*t6-2.0*t6*phi-1.0*t10-3.0*t12*t13-3.0*t12+6.0*
delta*t13+6.0*delta)/delta);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(0.5*(2.0*t1-2.0*t5-2.0*phi-3.0*delta*t8-3.0*delta)/delta);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441542E1-x*(0.15E1-f)-f)/(x*(a-3.0)+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t1;
  double t13;
  double t3;
  double t4;
  double t6;
  {
    t1 = log(phi);
    t3 = 1.0-phi;
    t4 = log(t3);
    t6 = phi*phi;
    t13 = log(rho);
    return(rho*(phi*t1+t3*t4-t6/2.0)/delta+0.15E1*phi*rho*t13+3.0*t3*rho*t13);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t1;
  double t3;
  double t8;
  {
    t1 = log(phi);
    t3 = log(1.0-phi);
    t8 = log(rho);
    return(rho*(t1-t3-phi)/delta-0.15E1*rho*t8);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho*(1/phi+1/(1.0-phi)-1.0)/delta);
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t12;
  double t13;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t12 = phi*delta;
    t13 = log(rho);
    return(0.5*(2.0*phi*t1+2.0*t6-2.0*t6*phi-1.0*t10-3.0*t12*t13-3.0*t12+6.0*
delta*t13+6.0*delta)/delta);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(-0.5*(-2.0*t1+2.0*t5+2.0*phi+3.0*delta*t8+3.0*delta)/delta);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441542E1-x*(0.15E1-f)-f)/(x*(a-3.0)+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t13;
  double t3;
  double t5;
  double t6;
  double t8;
  {
    t3 = log(phi);
    t5 = 1.0-phi;
    t6 = log(t5);
    t8 = phi*phi;
    t13 = log(rho);
    return(rho/delta*(phi*t3+t5*t6-t8/2.0)+0.15E1*phi*rho*t13+3.0*t5*rho*t13);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t3;
  double t5;
  double t8;
  {
    t3 = log(phi);
    t5 = log(1.0-phi);
    t8 = log(rho);
    return(rho/delta*(t3-t5-phi)-0.15E1*rho*t8);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho/delta*(1/phi+1/(1.0-phi)-1.0));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t12;
  double t13;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t12 = phi*delta;
    t13 = log(rho);
    return(0.5*(2.0*phi*t1+2.0*t6-2.0*t6*phi-1.0*t10-3.0*t12*t13-3.0*t12+6.0*
delta*t13+6.0*delta)/delta);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(0.5*(2.0*t1-2.0*t5-2.0*phi-3.0*t8*delta-3.0*delta)/delta);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441542E1-x*(0.15E1-f)-f)/(x*(a-3.0)+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t13;
  double t3;
  double t5;
  double t6;
  double t8;
  {
    t3 = log(phi);
    t5 = 1.0-phi;
    t6 = log(t5);
    t8 = phi*phi;
    t13 = log(rho);
    return(rho/delta*(phi*t3+t5*t6-t8/2.0)+0.15E1*phi*rho*t13+3.0*t5*rho*t13);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t3;
  double t5;
  double t8;
  {
    t3 = log(phi);
    t5 = log(1.0-phi);
    t8 = log(rho);
    return(rho/delta*(t3-t5-phi)-0.15E1*rho*t8);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho/delta*(1/phi+1/(1.0-phi)-1.0));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t12;
  double t13;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t12 = phi*delta;
    t13 = log(rho);
    return(0.5*(2.0*phi*t1+2.0*t6-2.0*t6*phi-1.0*t10-3.0*t12*t13-3.0*t12+6.0*
delta*t13+6.0*delta)/delta);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(0.5*(2.0*t1-2.0*t5-2.0*phi-3.0*t8*delta-3.0*delta)/delta);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441542E1-x*(0.15E1-f)-f)/(x*(a-3.0)+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t13;
  double t3;
  double t5;
  double t6;
  double t8;
  {
    t3 = log(phi);
    t5 = 1.0-phi;
    t6 = log(t5);
    t8 = phi*phi;
    t13 = log(rho);
    return(rho/delta*(phi*t3+t5*t6-t8/2.0)+0.15E1*phi*rho*t13+3.0*t5*rho*t13);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t3;
  double t5;
  double t8;
  {
    t3 = log(phi);
    t5 = log(1.0-phi);
    t8 = log(rho);
    return(rho/delta*(t3-t5-phi)-0.15E1*rho*t8);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho/delta*(1/phi+1/(1.0-phi)-1.0));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t12;
  double t13;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t12 = phi*delta;
    t13 = log(rho);
    return(0.5*(2.0*phi*t1+2.0*t6-2.0*t6*phi-1.0*t10-3.0*t12*t13-3.0*t12+6.0*
delta*t13+6.0*delta)/delta);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(0.5*(2.0*t1-2.0*t5-2.0*phi-3.0*t8*delta-3.0*delta)/delta);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441542E1-x*(0.15E1-f)-f)/(x*(a-3.0)+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t13;
  double t3;
  double t5;
  double t6;
  double t8;
  {
    t3 = log(phi);
    t5 = 1.0-phi;
    t6 = log(t5);
    t8 = phi*phi;
    t13 = log(rho);
    return(rho/delta*(phi*t3+t5*t6-t8/2.0)+0.15E1*phi*rho*t13+3.0*t5*rho*t13);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t3;
  double t5;
  double t8;
  {
    t3 = log(phi);
    t5 = log(1.0-phi);
    t8 = log(rho);
    return(rho/delta*(t3-t5-phi)-0.15E1*rho*t8);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho/delta*(1/phi+1/(1.0-phi)-1.0));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t11;
  double t12;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t11 = phi*delta;
    t12 = log(rho);
    return(-0.5*(-2.0*phi*t1-2.0*t6+2.0*t6*phi+t10+3.0*t11*t12+3.0*t11-6.0*
delta*t12-6.0*delta)/delta);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(-0.5*(-2.0*t1+2.0*t5+2.0*phi+3.0*delta*t8+3.0*delta)/delta);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441541679836E1-x*(0.15E1-f)-f)/(x*(a-3.0)+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t13;
  double t3;
  double t5;
  double t6;
  double t8;
  {
    t3 = log(phi);
    t5 = 1.0-phi;
    t6 = log(t5);
    t8 = phi*phi;
    t13 = log(rho);
    return(rho/delta*(phi*t3+t5*t6-t8/2.0)+0.15E1*phi*rho*t13+3.0*t5*rho*t13);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t3;
  double t5;
  double t8;
  {
    t3 = log(phi);
    t5 = log(1.0-phi);
    t8 = log(rho);
    return(rho/delta*(t3-t5-phi)-0.15E1*rho*t8);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho/delta*(1/phi+1/(1.0-phi)-1.0));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t11;
  double t12;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t11 = phi*delta;
    t12 = log(rho);
    return(-0.5*(-2.0*phi*t1-2.0*t6+2.0*t6*phi+t10+3.0*t11*t12+3.0*t11-6.0*
delta*t12-6.0*delta)/delta);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(-0.5*(-2.0*t1+2.0*t5+2.0*phi+3.0*delta*t8+3.0*delta)/delta);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441541679836E1-x*(0.15E1-f)-f)/(x*(a-3.0)+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t13;
  double t3;
  double t5;
  double t6;
  double t8;
  {
    t3 = log(phi);
    t5 = 1.0-phi;
    t6 = log(t5);
    t8 = phi*phi;
    t13 = log(rho);
    return(rho/delta*(phi*t3+t5*t6-t8/2.0)+0.15E1*phi*rho*t13+3.0*t5*rho*t13);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t3;
  double t5;
  double t8;
  {
    t3 = log(phi);
    t5 = log(1.0-phi);
    t8 = log(rho);
    return(rho/delta*(t3-t5-phi)-0.15E1*rho*t8);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho/delta*(1/phi+1/(1.0-phi)-1.0));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t12;
  double t13;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t12 = phi*delta;
    t13 = log(rho);
    return(0.5*(2.0*phi*t1+2.0*t6-2.0*t6*phi-1.0*t10-3.0*t12*t13-3.0*t12+6.0*
delta*t13+6.0*delta)/delta);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(0.5*(2.0*t1-2.0*t5-2.0*phi-3.0*t8*delta-3.0*delta)/delta);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441541679836E1-x*(0.15E1-f)-f)/(x*(a-3.0)+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t13;
  double t3;
  double t5;
  double t6;
  double t8;
  {
    t3 = log(phi);
    t5 = 1.0-phi;
    t6 = log(t5);
    t8 = phi*phi;
    t13 = log(rho);
    return(rho/delta*(phi*t3+t5*t6-t8/2.0)+0.15E1*phi*rho*t13+3.0*t5*rho*t13);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t3;
  double t5;
  double t8;
  {
    t3 = log(phi);
    t5 = log(1.0-phi);
    t8 = log(rho);
    return(rho/delta*(t3-t5-phi)-0.15E1*rho*t8);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho/delta*(1/phi+1/(1.0-phi)-1.0));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t12;
  double t13;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t12 = phi*delta;
    t13 = log(rho);
    return(0.5*(2.0*phi*t1+2.0*t6-2.0*t6*phi-1.0*t10-3.0*t12*t13-3.0*t12+6.0*
delta*t13+6.0*delta)/delta);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(0.5*(2.0*t1-2.0*t5-2.0*phi-3.0*t8*delta-3.0*delta)/delta);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441541679836E1-x*(0.15E1-f)-f)/(x*(a-3.0)+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t13;
  double t3;
  double t5;
  double t6;
  double t8;
  {
    t3 = log(phi);
    t5 = 1.0-phi;
    t6 = log(t5);
    t8 = phi*phi;
    t13 = log(rho);
    return(rho/delta*(phi*t3+t5*t6-t8/2.0)+0.15E1*phi*rho*t13+3.0*t5*rho*t13);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t3;
  double t5;
  double t8;
  {
    t3 = log(phi);
    t5 = log(1.0-phi);
    t8 = log(rho);
    return(rho/delta*(t3-t5-phi)-0.15E1*rho*t8);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho/delta*(1/phi+1/(1.0-phi)-1.0));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t12;
  double t13;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t12 = phi*delta;
    t13 = log(rho);
    return(0.5*(2.0*phi*t1+2.0*t6-2.0*t6*phi-1.0*t10-3.0*t12*t13-3.0*t12+6.0*
delta*t13+6.0*delta)/delta);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(0.5*(2.0*t1-2.0*t5-2.0*phi-3.0*delta*t8-3.0*delta)/delta);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441541679836E1-x*(0.15E1-f)-f)/(x*(a-3.0)+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t13;
  double t3;
  double t5;
  double t6;
  double t8;
  {
    t3 = log(phi);
    t5 = 1.0-phi;
    t6 = log(t5);
    t8 = phi*phi;
    t13 = log(rho);
    return(rho/delta*(phi*t3+t5*t6-t8/2.0)+0.15E1*phi*rho*t13+3.0*t5*rho*t13);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t3;
  double t5;
  double t8;
  {
    t3 = log(phi);
    t5 = log(1.0-phi);
    t8 = log(rho);
    return(rho/delta*(t3-t5-phi)-0.15E1*rho*t8);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho/delta*(1/phi+1/(1.0-phi)-1.0));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t12;
  double t13;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t12 = phi*delta;
    t13 = log(rho);
    return(0.5*(2.0*phi*t1+2.0*t6-2.0*t6*phi-1.0*t10-3.0*t12*t13-3.0*t12+6.0*
delta*t13+6.0*delta)/delta);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(0.5*(2.0*t1-2.0*t5-2.0*phi-3.0*delta*t8-3.0*delta)/delta);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441541679836E1-x*(0.15E1-f)-f)/(x*(a-3.0)+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t13;
  double t3;
  double t5;
  double t6;
  double t8;
  {
    t3 = log(phi);
    t5 = 1.0-phi;
    t6 = log(t5);
    t8 = phi*phi;
    t13 = log(rho);
    return(rho/delta*(phi*t3+t5*t6-t8/2.0)+0.15E1*phi*rho*t13+3.0*t5*rho*t13);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t3;
  double t5;
  double t8;
  {
    t3 = log(phi);
    t5 = log(1.0-phi);
    t8 = log(rho);
    return(rho/delta*(t3-t5-phi)-0.15E1*rho*t8);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho/delta*(1/phi+1/(1.0-phi)-1.0));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t12;
  double t13;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t12 = phi*delta;
    t13 = log(rho);
    return(0.5*(2.0*phi*t1+2.0*t6-2.0*t6*phi-1.0*t10-3.0*t12*t13-3.0*t12+6.0*
delta*t13+6.0*delta)/delta);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(0.5*(2.0*t1-2.0*t5-2.0*phi-3.0*t8*delta-3.0*delta)/delta);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441541679836E1-x*(0.15E1-f)-f)/(x*(a-3.0)+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t13;
  double t3;
  double t5;
  double t6;
  double t8;
  {
    t3 = log(phi);
    t5 = 1.0-phi;
    t6 = log(t5);
    t8 = phi*phi;
    t13 = log(rho);
    return(rho/delta*(phi*t3+t5*t6-t8/2.0)+0.15E1*phi*rho*t13+3.0*t5*rho*t13);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t3;
  double t5;
  double t8;
  {
    t3 = log(phi);
    t5 = log(1.0-phi);
    t8 = log(rho);
    return(rho/delta*(t3-t5-phi)-0.15E1*rho*t8);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho/delta*(1/phi+1/(1.0-phi)-1.0));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t11;
  double t12;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t11 = phi*delta;
    t12 = log(rho);
    return(-0.5*(-2.0*phi*t1-2.0*t6+2.0*t6*phi+t10+3.0*t11*t12+3.0*t11-6.0*
delta*t12-6.0*delta)/delta);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(-0.5*(-2.0*t1+2.0*t5+2.0*phi+3.0*delta*t8+3.0*delta)/delta);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441541679836E1-x*(0.15E1-f)-f)/(x*(a-3.0)+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t13;
  double t3;
  double t5;
  double t6;
  double t8;
  {
    t3 = log(phi);
    t5 = 1.0-phi;
    t6 = log(t5);
    t8 = phi*phi;
    t13 = log(rho);
    return(rho/delta*(phi*t3+t5*t6-t8/2.0)+0.15E1*phi*rho*t13+3.0*t5*rho*t13);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t3;
  double t5;
  double t8;
  {
    t3 = log(phi);
    t5 = log(1.0-phi);
    t8 = log(rho);
    return(rho/delta*(t3-t5-phi)-0.15E1*rho*t8);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho/delta*(1/phi+1/(1.0-phi)-1.0));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t12;
  double t13;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t12 = phi*delta;
    t13 = log(rho);
    return(0.5*(2.0*phi*t1+2.0*t6-2.0*t6*phi-1.0*t10-3.0*t12*t13-3.0*t12+6.0*
delta*t13+6.0*delta)/delta);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(-0.5*(-2.0*t1+2.0*t5+2.0*phi+3.0*t8*delta+3.0*delta)/delta);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441541679836E1-x*(0.15E1-f)-f)/(x*(a-3.0)+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t13;
  double t3;
  double t5;
  double t6;
  double t8;
  {
    t3 = log(phi);
    t5 = 1.0-phi;
    t6 = log(t5);
    t8 = phi*phi;
    t13 = log(rho);
    return(rho/delta*(phi*t3+t5*t6-t8/2.0)+0.15E1*phi*rho*t13+3.0*t5*rho*t13);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t3;
  double t5;
  double t8;
  {
    t3 = log(phi);
    t5 = log(1.0-phi);
    t8 = log(rho);
    return(rho/delta*(t3-t5-phi)-0.15E1*rho*t8);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho/delta*(1/phi+1/(1.0-phi)-1.0));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t12;
  double t13;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t12 = phi*delta;
    t13 = log(rho);
    return(0.5*(2.0*phi*t1+2.0*t6-2.0*t6*phi-1.0*t10-3.0*t12*t13-3.0*t12+6.0*
delta*t13+6.0*delta)/delta);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(-0.5*(-2.0*t1+2.0*t5+2.0*phi+3.0*t8*delta+3.0*delta)/delta);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441541679836E1-x*(0.15E1-f)-f)/(x*(a-3.0)+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t13;
  double t3;
  double t5;
  double t6;
  double t8;
  {
    t3 = log(phi);
    t5 = 1.0-phi;
    t6 = log(t5);
    t8 = phi*phi;
    t13 = log(rho);
    return(rho/delta*(phi*t3+t5*t6-t8/2.0)+0.15E1*phi*rho*t13+3.0*t5*rho*t13);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t3;
  double t5;
  double t8;
  {
    t3 = log(phi);
    t5 = log(1.0-phi);
    t8 = log(rho);
    return(rho/delta*(t3-t5-phi)-0.15E1*rho*t8);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho/delta*(1/phi+1/(1.0-phi)-1.0));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t12;
  double t13;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t12 = phi*delta;
    t13 = log(rho);
    return(0.5*(2.0*phi*t1+2.0*t6-2.0*t6*phi-1.0*t10-3.0*t12*t13-3.0*t12+6.0*
delta*t13+6.0*delta)/delta);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(-0.5*(-2.0*t1+2.0*t5+2.0*phi+3.0*t8*delta+3.0*delta)/delta);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441541679836E1-x*(0.15E1-f)-f)/(x*(a-3.0)+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t13;
  double t3;
  double t5;
  double t6;
  double t8;
  {
    t3 = log(phi);
    t5 = 1.0-phi;
    t6 = log(t5);
    t8 = phi*phi;
    t13 = log(rho);
    return(rho/delta*(phi*t3+t5*t6-t8/2.0)+0.15E1*phi*rho*t13+3.0*t5*rho*t13);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t3;
  double t5;
  double t8;
  {
    t3 = log(phi);
    t5 = log(1.0-phi);
    t8 = log(rho);
    return(rho/delta*(t3-t5-phi)-0.15E1*rho*t8);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho/delta*(1/phi+1/(1.0-phi)-1.0));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t12;
  double t13;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t12 = phi*delta;
    t13 = log(rho);
    return(0.5*(2.0*phi*t1+2.0*t6-2.0*t6*phi-1.0*t10-3.0*t12*t13-3.0*t12+6.0*
delta*t13+6.0*delta)/delta);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(0.5*(2.0*t1-2.0*t5-2.0*phi-3.0*t8*delta-3.0*delta)/delta);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441541679836E1-x*(0.15E1-f)-f)/(x*(a-3.0)+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t13;
  double t3;
  double t5;
  double t6;
  double t8;
  {
    t3 = log(phi);
    t5 = 1.0-phi;
    t6 = log(t5);
    t8 = phi*phi;
    t13 = log(rho);
    return(rho/delta*(phi*t3+t5*t6-t8/2.0)+0.15E1*phi*rho*t13+3.0*t5*rho*t13);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t3;
  double t5;
  double t8;
  {
    t3 = log(phi);
    t5 = log(1.0-phi);
    t8 = log(rho);
    return(rho/delta*(t3-t5-phi)-0.15E1*rho*t8);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho/delta*(1/phi+1/(1.0-phi)-1.0));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t12;
  double t13;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t12 = phi*delta;
    t13 = log(rho);
    return(0.5*(2.0*phi*t1+2.0*t6-2.0*t6*phi-1.0*t10-3.0*t12*t13-3.0*t12+6.0*
delta*t13+6.0*delta)/delta);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(0.5*(2.0*t1-2.0*t5-2.0*phi-3.0*t8*delta-3.0*delta)/delta);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441541679836E1-x*(0.15E1-f)-f)/(x*(a-3.0)+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t13;
  double t3;
  double t5;
  double t6;
  double t8;
  {
    t3 = log(phi);
    t5 = 1.0-phi;
    t6 = log(t5);
    t8 = phi*phi;
    t13 = log(rho);
    return(rho/delta*(phi*t3+t5*t6-t8/2.0)+0.15E1*phi*rho*t13+3.0*t5*rho*t13);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t3;
  double t5;
  double t8;
  {
    t3 = log(phi);
    t5 = log(1.0-phi);
    t8 = log(rho);
    return(rho/delta*(t3-t5-phi)-0.15E1*rho*t8);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho/delta*(1/phi+1/(1.0-phi)-1.0));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t12;
  double t13;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t12 = phi*delta;
    t13 = log(rho);
    return(0.5*(2.0*phi*t1+2.0*t6-2.0*t6*phi-1.0*t10-3.0*t12*t13-3.0*t12+6.0*
delta*t13+6.0*delta)/delta);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(0.5*(2.0*t1-2.0*t5-2.0*phi-3.0*t8*delta-3.0*delta)/delta);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441541679836E1-x*(0.15E1-f)-f)/(x*(a-3.0)+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t13;
  double t3;
  double t5;
  double t6;
  double t8;
  {
    t3 = log(phi);
    t5 = 1.0-phi;
    t6 = log(t5);
    t8 = phi*phi;
    t13 = log(rho);
    return(rho/delta*(phi*t3+t5*t6-t8/2.0)+0.15E1*phi*rho*t13+3.0*t5*rho*t13);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t3;
  double t5;
  double t8;
  {
    t3 = log(phi);
    t5 = log(1.0-phi);
    t8 = log(rho);
    return(rho/delta*(t3-t5-phi)-0.15E1*rho*t8);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho/delta*(1/phi+1/(1.0-phi)-1.0));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t12;
  double t13;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t12 = phi*delta;
    t13 = log(rho);
    return(0.5*(2.0*phi*t1+2.0*t6-2.0*t6*phi-1.0*t10-3.0*t12*t13-3.0*t12+6.0*
delta*t13+6.0*delta)/delta);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(0.5*(2.0*t1-2.0*t5-2.0*phi-3.0*t8*delta-3.0*delta)/delta);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441541679836E1-x*(0.15E1-f)-f)/(x*(a-3.0)+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t13;
  double t3;
  double t5;
  double t6;
  double t8;
  {
    t3 = log(phi);
    t5 = 1.0-phi;
    t6 = log(t5);
    t8 = phi*phi;
    t13 = log(rho);
    return(rho/delta*(phi*t3+t5*t6-t8/2.0)+0.15E1*phi*rho*t13+3.0*t5*rho*t13);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t3;
  double t5;
  double t8;
  {
    t3 = log(phi);
    t5 = log(1.0-phi);
    t8 = log(rho);
    return(rho/delta*(t3-t5-phi)-0.15E1*rho*t8);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho/delta*(1/phi+1/(1.0-phi)-1.0));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t12;
  double t13;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t12 = phi*delta;
    t13 = log(rho);
    return(0.5*(2.0*phi*t1+2.0*t6-2.0*t6*phi-1.0*t10-3.0*t12*t13-3.0*t12+6.0*
delta*t13+6.0*delta)/delta);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(0.5*(2.0*t1-2.0*t5-2.0*phi-3.0*t8*delta-3.0*delta)/delta);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441541679836E1-x*(0.15E1-f)-f)/(x*(a-3.0)+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t13;
  double t3;
  double t5;
  double t6;
  double t8;
  {
    t3 = log(phi);
    t5 = 1.0-phi;
    t6 = log(t5);
    t8 = phi*phi;
    t13 = log(rho);
    return(rho/delta*(phi*t3+t5*t6-t8/2.0)+0.15E1*phi*rho*t13+3.0*t5*rho*t13);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t3;
  double t5;
  double t8;
  {
    t3 = log(phi);
    t5 = log(1.0-phi);
    t8 = log(rho);
    return(rho/delta*(t3-t5-phi)-0.15E1*rho*t8);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho/delta*(1/phi+1/(1.0-phi)-1.0));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t12;
  double t13;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t12 = phi*delta;
    t13 = log(rho);
    return(0.5*(2.0*phi*t1+2.0*t6-2.0*t6*phi-1.0*t10-3.0*t12*t13-3.0*t12+6.0*
delta*t13+6.0*delta)/delta);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(-0.5*(-2.0*t1+2.0*t5+2.0*phi+3.0*t8*delta+3.0*delta)/delta);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441541679836E1-x*(0.15E1-f)-f)/(x*(a-3.0)+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t13;
  double t3;
  double t5;
  double t6;
  double t8;
  {
    t3 = log(phi);
    t5 = 1.0-phi;
    t6 = log(t5);
    t8 = phi*phi;
    t13 = log(rho);
    return(rho/delta*(phi*t3+t5*t6-t8/2.0)+0.15E1*phi*rho*t13+3.0*t5*rho*t13);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t3;
  double t5;
  double t8;
  {
    t3 = log(phi);
    t5 = log(1.0-phi);
    t8 = log(rho);
    return(rho/delta*(t3-t5-phi)-0.15E1*rho*t8);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho/delta*(1/phi+1/(1.0-phi)-1.0));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t12;
  double t13;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t12 = phi*delta;
    t13 = log(rho);
    return(0.5*(2.0*phi*t1+2.0*t6-2.0*t6*phi-1.0*t10-3.0*t12*t13-3.0*t12+6.0*
delta*t13+6.0*delta)/delta);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(0.5*(2.0*t1-2.0*t5-2.0*phi-3.0*delta*t8-3.0*delta)/delta);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441541679836E1-x*(0.15E1-f)-f)/(x*(a-3.0)+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t1;
  double t13;
  double t3;
  double t4;
  double t6;
  {
    t1 = log(phi);
    t3 = 1.0-phi;
    t4 = log(t3);
    t6 = phi*phi;
    t13 = log(rho);
    return(rho*(phi*t1+t3*t4-t6/2.0)/delta+0.15E1*phi*rho*t13+3.0*t3*rho*t13);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t1;
  double t3;
  double t8;
  {
    t1 = log(phi);
    t3 = log(1.0-phi);
    t8 = log(rho);
    return(rho*(t1-t3-phi)/delta-0.15E1*rho*t8);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho*(1/phi+1/(1.0-phi)-1.0)/delta);
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t12;
  double t13;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t12 = phi*delta;
    t13 = log(rho);
    return(0.5*(2.0*phi*t1+2.0*t6-2.0*t6*phi-1.0*t10-3.0*t12*t13-3.0*t12+6.0*
delta*t13+6.0*delta)/delta);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(-0.5*(-2.0*t1+2.0*t5+2.0*phi+3.0*t8*delta+3.0*delta)/delta);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441541679836E1-x*(0.15E1-f)-f)/(x*(a-3.0)+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t13;
  double t3;
  double t5;
  double t6;
  double t8;
  {
    t3 = log(phi);
    t5 = 1.0-phi;
    t6 = log(t5);
    t8 = phi*phi;
    t13 = log(rho);
    return(rho/delta*(phi*t3+t5*t6-t8/2.0)+0.15E1*phi*rho*t13+3.0*t5*rho*t13);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t3;
  double t5;
  double t8;
  {
    t3 = log(phi);
    t5 = log(1.0-phi);
    t8 = log(rho);
    return(rho/delta*(t3-t5-phi)-0.15E1*rho*t8);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho/delta*(1/phi+1/(1.0-phi)-1.0));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t12;
  double t13;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t12 = phi*delta;
    t13 = log(rho);
    return(0.5*(2.0*phi*t1+2.0*t6-2.0*t6*phi-1.0*t10-3.0*t12*t13-3.0*t12+6.0*
delta*t13+6.0*delta)/delta);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(0.5*(2.0*t1-2.0*t5-2.0*phi-3.0*delta*t8-3.0*delta)/delta);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441541679836E1-x*(0.15E1-f)-f)/(x*(a-3.0)+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t13;
  double t3;
  double t5;
  double t6;
  double t8;
  {
    t3 = log(phi);
    t5 = 1.0-phi;
    t6 = log(t5);
    t8 = phi*phi;
    t13 = log(rho);
    return(rho/delta*(phi*t3+t5*t6-t8/2.0)+0.15E1*phi*rho*t13+3.0*t5*rho*t13);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t3;
  double t5;
  double t8;
  {
    t3 = log(phi);
    t5 = log(1.0-phi);
    t8 = log(rho);
    return(rho/delta*(t3-t5-phi)-0.15E1*rho*t8);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho/delta*(1/phi+1/(1.0-phi)-1.0));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t12;
  double t13;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t12 = phi*delta;
    t13 = log(rho);
    return(0.5*(2.0*phi*t1+2.0*t6-2.0*t6*phi-1.0*t10-3.0*t12*t13-3.0*t12+6.0*
delta*t13+6.0*delta)/delta);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(0.5*(2.0*t1-2.0*t5-2.0*phi-3.0*delta*t8-3.0*delta)/delta);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441542E1-x*(b-f)-f)/(-0.15E1*x+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t13;
  double t3;
  double t5;
  double t6;
  double t8;
  {
    t3 = log(phi);
    t5 = 1.0-phi;
    t6 = log(t5);
    t8 = phi*phi;
    t13 = log(rho);
    return(rho/delta*(phi*t3+t5*t6-t8/2.0)+0.15E1*phi*rho*t13+3.0*t5*rho*t13);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t3;
  double t5;
  double t8;
  {
    t3 = log(phi);
    t5 = log(1.0-phi);
    t8 = log(rho);
    return(rho/delta*(t3-t5-phi)-0.15E1*rho*t8);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho/delta*(1/phi+1/(1.0-phi)-1.0));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t11;
  double t12;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t11 = phi*delta;
    t12 = log(rho);
    return(-0.5*(-2.0*phi*t1-2.0*t6+2.0*t6*phi+t10+3.0*t11*t12+3.0*t11-6.0*
delta*t12-6.0*delta)/delta);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(0.5*(2.0*t1-2.0*t5-2.0*phi-3.0*t8*delta-3.0*delta)/delta);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((Const-x*(b-f)-f)/(-0.115E2*x+13.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t13;
  double t3;
  double t5;
  double t6;
  double t8;
  {
    t3 = log(phi);
    t5 = 1.0-phi;
    t6 = log(t5);
    t8 = phi*phi;
    t13 = log(rho);
    return(rho/delta*(phi*t3+t5*t6-t8/2.0)+0.15E1*phi*rho*t13+13.0*t5*rho*t13);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t3;
  double t5;
  double t8;
  {
    t3 = log(phi);
    t5 = log(1.0-phi);
    t8 = log(rho);
    return(rho/delta*(t3-t5-phi)-0.115E2*rho*t8);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho/delta*(1/phi+1/(1.0-phi)-1.0));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t12;
  double t13;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t12 = phi*delta;
    t13 = log(rho);
    return(0.5*(2.0*phi*t1+2.0*t6-2.0*t6*phi-1.0*t10-23.0*t12*t13-23.0*t12+26.0
*delta*t13+26.0*delta)/delta);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.5*(23.0*phi-26.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(0.5*(2.0*t1-2.0*t5-2.0*phi-23.0*t8*delta-23.0*delta)/delta);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.115E2*phi*rho+13.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-46.0*phi+52.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((Const-x*(b-f)-f)/(-0.115E2*x+13.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t13;
  double t3;
  double t5;
  double t6;
  double t8;
  {
    t3 = log(phi);
    t5 = 1.0-phi;
    t6 = log(t5);
    t8 = phi*phi;
    t13 = log(rho);
    return(rho/delta*(phi*t3+t5*t6-t8/2.0)+0.15E1*phi*rho*t13+13.0*t5*rho*t13);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t3;
  double t5;
  double t8;
  {
    t3 = log(phi);
    t5 = log(1.0-phi);
    t8 = log(rho);
    return(rho/delta*(t3-t5-phi)-0.115E2*rho*t8);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho/delta*(1/phi+1/(1.0-phi)-1.0));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t12;
  double t13;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t12 = phi*delta;
    t13 = log(rho);
    return(0.5*(2.0*phi*t1+2.0*t6-2.0*t6*phi-1.0*t10-23.0*t12*t13-23.0*t12+26.0
*delta*t13+26.0*delta)/delta);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.5*(23.0*phi-26.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(0.5*(2.0*t1-2.0*t5-2.0*phi-23.0*t8*delta-23.0*delta)/delta);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.115E2*phi*rho+13.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-46.0*phi+52.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((Const-x*(b-f)-f)/(-0.115E2*x+13.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t13;
  double t3;
  double t5;
  double t6;
  double t8;
  {
    t3 = log(phi);
    t5 = 1.0-phi;
    t6 = log(t5);
    t8 = phi*phi;
    t13 = log(rho);
    return(rho/delta*(phi*t3+t5*t6-t8/2.0)+0.15E1*phi*rho*t13+13.0*t5*rho*t13);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t3;
  double t5;
  double t8;
  {
    t3 = log(phi);
    t5 = log(1.0-phi);
    t8 = log(rho);
    return(rho/delta*(t3-t5-phi)-0.115E2*rho*t8);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho/delta*(1/phi+1/(1.0-phi)-1.0));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t12;
  double t13;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t12 = phi*delta;
    t13 = log(rho);
    return(0.5*(2.0*phi*t1+2.0*t6-2.0*t6*phi-1.0*t10-23.0*t12*t13-23.0*t12+26.0
*delta*t13+26.0*delta)/delta);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.5*(23.0*phi-26.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(0.5*(2.0*t1-2.0*t5-2.0*phi-23.0*delta*t8-23.0*delta)/delta);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.115E2*phi*rho+13.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-46.0*phi+52.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441542E1-x*(b-f)-f)/(-0.15E1*x+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t13;
  double t3;
  double t5;
  double t6;
  double t8;
  {
    t3 = log(phi);
    t5 = 1.0-phi;
    t6 = log(t5);
    t8 = phi*phi;
    t13 = log(rho);
    return(rho/delta*(phi*t3+t5*t6-t8/2.0)+0.15E1*phi*rho*t13+3.0*t5*rho*t13);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t3;
  double t5;
  double t8;
  {
    t3 = log(phi);
    t5 = log(1.0-phi);
    t8 = log(rho);
    return(rho/delta*(t3-t5-phi)-0.15E1*rho*t8);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho/delta*(1/phi+1/(1.0-phi)-1.0));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t12;
  double t13;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t12 = phi*delta;
    t13 = log(rho);
    return(0.5*(2.0*phi*t1+2.0*t6-2.0*t6*phi-1.0*t10-3.0*t12*t13-3.0*t12+6.0*
delta*t13+6.0*delta)/delta);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(-0.5*(-2.0*t1+2.0*t5+2.0*phi+3.0*t8*delta+3.0*delta)/delta);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441542E1-x*(0.15E1-f)-f)/(x*(a-3.0)+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t13;
  double t3;
  double t5;
  double t6;
  double t8;
  {
    t3 = log(phi);
    t5 = 1.0-phi;
    t6 = log(t5);
    t8 = phi*phi;
    t13 = log(rho);
    return(rho/delta*(phi*t3+t5*t6-t8/2.0)+0.15E1*phi*rho*t13+3.0*t5*rho*t13);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t3;
  double t5;
  double t8;
  {
    t3 = log(phi);
    t5 = log(1.0-phi);
    t8 = log(rho);
    return(rho/delta*(t3-t5-phi)-0.15E1*rho*t8);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho/delta*(1/phi+1/(1.0-phi)-1.0));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t12;
  double t13;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t12 = phi*delta;
    t13 = log(rho);
    return(0.5*(2.0*phi*t1+2.0*t6-2.0*t6*phi-1.0*t10-3.0*t12*t13-3.0*t12+6.0*
delta*t13+6.0*delta)/delta);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(0.5*(2.0*t1-2.0*t5-2.0*phi-3.0*t8*delta-3.0*delta)/delta);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441542E1-x*(0.15E1-f)-f)/(x*(a-3.0)+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t13;
  double t3;
  double t5;
  double t6;
  double t8;
  {
    t3 = log(phi);
    t5 = 1.0-phi;
    t6 = log(t5);
    t8 = phi*phi;
    t13 = log(rho);
    return(rho/delta*(phi*t3+t5*t6-t8/2.0)+0.15E1*phi*rho*t13+3.0*t5*rho*t13);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t3;
  double t5;
  double t8;
  {
    t3 = log(phi);
    t5 = log(1.0-phi);
    t8 = log(rho);
    return(rho/delta*(t3-t5-phi)-0.15E1*rho*t8);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho/delta*(1/phi+1/(1.0-phi)-1.0));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t11;
  double t12;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t11 = phi*delta;
    t12 = log(rho);
    return(-0.5*(-2.0*phi*t1-2.0*t6+2.0*t6*phi+t10+3.0*t11*t12+3.0*t11-6.0*
delta*t12-6.0*delta)/delta);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(-0.5*(-2.0*t1+2.0*t5+2.0*phi+3.0*t8*delta+3.0*delta)/delta);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441542E1-x*(0.15E1-f)-f)/(x*(a-3.0)+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t13;
  double t3;
  double t5;
  double t6;
  double t8;
  {
    t3 = log(phi);
    t5 = 1.0-phi;
    t6 = log(t5);
    t8 = phi*phi;
    t13 = log(rho);
    return(rho/delta*(phi*t3+t5*t6-t8/2.0)+0.15E1*phi*rho*t13+3.0*t5*rho*t13);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t3;
  double t5;
  double t8;
  {
    t3 = log(phi);
    t5 = log(1.0-phi);
    t8 = log(rho);
    return(rho/delta*(t3-t5-phi)-0.15E1*rho*t8);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho/delta*(1/phi+1/(1.0-phi)-1.0));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t12;
  double t13;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t12 = phi*delta;
    t13 = log(rho);
    return(0.5*(2.0*phi*t1+2.0*t6-2.0*t6*phi-1.0*t10-3.0*t12*t13-3.0*t12+6.0*
delta*t13+6.0*delta)/delta);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(0.5*(2.0*t1-2.0*t5-2.0*phi-3.0*delta*t8-3.0*delta)/delta);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}

#include <math.h>
double evalRho(double x)
{
  {
    return(exp((0.2079441542E1-x*(0.15E1-f)-f)/(x*(a-3.0)+3.0)));
  }
}

#include <math.h>
double helmholtz(double rho,double phi)
{
  double t13;
  double t3;
  double t5;
  double t6;
  double t8;
  {
    t3 = log(phi);
    t5 = 1.0-phi;
    t6 = log(t5);
    t8 = phi*phi;
    t13 = log(rho);
    return(rho/delta*(phi*t3+t5*t6-t8/2.0)+0.15E1*phi*rho*t13+3.0*t5*rho*t13);
  }
}

#include <math.h>
double reactionSource(double rho,double phi)
{
  double t3;
  double t5;
  double t8;
  {
    t3 = log(phi);
    t5 = log(1.0-phi);
    t8 = log(rho);
    return(rho/delta*(t3-t5-phi)-0.15E1*rho*t8);
  }
}

#include <math.h>
double dphireactionSource(double rho,double phi)
{
  {
    return(rho/delta*(1/phi+1/(1.0-phi)-1.0));
  }
}

#include <math.h>
double chemicalPotential(double rho,double phi)
{
  double t1;
  double t10;
  double t11;
  double t12;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t11 = phi*delta;
    t12 = log(rho);
    return(-0.5*(-2.0*phi*t1-2.0*t6+2.0*t6*phi+t10+3.0*t11*t12+3.0*t11-6.0*
delta*t12-6.0*delta)/delta);
  }
}

double drhochemicalPotential(double rho,double phi)
{
  {
    return(-0.15E1*(phi-2.0)/rho);
  }
}

#include <math.h>
double dphichemicalPotential(double rho,double phi)
{
  double t1;
  double t5;
  double t8;
  {
    t1 = log(phi);
    t5 = log(1.0-1.0*phi);
    t8 = log(rho);
    return(-0.5*(-2.0*t1+2.0*t5+2.0*phi+3.0*delta*t8+3.0*delta)/delta);
  }
}

double pressure(double rho,double phi)
{
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }
}

#include <math.h>
double a(double rho,double phi)
{
  double t3;
  {
    t3 = sqrt(-6.0*phi+12.0);
    return(0.5*t3);
  }
}



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
    return(exp((0.4000000004E1-18.0*t3+45.0*t2-30.0*t6)/(-0.9E1*t3+0.225E2*t2
-0.15E2*t6+3.0)));
  }
}


inline double helmholtz ( double rho ,double phi ) const
{
  double t1;
  double t10;
  double t11;
  double t13;
  double t16;
  double t17;
  double t3;
  double t8;
  {
    t1 = phi*phi;
    t3 = pow(1.0-phi,2.0);
    t8 = t1*t1;
    t10 = 6.0*t8*phi;
    t11 = 15.0*t8;
    t13 = 10.0*t1*phi;
    t16 = log(rho);
    t17 = rho*t16;
    return(2.0*t1*t3/delta_+(t10-t11+t13)*(-0.25E1*rho+0.15E1*t17+0.15E1)+(1.0-
t10+t11-t13)*(-7.0*rho+3.0*t17+0.995940263E1));
  }
}


inline double reactionSource ( double rho ,double phi ) const
{
  double t1;
  double t11;
  double t16;
  double t18;
  double t19;
  double t2;
  double t4;
  double t7;
  {
    t1 = 1.0-phi;
    t2 = t1*t1;
    t4 = 1/delta_;
    t7 = phi*phi;
    t11 = t7*t7;
    t16 = 30.0*t11-60.0*t7*phi+30.0*t7;
    t18 = log(rho);
    t19 = rho*t18;
    return(4.0*phi*t2*t4-4.0*t7*t1*t4+t16*(-0.25E1*rho+0.15E1*t19+0.15E1)-t16*(
-7.0*rho+3.0*t19+0.995940263E1));
  }
}


inline double dphireactionSource ( double rho ,double phi ) const
{
  double t1;
  double t16;
  double t18;
  double t19;
  double t2;
  double t3;
  double t9;
  {
    t1 = 1.0-phi;
    t2 = t1*t1;
    t3 = 1/delta_;
    t9 = phi*phi;
    t16 = 120.0*t9*phi-180.0*t9+60.0*phi;
    t18 = log(rho);
    t19 = rho*t18;
    return(4.0*t2*t3-16.0*phi*t1*t3+4.0*t9*t3+t16*(-0.25E1*rho+0.15E1*t19+
0.15E1)-t16*(-7.0*rho+3.0*t19+0.995940263E1));
  }
}


inline double chemicalPotential ( double rho ,double phi ) const
{
  double t1;
  double t11;
  double t2;
  double t3;
  double t5;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = t2*phi;
    t5 = log(rho);
    t11 = t1*phi;
    return(18.0*t3-9.0*t3*t5-45.0*t2+0.225E2*t2*t5+30.0*t11-15.0*t11*t5-4.0+3.0
*t5);
  }
}

inline double drhochemicalPotential ( double rho ,double phi ) const
{
  double t1;
  double t2;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    return(-0.15E1*(6.0*t2*phi-15.0*t2+10.0*t1*phi-2.0)/rho);
  }
}


inline double dphichemicalPotential ( double rho ,double phi ) const
{
  double t1;
  double t2;
  double t4;
  double t7;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = log(rho);
    t7 = t1*phi;
    return(90.0*t2-45.0*t2*t4-180.0*t7+90.0*t7*t4+90.0*t1-45.0*t1*t4);
  }
}


inline double pressure ( double rho ,double phi ) const
{
  double t1;
  double t10;
  double t11;
  double t2;
  double t27;
  double t4;
  double t5;
  double t7;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t7 = 10.0*t1*phi;
    t10 = log(rho);
    t11 = rho*t10;
    t27 = pow(1.0-phi,2.0);
    return((t4-t5+t7)*(0.25E1*rho-0.15E1*t11-0.15E1+rho*(-1.0+0.15E1*t10))+(1.0
-t4+t5-t7)*(7.0*rho-3.0*t11-0.995940263E1+rho*(-4.0+3.0*t10))-2.0*t1*t27/delta_
);
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
    t9 = sqrt(-36.0*t2*phi+90.0*t2-60.0*t1*phi+12.0);
    return(0.5*t9);
  }
}


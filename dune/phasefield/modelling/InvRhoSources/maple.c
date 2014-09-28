
double solproc1(double x)
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
    return((0.113810037E2-0.5075641578E2*t3+0.1268910394E3*t2-0.845940263E2*t6)
/(-0.9E1*t3+0.225E2*t2-0.15E2*t6+3.0));
  }
}


inline double helmholtz ( double rho ,double phi ) const
{
  double t1;
  double t12;
  double t13;
  double t14;
  double t17;
  double t18;
  double t2;
  double t4;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = t1*phi;
    t12 = 6.0*t2*phi;
    t13 = 15.0*t2;
    t14 = 10.0*t4;
    t17 = log(rho);
    t18 = rho*t17;
    return(rho*(2.0*t2-4.0*t4+2.0*t1)/delta_+(t12-t13+t14)*(-0.25E1*rho+0.15E1*
t18+1.0)+(1.0-t12+t13-t14)*(-7.0*rho+3.0*t18+0.945940263E1));
  }
}


inline double reactionSource ( double rho ,double phi ) const
{
  double t1;
  double t10;
  double t14;
  double t16;
  double t17;
  double t2;
  {
    t1 = phi*phi;
    t2 = t1*phi;
    t10 = t1*t1;
    t14 = 30.0*t10-60.0*t2+30.0*t1;
    t16 = log(rho);
    t17 = rho*t16;
    return(rho*(8.0*t2-12.0*t1+4.0*phi)/delta_+t14*(1.0-0.25E1*rho+0.15E1*t17)-
t14*(0.945940263E1-7.0*rho+3.0*t17));
  }
}


inline double dphireactionSource ( double rho ,double phi ) const
{
  double t1;
  double t12;
  double t14;
  double t15;
  {
    t1 = phi*phi;
    t12 = 120.0*t1*phi-180.0*t1+60.0*phi;
    t14 = log(rho);
    t15 = rho*t14;
    return(rho*(24.0*t1-24.0*phi+4.0)/delta_+t12*(-0.25E1*rho+0.15E1*t15+1.0)-
t12*(-7.0*rho+3.0*t15+0.945940263E1));
  }
}


inline double chemicalPotential ( double rho ,double phi ) const
{
  double t1;
  double t11;
  double t12;
  double t13;
  double t15;
  double t2;
  double t4;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = t1*phi;
    t11 = 6.0*t2*phi;
    t12 = 15.0*t2;
    t13 = 10.0*t4;
    t15 = log(rho);
    return((2.0*t2-4.0*t4+2.0*t1)/delta_+(t11-t12+t13)*(-1.0+0.15E1*t15)+(1.0-
t11+t12-t13)*(-4.0+3.0*t15));
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
  double t10;
  double t14;
  double t5;
  double t7;
  {
    t1 = phi*phi;
    t5 = t1*phi*delta_;
    t7 = log(rho);
    t10 = t1*delta_;
    t14 = phi*delta_;
    return(-1.0*phi*(-8.0*t1+12.0*phi-4.0-90.0*t5+45.0*t5*t7+180.0*t10-90.0*t10
*t7-90.0*t14+45.0*t14*t7)/delta_);
  }
}

inline double pressure ( double rho ,double phi ) const
{
  double t1;
  double t10;
  double t2;
  double t3;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = t2*phi;
    t10 = t1*phi;
    return(-9.0*t3*rho+0.5075641578E2*t3+0.225E2*t2*rho-0.1268910394E3*t2-15.0*
t10*rho+0.845940263E2*t10+3.0*rho-0.945940263E1);
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


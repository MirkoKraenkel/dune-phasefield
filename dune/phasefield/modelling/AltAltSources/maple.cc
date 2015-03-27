double evalRho(double x)
{
  {
    return(0.18708E1*x+0.1292);
  }
}


inline double helmholtz ( double rho ,double phi ) const
{
  double t17;
  double t18;
  double t19;
  double t2;
  double t22;
  double t3;
  double t30;
  double t6;
  {
    t2 = phi*phi;
    t3 = t2*t2;
    t6 = t2*phi;
    t17 = 6.0*t3*phi;
    t18 = 15.0*t3;
    t19 = 10.0*t6;
    t22 = log(0.773993808E1*rho);
    t30 = log(0.5*rho);
    return(2.0*rho*A_*(t3+(2.0*alpha_-2.0)*t6+(-3.0*alpha_+1.0)*t2+alpha_)/delta_+
beta_*((t17-t18+t19)*(0.5*rho*(t22-1.0)+0.646E-1)+(1.0-t17+t18-t19)*(0.5*rho*(
t30-1.0)+0.1E1)));
  }
}


inline double reactionSource ( double rho ,double phi ) const
{
  double t18;
  double t2;
  double t22;
  double t24;
  double t3;
  double t32;
  {
    t2 = phi*phi;
    t3 = t2*phi;
    t18 = t2*t2;
    t22 = 30.0*t18-60.0*t3+30.0*t2;
    t24 = log(0.773993808E1*rho);
    t32 = log(0.5*rho);
    return(2.0*rho*A_*(4.0*t3+3.0*(2.0*alpha_-2.0)*t2+2.0*(-3.0*alpha_+1.0)*phi)/
delta_+beta_*(t22*(0.5*rho*(t24-1.0)+0.646E-1)-t22*(0.5*rho*(t32-1.0)+0.1E1)));
  }
}


inline double dphireactionSource ( double rho ,double phi ) const
{
  double t18;
  double t2;
  double t20;
  double t28;
  {
    t2 = phi*phi;
    t18 = 120.0*t2*phi-180.0*t2+60.0*phi;
    t20 = log(0.773993808E1*rho);
    t28 = log(0.5*rho);
    return(2.0*rho*A_*(12.0*t2+6.0*(2.0*alpha_-2.0)*phi-6.0*alpha_+2.0)/delta_+beta_
*(t18*(0.5*rho*(t20-1.0)+0.646E-1)-t18*(0.5*rho*(t28-1.0)+0.1E1)));
  }
}


inline double chemicalPotential ( double rho ,double phi ) const
{
  double t1;
  double t10;
  double t16;
  double t2;
  double t25;
  double t28;
  double t5;
  double t6;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t5 = t1*phi;
    t6 = A_*t5;
    t10 = A_*t1;
    t16 = beta_*delta_;
    t25 = log(rho);
    t28 = 0.2E11*A_*t2+0.4E11*t6*alpha_-0.4E11*t6-0.6E11*t10*alpha_+0.2E11*t10+
0.2E11*A_*alpha_+0.8218622606E11*t16*t2*phi-0.2054655651E12*t16*t2+
0.1369770434E12*t16*t5-3465735903.0*t16+5000000000.0*t16*t25;
    return(0.1E-9*t28/delta_);
  }
}

inline double drhochemicalPotential ( double rho ,double phi ) const
{
  {
    return(0.5*beta_/rho);
  }
}

inline double dphichemicalPotential ( double rho ,double phi ) const
{
  double t1;
  double t4;
  {
    t1 = phi*phi;
    t4 = A_*phi;
    return(0.1E-7*phi*(800000000.0*A_*t1+1200000000.0*t4*alpha_-1200000000.0*t4
-1200000000.0*A_*alpha_+400000000.0*A_+4109311303.0*beta_*t1*phi*delta_-8218622606.0
*beta_*t1*delta_+4109311303.0*beta_*phi*delta_)/delta_);
  }
}


inline double pressure ( double rho ,double phi ) const
{
  double t1;
  double t10;
  double t2;
  double t20;
  double t4;
  double t5;
  double t7;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t7 = 10.0*t1*phi;
    t10 = log(0.773993808E1*rho);
    t20 = log(0.5*rho);
    return(beta_*((t4-t5+t7)*(-0.5*rho*(t10-1.0)-0.646E-1+0.5*rho*t10)+(1.0-t4+
t5-t7)*(-0.5*rho*(t20-1.0)-0.1E1+0.5*rho*t20)));
  }
}


inline double a ( double rho ,double phi ) const
{
  double t1;
  {
    t1 = sqrt(beta_);
    return(0.7071067812*t1);
  }
}


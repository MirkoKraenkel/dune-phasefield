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
beta_*((t17-t18+t19)*(10.0*rho*(t22-1.0)+0.1292E1)+(1.0-t17+t18-t19)*(5.0*rho*(
t30-1.0)+10.0)));
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
delta_+beta_*(t22*(10.0*rho*(t24-1.0)+0.1292E1)-t22*(5.0*rho*(t32-1.0)+10.0)));
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
*(t18*(10.0*rho*(t20-1.0)+0.1292E1)-t18*(5.0*rho*(t28-1.0)+10.0)));
  }
}


inline double chemicalPotential ( double rho ,double phi ) const
{
  double t1;
  double t10;
  double t16;
  double t17;
  double t2;
  double t20;
  double t37;
  double t5;
  double t6;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t5 = t1*phi;
    t6 = A_*t5;
    t10 = A_*t1;
    t16 = beta_*delta_;
    t17 = t2*phi;
    t20 = log(rho);
    t37 = 2000000000.0*A_*t2+4000000000.0*t6*alpha_-4000000000.0*t6-6000000000.0*
t10*alpha_+2000000000.0*t10+2000000000.0*A_*alpha_+0.1435780367E12*t16*t17+0.3E11*
t16*t17*t20-0.3589450917E12*t16*t2-0.75E11*t16*t2*t20+0.2392967278E12*t16*t5+
0.5E11*t16*t5*t20-3465735903.0*t16+5000000000.0*t16*t20;
    return(0.1E-8*t37/delta_);
  }
}

inline double drhochemicalPotential ( double rho ,double phi ) const
{
  double t1;
  double t2;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    return(5.0*beta_*(6.0*t2*phi-15.0*t2+10.0*t1*phi+1.0)/rho);
  }
}


inline double dphichemicalPotential ( double rho ,double phi ) const
{
  double t1;
  double t12;
  double t15;
  double t16;
  double t19;
  double t24;
  double t29;
  double t4;
  {
    t1 = phi*phi;
    t4 = A_*phi;
    t12 = beta_*t1*phi;
    t15 = log(rho);
    t16 = delta_*t15;
    t19 = beta_*t1;
    t24 = beta_*phi;
    t29 = 16000000.0*A_*t1+24000000.0*t4*alpha_-24000000.0*t4-24000000.0*A_*alpha_+
8000000.0*A_+1435780367.0*t12*delta_+300000000.0*t12*t16-2871560734.0*t19*delta_
-600000000.0*t19*t16+1435780367.0*t24*delta_+300000000.0*t24*t16;
    return(0.5E-6*phi*t29/delta_);
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
    return(beta_*((t4-t5+t7)*(-10.0*rho*(t10-1.0)-0.1292E1+10.0*rho*t10)+(1.0-t4
+t5-t7)*(-5.0*rho*(t20-1.0)-10.0+5.0*rho*t20)));
  }
}


inline double a ( double rho ,double phi ) const
{
  double t1;
  double t10;
  double t2;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t10 = sqrt(beta_*(6.0*t2*phi-15.0*t2+10.0*t1*phi+1.0));
    return(0.2236067977E1*t10);
  }
}


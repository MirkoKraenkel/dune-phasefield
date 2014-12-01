
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
  double t17;
  double t18;
  double t19;
  double t2;
  double t22;
  double t23;
  double t3;
  double t6;
  {
    t2 = phi*phi;
    t3 = t2*t2;
    t6 = t2*phi;
    t17 = 6.0*t3*phi;
    t18 = 15.0*t3;
    t19 = 10.0*t6;
    t22 = log(rho);
    t23 = rho*t22;
    return(2.0*rho*A_*(t3+(2.0*alpha_-2.0)*t6+(-3.0*alpha_+1.0)*t2+alpha_)/delta_+
beta_*((t17-t18+t19)*(-0.25E1*rho+0.15E1*t23+1.0)+(1.0-t17+t18-t19)*(-7.0*rho+
3.0*t23+0.945940263E1)));
  }
}


inline double reactionSource ( double rho ,double phi ) const
{
  double t18;
  double t2;
  double t22;
  double t24;
  double t25;
  double t3;
  {
    t2 = phi*phi;
    t3 = t2*phi;
    t18 = t2*t2;
    t22 = 30.0*t18-60.0*t3+30.0*t2;
    t24 = log(rho);
    t25 = rho*t24;
    return(2.0*rho*A_*(4.0*t3+3.0*(2.0*alpha_-2.0)*t2+2.0*(-3.0*alpha_+1.0)*phi)/
delta_+beta_*(t22*(-0.25E1*rho+0.15E1*t25+1.0)-t22*(-7.0*rho+3.0*t25+
0.945940263E1)));
  }
}


inline double dphireactionSource ( double rho ,double phi ) const
{
  double t18;
  double t2;
  double t20;
  double t21;
  {
    t2 = phi*phi;
    t18 = 120.0*t2*phi-180.0*t2+60.0*phi;
    t20 = log(rho);
    t21 = rho*t20;
    return(2.0*rho*A_*(12.0*t2+6.0*(2.0*alpha_-2.0)*phi-6.0*alpha_+2.0)/delta_+beta_
*(t18*(-0.25E1*rho+0.15E1*t21+1.0)-t18*(-7.0*rho+3.0*t21+0.945940263E1)));
  }
}


inline double chemicalPotential ( double rho ,double phi ) const
{
  double t1;
  double t16;
  double t17;
  double t18;
  double t2;
  double t20;
  double t5;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t5 = t1*phi;
    t16 = 6.0*t2*phi;
    t17 = 15.0*t2;
    t18 = 10.0*t5;
    t20 = log(rho);
    return(2.0*A_*(t2+(2.0*alpha_-2.0)*t5+(-3.0*alpha_+1.0)*t1+alpha_)/delta_+beta_*(
(t16-t17+t18)*(-1.0+0.15E1*t20)+(1.0-t16+t17-t18)*(-4.0+3.0*t20)));
  }
}

inline double drhochemicalPotential ( double rho ,double phi ) const
{
  double t1;
  double t2;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    return(-0.15E1*beta_*(6.0*t2*phi-15.0*t2+10.0*t1*phi-2.0)/rho);
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
    t29 = 8.0*A_*t1+12.0*t4*alpha_-12.0*t4-12.0*A_*alpha_+4.0*A_+90.0*t12*delta_-45.0
*t12*t16-180.0*t19*delta_+90.0*t19*t16+90.0*t24*delta_-45.0*t24*t16;
    return(phi*t29/delta_);
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
    return(-0.1E-7*beta_*(900000000.0*t3*rho-5075641578.0*t3-2250000000.0*t2*rho
+0.1268910394E11*t2+1500000000.0*t10*rho-8459402630.0*t10-300000000.0*rho+
945940263.0));
  }
}


inline double a ( double rho ,double phi ) const
{
  double t1;
  double t11;
  double t2;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t11 = sqrt(-1.0*beta_*(6.0*t2*phi-15.0*t2+10.0*t1*phi-2.0));
    return(0.1224744871E1*t11);
  }
}


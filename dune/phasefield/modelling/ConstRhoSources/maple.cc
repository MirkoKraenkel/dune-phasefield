
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
  double t15;
  double t16;
  double t17;
  double t2;
  double t20;
  double t21;
  double t5;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t5 = t1*phi;
    t15 = 6.0*t2*phi;
    t16 = 15.0*t2;
    t17 = 10.0*t5;
    t20 = log(rho);
    t21 = rho*t20;
    return((2.0*t2+2.0*(2.0*alpha_-2.0)*t5+2.0*(-3.0*alpha_+1.0)*t1+2.0*alpha_)/
delta+beta_*((t15-t16+t17)*(-0.25E1*rho+0.15E1*t21+0.15E1)+(1.0-t15+t16-t17)*(
-7.0*rho+3.0*t21+0.995940263E1)));
  }
}


inline double reactionSource ( double rho ,double phi ) const
{
  double t1;
  double t15;
  double t19;
  double t2;
  double t21;
  double t22;
  {
    t1 = phi*phi;
    t2 = t1*phi;
    t15 = t1*t1;
    t19 = 30.0*t15-60.0*t2+30.0*t1;
    t21 = log(rho);
    t22 = rho*t21;
    return((8.0*t2+6.0*(2.0*alpha_-2.0)*t1+4.0*(-3.0*alpha_+1.0)*phi)/delta+beta_*
(t19*(-0.25E1*rho+0.15E1*t22+0.15E1)-t19*(-7.0*rho+3.0*t22+0.995940263E1)));
  }
}


inline double dphireactionSource ( double rho ,double phi ) const
{
  double t1;
  double t15;
  double t17;
  double t18;
  {
    t1 = phi*phi;
    t15 = 120.0*t1*phi-180.0*t1+60.0*phi;
    t17 = log(rho);
    t18 = rho*t17;
    return((24.0*t1+12.0*(2.0*alpha_-2.0)*phi-12.0*alpha_+4.0)/delta+beta_*(t15*(
-0.25E1*rho+0.15E1*t18+0.15E1)-t15*(-7.0*rho+3.0*t18+0.995940263E1)));
  }
}


inline double chemicalPotential ( double rho ,double phi ) const
{
  double t1;
  double t13;
  double t14;
  double t2;
  double t3;
  double t4;
  double t6;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = t2*phi;
    t4 = beta_*t3;
    t6 = log(rho);
    t9 = beta_*t2;
    t13 = t1*phi;
    t14 = beta_*t13;
    return(-6.0*t4+9.0*t4*t6+15.0*t9-0.225E2*t9*t6-10.0*t14+15.0*t14*t6-4.0+3.0
*t6+24.0*t3-18.0*t3*t6-60.0*t2+45.0*t2*t6+40.0*t13-30.0*t13*t6);
  }
}

inline double drhochemicalPotential ( double rho ,double phi ) const
{
  double t1;
  double t2;
  double t3;
  double t8;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = t2*phi;
    t8 = t1*phi;
    return(0.15E1*(6.0*beta_*t3-15.0*beta_*t2+10.0*beta_*t8+2.0-12.0*t3+30.0*t2
-20.0*t8)/rho);
  }
}


inline double dphichemicalPotential ( double rho ,double phi ) const
{
  double t1;
  double t13;
  double t2;
  double t3;
  double t5;
  double t8;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = beta_*t2;
    t5 = log(rho);
    t8 = t1*phi;
    t9 = beta_*t8;
    t13 = beta_*t1;
    return(-30.0*t3+45.0*t3*t5+60.0*t9-90.0*t9*t5-30.0*t13+45.0*t13*t5+120.0*t2
-90.0*t2*t5-240.0*t8+180.0*t8*t5+120.0*t1-90.0*t1*t5);
  }
}


inline double pressure ( double rho ,double phi ) const
{
  double t1;
  double t11;
  double t12;
  double t2;
  double t4;
  double t5;
  double t6;
  double t7;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t6 = t1*phi;
    t7 = 10.0*t6;
    t11 = log(rho);
    t12 = rho*t11;
    return(beta_*(t4-t5+t7)*(0.25E1*rho-0.15E1*t12-0.15E1+rho*(-1.0+0.15E1*t11))
+(1.0-t4+t5-t7)*(7.0*rho-3.0*t12-0.995940263E1+rho*(-4.0+3.0*t11))-(2.0*t2+2.0*
(2.0*alpha_-2.0)*t6+2.0*(-3.0*alpha_+1.0)*t1+2.0*alpha_)/delta);
  }
}


inline double a ( double rho ,double phi ) const
{
  double t1;
  double t15;
  double t2;
  double t3;
  double t8;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = t2*phi;
    t8 = t1*phi;
    t15 = sqrt(36.0*beta_*t3-90.0*beta_*t2+60.0*beta_*t8+12.0-72.0*t3+180.0*t2
-120.0*t8);
    return(0.5*t15);
  }
}


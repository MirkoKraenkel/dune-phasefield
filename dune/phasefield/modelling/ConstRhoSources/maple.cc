
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
  double t16;
  double t17;
  double t18;
  double t2;
  double t21;
  double t22;
  double t5;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t5 = t1*phi;
    t16 = 6.0*t2*phi;
    t17 = 15.0*t2;
    t18 = 10.0*t5;
    t21 = log(rho);
    t22 = rho*t21;
    return(2.0*A_*(t2+(2.0*alpha_-2.0)*t5+(-3.0*alpha_+1.0)*t1+alpha_)/delta_+beta_*(
(t16-t17+t18)*(-0.25E1*rho+0.15E1*t22+0.2921601062E1)+(1.0-t16+t17-t18)*(-7.0*
rho+3.0*t22+0.1138100368E2)));
  }
}


inline double reactionSource ( double rho ,double phi ) const
{
  double t1;
  double t17;
  double t2;
  double t21;
  double t23;
  double t24;
  {
    t1 = phi*phi;
    t2 = t1*phi;
    t17 = t1*t1;
    t21 = 30.0*t17-60.0*t2+30.0*t1;
    t23 = log(rho);
    t24 = rho*t23;
    return(2.0*A_*(4.0*t2+3.0*(2.0*alpha_-2.0)*t1+2.0*(-3.0*alpha_+1.0)*phi)/delta_
+beta_*(t21*(-0.25E1*rho+0.15E1*t24+0.2921601062E1)-t21*(-7.0*rho+3.0*t24+
0.1138100368E2)));
  }
}


inline double dphireactionSource ( double rho ,double phi ) const
{
  double t1;
  double t17;
  double t19;
  double t20;
  {
    t1 = phi*phi;
    t17 = 120.0*t1*phi-180.0*t1+60.0*phi;
    t19 = log(rho);
    t20 = rho*t19;
    return(2.0*A_*(12.0*t1+6.0*(2.0*alpha_-2.0)*phi-6.0*alpha_+2.0)/delta_+beta_*(
t17*(-0.25E1*rho+0.15E1*t20+0.2921601062E1)-t17*(-7.0*rho+3.0*t20+
0.1138100368E2)));
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
    return(beta_*(t4-t5+t7)*(0.25E1*rho-0.15E1*t12-0.2921601062E1+rho*(-1.0+
0.15E1*t11))+(1.0-t4+t5-t7)*(7.0*rho-3.0*t12-0.1138100368E2+rho*(-4.0+3.0*t11))
-2.0*A_*(t2+(2.0*alpha_-2.0)*t6+(-3.0*alpha_+1.0)*t1+alpha_)/delta_);
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


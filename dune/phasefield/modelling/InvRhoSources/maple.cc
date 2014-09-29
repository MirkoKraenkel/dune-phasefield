
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
    return(rho*(2.0*t2+2.0*(2.0*alpha_-2.0)*t5+2.0*(-3.0*alpha_+1.0)*t1+2.0*alpha_
)/delta_+beta_*((t16-t17+t18)*(-0.25E1*rho+0.15E1*t22+1.0)+(1.0-t16+t17-t18)*(
-7.0*rho+3.0*t22+0.945940263E1)));
  }
}


inline double reactionSource ( double rho ,double phi ) const
{
  double t1;
  double t16;
  double t2;
  double t20;
  double t22;
  double t23;
  {
    t1 = phi*phi;
    t2 = t1*phi;
    t16 = t1*t1;
    t20 = 30.0*t16-60.0*t2+30.0*t1;
    t22 = log(rho);
    t23 = rho*t22;
    return(rho*(8.0*t2+6.0*(2.0*alpha_-2.0)*t1+4.0*(-3.0*alpha_+1.0)*phi)/delta_+
beta_*(t20*(-0.25E1*rho+0.15E1*t23+1.0)-t20*(-7.0*rho+3.0*t23+0.945940263E1)));
  }
}


inline double dphireactionSource ( double rho ,double phi ) const
{
  double t1;
  double t16;
  double t18;
  double t19;
  {
    t1 = phi*phi;
    t16 = 120.0*t1*phi-180.0*t1+60.0*phi;
    t18 = log(rho);
    t19 = rho*t18;
    return(rho*(24.0*t1+12.0*(2.0*alpha_-2.0)*phi-12.0*alpha_+4.0)/delta_+beta_*(
t16*(-0.25E1*rho+0.15E1*t19+1.0)-t16*(-7.0*rho+3.0*t19+0.945940263E1)));
  }
}


inline double chemicalPotential ( double rho ,double phi ) const
{
  double t1;
  double t15;
  double t16;
  double t17;
  double t2;
  double t20;
  double t5;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t5 = t1*phi;
    t15 = 6.0*t2*phi;
    t16 = 15.0*t2;
    t17 = 10.0*t5;
    t20 = log(rho);
    return((2.0*t2+2.0*(2.0*alpha_-2.0)*t5+2.0*(-3.0*alpha_+1.0)*t1+2.0*alpha_)/
delta_+beta_*(t15-t16+t17)*(-1.0+0.15E1*t20)+(1.0-t15+t16-t17)*(-4.0+3.0*t20));
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
  double t11;
  double t12;
  double t15;
  double t20;
  double t25;
  double t29;
  double t33;
  double t37;
  double t7;
  double t8;
  {
    t1 = phi*phi;
    t7 = t1*phi;
    t8 = beta_*t7;
    t11 = log(rho);
    t12 = delta_*t11;
    t15 = beta_*t1;
    t20 = beta_*phi;
    t25 = t7*delta_;
    t29 = t1*delta_;
    t33 = phi*delta_;
    t37 = 8.0*t1+12.0*phi*alpha_-12.0*phi-12.0*alpha_+4.0-30.0*t8*delta_+45.0*t8*
t12+60.0*t15*delta_-90.0*t15*t12-30.0*t20*delta_+45.0*t20*t12+120.0*t25-90.0*t25*
t11-240.0*t29+180.0*t29*t11+120.0*t33-90.0*t33*t11;
    return(phi*t37/delta_);
  }
}

inline double pressure ( double rho ,double phi ) const
{
  double t1;
  double t12;
  double t13;
  double t2;
  double t3;
  double t4;
  double t8;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = t2*phi;
    t4 = beta_*t3;
    t8 = beta_*t2;
    t12 = t1*phi;
    t13 = beta_*t12;
    return(9.0*t4*rho-6.0*t4-0.225E2*t8*rho+15.0*t8+15.0*t13*rho-10.0*t13+3.0*
rho-0.945940263E1-18.0*t3*rho+0.5675641578E2*t3+45.0*t2*rho-0.1418910394E3*t2
-30.0*t12*rho+0.945940263E2*t12);
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


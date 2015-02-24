
double solproc1(double x)
{
  double t1;
  double t2;
  double t4;
  double t5;
  double t7;
  {
    t1 = x*x;
    t2 = t1*t1;
    t4 = 0.3E1*t2*x;
    t5 = 0.75E1*t2;
    t7 = 0.5E1*t1*x;
    return((0.2E1-t4+t5-t7)/(1.0+t4-t5+t7));
  }
}


inline double helmholtz ( double rho ,double phi ) const
{
  double t17;
  double t18;
  double t19;
  double t2;
  double t21;
  double t24;
  double t25;
  double t28;
  double t3;
  double t6;
  {
    t2 = phi*phi;
    t3 = t2*t2;
    t6 = t2*phi;
    t17 = 6.0*t3*phi;
    t18 = 15.0*t3;
    t19 = 10.0*t6;
    t21 = log(2.0);
    t24 = log(rho);
    t25 = rho*t24;
    t28 = (t21-0.15E-9)*rho;
    return(2.0*rho*A_*(t3+(2.0*alpha_-2.0)*t6+(-3.0*alpha_+1.0)*t2+alpha_)/delta_+
beta_*((t17-t18+t19)*((t21-0.15E1)*rho+0.15E1*t25-t28+0.15E1)+(1.0-t17+t18-t19)*
(-rho+t25+0.2E1-t28)));
  }
}


inline double reactionSource ( double rho ,double phi ) const
{
  double t18;
  double t2;
  double t22;
  double t23;
  double t26;
  double t27;
  double t3;
  double t30;
  {
    t2 = phi*phi;
    t3 = t2*phi;
    t18 = t2*t2;
    t22 = 30.0*t18-60.0*t3+30.0*t2;
    t23 = log(2.0);
    t26 = log(rho);
    t27 = rho*t26;
    t30 = (t23-0.15E-9)*rho;
    return(2.0*rho*A_*(4.0*t3+3.0*(2.0*alpha_-2.0)*t2+2.0*(-3.0*alpha_+1.0)*phi)/
delta_+beta_*(t22*((t23-0.15E1)*rho+0.15E1*t27-t30+0.15E1)-t22*(-rho+t27+0.2E1-
t30)));
  }
}


inline double dphireactionSource ( double rho ,double phi ) const
{
  double t18;
  double t19;
  double t2;
  double t22;
  double t23;
  double t26;
  {
    t2 = phi*phi;
    t18 = 120.0*t2*phi-180.0*t2+60.0*phi;
    t19 = log(2.0);
    t22 = log(rho);
    t23 = rho*t22;
    t26 = (t19-0.15E-9)*rho;
    return(2.0*rho*A_*(12.0*t2+6.0*(2.0*alpha_-2.0)*phi-6.0*alpha_+2.0)/delta_+beta_
*(t18*((t19-0.15E1)*rho+0.15E1*t23-t26+0.15E1)-t18*(-rho+t23+0.2E1-t26)));
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
  double t21;
  double t5;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t5 = t1*phi;
    t16 = 6.0*t2*phi;
    t17 = 15.0*t2;
    t18 = 10.0*t5;
    t20 = log(2.0);
    t21 = log(rho);
    return(2.0*A_*(t2+(2.0*alpha_-2.0)*t5+(-3.0*alpha_+1.0)*t1+alpha_)/delta_+beta_*(
(t16-t17+t18)*(t20+0.15E1*t21)+(1.0-t16+t17-t18)*t21));
  }
}

inline double drhochemicalPotential ( double rho ,double phi ) const
{
  double t1;
  double t2;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    return(0.5*beta_*(6.0*t2*phi-15.0*t2+10.0*t1*phi+2.0)/rho);
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
    t29 = 400000000.0*A_*t1+600000000.0*t4*alpha_-600000000.0*t4-600000000.0*A_*
alpha_+200000000.0*A_+1039720771.0*t12*delta_+750000000.0*t12*t16-2079441542.0*t19
*delta_-1500000000.0*t19*t16+1039720771.0*t24*delta_+750000000.0*t24*t16;
    return(0.2E-7*phi*t29/delta_);
  }
}

inline double pressure ( double rho ,double phi ) const
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
    return(0.5*beta_*(6.0*t3*rho-15.0*t2*rho+10.0*t8*rho+2.0*rho-1.0+6.0*t3-15.0
*t2+10.0*t8));
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
    t10 = sqrt(beta_*(6.0*t2*phi-15.0*t2+10.0*t1*phi+2.0));
    return(0.7071067812*t10);
  }
}


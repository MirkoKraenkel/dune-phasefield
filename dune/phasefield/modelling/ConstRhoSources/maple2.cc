
double evalRho(double x)
{
  double t1;
  double t2;
  double t3;
  double t4;
  double t7;
  {
    t1 = log(2.0);
    t2 = x*x;
    t3 = t2*t2;
    t4 = t3*x;
    t7 = t2*x;
    return(exp((-0.15E-9+t1-(6.0*t4-15.0*t3+10.0*t7)*t1)/(1.0+0.3E1*t4-0.75E1*
t3+0.5E1*t7)));
  }
}


inline double helmholtz ( double rho ,double phi ) const
{
  double t1;
  double t16;
  double t17;
  double t18;
  double t2;
  double t20;
  double t23;
  double t24;
  double t27;
  double t5;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t5 = t1*phi;
    t16 = 6.0*t2*phi;
    t17 = 15.0*t2;
    t18 = 10.0*t5;
    t20 = log(2.0);
    t23 = log(rho);
    t24 = rho*t23;
    t27 = (t20-0.15E-9)*rho;
    return(2.0*A_*(t2+(2.0*alpha_-2.0)*t5+(-3.0*alpha_+1.0)*t1+alpha_)/delta_+beta_*(
(t16-t17+t18)*((t20-0.15E1)*rho+0.15E1*t24-t27+0.15E1)+(1.0-t16+t17-t18)*(-rho+
t24+0.2E1-t27)));
  }
}


inline double reactionSource ( double rho ,double phi,double phiMid, double phiOld) const
{
  double t1;
  double t14;
  double t17;
  double t2;
  double t21;
  double t22;
  double t25;
  double t26;
  double t29;
  double t30;
  double t33;
  double t45;
  {
    t1 = phiMid*phiMid;
    t2 = t1*phiMid;
    t14 = 1/delta_;
    t17 = t1*t1;
    t21 = 30.0*t17-60.0*t2+30.0*t1;
    t22 = log(2.0);
    t25 = log(rho);
    t26 = rho*t25;
    t29 = (t22-0.15E-9)*rho;
    t30 = (t22-0.15E1)*rho+0.15E1*t26-t29+0.15E1;
    t33 = -rho+t26+0.2E1-t29;
    t45 = 360.0*t1-360.0*phiMid+60.0;
    return(2.0*A_*(4.0*t2+3.0*(2.0*alpha_-2.0)*t1+2.0*(-3.0*alpha_+1.0)*phiMid)*
t14+beta_*(t21*t30-t21*t33)+(2.0*A_*(24.0*phiMid+12.0*alpha_-12.0)*t14+beta_*(t45*
t30-t45*t33))*(phi-phiOld)/24.0);
  }
}

inline double dphireactionSource ( double rho ,double phi ) const
{
  {
    return(0);
  }
}


inline double chemicalPotential ( double rho , double rhoMid, double rhoOld, double phi ) const
{
  double t1;
  double t13;
  double t14;
  double t19;
  double t2;
  double t20;
  double t4;
  double t5;
  double t7;
  double t8;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t7 = 10.0*t1*phi;
    t8 = t4-t5+t7;
    t9 = log(rhoMid);
    t13 = 1.0-t4+t5-t7;
    t14 = log(2.0);
    t19 = rhoMid*rhoMid;
    t20 = 1/t19;
    return(beta_*(t8*(0.15E-9+0.15E1*t9)+t13*(0.15E-9+t9-t14))+beta_*(-0.15E1*t8*
t20-t13*t20)*(rho-rhoOld)/24.0);
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
  double t4;
  {
    t1 = phi*phi;
    t4 = log(rho);
    return(0.15E-7*beta_*t1*(1386294361.0*t1+1000000000.0*t1*t4-2772588722.0*phi
-2000000000.0*phi*t4+1386294361.0+1000000000.0*t4));
  }
}


inline double pressure ( double rho ,double phi ) const
{
  double t1;
  double t12;
  double t2;
  double t4;
  double t5;
  double t6;
  double t7;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t6 = t1*phi;
    t7 = 10.0*t6;
    t9 = log(2.0);
    t12 = log(rho);
    return(beta_*((t4-t5+t7)*(-(t9-0.15E1)*rho-0.15E1*rho*t12+rho*(t9+0.15E1*t12
))+(1.0-t4+t5-t7)*(rho-0.5))-2.0*A_*(t2+(2.0*alpha_-2.0)*t6+(-3.0*alpha_+1.0)*t1+
alpha_)/delta_);
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


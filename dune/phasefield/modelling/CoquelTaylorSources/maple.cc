
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


inline double reactionSource ( double rho ,double phi ,double old ) const
{
  double t17;
  double t20;
  double t24;
  double t25;
  double t28;
  double t29;
  double t3;
  double t32;
  double t33;
  double t36;
  double t4;
  double t5;
  double t50;
  {
    t3 = 0.5*phi+0.5*old;
    t4 = t3*t3;
    t5 = t4*t3;
    t17 = 1/delta_;
    t20 = t4*t4;
    t24 = 30.0*t20-60.0*t5+30.0*t4;
    t25 = log(2.0);
    t28 = log(rho);
    t29 = rho*t28;
    t32 = (t25-0.15E-9)*rho;
    t33 = (t25-0.15E1)*rho+0.15E1*t29-t32+0.15E1;
    t36 = -rho+t29+0.2E1-t32;
    t50 = 360.0*t4-0.18E3*phi-0.18E3*old+60.0;
    return(2.0*A_*(4.0*t5+3.0*(2.0*alpha_-2.0)*t4+2.0*(-3.0*alpha_+1.0)*t3)*t17+
beta_*(t24*t33-t24*t36)+(2.0*A_*(0.12E2*phi+0.12E2*old+12.0*alpha_-12.0)*t17+beta_*
(t50*t33-t50*t36))*(phi-old)/24.0);
  }
}


inline double dphireactionSource ( double rho ,double phi ,double old ) const
{
  double t13;
  double t21;
  double t22;
  double t25;
  double t26;
  double t29;
  double t3;
  double t30;
  double t33;
  double t4;
  double t41;
  double t61;
  {
    t3 = 0.5*phi+0.5*old;
    t4 = t3*t3;
    t13 = 1/delta_;
    t21 = 0.6E2*t4*t3-0.9E2*t4+0.15E2*phi+0.15E2*old;
    t22 = log(2.0);
    t25 = log(rho);
    t26 = rho*t25;
    t29 = (t22-0.15E-9)*rho;
    t30 = (t22-0.15E1)*rho+0.15E1*t26-t29+0.15E1;
    t33 = -rho+t26+0.2E1-t29;
    t41 = 0.18E3*phi+0.18E3*old-0.18E3;
    t61 = 360.0*t4-0.18E3*phi-0.18E3*old+60.0;
    return(2.0*A_*(0.6E1*t4+0.3E1*(2.0*alpha_-2.0)*t3-0.3E1*alpha_+0.1E1)*t13+beta_
*(t21*t30-t21*t33)+(0.24E2*A_*t13+beta_*(t41*t30-t41*t33))*(phi-old)/24.0+A_*(
0.12E2*phi+0.12E2*old+12.0*alpha_-12.0)*t13/12.0+beta_*(t61*t30-t61*t33)/24.0);
  }
}


inline double chemicalPotential ( double rho ,double phi ,double old ) const
{
  double t1;
  double t11;
  double t12;
  double t16;
  double t17;
  double t2;
  double t22;
  double t23;
  double t4;
  double t5;
  double t7;
  double t8;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t7 = 10.0*t1*phi;
    t8 = t4-t5+t7;
    t11 = 0.5*rho+0.5*old;
    t12 = log(t11);
    t16 = 1.0-t4+t5-t7;
    t17 = log(2.0);
    t22 = t11*t11;
    t23 = 1/t22;
    return(beta_*(t8*(0.15E-9+0.15E1*t12)+t16*(0.15E-9+t12-t17))+beta_*(-0.15E1*
t8*t23-t16*t23)*(rho-old)/24.0);
  }
}

inline double drhochemicalPotential ( double rho ,double phi ,double old ) const
{
  double t1;
  double t11;
  double t12;
  double t15;
  double t2;
  double t20;
  double t22;
  double t32;
  double t4;
  double t5;
  double t7;
  double t8;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t7 = 10.0*t1*phi;
    t8 = t4-t5+t7;
    t11 = 0.5*rho+0.5*old;
    t12 = 1/t11;
    t15 = 1.0-t4+t5-t7;
    t20 = t11*t11;
    t22 = 1/t20/t11;
    t32 = 1/t20;
    return(beta_*(0.75*t8*t12+0.5*t15*t12)+beta_*(0.15E1*t8*t22+0.1E1*t15*t22)*(
rho-old)/24.0+beta_*(-0.15E1*t8*t32-t15*t32)/24.0);
  }
}


inline double dphichemicalPotential ( double rho ,double phi ,double old ) const
{
  double t1;
  double t10;
  double t11;
  double t15;
  double t16;
  double t2;
  double t21;
  double t22;
  double t7;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t7 = 30.0*t2-60.0*t1*phi+30.0*t1;
    t10 = 0.5*rho+0.5*old;
    t11 = log(t10);
    t15 = -t7;
    t16 = log(2.0);
    t21 = t10*t10;
    t22 = 1/t21;
    return(beta_*(t7*(0.15E-9+0.15E1*t11)+t15*(0.15E-9+t11-t16))+beta_*(-0.15E1*
t7*t22-t15*t22)*(rho-old)/24.0);
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



inline double helmholtz ( double rho ,double phi ) const
{
  double t1;
  double t11;
  double t12;
  double t13;
  double t15;
  double t16;
  double t2;
  double t3;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = t1*phi;
    t11 = 6.0*t2*phi;
    t12 = 15.0*t2;
    t13 = 10.0*t3;
    t15 = log(rho);
    t16 = rho*t15;
    return(0.2*A_*(t2-2.0*t3+t1)/delta_+0.1*(t11-t12+t13)*(-rho+t16+0.5)+0.1*(1.0
-t11+t12-t13)*(-0.8068528194400547*rho+0.15E1*t16));
  }
}


inline double reactionSource ( double rho ,double phi ,double old ) const
{
  double t12;
  double t15;
  double t19;
  double t20;
  double t21;
  double t22;
  double t28;
  double t3;
  double t4;
  double t40;
  double t48;
  double t5;
  {
    t3 = 0.5*phi+0.5*old;
    t4 = t3*t3;
    t5 = t4*t3;
    t12 = 1/delta_;
    t15 = t4*t4;
    t19 = 30.0*t15-60.0*t5+30.0*t4;
    t20 = log(rho);
    t21 = rho*t20;
    t22 = -rho+t21+0.5;
    t28 = -0.8068528194400547*rho+0.15E1*t21;
    t40 = 360.0*t4-0.18E3*phi-0.18E3*old+60.0;
    t48 = pow(phi-old,2.0);
    return(0.2*A_*(4.0*t5-6.0*t4+0.1E1*phi+0.1E1*old)*t12+0.1*t19*t22-0.1*t19*
t28+(0.2*A_*(0.12E2*phi+0.12E2*old-12.0)*t12+0.1*t40*t22-0.1*t40*t28)*t48/24.0);
  }
}


inline double dphireactionSource ( double rho ,double phi ,double old ) const
{
  double t10;
  double t18;
  double t19;
  double t20;
  double t21;
  double t27;
  double t3;
  double t34;
  double t4;
  double t41;
  double t42;
  double t54;
  {
    t3 = 0.5*phi+0.5*old;
    t4 = t3*t3;
    t10 = 1/delta_;
    t18 = 0.6E2*t4*t3-0.9E2*t4+0.15E2*phi+0.15E2*old;
    t19 = log(rho);
    t20 = rho*t19;
    t21 = -rho+t20+0.5;
    t27 = -0.8068528194400547*rho+0.15E1*t20;
    t34 = 0.18E3*phi+0.18E3*old-0.18E3;
    t41 = phi-old;
    t42 = t41*t41;
    t54 = 360.0*t4-0.18E3*phi-0.18E3*old+60.0;
    return(0.2*A_*(0.6E1*t4-0.3E1*phi-0.3E1*old+0.1E1)*t10+0.1*t18*t21-0.1*t18*
t27+(0.24E1*A_*t10+0.1*t34*t21-0.1*t34*t27)*t42/24.0+(0.2*A_*(0.12E2*phi+0.12E2*
old-12.0)*t10+0.1*t54*t21-0.1*t54*t27)*t41/12.0);
  }
}


inline double drhoreactionSource ( double rho ,double phi ,double old ) const
{
  double t10;
  double t11;
  double t16;
  double t22;
  double t3;
  double t30;
  double t4;
  double t5;
  {
    t3 = 0.5*phi+0.5*old;
    t4 = t3*t3;
    t5 = t4*t4;
    t10 = 30.0*t5-60.0*t4*t3+30.0*t4;
    t11 = log(rho);
    t16 = 0.6931471805599453+0.15E1*t11;
    t22 = 360.0*t4-0.18E3*phi-0.18E3*old+60.0;
    t30 = pow(phi-old,2.0);
    return(0.1*t10*t11-0.1*t10*t16+(0.1*t22*t11-0.1*t22*t16)*t30/24.0);
  }
}


inline double chemicalPotential ( double rho ,double phi ,double old ) const
{
  double t1;
  double t11;
  double t12;
  double t15;
  double t2;
  double t20;
  double t21;
  double t28;
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
    t15 = 1.0-t4+t5-t7;
    t20 = t11*t11;
    t21 = 1/t20;
    t28 = pow(rho-old,2.0);
    return(0.1*t8*t12+0.1*t15*(0.6931471805599453+0.15E1*t12)+(-0.1*t8*t21-0.15
*t15*t21)*t28/24.0);
  }
}

inline double drhochemicalPotential ( double rho ,double phi ,double old ) const
{
  double t1;
  double t11;
  double t12;
  double t15;
  double t18;
  double t2;
  double t20;
  double t26;
  double t27;
  double t30;
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
    t18 = t11*t11;
    t20 = 1/t18/t11;
    t26 = rho-old;
    t27 = t26*t26;
    t30 = 1/t18;
    return(0.5E-1*t8*t12+0.75E-1*t15*t12+(0.1*t8*t20+0.15*t15*t20)*t27/24.0+(
-0.1*t8*t30-0.15*t15*t30)*t26/12.0);
  }
}


inline double dphichemicalPotential ( double rho ,double phi ,double old ) const
{
  double t1;
  double t10;
  double t11;
  double t14;
  double t19;
  double t2;
  double t20;
  double t27;
  double t7;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t7 = 30.0*t2-60.0*t1*phi+30.0*t1;
    t10 = 0.5*rho+0.5*old;
    t11 = log(t10);
    t14 = -t7;
    t19 = t10*t10;
    t20 = 1/t19;
    t27 = pow(rho-old,2.0);
    return(0.1*t7*t11+0.1*t14*(0.6931471805599453+0.15E1*t11)+(-0.1*t7*t20-0.15
*t14*t20)*t27/24.0);
  }
}


inline double pressure ( double rho ,double phi ) const
{
  double t1;
  double t13;
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
    t13 = log(rho);
    return((t4-t5+t7)*(rho-0.5)+(1.0-t4+t5-t7)*(0.8068528194400547*rho-0.15E1*
rho*t13+rho*(0.6931471805599453+0.15E1*t13))-2.0*A_*(t2-2.0*t6+t1)/delta_);
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
    t9 = sqrt(-12.0*t2*phi+30.0*t2-20.0*t1*phi+6.0);
    return(0.5*t9);
  }
}


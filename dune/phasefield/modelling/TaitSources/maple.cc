double evalRho(double x)
{
  {
    return(0.4958034541*x+0.1065766549);
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
  double t22;
  double t29;
  double t5;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t5 = t1*phi;
    t16 = 6.0*t2*phi;
    t17 = 15.0*t2;
    t18 = 10.0*t5;
    t20 = rho*rho;
    t22 = t20*t20;
    t29 = log(rho);
    return(2.0*A_*(t2+(2.0*alpha_-2.0)*t5+(-3.0*alpha_+1.0)*t1+alpha_)/delta_+(t16-
t17+t18)*(0.35*t22*t20*rho-0.1170549602*rho+0.6043849679E-1)+(1.0-t16+t17-t18)*
(0.1753186565*rho*t29+0.2172006685*rho+0.186848759E-1));
  }
}


inline double reactionSource ( double rho ,double phi ) const
{
  double t1;
  double t17;
  double t2;
  double t21;
  double t22;
  double t24;
  double t31;
  {
    t1 = phi*phi;
    t2 = t1*phi;
    t17 = t1*t1;
    t21 = 30.0*t17-60.0*t2+30.0*t1;
    t22 = rho*rho;
    t24 = t22*t22;
    t31 = log(rho);
    return(2.0*A_*(4.0*t2+3.0*(2.0*alpha_-2.0)*t1+2.0*(-3.0*alpha_+1.0)*phi)/delta_
+t21*(0.35*t24*t22*rho-0.1170549602*rho+0.6043849679E-1)-t21*(0.1753186565*rho*
t31+0.2172006685*rho+0.186848759E-1));
  }
}


inline double dphireactionSource ( double rho ,double phi ) const
{
  double t1;
  double t17;
  double t18;
  double t20;
  double t27;
  {
    t1 = phi*phi;
    t17 = 120.0*t1*phi-180.0*t1+60.0*phi;
    t18 = rho*rho;
    t20 = t18*t18;
    t27 = log(rho);
    return(2.0*A_*(12.0*t1+6.0*(2.0*alpha_-2.0)*phi-6.0*alpha_+2.0)/delta_+t17*(
0.35*t20*t18*rho-0.1170549602*rho+0.6043849679E-1)-t17*(0.1753186565*rho*t27+
0.2172006685*rho+0.186848759E-1));
  }
}


inline double chemicalPotential ( double rho ,double phi ) const
{
  double t1;
  double t13;
  double t17;
  double t2;
  double t3;
  double t4;
  double t5;
  double t6;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = t2*phi;
    t4 = rho*rho;
    t5 = t4*t4;
    t6 = t5*t4;
    t13 = t1*phi;
    t17 = log(rho);
    return(0.147E2*t3*t6-0.3057445711E1*t3-0.3675E2*t2*t6+0.7643614278E1*t2+
0.245E2*t13*t6-0.5095742852E1*t13+0.1753186565*t17+0.392519325-0.1051911939E1*
t3*t17+0.2629779848E1*t2*t17-0.1753186565E1*t13*t17);
  }
}

inline double drhochemicalPotential ( double rho ,double phi ) const
{
  double t1;
  double t11;
  double t2;
  double t3;
  double t4;
  double t5;
  double t6;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = t2*phi;
    t4 = rho*rho;
    t5 = t4*t4;
    t6 = t5*t4;
    t11 = t1*phi;
    return(0.5E-9*(350637313.0+0.1764E12*t3*t6-0.441E12*t2*t6+0.294E12*t11*t6
-2103823878.0*t3+5259559695.0*t2-3506373130.0*t11)/rho);
  }
}


inline double dphichemicalPotential ( double rho ,double phi ) const
{
  double t1;
  double t16;
  double t2;
  double t3;
  double t4;
  double t5;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = rho*rho;
    t4 = t3*t3;
    t5 = t4*t3;
    t9 = t1*phi;
    t16 = log(rho);
    return(0.735E2*t2*t5-0.1528722856E2*t2-147.0*t9*t5+0.3057445711E2*t9+
0.735E2*t1*t5-0.1528722856E2*t1-0.5259559695E1*t2*t16+0.1051911939E2*t9*t16
-0.5259559695E1*t1*t16);
  }
}


inline double pressure ( double rho ,double phi ) const
{
  double t1;
  double t11;
  double t2;
  double t22;
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
    t9 = rho*rho;
    t11 = t9*t9;
    t22 = log(rho);
    return((t4-t5+t7)*(-0.35*t11*t9*rho+0.1170549602*rho-0.6043849679E-1+rho*(
0.245E1*t11*t9-0.1170549602))+(1.0-t4+t5-t7)*(-0.1753186565*rho*t22
-0.2172006685*rho-0.186848759E-1+rho*(0.1753186565*t22+0.392519325))-2.0*A_*(t2+
(2.0*alpha_-2.0)*t6+(-3.0*alpha_+1.0)*t1+alpha_)/delta_);
  }
}


inline double a ( double rho ,double phi ) const
{
  double t1;
  double t11;
  double t18;
  double t2;
  double t3;
  double t4;
  double t5;
  double t6;
  {
    t1 = rho*rho;
    t2 = t1*t1;
    t3 = t2*t1;
    t4 = phi*phi;
    t5 = t4*t4;
    t6 = t5*phi;
    t11 = t4*phi;
    t18 = sqrt(0.882E12*t3*t6-0.2205E13*t5*t3+0.147E13*t11*t3+1753186565.0
-0.1051911939E11*t6+0.2629779848E11*t5-0.1753186565E11*t11);
    return(0.1E-4*t18);
  }
}


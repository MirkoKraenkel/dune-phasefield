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
  double t21;
  double t5;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t5 = t1*phi;
    t16 = 6.0*t2*phi;
    t17 = 15.0*t2;
    t18 = 10.0*t5;
    t20 = log(rho);
    t21 = rho*t20;
    return(2.0*A_*(t2+(2.0*alpha_-2.0)*t5+(-3.0*alpha_+1.0)*t1+alpha_)/delta_+(t16-
t17+t18)*(0.7609303949E-1*t21-0.3752401776E-1*rho+0.4583693344E-1)+(1.0-t16+t17
-t18)*(0.3001504024E-1*t21+0.3718535686E-1*rho+0.3198902536E-2));
  }
}


inline double reactionSource ( double rho ,double phi ) const
{
  double t1;
  double t17;
  double t2;
  double t21;
  double t22;
  double t23;
  {
    t1 = phi*phi;
    t2 = t1*phi;
    t17 = t1*t1;
    t21 = 30.0*t17-60.0*t2+30.0*t1;
    t22 = log(rho);
    t23 = rho*t22;
    return(2.0*A_*(4.0*t2+3.0*(2.0*alpha_-2.0)*t1+2.0*(-3.0*alpha_+1.0)*phi)/delta_
+t21*(0.7609303949E-1*t23-0.3752401776E-1*rho+0.4583693344E-1)-t21*(
0.3001504024E-1*t23+0.3718535686E-1*rho+0.3198902536E-2));
  }
}


inline double dphireactionSource ( double rho ,double phi ) const
{
  double t1;
  double t17;
  double t18;
  double t19;
  {
    t1 = phi*phi;
    t17 = 120.0*t1*phi-180.0*t1+60.0*phi;
    t18 = log(rho);
    t19 = rho*t18;
    return(2.0*A_*(12.0*t1+6.0*(2.0*alpha_-2.0)*phi-6.0*alpha_+2.0)/delta_+t17*(
0.7609303949E-1*t19-0.3752401776E-1*rho+0.4583693344E-1)-t17*(0.3001504024E-1*
t19+0.3718535686E-1*rho+0.3198902536E-2));
  }
}


inline double chemicalPotential ( double rho ,double phi ) const
{
  double t1;
  double t11;
  double t2;
  double t3;
  double t4;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = t2*phi;
    t4 = log(rho);
    t11 = t1*phi;
    return(0.672003971E-1+0.2764679955*t3*t4-0.1717882522*t3-0.6911699888*t2*t4
+0.4294706306*t2+0.4607799925*t11*t4-0.2863137537*t11+0.3001504024E-1*t4);
  }
}

inline double drhochemicalPotential ( double rho ,double phi ) const
{
  double t1;
  double t2;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    return(0.1E-10*(0.2764679955E11*t2*phi-0.6911699888E11*t2+0.4607799925E11*
t1*phi+3001504024.0)/rho);
  }
}


inline double dphichemicalPotential ( double rho ,double phi ) const
{
  double t1;
  double t2;
  double t3;
  double t7;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = log(rho);
    t7 = t1*phi;
    return(0.1382339978E1*t2*t3-0.8589412611*t2-0.2764679955E1*t7*t3+
0.1717882522E1*t7+0.1382339978E1*t1*t3-0.8589412611*t1);
  }
}


inline double pressure ( double rho ,double phi ) const
{
  double t1;
  double t10;
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
    t9 = log(rho);
    t10 = rho*t9;
    return((t4-t5+t7)*(-0.7609303949E-1*t10+0.3752401776E-1*rho-0.4583693344E-1
+rho*(0.7609303949E-1*t9+0.3856902173E-1))+(1.0-t4+t5-t7)*(-0.3001504024E-1*t10
-0.3718535686E-1*rho-0.3198902536E-2+rho*(0.3001504024E-1*t9+0.672003971E-1))
-2.0*A_*(t2+(2.0*alpha_-2.0)*t6+(-3.0*alpha_+1.0)*t1+alpha_)/delta_);
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
    t9 = sqrt(0.6911699885E11*t2*phi-0.1727924971E12*t2+0.1151949981E12*t1*phi+
7503760060.0);
    return(0.2E-5*t9);
  }
}


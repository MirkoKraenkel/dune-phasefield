int evalRho(double x)
{
  {
    return(0);
  }
}


inline double helmholtz ( double rho ,double phi ) const
{
  double t1;
  double t10;
  double t11;
  double t13;
  double t15;
  double t16;
  double t3;
  double t8;
  {
    t1 = phi*phi;
    t3 = pow(1.0-phi,2.0);
    t8 = t1*t1;
    t10 = 6.0*t8*phi;
    t11 = 15.0*t8;
    t13 = 10.0*t1*phi;
    t15 = log(rho);
    t16 = rho*t15;
    return(2.0*t1*t3/delta_+(t10-t11+t13)*(0.702782097E1*t16-0.1067506963E2*rho+
11.0)+(1.0-t10+t11-t13)*(0.1519697146E1*t16-0.318985105*rho));
  }
}


inline double reactionSource ( double rho ,double phi ) const
{
  double t1;
  double t11;
  double t16;
  double t17;
  double t18;
  double t2;
  double t4;
  double t7;
  {
    t1 = 1.0-phi;
    t2 = t1*t1;
    t4 = 1/delta_;
    t7 = phi*phi;
    t11 = t7*t7;
    t16 = 30.0*t11-60.0*t7*phi+30.0*t7;
    t17 = log(rho);
    t18 = rho*t17;
    return(4.0*phi*t2*t4-4.0*t7*t1*t4+t16*(11.0+0.702782097E1*t18
-0.1067506963E2*rho)-t16*(0.1519697146E1*t18-0.318985105*rho));
  }
}


inline double dphireactionSource ( double rho ,double phi ) const
{
  double t1;
  double t16;
  double t17;
  double t18;
  double t2;
  double t3;
  double t9;
  {
    t1 = 1.0-phi;
    t2 = t1*t1;
    t3 = 1/delta_;
    t9 = phi*phi;
    t16 = 120.0*t9*phi-180.0*t9+60.0*phi;
    t17 = log(rho);
    t18 = rho*t17;
    return(4.0*t2*t3-16.0*phi*t1*t3+4.0*t9*t3+t16*(11.0+0.702782097E1*t18
-0.1067506963E2*rho)-t16*(0.1519697146E1*t18-0.318985105*rho));
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
    return(0.3304874294E2*t3*t4-0.2908776421E2*t3-0.8262185736E2*t2*t4+
0.7271941052E2*t2+0.5508123824E2*t11*t4-0.4847960701E2*t11+0.1519697146E1*t4+
0.1200712041E1);
  }
}

inline double drhochemicalPotential ( double rho ,double phi ) const
{
  double t1;
  double t2;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    return(0.2E-8*(0.1652437147E11*t2*phi-0.4131092868E11*t2+0.2754061912E11*t1
*phi+759848573.0)/rho);
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
    return(0.1652437147E3*t2*t3-0.145438821E3*t2-0.3304874294E3*t7*t3+
0.2908776421E3*t7+0.1652437147E3*t1*t3-0.145438821E3*t1);
  }
}


inline double pressure ( double rho ,double phi ) const
{
  double t1;
  double t16;
  double t2;
  double t21;
  double t4;
  double t5;
  double t7;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t7 = 10.0*t1*phi;
    t9 = log(rho);
    t16 = -0.702782097E1*rho*t9+0.1067506963E2*rho-11.0+rho*(0.702782097E1*t9
-0.3647248658E1);
    t21 = pow(1.0-phi,2.0);
    return((t4-t5+t7)*t16+(1.0-t4+t5-t7)*t16-2.0*t1*t21/delta_);
  }
}

inline double a ( double rho ,double phi ) const
{
  {
    return(0.2651003767E1);
  }
}


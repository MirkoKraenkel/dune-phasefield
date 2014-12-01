int evalRho(double x)
{
  {
    return(0);
  }
}


inline double helmholtz ( double rho ,double phi ) const
{
  double t10;
  double t2;
  double t4;
  double t5;
  {
    t2 = log(phi);
    t4 = 1.0-phi;
    t5 = log(t4);
    t10 = log(rho);
    return(theta_*rho*(phi*t2+t4*t5)+phi*rho*t10+2.0*t4*rho*t10);
  }
}


inline double reactionSource ( double rho ,double phi ) const
{
  double t2;
  double t4;
  double t7;
  {
    t2 = log(phi);
    t4 = log(1.0-phi);
    t7 = log(rho);
    return(theta_*rho*(t2-t4)-rho*t7);
  }
}


inline double dphireactionSource ( double rho ,double phi ) const
{
  {
    return(theta_*rho*(1/phi+1/(1.0-phi)));
  }
}


inline double chemicalPotential ( double rho ,double phi ) const
{
  double t1;
  double t4;
  double t6;
  {
    t1 = log(phi);
    t4 = log(1.0-phi);
    t6 = log(rho);
    return(phi*t1+t4-t4*phi-phi*t6-phi+2.0*t6+2.0);
  }
}

inline double drhochemicalPotential ( double rho ,double phi ) const
{
  {
    return(-(phi-2.0)/rho);
  }
}


inline double dphichemicalPotential ( double rho ,double phi ) const
{
  double t1;
  double t3;
  double t4;
  {
    t1 = log(phi);
    t3 = log(1.0-phi);
    t4 = log(rho);
    return(t1-1.0-t3-t4);
  }
}

inline double pressure ( double rho ,double phi ) const
{
  {
    return(-rho*(phi-2.0));
  }
}


inline double a ( double rho ,double phi ) const
{
  {
    return(sqrt(-phi+2.0));
  }
}


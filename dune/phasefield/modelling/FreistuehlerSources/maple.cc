
double evalRho(double x)
{
  {
    return(exp((0.1386294361E1-x*(1.0-f)-f)/(x*(a-2.0)+2.0)));
  }
}


inline double helmholtz ( double rho ,double phi ) const
{
  double t12;
  double t2;
  double t4;
  double t5;
  double t7;
  {
    t2 = log(phi);
    t4 = 1.0-phi;
    t5 = log(t4);
    t7 = phi*phi;
    t12 = log(rho);
    return(theta_*rho*(phi*t2+t4*t5-t7/2.0)+phi*rho*t12+2.0*t4*rho*t12);
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
    return((theta_*rho*(t2-t4-phi)-rho*t7)/delta_);
  }
}


inline double dphireactionSource ( double rho ,double phi ) const
{
  {
    return(theta_*rho*(1/phi+1/(1.0-phi)-1.0)/delta_);
  }
}


inline double chemicalPotential ( double rho ,double phi ) const
{
  double t1;
  double t10;
  double t11;
  double t12;
  double t6;
  {
    t1 = log(phi);
    t6 = log(1.0-1.0*phi);
    t10 = phi*phi;
    t11 = phi*delta_;
    t12 = log(rho);
    return(-0.5*(-2.0*phi*t1-2.0*t6+2.0*t6*phi+t10+2.0*t11*t12+2.0*t11-4.0*
delta_*t12-4.0*delta_)/delta_);
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
  double t2;
  double t4;
  double t6;
  {
    t1 = log(phi);
    t2 = 1.0*phi;
    t4 = log(1.0-t2);
    t6 = log(rho);
    return((t1-1.0*t4-t2-1.0*delta_*t6-1.0*delta_)/delta_);
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


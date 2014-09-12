
inline double helmholtz ( double rho ) const
{
  double t5;
  {
    t5 = log(rho/(3.0-rho));
    return(rho*(-3.0*rho+0.24E1*t5+0.3979297858E1));
  }
}


inline double chemicalPotential ( double rho ) const
{
  double t2;
  double t3;
  double t5;
  double t7;
  {
    t2 = 3.0-rho;
    t3 = 1/t2;
    t5 = log(rho*t3);
    t7 = t2*t2;
    return(-3.0*rho+0.24E1*t5+0.3979297858E1+rho*(-3.0+0.24E1*(t3+rho/t7)/rho*
t2));
  }
}


inline double drhochemicalPotential ( double rho ) const
{
  double t10;
  double t2;
  {
    t2 = rho*rho;
    t10 = pow(-3.0+rho,2.0);
    return(-0.12E1*(45.0*rho-30.0*t2+5.0*t2*rho-18.0)/rho/t10);
  }
}


inline double pressure ( double rho ) const
{
  double t1;
  double t2;
  double t3;
  double t5;
  double t6;
  double t9;
  {
    t1 = 3.0*rho;
    t2 = 3.0-rho;
    t3 = 1/t2;
    t5 = log(rho*t3);
    t6 = 0.24E1*t5;
    t9 = t2*t2;
    return(-rho*(-t1+t6+0.3979297858E1)+rho*(-t1+t6+0.3979297858E1+rho*(-3.0+
0.24E1*(t3+rho/t9)/rho*t2)));
  }
}


inline double a ( double rho ) const
{
  double t2;
  double t8;
  {
    t2 = rho*rho;
    t8 = pow(-3.0+rho,2.0);
    return(-0.12E1*(45.0*rho-30.0*t2+5.0*t2*rho-18.0)/t8);
  }
}


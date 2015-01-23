inline double helmholtz ( double rho ,double phi ) const
{
  double t1;
  double t16;
  double t17;
  double t18;
  double t2;
  double t5;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t5 = t1*phi;
    t16 = 6.0*t2*phi;
    t17 = 15.0*t2;
    t18 = 10.0*t5;
    return(2.0*A_*(t2+(2.0*alpha_-2.0)*t5+(-3.0*alpha_+1.0)*t1+alpha_)/delta_+(t16-
t17+t18)*psi1( rho )+(1.0-t16+t17-t18)*psi2( rho ));
  }
}

inline double reactionSource ( double rho ,double phi ) const
{
  double t1;
  double t17;
  double t2;
  double t21;
  {
    t1 = phi*phi;
    t2 = t1*phi;
    t17 = t1*t1;
    t21 = 30.0*t17-60.0*t2+30.0*t1;
    return(2.0*A_*(4.0*t2+3.0*(2.0*alpha_-2.0)*t1+2.0*(-3.0*alpha_+1.0)*phi)/delta_
+t21*psi1( rho )-t21*psi2( rho ));
  }
}

inline double dphireactionSource ( double rho ,double phi ) const
{
  double t1;
  double t17;
  {
    t1 = phi*phi;
    t17 = 120.0*t1*phi-180.0*t1+60.0*phi;
    return(2.0*A_*(12.0*t1+6.0*(2.0*alpha_-2.0)*phi-6.0*alpha_+2.0)/delta_+t17*psi1( rho )-
t17*psi2( rho ));
  }
}

inline double chemicalPotential ( double rho ,double phi ) const
{
  double t1;
  double t2;
  double t4;
  double t5;
  double t7;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t7 = 10.0*t1*phi;
    return((t4-t5+t7)*mu1( rho )+(1.0-t4+t5-t7)*mu2( rho ));
  }
}

inline double dphichemicalPotential ( double rho ,double phi ) const
{
  double t1;
  double t2;
  double t7;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t7 = 30.0*t2-60.0*t1*phi+30.0*t1;
    return(t7*mu1( rho )-t7*mu2( rho ));
  }
}

inline double drhochemicalPotential ( double rho ,double phi ) const
{
  double t1;
  double t2;
  double t4;
  double t5;
  double t7;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t7 = 10.0*t1*phi;
    return((t4-t5+t7)*dmu1( rho )+(1.0-t4+t5-t7)*dmu2( rho ));
  }
}

inline double pressure ( double rho ,double phi ) const
{
  double t1;
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
    return((t4-t5+t7)*p1( rho )+(1.0-t4+t5-t7)*p2( rho )-2.0*A_*(t2+(2.0*alpha_-2.0)*t6+(
-3.0*alpha_+1.0)*t1+alpha_)/delta_);
  }
}


inline double a ( double rho ,double phi ) const
{
  double t1;
  double t2;
  double t4;
  double t5;
  double t7;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t7 = 10.0*t1*phi;
    return(sqrt((t4-t5+t7)*dp1+(1.0-t4+t5-t7)*dp0));
  }
}


inline double psi1 ( double rho ) const
{
  double t10;
  double t6;
  {
    t6 = log(rho/(3.0-rho));
    t10 = log(rho);
    if( rho <= 0.1280223018 )
      return(rho*(-3.0*rho+0.1866666667E1*t6));
    else 
      return(-0.4586570656+0.1268660419E1*rho*t10);
  }
}


inline double psi2 ( double rho ) const
{
  double t10;
  double t6;
  {
    t6 = log(rho/(3.0-rho));
    t10 = log(rho);
    if( rho <= 0.1280223018 )
      return(rho*(-3.0*rho+0.1866666667E1*t6));
    else 
      return(-0.4586570656+0.1268660419E1*rho*t10);
  }
}


inline double mu1 ( double rho ) const
{
  double t19;
  double t3;
  double t4;
  double t6;
  double t8;
  {
    t3 = 3.0-rho;
    t4 = 1/t3;
    t6 = log(rho*t4);
    t8 = t3*t3;
    t19 = log(rho);
    if( rho <= 0.1280223018 )
      return(-3.0*rho+0.1866666667E1*t6+rho*(-3.0+0.1866666667E1*(t4+rho/t8)/
rho*t3));
    else 
      return(0.1268660419E1*t19+0.1268660419E1);
  }
}


inline double mu2 ( double rho ) const
{
  double t11;
  double t2;
  double t6;
  double t7;
  double t9;
  {
    t2 = log(rho);
    t6 = 3.0-rho;
    t7 = 1/t6;
    t9 = log(rho*t7);
    t11 = t6*t6;
    if( rho < 0.2140442548E1 )
      return(0.9895721682E1*t2-0.1215539003E2);
    else 
      return(-3.0*rho+0.1866666667E1*t9+rho*(-3.0+0.1866666667E1*(t7+rho/t11)/
rho*t6));
  }
}


inline double p1 ( double rho ) const
{
  double t19;
  double t3;
  double t4;
  double t6;
  double t8;
  {
    t3 = 3.0-rho;
    t4 = 1/t3;
    t6 = log(rho*t4);
    t8 = t3*t3;
    t19 = log(rho);
    if( rho <= 0.1280223018 )
      return(-3.0*rho+0.1866666667E1*t6+rho*(-3.0+0.1866666667E1*(t4+rho/t8)/
rho*t3));
    else 
      return(0.1268660419E1*t19+0.1268660419E1);
  }
}


inline double p2 ( double rho ) const
{
  double t19;
  double t3;
  double t4;
  double t6;
  double t8;
  {
    t3 = 3.0-rho;
    t4 = 1/t3;
    t6 = log(rho*t4);
    t8 = t3*t3;
    t19 = log(rho);
    if( rho <= 0.1280223018 )
      return(-3.0*rho+0.1866666667E1*t6+rho*(-3.0+0.1866666667E1*(t4+rho/t8)/
rho*t3));
    else 
      return(0.1268660419E1*t19+0.1268660419E1);
  }
}

inline double dmu1 ( double rho ) const
{
  double t2;
  double t20;
  double t4;
  double t5;
  double t7;
  double t8;
  double t9;
  {
    t2 = 3.0-rho;
    t4 = t2*t2;
    t5 = 1/t4;
    t7 = 1/t2+rho*t5;
    t8 = 1/rho;
    t9 = t7*t8;
    t20 = rho*rho;
    if( rho < 0.1280223018 )
      return(-6.0+0.3733333334E1*t9*t2+rho*(0.1866666667E1*(2.0*t5+2.0*rho/t4/
t2)*t8*t2-0.1866666667E1*t7/t20*t2-0.1866666667E1*t9));
    else 
      return(0.1268660419E1*t8);
  }
}

inline double dmu2 ( double rho ) const
{
  double t10;
  double t2;
  double t21;
  double t4;
  double t6;
  double t7;
  double t9;
  {
    t2 = 1/rho;
    t4 = 3.0-rho;
    t6 = t4*t4;
    t7 = 1/t6;
    t9 = 1/t4+rho*t7;
    t10 = t9*t2;
    t21 = rho*rho;
    if( rho < 0.2140442548E1 )
      return(0.1268660419E1*t2);
    else 
      return(-6.0+0.3733333334E1*t10*t4+rho*(0.1866666667E1*(2.0*t7+2.0*rho/t6/
t4)*t2*t4-0.1866666667E1*t9/t21*t4-0.1866666667E1*t10));
  }
}


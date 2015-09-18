inline double mwpliq ( double x ) const
{
  {
    return(0);
  }
}

inline double mwpvap ( double x ) const
{
  {
    return(2);
  }
}


inline double evalRho ( double x ) const
{
  double t1;
  double t2;
  double t8;
  {
    t1 = x*x;
    t2 = t1*t1;
    t8 = 6.0*t2*x-15.0*t2+10.0*t1*x;
    return(exp((-t8*(bv-bl)-bl)/(t8*(cv-cl)+cl)));
  }
}


inline double mwpliq ( double x ) const
{
  {
    return(0.7367096199);
  }
}

inline double mwpvap ( double x ) const
{
  {
    return(0.4768726654E-1);
  }
}


inline double evalRho ( double x ) const
{
  double t1;
  double t2;
  double t3;
  double t6;
  {
    t1 = x*x;
    t2 = t1*t1;
    t3 = t2*x;
    t6 = t1*x;
    return(exp((-0.2891561136E2+0.1204277045E3*t3-0.3010692614E3*t2+
0.2007128409E3*t6)/(0.9502053219E1+0.116654747E3*t3-0.2916368674E3*t2+
0.1944245783E3*t6)));
  }
}


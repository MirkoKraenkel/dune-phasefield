inline double exactrho (double t, double x, double y ) const
{
  {
    return(0.15E1);
  }
}


inline double exactv1 (double t, double x, double y ) const
{
  double t3;
  double t6;
  {
    t3 = cos(2.0*0.3141592653589793E1*t);
    t6 = sin(2.0*0.3141592653589793E1*x);
    return(t3*t6);
  }
}

inline double exactv2 (double t, double x, double y ) const
{
  {
    return(0);
  }
}


inline double exactphi (double t, double x, double y ) const
{
  double t3;
  double t6;
  {
    t3 = cos(2.0*0.3141592653589793E1*t);
    t6 = cos(2.0*0.3141592653589793E1*x);
    return(0.5*t3*t6+0.5);
  }
}

inline double exactsigma2 (double t, double x, double y ) const
{
  {
    return(0);
  }
}


inline double exactsigma1 (double t, double x, double y ) const
{
  {
    return(-0.1E1*cos(2.0*0.3141592653589793E1*t)*sin(2.0*0.3141592653589793E1*
x)*0.3141592653589793E1);
  }
}


inline double exacttau (double t, double x, double y ) const
{
  {
    return(0.2*A_*(4.0*pow(0.5*cos(2.0*0.3141592653589793E1*t)*cos(2.0*
0.3141592653589793E1*x)+0.5,3.0)-6.0*pow(0.5*cos(2.0*0.3141592653589793E1*t)*
cos(2.0*0.3141592653589793E1*x)+0.5,2.0)+0.1E1*cos(2.0*0.3141592653589793E1*t)*
cos(2.0*0.3141592653589793E1*x)+0.1E1)/delta_-0.2814588057631238*pow(0.5*cos(2.0
*0.3141592653589793E1*t)*cos(2.0*0.3141592653589793E1*x)+0.5,4.0)+
0.5629176115262475*pow(0.5*cos(2.0*0.3141592653589793E1*t)*cos(2.0*
0.3141592653589793E1*x)+0.5,3.0)-0.2814588057631238*pow(0.5*cos(2.0*
0.3141592653589793E1*t)*cos(2.0*0.3141592653589793E1*x)+0.5,2.0)+0.2E1*delta_*A_*
cos(2.0*0.3141592653589793E1*t)*cos(2.0*0.3141592653589793E1*x)*
0.3141592653589793E1*0.3141592653589793E1);
  }
}


inline double exactmu (double t, double x, double y ) const
{
  {
    return(0.1301344842722192-0.5375278407684165*pow(0.5*cos(2.0*
0.3141592653589793E1*t)*cos(2.0*0.3141592653589793E1*x)+0.5,5.0)+
0.1343819601921041E1*pow(0.5*cos(2.0*0.3141592653589793E1*t)*cos(2.0*
0.3141592653589793E1*x)+0.5,4.0)-0.8958797346140275*pow(0.5*cos(2.0*
0.3141592653589793E1*t)*cos(2.0*0.3141592653589793E1*x)+0.5,3.0)+0.5*pow(cos(
2.0*0.3141592653589793E1*t),2.0)*pow(sin(2.0*0.3141592653589793E1*x),2.0));
  }
}


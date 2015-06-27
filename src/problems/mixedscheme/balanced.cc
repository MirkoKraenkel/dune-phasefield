inline double evalRho(double x) const
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
    return(exp((0.6931471804-0.4158883084E1*t3+0.1039720771E2*t2-0.6931471806E1
*t6)/(0.3E1*t3-0.75E1*t2+0.5E1*t6+1.0)));
  }
}


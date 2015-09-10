
inline double rhsRho (double t, double x, double y ) const
{
  double t3;
  double t6;
  {
    t3 = cos(2.0*0.3141592653589793E1*t);
    t6 = cos(2.0*0.3141592653589793E1*x);
    return(0.3E1*t3*t6*0.3141592653589793E1);
  }
}


inline double rhsV1 (double t, double x, double y ) const
{
  double t10;
  double t11;
  double t13;
  double t17;
  double t19;
  double t2;
  double t20;
  double t21;
  double t23;
  double t26;
  double t3;
  double t46;
  double t6;
  double t7;
  {
    t2 = 2.0*0.3141592653589793E1*t;
    t3 = sin(t2);
    t6 = 2.0*0.3141592653589793E1*x;
    t7 = sin(t6);
    t10 = cos(t2);
    t11 = t10*t10;
    t13 = cos(t6);
    t17 = t10*t13;
    t19 = 0.5*t17+0.5;
    t20 = t19*t19;
    t21 = t20*t20;
    t23 = t7*0.3141592653589793E1;
    t26 = t20*t19;
    t46 = 0.3141592653589793E1*0.3141592653589793E1;
    return(-0.3E1*t3*0.3141592653589793E1*t7+0.3E1*t11*t7*t13*
0.3141592653589793E1+0.4031458805763124E1*t21*t10*t23-0.8062917611526248E1*t26*
t10*t23+0.4031458805763124E1*t20*t10*t23+0.1E1*t10*t7*0.3141592653589793E1*(0.2
*A_*(4.0*t26-6.0*t20+0.1E1*t17+0.1E1)/delta_-0.2814588057631238*t21+
0.5629176115262475*t26-0.2814588057631238*t20+0.2E1*delta_*A_*t17*t46)+4.0*mu1Liq_*
t10*t7*t46);
  }
}

inline double rhsV2 (double t, double x, double y ) const
{
  {
    return(0.0);
  }
}


inline double rhsPhi (double t, double x, double y ) const
{
  double t10;
  double t11;
  double t12;
  double t13;
  double t17;
  double t19;
  double t2;
  double t20;
  double t21;
  double t3;
  double t30;
  double t35;
  double t6;
  double t7;
  {
    t2 = 2.0*0.3141592653589793E1*t;
    t3 = sin(t2);
    t6 = 2.0*0.3141592653589793E1*x;
    t7 = cos(t6);
    t10 = cos(t2);
    t11 = t10*t10;
    t12 = sin(t6);
    t13 = t12*t12;
    t17 = t10*t7;
    t19 = 0.5*t17+0.5;
    t20 = t19*t19;
    t21 = t20*t19;
    t30 = t20*t20;
    t35 = 0.3141592653589793E1*0.3141592653589793E1;
    return(-0.1E1*t3*0.3141592653589793E1*t7-0.1E1*t11*t13*0.3141592653589793E1
+0.1333333333333333*A_*(4.0*t21-6.0*t20+0.1E1*t17+0.1E1)/delta_
-0.1876392038420825*t30+0.375278407684165*t21-0.1876392038420825*t20+
0.1333333333333333E1*delta_*A_*t17*t35);
  }
}


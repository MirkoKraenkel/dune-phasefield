
inline double rhsRho (double t, double x, double y ) const
{
  double t10;
  double t11;
  double t12;
  double t13;
  double t2;
  double t3;
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
    return(-0.1E1*t3*0.3141592653589793E1*t7-0.1E1*t11*t13*0.3141592653589793E1
+2.0*(0.5*t10*t7+0.15E1)*t10*t7*0.3141592653589793E1);
  }
}


inline double rhsV1 (double t, double x, double y ) const
{
  double t10;
  double t103;
  double t12;
  double t13;
  double t16;
  double t18;
  double t2;
  double t3;
  double t32;
  double t38;
  double t39;
  double t40;
  double t44;
  double t5;
  double t51;
  double t52;
  double t56;
  double t57;
  double t58;
  double t6;
  double t62;
  double t7;
  double t8;
  double t9;
  double t91;
  double t92;
  {
    t2 = 2.0*0.3141592653589793E1*t;
    t3 = cos(t2);
    t5 = 2.0*0.3141592653589793E1*x;
    t6 = cos(t5);
    t7 = t3*t6;
    t8 = 0.5*t7;
    t9 = t8+0.15E1;
    t10 = sin(t2);
    t12 = sin(t5);
    t13 = 0.3141592653589793E1*t12;
    t16 = t3*t3;
    t18 = t12*t12;
    t32 = t6*0.3141592653589793E1;
    t38 = t8+0.5;
    t39 = t38*t38;
    t40 = t39*t39;
    t44 = t39*t38;
    t51 = -0.3E2*t40*t3*t13+0.6E2*t44*t3*t13-0.3E2*t39*t3*t13;
    t52 = log(t9);
    t56 = 6.0*t40*t38;
    t57 = 15.0*t40;
    t58 = 10.0*t44;
    t62 = t13/t9;
    t91 = 30.0*t40-60.0*t44+30.0*t39;
    t92 = t9*t52;
    t103 = 0.3141592653589793E1*0.3141592653589793E1;
    return(-2.0*t9*t10*t13-0.1E1*t16*t3*t18*t12*0.3141592653589793E1+0.2E1*t9*
t16*t12*t6*0.3141592653589793E1-(-0.1E1*t16*t18*0.3141592653589793E1+2.0*t9*t3*
t32)*t3*t12+t9*(0.1*t51*t52-0.1*(t56-t57+t58)*t3*t62-0.1*t51*(
0.6931471805599453+0.15E1*t52)-0.15*(1.0-t56+t57-t58)*t3*t62+0.2E1*t16*t12*t32)
+0.1E1*t3*t12*0.3141592653589793E1*(0.2*A_*(4.0*t44-6.0*t39+0.1E1*t7+0.1E1)/
delta_+0.1*t91*(-t8-0.1E1+t92)-0.1*t91*(-0.1210279229160082E1-0.4034264097200273
*t7+0.15E1*t92)+0.2E1*delta_*A_*t7*t103)+4.0*mu1Liq_*t3*t12*t103);
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
  double t18;
  double t19;
  double t2;
  double t20;
  double t21;
  double t3;
  double t30;
  double t34;
  double t35;
  double t36;
  double t37;
  double t48;
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
    t18 = 0.5*t17;
    t19 = t18+0.5;
    t20 = t19*t19;
    t21 = t20*t19;
    t30 = t20*t20;
    t34 = 30.0*t30-60.0*t21+30.0*t20;
    t35 = t18+0.15E1;
    t36 = log(t35);
    t37 = t35*t36;
    t48 = 0.3141592653589793E1*0.3141592653589793E1;
    return(-0.1E1*t3*0.3141592653589793E1*t7-0.1E1*t11*t13*0.3141592653589793E1
+(0.2*A_*(4.0*t21-6.0*t20+0.1E1*t17+0.1E1)/delta_+0.1*t34*(-t18-0.1E1+t37)-0.1*
t34*(-0.1210279229160082E1-0.4034264097200273*t17+0.15E1*t37)+0.2E1*delta_*A_*t17
*t48)/t35);
  }
}


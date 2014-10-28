
inline double rhosol ( double x ) const
{
  double t11;
  double t3;
  double t5;
  double t6;
  double t7;
  double t8;
  {
    t3 = tanh(x/delta_);
    t5 = 0.5*t3+0.5;
    t6 = t5*t5;
    t7 = t6*t6;
    t8 = t7*t5;
    t11 = t6*t5;
    return(exp((0.4000000004E1-18.0*t8+45.0*t7-30.0*t11)/(-0.9E1*t8+0.225E2*t7
-0.15E2*t11+3.0)));
  }
}


inline double gradrho ( double x ) const
{
  double t1;
  double t11;
  double t13;
  double t15;
  double t18;
  double t21;
  double t25;
  double t26;
  double t3;
  double t31;
  double t32;
  double t42;
  double t5;
  double t6;
  double t7;
  double t8;
  double t9;
  {
    t1 = 1/delta_;
    t3 = tanh(t1*x);
    t5 = 0.5*t3+0.5;
    t6 = t5*t5;
    t7 = t6*t6;
    t8 = t3*t3;
    t9 = 1.0-t8;
    t11 = t7*t9*t1;
    t13 = t6*t5;
    t15 = t13*t9*t1;
    t18 = t6*t9*t1;
    t21 = t7*t5;
    t25 = 3.0-0.9E1*t21+0.225E2*t7-0.15E2*t13;
    t26 = 1/t25;
    t31 = 0.4000000004E1-18.0*t21+45.0*t7-30.0*t13;
    t32 = t25*t25;
    t42 = exp(t31*t26);
    return(((-0.45E2*t11+0.9E2*t15-0.45E2*t18)*t26-t31/t32*(-0.225E2*t11+0.45E2
*t15-0.225E2*t18))*t42);
  }
}


inline double gradphi ( double x ) const
{
  {
    return(0.5*(1.0-pow(tanh(x/delta_),2.0))/delta_);
  }
}


inline double thetasol ( double x ) const
{
  double t1;
  double t21;
  double t25;
  double t26;
  double t3;
  double t37;
  double t39;
  double t40;
  double t5;
  double t51;
  double t6;
  double t7;
  {
    t1 = 1/delta_;
    t3 = tanh(t1*x);
    t5 = 0.5*t3+0.5;
    t6 = t5*t5;
    t7 = t6*t5;
    t21 = t6*t6;
    t25 = 30.0*t21-60.0*t7+30.0*t6;
    t26 = t21*t5;
    t37 = exp((0.4000000004E1-18.0*t26+45.0*t21-30.0*t7)/(-0.9E1*t26+0.225E2*
t21-0.15E2*t7+3.0));
    t39 = log(t37);
    t40 = t37*t39;
    t51 = t3*t3;
    return(2.0*A_*(4.0*t7+3.0*(2.0*alpha_-2.0)*t6+2.0*(-3.0*alpha_+1.0)*t5)*t1+
beta_*(t25*(-0.25E1*t37+0.15E1*t40+0.2921601062E1)-t25*(-7.0*t37+3.0*t40+
0.1138100368E2))+0.1E1*t1*t3*(1.0-t51));
  }
}


inline double phiSource ( double x ) const
{
  double t1;
  double t11;
  double t20;
  double t22;
  double t24;
  double t25;
  double t26;
  double t3;
  double t40;
  double t44;
  double t45;
  double t5;
  double t56;
  double t58;
  double t59;
  double t6;
  double t7;
  double t70;
  double t8;
  {
    t1 = 1/delta_;
    t3 = tanh(t1*x);
    t5 = 0.5*t3+0.5;
    t6 = t5*t5;
    t7 = t6*t6;
    t8 = t7*t5;
    t11 = t6*t5;
    t20 = exp((0.4000000004E1-18.0*t8+45.0*t7-30.0*t11)/(-0.9E1*t8+0.225E2*t7
-0.15E2*t11+3.0));
    t22 = tanh(t1*t20);
    t24 = 0.5*t22+0.5;
    t25 = t24*t24;
    t26 = t25*t24;
    t40 = t25*t25;
    t44 = 30.0*t40-60.0*t26+30.0*t25;
    t45 = t40*t24;
    t56 = exp((0.4000000004E1-18.0*t45+45.0*t40-30.0*t26)/(-0.9E1*t45+0.225E2*
t40-0.15E2*t26+3.0));
    t58 = log(t56);
    t59 = t58*t56;
    t70 = t22*t22;
    return((2.0*A_*(4.0*t26+3.0*(2.0*alpha_-2.0)*t25+2.0*(-3.0*alpha_+1.0)*t24)*t1
+beta_*(t44*(-0.25E1*t56+0.15E1*t59+0.2921601062E1)-t44*(-7.0*t56+3.0*t59+
0.1138100368E2))+0.1E1*t1*t22*(1.0-t70))/t20);
  }
}


inline double musol ( double x ) const
{
  double t13;
  double t22;
  double t23;
  double t27;
  double t3;
  double t31;
  double t5;
  double t6;
  double t7;
  double t8;
  {
    t3 = tanh(x/delta_);
    t5 = 0.5*t3+0.5;
    t6 = t5*t5;
    t7 = t6*t6;
    t8 = t7*t5;
    t13 = t6*t5;
    t22 = exp((0.4000000004E1-18.0*t8+45.0*t7-30.0*t13)/(-0.9E1*t8+0.225E2*t7
-0.15E2*t13+3.0));
    t23 = log(t22);
    t27 = beta_*t7;
    t31 = beta_*t13;
    return(-6.0*beta_*t8+9.0*t23*beta_*t8+15.0*t27-0.225E2*t27*t23-10.0*t31+15.0*
t31*t23-4.0+3.0*t23+24.0*t8-18.0*t8*t23-60.0*t7+45.0*t7*t23+40.0*t13-30.0*t13*
t23);
  }
}


inline double veloSource ( double x ) const
{
  double t1;
  double t10;
  double t124;
  double t24;
  double t28;
  double t29;
  double t3;
  double t33;
  double t37;
  double t38;
  double t4;
  double t40;
  double t42;
  double t43;
  double t5;
  double t59;
  double t6;
  double t63;
  double t66;
  double t69;
  double t72;
  double t73;
  double t8;
  double t80;
  double t81;
  double t88;
  double t9;
  double t92;
  double t97;
  {
    t1 = 1/delta_;
    t3 = tanh(x*t1);
    t4 = t3*t3;
    t5 = 1.0-t4;
    t6 = t5*t1;
    t8 = 0.5*t3+0.5;
    t9 = t8*t8;
    t10 = t9*t8;
    t24 = t9*t9;
    t28 = 30.0*t24-60.0*t10+30.0*t9;
    t29 = t24*t8;
    t33 = 0.4000000004E1-18.0*t29+45.0*t24-30.0*t10;
    t37 = -0.9E1*t29+0.225E2*t24-0.15E2*t10+3.0;
    t38 = 1/t37;
    t40 = exp(t33*t38);
    t42 = log(t40);
    t43 = t42*t40;
    t59 = beta_*t24;
    t63 = t24*t5*t1;
    t66 = t10*t5*t1;
    t69 = t9*t5*t1;
    t72 = (-0.45E2*t63+0.9E2*t66-0.45E2*t69)*t38;
    t73 = t37*t37;
    t80 = t33/t73*(-0.225E2*t63+0.45E2*t66-0.225E2*t69);
    t81 = t72-t80;
    t88 = beta_*t10;
    t92 = t42*t5*t1;
    t97 = beta_*t9;
    t124 = -0.15E2*t59*t6+9.0*t81*beta_*t29+0.225E2*t42*beta_*t63+0.3E2*t88*t6
-0.45E2*t88*t92-0.225E2*t59*t81-0.15E2*t97*t6+0.225E2*t97*t92+15.0*t88*t81+3.0*
t72-3.0*t80+0.6E2*t63-0.45E2*t24*t42*t6-18.0*t29*t81-0.12E3*t66+0.9E2*t10*t42*
t6+45.0*t24*t81+0.6E2*t69-0.45E2*t9*t42*t6-30.0*t10*t81;
    return(-0.5*t6*(2.0*A_*(4.0*t10+3.0*(2.0*alpha_-2.0)*t9+2.0*(-3.0*alpha_+1.0)*
t8)*t1+beta_*(t28*(-0.25E1*t40+0.15E1*t43+0.2921601062E1)-t28*(-7.0*t40+3.0*t43+
0.1138100368E2))+0.1E1*t1*t3*t5)+t40*t124);
  }
}


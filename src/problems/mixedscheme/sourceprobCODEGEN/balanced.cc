
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
  double t19;
  double t23;
  double t24;
  double t3;
  double t35;
  double t37;
  double t38;
  double t49;
  double t5;
  double t6;
  double t7;
  {
    t1 = 1/delta_;
    t3 = tanh(t1*x);
    t5 = 0.5*t3+0.5;
    t6 = t5*t5;
    t7 = t6*t5;
    t19 = t6*t6;
    t23 = 30.0*t19-60.0*t7+30.0*t6;
    t24 = t19*t5;
    t35 = exp((0.4000000004E1-18.0*t24+45.0*t19-30.0*t7)/(-0.9E1*t24+0.225E2*
t19-0.15E2*t7+3.0));
    t37 = log(t35);
    t38 = t35*t37;
    t49 = t3*t3;
    return((8.0*t7+6.0*(2.0*alpha_-2.0)*t6+4.0*(-3.0*alpha_+1.0)*t5)*t1+beta_*(t23
*(-0.25E1*t35+0.15E1*t38+0.15E1)-t23*(-7.0*t35+3.0*t38+0.995940263E1))+0.1E1*t1
*t3*(1.0-t49));
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
  double t38;
  double t42;
  double t43;
  double t5;
  double t54;
  double t56;
  double t57;
  double t6;
  double t68;
  double t7;
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
    t38 = t25*t25;
    t42 = 30.0*t38-60.0*t26+30.0*t25;
    t43 = t38*t24;
    t54 = exp((0.4000000004E1-18.0*t43+45.0*t38-30.0*t26)/(-0.9E1*t43+0.225E2*
t38-0.15E2*t26+3.0));
    t56 = log(t54);
    t57 = t56*t54;
    t68 = t22*t22;
    return(((8.0*t26+6.0*(2.0*alpha_-2.0)*t25+4.0*(-3.0*alpha_+1.0)*t24)*t1+beta_*
(t42*(-0.25E1*t54+0.15E1*t57+0.15E1)-t42*(-7.0*t54+3.0*t57+0.995940263E1))+
0.1E1*t1*t22*(1.0-t68))/t20);
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
  double t122;
  double t22;
  double t26;
  double t27;
  double t3;
  double t31;
  double t35;
  double t36;
  double t38;
  double t4;
  double t40;
  double t41;
  double t5;
  double t57;
  double t6;
  double t61;
  double t64;
  double t67;
  double t70;
  double t71;
  double t78;
  double t79;
  double t8;
  double t86;
  double t9;
  double t90;
  double t95;
  {
    t1 = 1/delta_;
    t3 = tanh(x*t1);
    t4 = t3*t3;
    t5 = 1.0-t4;
    t6 = t5*t1;
    t8 = 0.5*t3+0.5;
    t9 = t8*t8;
    t10 = t9*t8;
    t22 = t9*t9;
    t26 = 30.0*t22-60.0*t10+30.0*t9;
    t27 = t22*t8;
    t31 = 0.4000000004E1-18.0*t27+45.0*t22-30.0*t10;
    t35 = -0.9E1*t27+0.225E2*t22-0.15E2*t10+3.0;
    t36 = 1/t35;
    t38 = exp(t31*t36);
    t40 = log(t38);
    t41 = t40*t38;
    t57 = beta_*t22;
    t61 = t22*t5*t1;
    t64 = t10*t5*t1;
    t67 = t9*t5*t1;
    t70 = (-0.45E2*t61+0.9E2*t64-0.45E2*t67)*t36;
    t71 = t35*t35;
    t78 = t31/t71*(-0.225E2*t61+0.45E2*t64-0.225E2*t67);
    t79 = t70-t78;
    t86 = beta_*t10;
    t90 = t40*t5*t1;
    t95 = beta_*t9;
    t122 = -0.15E2*t57*t6+9.0*t79*beta_*t27+0.225E2*t40*beta_*t61+0.3E2*t86*t6
-0.45E2*t86*t90-0.225E2*t57*t79-0.15E2*t95*t6+0.225E2*t95*t90+15.0*t86*t79+3.0*
t70-3.0*t78+0.6E2*t61-0.45E2*t22*t40*t6-18.0*t27*t79-0.12E3*t64+0.9E2*t10*t40*
t6+45.0*t22*t79+0.6E2*t67-0.45E2*t9*t40*t6-30.0*t10*t79;
    return(-0.5*t6*((8.0*t10+6.0*(2.0*alpha_-2.0)*t9+4.0*(-3.0*alpha_+1.0)*t8)*t1
+beta_*(t26*(-0.25E1*t38+0.15E1*t41+0.15E1)-t26*(-7.0*t38+3.0*t41+0.995940263E1)
)+0.1E1*t1*t3*t5)+t38*t122);
  }
}


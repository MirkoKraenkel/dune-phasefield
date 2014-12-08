
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
    return((0.113810037E2-0.5075641578E2*t8+0.1268910394E3*t7-0.845940263E2*t11
)/(-0.9E1*t8+0.225E2*t7-0.15E2*t11+3.0));
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
  double t3;
  double t32;
  double t5;
  double t6;
  double t7;
  double t8;
  double t9;
  {
    t1 = 1/delta_;
    t3 = tanh(x*t1);
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
    t25 = -0.9E1*t21+0.225E2*t7-0.15E2*t13+3.0;
    t32 = t25*t25;
    return((-0.1268910394E3*t11+0.2537820788E3*t15-0.1268910394E3*t18)/t25-(
0.113810037E2-0.5075641578E2*t21+0.1268910394E3*t7-0.845940263E2*t13)/t32*(
-0.225E2*t11+0.45E2*t15-0.225E2*t18));
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
  double t11;
  double t13;
  double t17;
  double t18;
  double t19;
  double t3;
  double t37;
  double t39;
  double t40;
  double t5;
  double t51;
  double t52;
  double t59;
  double t6;
  double t62;
  double t65;
  double t69;
  double t7;
  double t79;
  double t8;
  {
    t1 = 1/delta_;
    t3 = tanh(x*t1);
    t5 = 0.5*t3+0.5;
    t6 = t5*t5;
    t7 = t6*t6;
    t8 = t7*t5;
    t11 = t6*t5;
    t13 = 0.113810037E2-0.5075641578E2*t8+0.1268910394E3*t7-0.845940263E2*t11;
    t17 = -0.9E1*t8+0.225E2*t7-0.15E2*t11+3.0;
    t18 = 1/t17;
    t19 = t13*t18;
    t37 = 30.0*t7-60.0*t11+30.0*t6;
    t39 = log(t19);
    t40 = t19*t39;
    t51 = t3*t3;
    t52 = 1.0-t51;
    t59 = t7*t52*t1;
    t62 = t11*t52*t1;
    t65 = t6*t52*t1;
    t69 = t17*t17;
    t79 = t13*t13;
    return(2.0*t19*A_*(4.0*t11+3.0*(2.0*alpha_-2.0)*t6+2.0*(-3.0*alpha_+1.0)*t5)*
t1+beta_*(t37*(-0.25E1*t19+0.15E1*t40+1.0)-t37*(-7.0*t19+3.0*t40+0.945940263E1))
+0.1E1*t1*t3*t52/t13*t17+0.5*((-0.1268910394E3*t59+0.2537820788E3*t62
-0.1268910394E3*t65)*t18-t13/t69*(-0.225E2*t59+0.45E2*t62-0.225E2*t65))*t52/t79
*t69);
  }
}


inline double phiSource ( double x ) const
{
  double t1;
  double t11;
  double t13;
  double t17;
  double t18;
  double t19;
  double t3;
  double t37;
  double t39;
  double t40;
  double t5;
  double t51;
  double t52;
  double t53;
  double t59;
  double t6;
  double t62;
  double t65;
  double t69;
  double t7;
  double t79;
  double t8;
  {
    t1 = 1/delta_;
    t3 = tanh(x*t1);
    t5 = 0.5*t3+0.5;
    t6 = t5*t5;
    t7 = t6*t6;
    t8 = t7*t5;
    t11 = t6*t5;
    t13 = 0.113810037E2-0.5075641578E2*t8+0.1268910394E3*t7-0.845940263E2*t11;
    t17 = -0.9E1*t8+0.225E2*t7-0.15E2*t11+3.0;
    t18 = 1/t17;
    t19 = t13*t18;
    t37 = 30.0*t7-60.0*t11+30.0*t6;
    t39 = log(t19);
    t40 = t19*t39;
    t51 = t3*t3;
    t52 = 1.0-t51;
    t53 = 1/t13;
    t59 = t7*t52*t1;
    t62 = t11*t52*t1;
    t65 = t6*t52*t1;
    t69 = t17*t17;
    t79 = t13*t13;
    return((2.0*t19*A_*(4.0*t11+3.0*(2.0*alpha_-2.0)*t6+2.0*(-3.0*alpha_+1.0)*t5)*
t1+beta_*(t37*(-0.25E1*t19+0.15E1*t40+1.0)-t37*(-7.0*t19+3.0*t40+0.945940263E1))
+0.1E1*t1*t3*t52*t53*t17+0.5*((-0.1268910394E3*t59+0.2537820788E3*t62
-0.1268910394E3*t65)*t18-t13/t69*(-0.225E2*t59+0.45E2*t62-0.225E2*t65))*t52/t79
*t69)*t53*t17);
  }
}


inline double musol ( double x ) const
{
  double t1;
  double t10;
  double t19;
  double t20;
  double t21;
  double t22;
  double t27;
  double t3;
  double t31;
  double t34;
  double t44;
  double t46;
  double t48;
  double t5;
  double t50;
  double t6;
  double t7;
  {
    t1 = 1/delta_;
    t3 = tanh(x*t1);
    t5 = 0.5*t3+0.5;
    t6 = t5*t5;
    t7 = t6*t6;
    t10 = t6*t5;
    t19 = t7*t5;
    t20 = 6.0*t19;
    t21 = 15.0*t7;
    t22 = 10.0*t10;
    t27 = 0.113810037E2-0.5075641578E2*t19+0.1268910394E3*t7-0.845940263E2*t10;
    t31 = -0.9E1*t19+0.225E2*t7-0.15E2*t10+3.0;
    t34 = log(t27/t31);
    t44 = t3*t3;
    t46 = pow(1.0-t44,2.0);
    t48 = t27*t27;
    t50 = t31*t31;
    return(2.0*A_*(t7+(2.0*alpha_-2.0)*t10+(-3.0*alpha_+1.0)*t6+alpha_)*t1+beta_*((
t20-t21+t22)*(-1.0+0.15E1*t34)+(1.0-t20+t21-t22)*(-4.0+3.0*t34))-0.125*t1*t46/
t48*t50);
  }
}


inline double veloSource ( double x ) const
{
  double t1;
  double t10;
  double t100;
  double t104;
  double t105;
  double t106;
  double t109;
  double t11;
  double t121;
  double t123;
  double t128;
  double t14;
  double t16;
  double t20;
  double t21;
  double t22;
  double t26;
  double t3;
  double t30;
  double t4;
  double t40;
  double t42;
  double t43;
  double t5;
  double t54;
  double t6;
  double t60;
  double t63;
  double t66;
  double t68;
  double t70;
  double t76;
  double t78;
  double t8;
  double t80;
  double t81;
  double t82;
  double t9;
  {
    t1 = 1/delta_;
    t3 = tanh(x*t1);
    t4 = t3*t3;
    t5 = 1.0-t4;
    t6 = t5*t1;
    t8 = 0.5*t3+0.5;
    t9 = t8*t8;
    t10 = t9*t9;
    t11 = t10*t8;
    t14 = t9*t8;
    t16 = 0.113810037E2-0.5075641578E2*t11+0.1268910394E3*t10-0.845940263E2*t14
;
    t20 = -0.9E1*t11+0.225E2*t10-0.15E2*t14+3.0;
    t21 = 1/t20;
    t22 = t16*t21;
    t26 = (2.0*alpha_-2.0)*t9;
    t30 = (-3.0*alpha_+1.0)*t8;
    t40 = 30.0*t10-60.0*t14+30.0*t9;
    t42 = log(t22);
    t43 = t22*t42;
    t54 = 1/t16;
    t60 = t10*t5*t1;
    t63 = t14*t5*t1;
    t66 = t9*t5*t1;
    t68 = -0.1268910394E3*t60+0.2537820788E3*t63-0.1268910394E3*t66;
    t70 = t20*t20;
    t76 = -0.225E2*t60+0.45E2*t63-0.225E2*t66;
    t78 = t68*t21-t16/t70*t76;
    t80 = t16*t16;
    t81 = 1/t80;
    t82 = t81*t70;
    t100 = 0.15E2*t60-0.3E2*t63+0.15E2*t66;
    t104 = 6.0*t11;
    t105 = 15.0*t10;
    t106 = 10.0*t14;
    t109 = t54*t20;
    t121 = delta_*delta_;
    t123 = t5*t5;
    t128 = t1*t123;
    return(-0.5*t6*(2.0*t22*A_*(4.0*t14+3.0*t26+2.0*t30)*t1+beta_*(t40*(-0.25E1*
t22+0.15E1*t43+1.0)-t40*(-7.0*t22+3.0*t43+0.945940263E1))+0.1E1*t1*t3*t5*t54*
t20+0.5*t78*t5*t82)+t22*(2.0*A_*(0.2E1*t63+0.15E1*t26*t6+0.1E1*t30*t6)*t1+beta_*(
t100*(-1.0+0.15E1*t42)+0.15E1*(t104-t105+t106)*t78*t109-t100*(-4.0+3.0*t42)+3.0
*(1.0-t104+t105-t106)*t78*t109)+0.5/t121*t123*t82*t3+0.25*t128/t80/t16*t70*t68
-0.25*t128*t81*t20*t76));
  }
}


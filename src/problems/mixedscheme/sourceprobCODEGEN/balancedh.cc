
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
    t25 = 3.0-0.9E1*t21+0.225E2*t7-0.15E2*t13;
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
  double t35;
  double t37;
  double t38;
  double t49;
  double t5;
  double t50;
  double t57;
  double t6;
  double t60;
  double t63;
  double t67;
  double t7;
  double t77;
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
    t35 = 30.0*t7-60.0*t11+30.0*t6;
    t37 = log(t19);
    t38 = t19*t37;
    t49 = t3*t3;
    t50 = 1.0-t49;
    t57 = t7*t50*t1;
    t60 = t11*t50*t1;
    t63 = t6*t50*t1;
    t67 = t17*t17;
    t77 = t13*t13;
    return(t19*(8.0*t11+6.0*(2.0*alpha_-2.0)*t6+4.0*(-3.0*alpha_+1.0)*t5)*t1+beta_
*(t35*(-0.25E1*t19+0.15E1*t38+1.0)-t35*(-7.0*t19+3.0*t38+0.945940263E1))+0.1E1*
t1*t3*t50/t13*t17+0.5*((-0.1268910394E3*t57+0.2537820788E3*t60-0.1268910394E3*
t63)*t18-t13/t67*(-0.225E2*t57+0.45E2*t60-0.225E2*t63))*t50/t77*t67);
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
  double t35;
  double t37;
  double t38;
  double t49;
  double t5;
  double t50;
  double t51;
  double t57;
  double t6;
  double t60;
  double t63;
  double t67;
  double t7;
  double t77;
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
    t35 = 30.0*t7-60.0*t11+30.0*t6;
    t37 = log(t19);
    t38 = t19*t37;
    t49 = t3*t3;
    t50 = 1.0-t49;
    t51 = 1/t13;
    t57 = t7*t50*t1;
    t60 = t11*t50*t1;
    t63 = t6*t50*t1;
    t67 = t17*t17;
    t77 = t13*t13;
    return((t19*(8.0*t11+6.0*(2.0*alpha_-2.0)*t6+4.0*(-3.0*alpha_+1.0)*t5)*t1+
beta_*(t35*(-0.25E1*t19+0.15E1*t38+1.0)-t35*(-7.0*t19+3.0*t38+0.945940263E1))+
0.1E1*t1*t3*t50*t51*t17+0.5*((-0.1268910394E3*t57+0.2537820788E3*t60
-0.1268910394E3*t63)*t18-t13/t67*(-0.225E2*t57+0.45E2*t60-0.225E2*t63))*t50/t77
*t67)*t51*t17);
  }
}


inline double musol ( double x ) const
{
  double t1;
  double t10;
  double t18;
  double t19;
  double t20;
  double t21;
  double t27;
  double t3;
  double t31;
  double t34;
  double t42;
  double t44;
  double t46;
  double t48;
  double t5;
  double t6;
  double t7;
  {
    t1 = 1/delta_;
    t3 = tanh(x*t1);
    t5 = 0.5*t3+0.5;
    t6 = t5*t5;
    t7 = t6*t6;
    t10 = t6*t5;
    t18 = t7*t5;
    t19 = 6.0*t18;
    t20 = 15.0*t7;
    t21 = 10.0*t10;
    t27 = 0.113810037E2-0.5075641578E2*t18+0.1268910394E3*t7-0.845940263E2*t10;
    t31 = -0.9E1*t18+0.225E2*t7-0.15E2*t10+3.0;
    t34 = log(t27/t31);
    t42 = t3*t3;
    t44 = pow(1.0-t42,2.0);
    t46 = t27*t27;
    t48 = t31*t31;
    return((2.0*t7+2.0*(2.0*alpha_-2.0)*t10+2.0*(-3.0*alpha_+1.0)*t6+2.0*alpha_)*
t1+beta_*(t19-t20+t21)*(-1.0+0.15E1*t34)+(1.0-t19+t20-t21)*(-4.0+3.0*t34)-0.125*
t1*t44/t46*t48);
  }
}


inline double veloSource ( double x ) const
{
  double t1;
  double t10;
  double t101;
  double t102;
  double t103;
  double t11;
  double t118;
  double t120;
  double t125;
  double t14;
  double t16;
  double t20;
  double t21;
  double t22;
  double t26;
  double t3;
  double t30;
  double t38;
  double t4;
  double t40;
  double t41;
  double t5;
  double t52;
  double t58;
  double t6;
  double t61;
  double t64;
  double t66;
  double t68;
  double t74;
  double t76;
  double t78;
  double t79;
  double t8;
  double t80;
  double t9;
  double t96;
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
    t38 = 30.0*t10-60.0*t14+30.0*t9;
    t40 = log(t22);
    t41 = t22*t40;
    t52 = 1/t16;
    t58 = t10*t5*t1;
    t61 = t14*t5*t1;
    t64 = t9*t5*t1;
    t66 = -0.1268910394E3*t58+0.2537820788E3*t61-0.1268910394E3*t64;
    t68 = t20*t20;
    t74 = -0.225E2*t58+0.45E2*t61-0.225E2*t64;
    t76 = t66*t21-t16/t68*t74;
    t78 = t16*t16;
    t79 = 1/t78;
    t80 = t79*t68;
    t96 = 0.15E2*t58-0.3E2*t61+0.15E2*t64;
    t101 = 6.0*t11;
    t102 = 15.0*t10;
    t103 = 10.0*t14;
    t118 = delta_*delta_;
    t120 = t5*t5;
    t125 = t1*t120;
    return(-0.5*t6*(t22*(8.0*t14+6.0*t26+4.0*t30)*t1+beta_*(t38*(1.0-0.25E1*t22+
0.15E1*t41)-t38*(0.945940263E1-7.0*t22+3.0*t41))+0.1E1*t1*t3*t5*t52*t20+0.5*t76
*t5*t80)+t22*((0.4E1*t61+0.3E1*t26*t6+0.2E1*t30*t6)*t1+beta_*t96*(-1.0+0.15E1*
t40)+0.15E1*beta_*(t101-t102+t103)*t76*t52*t20-t96*(-4.0+3.0*t40)+3.0*(1.0-t101+
t102-t103)*t76*t52*t20+0.5/t118*t120*t80*t3+0.25*t125/t78/t16*t68*t66-0.25*t125
*t79*t20*t74));
  }
}


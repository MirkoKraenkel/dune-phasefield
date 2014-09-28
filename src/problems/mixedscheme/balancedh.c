
inline double rhosol ( double x) const
{
  double t10;
  double t13;
  double t5;
  double t7;
  double t8;
  double t9;
  {
    t5 = tanh(x/lambda_/delta_);
    t7 = 0.5*t5+0.5;
    t8 = t7*t7;
    t9 = t8*t8;
    t10 = t9*t7;
    t13 = t8*t7;
    return((0.113810037E2-0.5075641578E2*t10+0.1268910394E3*t9-0.845940263E2*
t13)/(-0.9E1*t10+0.225E2*t9-0.15E2*t13+3.0));
  }
}


inline double gradrho ( double x) const
{
  double t1;
  double t10;
  double t11;
  double t13;
  double t14;
  double t16;
  double t18;
  double t21;
  double t24;
  double t28;
  double t3;
  double t35;
  double t5;
  double t7;
  double t8;
  double t9;
  {
    t1 = 1/lambda_;
    t3 = 1/delta_;
    t5 = tanh(x*t1*t3);
    t7 = 0.5*t5+0.5;
    t8 = t7*t7;
    t9 = t8*t8;
    t10 = t5*t5;
    t11 = 1.0-t10;
    t13 = t1*t3;
    t14 = t9*t11*t13;
    t16 = t8*t7;
    t18 = t16*t11*t13;
    t21 = t8*t11*t13;
    t24 = t9*t7;
    t28 = -0.9E1*t24+0.225E2*t9-0.15E2*t16+3.0;
    t35 = t28*t28;
    return((-0.1268910394E3*t14+0.2537820788E3*t18-0.1268910394E3*t21)/t28-(
0.113810037E2-0.5075641578E2*t24+0.1268910394E3*t9-0.845940263E2*t16)/t35*(
-0.225E2*t14+0.45E2*t18-0.225E2*t21));
  }
}


inline double gradphi ( double x) const
{
  {
    return(0.5*(1.0-pow(tanh(x/lambda_/delta_),2.0))/lambda_/delta_);
  }
}


inline double thetasol ( double x) const
{
  double t1;
  double t10;
  double t13;
  double t15;
  double t16;
  double t2;
  double t20;
  double t21;
  double t23;
  double t24;
  double t35;
  double t36;
  double t38;
  double t39;
  double t5;
  double t50;
  double t51;
  double t53;
  double t59;
  double t6;
  double t60;
  double t63;
  double t66;
  double t7;
  double t70;
  double t8;
  double t80;
  double t9;
  {
    t1 = 1/delta_;
    t2 = 1/lambda_;
    t5 = tanh(x*t2*t1);
    t6 = 0.5*t5;
    t7 = t6+0.5;
    t8 = t7*t7;
    t9 = t8*t8;
    t10 = t9*t7;
    t13 = t8*t7;
    t15 = 0.113810037E2-0.5075641578E2*t10+0.1268910394E3*t9-0.845940263E2*t13;
    t16 = t1*t15;
    t20 = -0.9E1*t10+0.225E2*t9-0.15E2*t13+3.0;
    t21 = 1/t20;
    t23 = 0.5-t6;
    t24 = t23*t23;
    t35 = 30.0*t9-60.0*t13+30.0*t8;
    t36 = t15*t21;
    t38 = log(t36);
    t39 = t36*t38;
    t50 = t5*t5;
    t51 = 1.0-t50;
    t53 = lambda_*lambda_;
    t59 = t2*t1;
    t60 = t9*t51*t59;
    t63 = t13*t51*t59;
    t66 = t8*t51*t59;
    t70 = t20*t20;
    t80 = t15*t15;
    return(4.0*t16*t21*t7*t24-4.0*t16*t21*t8*t23+t35*(1.0-0.25E1*t36+0.15E1*t39
)-t35*(0.945940263E1-7.0*t36+3.0*t39)+0.1E1*t1/t15*t20*t5*t51/t53+0.5*((
-0.1268910394E3*t60+0.2537820788E3*t63-0.1268910394E3*t66)*t21-t15/t70*(
-0.225E2*t60+0.45E2*t63-0.225E2*t66))*t51*t2/t80*t70);
  }
}


inline double phiSource ( double x) const
{
  double t1;
  double t10;
  double t13;
  double t15;
  double t16;
  double t2;
  double t20;
  double t21;
  double t23;
  double t24;
  double t35;
  double t36;
  double t38;
  double t39;
  double t47;
  double t5;
  double t50;
  double t51;
  double t53;
  double t59;
  double t6;
  double t60;
  double t63;
  double t66;
  double t7;
  double t70;
  double t8;
  double t80;
  double t9;
  {
    t1 = 1/delta_;
    t2 = 1/lambda_;
    t5 = tanh(x*t2*t1);
    t6 = 0.5*t5;
    t7 = t6+0.5;
    t8 = t7*t7;
    t9 = t8*t8;
    t10 = t9*t7;
    t13 = t8*t7;
    t15 = 0.113810037E2-0.5075641578E2*t10+0.1268910394E3*t9-0.845940263E2*t13;
    t16 = t1*t15;
    t20 = -0.9E1*t10+0.225E2*t9-0.15E2*t13+3.0;
    t21 = 1/t20;
    t23 = 0.5-t6;
    t24 = t23*t23;
    t35 = 30.0*t9-60.0*t13+30.0*t8;
    t36 = t15*t21;
    t38 = log(t36);
    t39 = t36*t38;
    t47 = 1/t15;
    t50 = t5*t5;
    t51 = 1.0-t50;
    t53 = lambda_*lambda_;
    t59 = t2*t1;
    t60 = t9*t51*t59;
    t63 = t13*t51*t59;
    t66 = t8*t51*t59;
    t70 = t20*t20;
    t80 = t15*t15;
    return((4.0*t16*t21*t7*t24-4.0*t16*t21*t8*t23+t35*(1.0-0.25E1*t36+0.15E1*
t39)-t35*(0.945940263E1-7.0*t36+3.0*t39)+0.1E1*t1*t47*t20*t5*t51/t53+0.5*((
-0.1268910394E3*t60+0.2537820788E3*t63-0.1268910394E3*t66)*t21-t15/t70*(
-0.225E2*t60+0.45E2*t63-0.225E2*t66))*t51*t2/t80*t70)*t47*t20);
  }
}


inline double musol ( double x) const
{
  double t1;
  double t11;
  double t14;
  double t15;
  double t16;
  double t17;
  double t18;
  double t19;
  double t24;
  double t28;
  double t31;
  double t39;
  double t41;
  double t43;
  double t45;
  double t48;
  double t5;
  double t6;
  double t7;
  double t8;
  {
    t1 = 1/delta_;
    t5 = tanh(x/lambda_*t1);
    t6 = 0.5*t5;
    t7 = t6+0.5;
    t8 = t7*t7;
    t11 = pow(0.5-t6,2.0);
    t14 = t8*t8;
    t15 = t14*t7;
    t16 = 6.0*t15;
    t17 = 15.0*t14;
    t18 = t8*t7;
    t19 = 10.0*t18;
    t24 = 0.113810037E2-0.5075641578E2*t15+0.1268910394E3*t14-0.845940263E2*t18
;
    t28 = -0.9E1*t15+0.225E2*t14-0.15E2*t18+3.0;
    t31 = log(t24/t28);
    t39 = t5*t5;
    t41 = pow(1.0-t39,2.0);
    t43 = lambda_*lambda_;
    t45 = t24*t24;
    t48 = t28*t28;
    return(2.0*t1*t8*t11+(t16-t17+t19)*(-1.0+0.15E1*t31)+(1.0-t16+t17-t19)*(
-4.0+3.0*t31)-0.125*t1*t41/t43/t45*t48);
  }
}


inline double veloSource ( double x) const
{
  double t1;
  double t10;
  double t13;
  double t15;
  double t19;
  double t20;
  double t21;
  double t22;
  double t23;
  double t25;
  double t26;
  double t27;
  double t28;
  double t3;
  double t39;
  double t40;
  double t43;
  double t46;
  double t48;
  double t49;
  double t5;
  double t53;
  double t54;
  double t55;
  double t6;
  double t60;
  double t62;
  double t68;
  double t7;
  double t70;
  double t73;
  double t8;
  double t83;
  double t85;
  double t89;
  double t9;
  double t90;
  double t97;
  {
    t1 = 1/lambda_;
    t3 = 1/delta_;
    t5 = tanh(x*t1*t3);
    t6 = 0.5*t5;
    t7 = t6+0.5;
    t8 = t7*t7;
    t9 = t8*t8;
    t10 = t9*t7;
    t13 = t8*t7;
    t15 = 0.113810037E2-0.5075641578E2*t10+0.1268910394E3*t9-0.845940263E2*t13;
    t19 = -0.9E1*t10+0.225E2*t9-0.15E2*t13+3.0;
    t20 = 1/t19;
    t21 = t15*t20;
    t22 = delta_*delta_;
    t23 = 1/t22;
    t25 = 0.5-t6;
    t26 = t25*t25;
    t27 = t5*t5;
    t28 = 1.0-t27;
    t39 = t1*t3;
    t40 = t9*t28*t39;
    t43 = t13*t28*t39;
    t46 = t8*t28*t39;
    t48 = 0.15E2*t40-0.3E2*t43+0.15E2*t46;
    t49 = log(t21);
    t53 = 6.0*t10;
    t54 = 15.0*t9;
    t55 = 10.0*t13;
    t60 = -0.1268910394E3*t40+0.2537820788E3*t43-0.1268910394E3*t46;
    t62 = t19*t19;
    t68 = -0.225E2*t40+0.45E2*t43-0.225E2*t46;
    t70 = t60*t20-t15/t62*t68;
    t73 = 1/t15*t19;
    t83 = t28*t28;
    t85 = lambda_*lambda_;
    t89 = t15*t15;
    t90 = 1/t89;
    t97 = t3*t83/t85;
    return(t21*(0.2E1*t23*t7*t26*t28*t1-0.2E1*t23*t8*t25*t28*t1+t48*(-1.0+
0.15E1*t49)+0.15E1*(t53-t54+t55)*t70*t73-t48*(-4.0+3.0*t49)+3.0*(1.0-t53+t54-
t55)*t70*t73+0.5*t23*t83/t85/lambda_*t90*t62*t5+0.25*t97/t89/t15*t62*t60-0.25*
t97*t90*t19*t68));
  }
}


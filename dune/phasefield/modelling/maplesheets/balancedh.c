
inline double rhosol ( double x) const
{
  double t11;
  double t3;
  double t5;
  double t6;
  double t7;
  double t8;
  {
    t3 = tanh(x/delta__);
    t5 = 0.5*t3+0.5;
    t6 = t5*t5;
    t7 = t6*t6;
    t8 = t7*t5;
    t11 = t6*t5;
    return((0.113810037E2-0.5075641578E2*t8+0.1268910394E3*t7-0.845940263E2*t11
)/(-0.9E1*t8+0.225E2*t7-0.15E2*t11+3.0));
  }
}


inline double gradrho ( double x) const
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
    t1 = 1/delta__;
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


inline double gradphi ( double x) const
{
  {
    return(0.5*(1.0-pow(tanh(x/delta__),2.0))/delta__);
  }
}


double thetasol1(double x)
{
  double t1;
  double t11;
  double t13;
  double t19;
  double t29;
  double t3;
  double t30;
  double t32;
  double t33;
  double t5;
  double t6;
  double t7;
  double t8;
  {
    t1 = 1/delta__;
    t3 = tanh(x*t1);
    t5 = 0.5*t3+0.5;
    t6 = t5*t5;
    t7 = t6*t6;
    t8 = t7*t5;
    t11 = t6*t5;
    t13 = 0.113810037E2-0.5075641578E2*t8+0.1268910394E3*t7-0.845940263E2*t11;
    t19 = 1/(-0.9E1*t8+0.225E2*t7-0.15E2*t11+3.0);
    t29 = 30.0*t7-60.0*t11+30.0*t6;
    t30 = t13*t19;
    t32 = log(t30);
    t33 = t30*t32;
    return(t1*t13*t19*(8.0*t11-12.0*t6+0.2E1*t3+0.2E1)+t29*(-0.25E1*t30+0.15E1*
t33+1.0)-t29*(-7.0*t30+3.0*t33+0.945940263E1));
  }
}


inline double thetasol2 ( double x) const
{
  double t1;
  double t11;
  double t13;
  double t19;
  double t21;
  double t22;
  double t27;
  double t3;
  double t30;
  double t33;
  double t38;
  double t48;
  double t5;
  double t6;
  double t7;
  double t8;
  {
    t1 = 1/delta__;
    t3 = tanh(x*t1);
    t5 = 0.5*t3+0.5;
    t6 = t5*t5;
    t7 = t6*t6;
    t8 = t7*t5;
    t11 = t6*t5;
    t13 = 0.113810037E2-0.5075641578E2*t8+0.1268910394E3*t7-0.845940263E2*t11;
    t19 = -0.9E1*t8+0.225E2*t7-0.15E2*t11+3.0;
    t21 = t3*t3;
    t22 = 1.0-t21;
    t27 = t7*t22*t1;
    t30 = t11*t22*t1;
    t33 = t6*t22*t1;
    t38 = t19*t19;
    t48 = t13*t13;
    return(0.1E1*t1/t13*t19*t3*t22+0.5*((-0.1268910394E3*t27+0.2537820788E3*t30
-0.1268910394E3*t33)/t19-t13/t38*(-0.225E2*t27+0.45E2*t30-0.225E2*t33))*t22/t48
*t38);
  }
}


inline double phiSource ( double x) const
{
  double t1;
  double t11;
  double t13;
  double t18;
  double t19;
  double t29;
  double t3;
  double t30;
  double t32;
  double t33;
  double t41;
  double t44;
  double t45;
  double t5;
  double t50;
  double t53;
  double t56;
  double t6;
  double t60;
  double t7;
  double t70;
  double t8;
  {
    t1 = 1/delta__;
    t3 = tanh(x*t1);
    t5 = 0.5*t3+0.5;
    t6 = t5*t5;
    t7 = t6*t6;
    t8 = t7*t5;
    t11 = t6*t5;
    t13 = 0.113810037E2-0.5075641578E2*t8+0.1268910394E3*t7-0.845940263E2*t11;
    t18 = -0.9E1*t8+0.225E2*t7-0.15E2*t11+3.0;
    t19 = 1/t18;
    t29 = 30.0*t7-60.0*t11+30.0*t6;
    t30 = t13*t19;
    t32 = log(t30);
    t33 = t30*t32;
    t41 = 1/t13;
    t44 = t3*t3;
    t45 = 1.0-t44;
    t50 = t7*t45*t1;
    t53 = t11*t45*t1;
    t56 = t6*t45*t1;
    t60 = t18*t18;
    t70 = t13*t13;
    return((t1*t13*t19*(8.0*t11-12.0*t6+0.2E1*t3+0.2E1)+t29*(-0.25E1*t30+0.15E1
*t33+1.0)-t29*(-7.0*t30+3.0*t33+0.945940263E1)+0.1E1*t1*t41*t18*t3*t45+0.5*((
-0.1268910394E3*t50+0.2537820788E3*t53-0.1268910394E3*t56)*t19-t13/t60*(
-0.225E2*t50+0.45E2*t53-0.225E2*t56))*t45/t70*t60)*t41*t18);
  }
}


inline double musol ( double x) const
{
  double t1;
  double t14;
  double t15;
  double t16;
  double t17;
  double t22;
  double t26;
  double t29;
  double t3;
  double t37;
  double t39;
  double t41;
  double t43;
  double t5;
  double t6;
  double t7;
  double t9;
  {
    t1 = 1/delta__;
    t3 = tanh(x*t1);
    t5 = 0.5*t3+0.5;
    t6 = t5*t5;
    t7 = t6*t6;
    t9 = t6*t5;
    t14 = t7*t5;
    t15 = 6.0*t14;
    t16 = 15.0*t7;
    t17 = 10.0*t9;
    t22 = 0.113810037E2-0.5075641578E2*t14+0.1268910394E3*t7-0.845940263E2*t9;
    t26 = 3.0-0.9E1*t14+0.225E2*t7-0.15E2*t9;
    t29 = log(t22/t26);
    t37 = t3*t3;
    t39 = pow(1.0-t37,2.0);
    t41 = t22*t22;
    t43 = t26*t26;
    return(t1*(2.0*t7-4.0*t9+2.0*t6)+(t15-t16+t17)*(-1.0+0.15E1*t29)+(1.0-t15+
t16-t17)*(-4.0+3.0*t29)-0.125*t1*t39/t41*t43);
  }
}


inline double veloSource ( double x) const
{
  double t1;
  double t10;
  double t108;
  double t11;
  double t110;
  double t115;
  double t14;
  double t16;
  double t21;
  double t22;
  double t3;
  double t32;
  double t33;
  double t35;
  double t36;
  double t4;
  double t44;
  double t5;
  double t51;
  double t54;
  double t57;
  double t59;
  double t61;
  double t67;
  double t69;
  double t71;
  double t72;
  double t73;
  double t8;
  double t89;
  double t9;
  double t93;
  double t94;
  double t95;
  double t98;
  {
    t1 = 1/delta__;
    t3 = tanh(x*t1);
    t4 = t3*t3;
    t5 = 1.0-t4;
    t8 = 0.5*t3+0.5;
    t9 = t8*t8;
    t10 = t9*t9;
    t11 = t10*t8;
    t14 = t9*t8;
    t16 = 0.113810037E2-0.5075641578E2*t11+0.1268910394E3*t10-0.845940263E2*t14
;
    t21 = -0.9E1*t11+0.225E2*t10-0.15E2*t14+3.0;
    t22 = 1/t21;
    t32 = 30.0*t10-60.0*t14+30.0*t9;
    t33 = t16*t22;
    t35 = log(t33);
    t36 = t33*t35;
    t44 = 1/t16;
    t51 = t10*t5*t1;
    t54 = t14*t5*t1;
    t57 = t5*t9*t1;
    t59 = -0.1268910394E3*t51+0.2537820788E3*t54-0.1268910394E3*t57;
    t61 = t21*t21;
    t67 = -0.225E2*t51+0.45E2*t54-0.225E2*t57;
    t69 = t59*t22-t16/t61*t67;
    t71 = t16*t16;
    t72 = 1/t71;
    t73 = t72*t61;
    t89 = 0.15E2*t51-0.3E2*t54+0.15E2*t57;
    t93 = 6.0*t11;
    t94 = 15.0*t10;
    t95 = 10.0*t14;
    t98 = t44*t21;
    t108 = delta__*delta__;
    t110 = t5*t5;
    t115 = t1*t110;
    return(-0.5*t5*t1*(t1*t16*t22*(8.0*t14-12.0*t9+0.2E1*t3+0.2E1)+t32*(-0.25E1
*t33+0.15E1*t36+1.0)-t32*(-7.0*t33+3.0*t36+0.945940263E1)+0.1E1*t1*t44*t21*t3*
t5+0.5*t69*t5*t73)+t33*(t1*(0.4E1*t54-0.6E1*t57+0.2E1*t8*t5*t1)+t89*(-1.0+
0.15E1*t35)+0.15E1*(t93-t94+t95)*t69*t98-t89*(-4.0+3.0*t35)+3.0*(1.0-t93+t94-
t95)*t69*t98+0.5/t108*t110*t73*t3+0.25*t115/t71/t16*t61*t59-0.25*t115*t72*t21*
t67));
  }
}


inline double rhosol ( double x) const
{
  double t10;
  double t2;
  double t4;
  double t5;
  double t6;
  double t7;
  {
    t2 = tanh(0.1E2*x);
    t4 = 0.5*t2+0.5;
    t5 = t4*t4;
    t6 = t5*t5;
    t7 = t6*t4;
    t10 = t5*t4;
    return((0.113810037E2-0.5075641578E2*t7+0.1268910394E3*t6-0.845940263E2*t10
)/(-0.9E1*t7+0.225E2*t6-0.15E2*t10+3.0));
  }
}


inline double gradrho ( double x) const
{
  double t10;
  double t12;
  double t13;
  double t15;
  double t18;
  double t2;
  double t22;
  double t29;
  double t4;
  double t5;
  double t6;
  double t7;
  double t9;
  {
    t2 = tanh(0.1E2*x);
    t4 = 0.5*t2+0.5;
    t5 = t4*t4;
    t6 = t5*t5;
    t7 = t2*t2;
    t9 = 5.0-5.0*t7;
    t10 = t6*t9;
    t12 = t5*t4;
    t13 = t12*t9;
    t15 = t5*t9;
    t18 = t6*t4;
    t22 = -0.9E1*t18+0.225E2*t6-0.15E2*t12+3.0;
    t29 = t22*t22;
    return((-0.2537820789E3*t10+0.5075641576E3*t13-0.2537820789E3*t15)/t22-(
0.113810037E2-0.5075641578E2*t18+0.1268910394E3*t6-0.845940263E2*t12)/t29*(
-0.45E2*t10+0.9E2*t13-0.45E2*t15));
  }
}


inline double gradphi ( double x) const
{
  {
    return(5.0-5.0*pow(tanh(0.1E2*x),2.0));
  }
}


double thetasol1(double x)
{
  double t10;
  double t18;
  double t2;
  double t28;
  double t30;
  double t31;
  double t4;
  double t5;
  double t6;
  double t7;
  {
    t2 = tanh(0.1E2*x);
    t4 = 0.5*t2+0.5;
    t5 = t4*t4;
    t6 = t5*t5;
    t7 = t6*t4;
    t10 = t5*t4;
    t18 = (0.113810037E2-0.5075641578E2*t7+0.1268910394E3*t6-0.845940263E2*t10)
/(-0.9E1*t7+0.225E2*t6-0.15E2*t10+3.0);
    t28 = 30.0*t6-60.0*t10+30.0*t5;
    t30 = log(t18);
    t31 = t18*t30;
    return(0.1E2*t18*(8.0*t10-12.0*t5+0.2E1*t2+0.2E1)+t28*(-0.25E1*t18+0.15E1*
t31+1.0)-t28*(-7.0*t18+3.0*t31+0.945940263E1));
  }
}


inline double thetasol2 ( double x) const
{
  double t10;
  double t12;
  double t17;
  double t19;
  double t2;
  double t26;
  double t27;
  double t29;
  double t31;
  double t36;
  double t4;
  double t46;
  double t5;
  double t6;
  double t7;
  {
    t2 = tanh(0.1E2*x);
    t4 = 0.5*t2+0.5;
    t5 = t4*t4;
    t6 = t5*t5;
    t7 = t6*t4;
    t10 = t5*t4;
    t12 = 0.113810037E2-0.5075641578E2*t7+0.1268910394E3*t6-0.845940263E2*t10;
    t17 = -0.9E1*t7+0.225E2*t6-0.15E2*t10+3.0;
    t19 = t2*t2;
    t26 = 5.0-5.0*t19;
    t27 = t6*t26;
    t29 = t10*t26;
    t31 = t5*t26;
    t36 = t17*t17;
    t46 = t12*t12;
    return(0.1E1/t12*t17*t2*(0.1E2-0.1E2*t19)+0.1*((-0.2537820789E3*t27+
0.5075641576E3*t29-0.2537820789E3*t31)/t17-t12/t36*(-0.45E2*t27+0.9E2*t29
-0.45E2*t31))*t26/t46*t36);
  }
}


inline double phiSource ( double x) const
{
  double t10;
  double t12;
  double t16;
  double t17;
  double t18;
  double t2;
  double t28;
  double t30;
  double t31;
  double t39;
  double t4;
  double t41;
  double t48;
  double t49;
  double t5;
  double t51;
  double t53;
  double t57;
  double t6;
  double t67;
  double t7;
  {
    t2 = tanh(0.1E2*x);
    t4 = 0.5*t2+0.5;
    t5 = t4*t4;
    t6 = t5*t5;
    t7 = t6*t4;
    t10 = t5*t4;
    t12 = 0.113810037E2-0.5075641578E2*t7+0.1268910394E3*t6-0.845940263E2*t10;
    t16 = -0.9E1*t7+0.225E2*t6-0.15E2*t10+3.0;
    t17 = 1/t16;
    t18 = t12*t17;
    t28 = 30.0*t6-60.0*t10+30.0*t5;
    t30 = log(t18);
    t31 = t18*t30;
    t39 = 1/t12;
    t41 = t2*t2;
    t48 = 5.0-5.0*t41;
    t49 = t6*t48;
    t51 = t10*t48;
    t53 = t5*t48;
    t57 = t16*t16;
    t67 = t12*t12;
    return((0.1E2*t18*(8.0*t10-12.0*t5+0.2E1*t2+0.2E1)+t28*(-0.25E1*t18+0.15E1*
t31+1.0)-t28*(-7.0*t18+3.0*t31+0.945940263E1)+0.1E1*t39*t16*t2*(0.1E2-0.1E2*t41
)+0.1*((-0.2537820789E3*t49+0.5075641576E3*t51-0.2537820789E3*t53)*t17-t12/t57*
(-0.45E2*t49+0.9E2*t51-0.45E2*t53))*t48/t67*t57)*t39*t16);
  }
}


inline double musol ( double x) const
{
  double t11;
  double t12;
  double t13;
  double t14;
  double t19;
  double t2;
  double t23;
  double t26;
  double t34;
  double t37;
  double t38;
  double t4;
  double t41;
  double t5;
  double t6;
  double t8;
  {
    t2 = tanh(0.1E2*x);
    t4 = 0.5*t2+0.5;
    t5 = t4*t4;
    t6 = t5*t5;
    t8 = t5*t4;
    t11 = t6*t4;
    t12 = 6.0*t11;
    t13 = 15.0*t6;
    t14 = 10.0*t8;
    t19 = 0.113810037E2-0.5075641578E2*t11+0.1268910394E3*t6-0.845940263E2*t8;
    t23 = -0.9E1*t11+0.225E2*t6-0.15E2*t8+3.0;
    t26 = log(t19/t23);
    t34 = t2*t2;
    t37 = pow(5.0-5.0*t34,2.0);
    t38 = t19*t19;
    t41 = t23*t23;
    return(0.2E2*t6-0.4E2*t8+0.2E2*t5+(t12-t13+t14)*(-1.0+0.15E1*t26)+(1.0-t12+
t13-t14)*(-4.0+3.0*t26)-0.5E-1*t37/t38*t41);
  }
}


inline double veloSource ( double x) const
{
  double t10;
  double t104;
  double t13;
  double t15;
  double t19;
  double t2;
  double t20;
  double t21;
  double t3;
  double t31;
  double t33;
  double t34;
  double t43;
  double t45;
  double t49;
  double t5;
  double t51;
  double t53;
  double t55;
  double t57;
  double t63;
  double t65;
  double t67;
  double t68;
  double t7;
  double t8;
  double t81;
  double t85;
  double t86;
  double t87;
  double t9;
  {
    t2 = tanh(0.1E2*x);
    t3 = t2*t2;
    t5 = 5.0-5.0*t3;
    t7 = 0.5*t2+0.5;
    t8 = t7*t7;
    t9 = t8*t8;
    t10 = t9*t7;
    t13 = t8*t7;
    t15 = 0.113810037E2-0.5075641578E2*t10+0.1268910394E3*t9-0.845940263E2*t13;
    t19 = -0.9E1*t10+0.225E2*t9-0.15E2*t13+3.0;
    t20 = 1/t19;
    t21 = t15*t20;
    t31 = 30.0*t9-60.0*t13+30.0*t8;
    t33 = log(t21);
    t34 = t21*t33;
    t43 = 1/t15*t19;
    t45 = 0.1E2-0.1E2*t3;
    t49 = t5*t9;
    t51 = t13*t5;
    t53 = t8*t5;
    t55 = -0.2537820789E3*t49+0.5075641576E3*t51-0.2537820789E3*t53;
    t57 = t19*t19;
    t63 = -0.45E2*t49+0.9E2*t51-0.45E2*t53;
    t65 = t55*t20-t15/t57*t63;
    t67 = t15*t15;
    t68 = 1/t67;
    t81 = 30.0*t49-60.0*t51+30.0*t53;
    t85 = 6.0*t10;
    t86 = 15.0*t9;
    t87 = 10.0*t13;
    t104 = t5*t5;
    return(-t5*(0.1E2*t21*(8.0*t13-12.0*t8+0.2E1*t2+0.2E1)+t31*(-0.25E1*t21+
0.15E1*t34+1.0)-t31*(-7.0*t21+3.0*t34+0.945940263E1)+0.1E1*t43*t2*t45+0.1*t65*
t5*t68*t57)+t21*(0.8E2*t51-0.12E3*t53+0.4E2*t7*t5+t81*(-1.0+0.15E1*t33)+0.15E1*
(t85-t86+t87)*t65*t43-t81*(-4.0+3.0*t33)+3.0*(1.0-t85+t86-t87)*t65*t43+0.1E1*t5
*t68*t57*t2*t45+0.1*t104/t67/t15*t57*t55-0.1*t104*t68*t19*t63));
  }
}


inline double rhosol ( double x) const
{
  double t10;
  double t2;
  double t4;
  double t5;
  double t6;
  double t7;
  {
    t2 = tanh(0.1E2*x);
    t4 = 0.5*t2+0.5;
    t5 = t4*t4;
    t6 = t5*t5;
    t7 = t6*t4;
    t10 = t5*t4;
    return((0.113810037E2-0.5075641578E2*t7+0.1268910394E3*t6-0.845940263E2*t10
)/(-0.9E1*t7+0.225E2*t6-0.15E2*t10+3.0));
  }
}


inline double gradrho ( double x) const
{
  double t10;
  double t12;
  double t13;
  double t15;
  double t18;
  double t2;
  double t22;
  double t29;
  double t4;
  double t5;
  double t6;
  double t7;
  double t9;
  {
    t2 = tanh(0.1E2*x);
    t4 = 0.5*t2+0.5;
    t5 = t4*t4;
    t6 = t5*t5;
    t7 = t2*t2;
    t9 = 5.0-5.0*t7;
    t10 = t6*t9;
    t12 = t5*t4;
    t13 = t12*t9;
    t15 = t5*t9;
    t18 = t6*t4;
    t22 = -0.9E1*t18+0.225E2*t6-0.15E2*t12+3.0;
    t29 = t22*t22;
    return((-0.2537820789E3*t10+0.5075641576E3*t13-0.2537820789E3*t15)/t22-(
0.113810037E2-0.5075641578E2*t18+0.1268910394E3*t6-0.845940263E2*t12)/t29*(
-0.45E2*t10+0.9E2*t13-0.45E2*t15));
  }
}


inline double gradphi ( double x) const
{
  {
    return(5.0-5.0*pow(tanh(0.1E2*x),2.0));
  }
}


double thetasol1(double x)
{
  double t10;
  double t18;
  double t2;
  double t28;
  double t30;
  double t31;
  double t4;
  double t5;
  double t6;
  double t7;
  {
    t2 = tanh(0.1E2*x);
    t4 = 0.5*t2+0.5;
    t5 = t4*t4;
    t6 = t5*t5;
    t7 = t6*t4;
    t10 = t5*t4;
    t18 = (0.113810037E2-0.5075641578E2*t7+0.1268910394E3*t6-0.845940263E2*t10)
/(-0.9E1*t7+0.225E2*t6-0.15E2*t10+3.0);
    t28 = 30.0*t6-60.0*t10+30.0*t5;
    t30 = log(t18);
    t31 = t18*t30;
    return(0.1E2*t18*(8.0*t10-12.0*t5+0.2E1*t2+0.2E1)+t28*(-0.25E1*t18+0.15E1*
t31+1.0)-t28*(-7.0*t18+3.0*t31+0.945940263E1));
  }
}


inline double thetasol2 ( double x) const
{
  double t10;
  double t12;
  double t17;
  double t19;
  double t2;
  double t26;
  double t27;
  double t29;
  double t31;
  double t36;
  double t4;
  double t46;
  double t5;
  double t6;
  double t7;
  {
    t2 = tanh(0.1E2*x);
    t4 = 0.5*t2+0.5;
    t5 = t4*t4;
    t6 = t5*t5;
    t7 = t6*t4;
    t10 = t5*t4;
    t12 = 0.113810037E2-0.5075641578E2*t7+0.1268910394E3*t6-0.845940263E2*t10;
    t17 = -0.9E1*t7+0.225E2*t6-0.15E2*t10+3.0;
    t19 = t2*t2;
    t26 = 5.0-5.0*t19;
    t27 = t6*t26;
    t29 = t10*t26;
    t31 = t5*t26;
    t36 = t17*t17;
    t46 = t12*t12;
    return(0.1E1/t12*t17*t2*(0.1E2-0.1E2*t19)+0.1*((-0.2537820789E3*t27+
0.5075641576E3*t29-0.2537820789E3*t31)/t17-t12/t36*(-0.45E2*t27+0.9E2*t29
-0.45E2*t31))*t26/t46*t36);
  }
}


inline double phiSource ( double x) const
{
  double t10;
  double t12;
  double t16;
  double t17;
  double t18;
  double t2;
  double t28;
  double t30;
  double t31;
  double t39;
  double t4;
  double t41;
  double t48;
  double t49;
  double t5;
  double t51;
  double t53;
  double t57;
  double t6;
  double t67;
  double t7;
  {
    t2 = tanh(0.1E2*x);
    t4 = 0.5*t2+0.5;
    t5 = t4*t4;
    t6 = t5*t5;
    t7 = t6*t4;
    t10 = t5*t4;
    t12 = 0.113810037E2-0.5075641578E2*t7+0.1268910394E3*t6-0.845940263E2*t10;
    t16 = -0.9E1*t7+0.225E2*t6-0.15E2*t10+3.0;
    t17 = 1/t16;
    t18 = t12*t17;
    t28 = 30.0*t6-60.0*t10+30.0*t5;
    t30 = log(t18);
    t31 = t18*t30;
    t39 = 1/t12;
    t41 = t2*t2;
    t48 = 5.0-5.0*t41;
    t49 = t6*t48;
    t51 = t10*t48;
    t53 = t5*t48;
    t57 = t16*t16;
    t67 = t12*t12;
    return((0.1E2*t18*(8.0*t10-12.0*t5+0.2E1*t2+0.2E1)+t28*(-0.25E1*t18+0.15E1*
t31+1.0)-t28*(-7.0*t18+3.0*t31+0.945940263E1)+0.1E1*t39*t16*t2*(0.1E2-0.1E2*t41
)+0.1*((-0.2537820789E3*t49+0.5075641576E3*t51-0.2537820789E3*t53)*t17-t12/t57*
(-0.45E2*t49+0.9E2*t51-0.45E2*t53))*t48/t67*t57)*t39*t16);
  }
}


inline double musol ( double x) const
{
  double t11;
  double t12;
  double t13;
  double t14;
  double t19;
  double t2;
  double t23;
  double t26;
  double t34;
  double t37;
  double t38;
  double t4;
  double t41;
  double t5;
  double t6;
  double t8;
  {
    t2 = tanh(0.1E2*x);
    t4 = 0.5*t2+0.5;
    t5 = t4*t4;
    t6 = t5*t5;
    t8 = t5*t4;
    t11 = t6*t4;
    t12 = 6.0*t11;
    t13 = 15.0*t6;
    t14 = 10.0*t8;
    t19 = 0.113810037E2-0.5075641578E2*t11+0.1268910394E3*t6-0.845940263E2*t8;
    t23 = -0.9E1*t11+0.225E2*t6-0.15E2*t8+3.0;
    t26 = log(t19/t23);
    t34 = t2*t2;
    t37 = pow(5.0-5.0*t34,2.0);
    t38 = t19*t19;
    t41 = t23*t23;
    return(0.2E2*t6-0.4E2*t8+0.2E2*t5+(t12-t13+t14)*(-1.0+0.15E1*t26)+(1.0-t12+
t13-t14)*(-4.0+3.0*t26)-0.5E-1*t37/t38*t41);
  }
}


inline double veloSource ( double x) const
{
  double t10;
  double t104;
  double t13;
  double t15;
  double t19;
  double t2;
  double t20;
  double t21;
  double t3;
  double t31;
  double t33;
  double t34;
  double t43;
  double t45;
  double t49;
  double t5;
  double t51;
  double t53;
  double t55;
  double t57;
  double t63;
  double t65;
  double t67;
  double t68;
  double t7;
  double t8;
  double t81;
  double t85;
  double t86;
  double t87;
  double t9;
  {
    t2 = tanh(0.1E2*x);
    t3 = t2*t2;
    t5 = 5.0-5.0*t3;
    t7 = 0.5*t2+0.5;
    t8 = t7*t7;
    t9 = t8*t8;
    t10 = t9*t7;
    t13 = t8*t7;
    t15 = 0.113810037E2-0.5075641578E2*t10+0.1268910394E3*t9-0.845940263E2*t13;
    t19 = -0.9E1*t10+0.225E2*t9-0.15E2*t13+3.0;
    t20 = 1/t19;
    t21 = t15*t20;
    t31 = 30.0*t9-60.0*t13+30.0*t8;
    t33 = log(t21);
    t34 = t21*t33;
    t43 = 1/t15*t19;
    t45 = 0.1E2-0.1E2*t3;
    t49 = t5*t9;
    t51 = t13*t5;
    t53 = t8*t5;
    t55 = -0.2537820789E3*t49+0.5075641576E3*t51-0.2537820789E3*t53;
    t57 = t19*t19;
    t63 = -0.45E2*t49+0.9E2*t51-0.45E2*t53;
    t65 = t55*t20-t15/t57*t63;
    t67 = t15*t15;
    t68 = 1/t67;
    t81 = 30.0*t49-60.0*t51+30.0*t53;
    t85 = 6.0*t10;
    t86 = 15.0*t9;
    t87 = 10.0*t13;
    t104 = t5*t5;
    return(-t5*(0.1E2*t21*(0.2E1+8.0*t13-12.0*t8+0.2E1*t2)+t31*(1.0-0.25E1*t21+
0.15E1*t34)-t31*(0.945940263E1-7.0*t21+3.0*t34)+0.1E1*t43*t2*t45+0.1*t65*t5*t68
*t57)+t21*(0.8E2*t51-0.12E3*t53+0.4E2*t7*t5+t81*(-1.0+0.15E1*t33)+0.15E1*(t85-
t86+t87)*t65*t43-t81*(-4.0+3.0*t33)+3.0*(1.0-t85+t86-t87)*t65*t43+0.1E1*t5*t68*
t57*t2*t45+0.1*t104/t67/t15*t57*t55-0.1*t104*t68*t19*t63));
  }
}


inline double rhosol ( double x) const
{
  double t10;
  double t2;
  double t4;
  double t5;
  double t6;
  double t7;
  {
    t2 = tanh(0.1E2*x);
    t4 = 0.5*t2+0.5;
    t5 = t4*t4;
    t6 = t5*t5;
    t7 = t6*t4;
    t10 = t5*t4;
    return((0.113810037E2-0.5075641578E2*t7+0.1268910394E3*t6-0.845940263E2*t10
)/(-0.9E1*t7+0.225E2*t6-0.15E2*t10+3.0));
  }
}


inline double gradrho ( double x) const
{
  double t10;
  double t12;
  double t13;
  double t15;
  double t18;
  double t2;
  double t22;
  double t29;
  double t4;
  double t5;
  double t6;
  double t7;
  double t9;
  {
    t2 = tanh(0.1E2*x);
    t4 = 0.5*t2+0.5;
    t5 = t4*t4;
    t6 = t5*t5;
    t7 = t2*t2;
    t9 = 5.0-5.0*t7;
    t10 = t6*t9;
    t12 = t5*t4;
    t13 = t12*t9;
    t15 = t5*t9;
    t18 = t6*t4;
    t22 = -0.9E1*t18+0.225E2*t6-0.15E2*t12+3.0;
    t29 = t22*t22;
    return((-0.2537820789E3*t10+0.5075641576E3*t13-0.2537820789E3*t15)/t22-(
0.113810037E2-0.5075641578E2*t18+0.1268910394E3*t6-0.845940263E2*t12)/t29*(
-0.45E2*t10+0.9E2*t13-0.45E2*t15));
  }
}


inline double gradphi ( double x) const
{
  {
    return(5.0-5.0*pow(tanh(0.1E2*x),2.0));
  }
}


double thetasol1(double x)
{
  double t10;
  double t18;
  double t2;
  double t28;
  double t30;
  double t31;
  double t4;
  double t5;
  double t6;
  double t7;
  {
    t2 = tanh(0.1E2*x);
    t4 = 0.5*t2+0.5;
    t5 = t4*t4;
    t6 = t5*t5;
    t7 = t6*t4;
    t10 = t5*t4;
    t18 = (0.113810037E2-0.5075641578E2*t7+0.1268910394E3*t6-0.845940263E2*t10)
/(-0.9E1*t7+0.225E2*t6-0.15E2*t10+3.0);
    t28 = 30.0*t6-60.0*t10+30.0*t5;
    t30 = log(t18);
    t31 = t18*t30;
    return(0.1E2*t18*(8.0*t10-12.0*t5+0.2E1*t2+0.2E1)+t28*(-0.25E1*t18+0.15E1*
t31+1.0)-t28*(-7.0*t18+3.0*t31+0.945940263E1));
  }
}


inline double thetasol2 ( double x) const
{
  double t10;
  double t12;
  double t17;
  double t19;
  double t2;
  double t26;
  double t27;
  double t29;
  double t31;
  double t36;
  double t4;
  double t46;
  double t5;
  double t6;
  double t7;
  {
    t2 = tanh(0.1E2*x);
    t4 = 0.5*t2+0.5;
    t5 = t4*t4;
    t6 = t5*t5;
    t7 = t6*t4;
    t10 = t5*t4;
    t12 = 0.113810037E2-0.5075641578E2*t7+0.1268910394E3*t6-0.845940263E2*t10;
    t17 = -0.9E1*t7+0.225E2*t6-0.15E2*t10+3.0;
    t19 = t2*t2;
    t26 = 5.0-5.0*t19;
    t27 = t6*t26;
    t29 = t10*t26;
    t31 = t5*t26;
    t36 = t17*t17;
    t46 = t12*t12;
    return(0.1E1/t12*t17*t2*(0.1E2-0.1E2*t19)+0.1*((-0.2537820789E3*t27+
0.5075641576E3*t29-0.2537820789E3*t31)/t17-t12/t36*(-0.45E2*t27+0.9E2*t29
-0.45E2*t31))*t26/t46*t36);
  }
}


inline double phiSource ( double x) const
{
  double t10;
  double t12;
  double t16;
  double t17;
  double t18;
  double t2;
  double t28;
  double t30;
  double t31;
  double t39;
  double t4;
  double t41;
  double t48;
  double t49;
  double t5;
  double t51;
  double t53;
  double t57;
  double t6;
  double t67;
  double t7;
  {
    t2 = tanh(0.1E2*x);
    t4 = 0.5*t2+0.5;
    t5 = t4*t4;
    t6 = t5*t5;
    t7 = t6*t4;
    t10 = t5*t4;
    t12 = 0.113810037E2-0.5075641578E2*t7+0.1268910394E3*t6-0.845940263E2*t10;
    t16 = -0.9E1*t7+0.225E2*t6-0.15E2*t10+3.0;
    t17 = 1/t16;
    t18 = t12*t17;
    t28 = 30.0*t6-60.0*t10+30.0*t5;
    t30 = log(t18);
    t31 = t18*t30;
    t39 = 1/t12;
    t41 = t2*t2;
    t48 = 5.0-5.0*t41;
    t49 = t6*t48;
    t51 = t10*t48;
    t53 = t5*t48;
    t57 = t16*t16;
    t67 = t12*t12;
    return((0.1E2*t18*(8.0*t10-12.0*t5+0.2E1*t2+0.2E1)+t28*(-0.25E1*t18+0.15E1*
t31+1.0)-t28*(-7.0*t18+3.0*t31+0.945940263E1)+0.1E1*t39*t16*t2*(0.1E2-0.1E2*t41
)+0.1*((-0.2537820789E3*t49+0.5075641576E3*t51-0.2537820789E3*t53)*t17-t12/t57*
(-0.45E2*t49+0.9E2*t51-0.45E2*t53))*t48/t67*t57)*t39*t16);
  }
}


inline double musol ( double x) const
{
  double t11;
  double t12;
  double t13;
  double t14;
  double t19;
  double t2;
  double t23;
  double t26;
  double t34;
  double t37;
  double t38;
  double t4;
  double t41;
  double t5;
  double t6;
  double t8;
  {
    t2 = tanh(0.1E2*x);
    t4 = 0.5*t2+0.5;
    t5 = t4*t4;
    t6 = t5*t5;
    t8 = t5*t4;
    t11 = t6*t4;
    t12 = 6.0*t11;
    t13 = 15.0*t6;
    t14 = 10.0*t8;
    t19 = 0.113810037E2-0.5075641578E2*t11+0.1268910394E3*t6-0.845940263E2*t8;
    t23 = -0.9E1*t11+0.225E2*t6-0.15E2*t8+3.0;
    t26 = log(t19/t23);
    t34 = t2*t2;
    t37 = pow(5.0-5.0*t34,2.0);
    t38 = t19*t19;
    t41 = t23*t23;
    return(0.2E2*t6-0.4E2*t8+0.2E2*t5+(t12-t13+t14)*(-1.0+0.15E1*t26)+(1.0-t12+
t13-t14)*(-4.0+3.0*t26)-0.5E-1*t37/t38*t41);
  }
}


inline double veloSource ( double x) const
{
  double t10;
  double t104;
  double t13;
  double t15;
  double t19;
  double t2;
  double t20;
  double t21;
  double t3;
  double t31;
  double t33;
  double t34;
  double t43;
  double t45;
  double t49;
  double t5;
  double t51;
  double t53;
  double t55;
  double t57;
  double t63;
  double t65;
  double t67;
  double t68;
  double t7;
  double t8;
  double t81;
  double t85;
  double t86;
  double t87;
  double t9;
  {
    t2 = tanh(0.1E2*x);
    t3 = t2*t2;
    t5 = 5.0-5.0*t3;
    t7 = 0.5*t2+0.5;
    t8 = t7*t7;
    t9 = t8*t8;
    t10 = t9*t7;
    t13 = t8*t7;
    t15 = 0.113810037E2-0.5075641578E2*t10+0.1268910394E3*t9-0.845940263E2*t13;
    t19 = -0.9E1*t10+0.225E2*t9-0.15E2*t13+3.0;
    t20 = 1/t19;
    t21 = t15*t20;
    t31 = 30.0*t9-60.0*t13+30.0*t8;
    t33 = log(t21);
    t34 = t21*t33;
    t43 = 1/t15*t19;
    t45 = 0.1E2-0.1E2*t3;
    t49 = t9*t5;
    t51 = t13*t5;
    t53 = t8*t5;
    t55 = -0.2537820789E3*t49+0.5075641576E3*t51-0.2537820789E3*t53;
    t57 = t19*t19;
    t63 = -0.45E2*t49+0.9E2*t51-0.45E2*t53;
    t65 = t55*t20-t15/t57*t63;
    t67 = t15*t15;
    t68 = 1/t67;
    t81 = 30.0*t49-60.0*t51+30.0*t53;
    t85 = 6.0*t10;
    t86 = 15.0*t9;
    t87 = 10.0*t13;
    t104 = t5*t5;
    return(-t5*(0.1E2*t21*(8.0*t13-12.0*t8+0.2E1*t2+0.2E1)+t31*(-0.25E1*t21+
0.15E1*t34+1.0)-t31*(-7.0*t21+3.0*t34+0.945940263E1)+0.1E1*t43*t2*t45+0.1*t65*
t5*t68*t57)+t21*(0.8E2*t51-0.12E3*t53+0.4E2*t7*t5+t81*(-1.0+0.15E1*t33)+0.15E1*
(t85-t86+t87)*t65*t43-t81*(-4.0+3.0*t33)+3.0*(1.0-t85+t86-t87)*t65*t43+0.1E1*t5
*t68*t57*t2*t45+0.1*t104/t67/t15*t57*t55-0.1*t104*t68*t19*t63));
  }
}


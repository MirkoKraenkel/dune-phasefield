
inline double rhosol ( double x) const
{
  double t11;
  double t3;
  double t5;
  double t6;
  double t7;
  double t8;
  {
    t3 = tanh(1/delta_*x);
    t5 = 0.5*t3+0.5;
    t6 = t5*t5;
    t7 = t6*t6;
    t8 = t5*t7;
    t11 = t6*t5;
    return((0.113810037E2-0.5075641578E2*t8+0.1268910394E3*t7-0.845940263E2*t11
)/(3.0-0.9E1*t8+0.225E2*t7-0.15E2*t11));
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
    t21 = t5*t7;
    t25 = 3.0-0.9E1*t21+0.225E2*t7-0.15E2*t13;
    t32 = t25*t25;
    return((-0.1268910394E3*t11+0.2537820788E3*t15-0.1268910394E3*t18)/t25-(
0.113810037E2-0.5075641578E2*t21+0.1268910394E3*t7-0.845940263E2*t13)/t32*(
-0.225E2*t11+0.45E2*t15-0.225E2*t18));
  }
}


inline double gradphi ( double x) const
{
  {
    return(0.5*(1.0-pow(tanh(1/delta_*x),2.0))/delta_);
  }
}


inline double thetasol ( double x) const
{
  double t1;
  double t11;
  double t13;
  double t14;
  double t18;
  double t19;
  double t21;
  double t22;
  double t3;
  double t33;
  double t34;
  double t36;
  double t37;
  double t4;
  double t48;
  double t49;
  double t5;
  double t54;
  double t57;
  double t6;
  double t60;
  double t64;
  double t7;
  double t74;
  double t8;
  {
    t1 = 1/delta_;
    t3 = tanh(t1*x);
    t4 = 0.5*t3;
    t5 = t4+0.5;
    t6 = t5*t5;
    t7 = t6*t6;
    t8 = t5*t7;
    t11 = t6*t5;
    t13 = 0.113810037E2-0.5075641578E2*t8+0.1268910394E3*t7-0.845940263E2*t11;
    t14 = t1*t13;
    t18 = 3.0-0.9E1*t8+0.225E2*t7-0.15E2*t11;
    t19 = 1/t18;
    t21 = 0.5-t4;
    t22 = t21*t21;
    t33 = 30.0*t7-60.0*t11+30.0*t6;
    t34 = t13*t19;
    t36 = log(t34);
    t37 = t34*t36;
    t48 = t3*t3;
    t49 = 1.0-t48;
    t54 = t7*t49*t1;
    t57 = t11*t49*t1;
    t60 = t6*t49*t1;
    t64 = t18*t18;
    t74 = t13*t13;
    return(4.0*t14*t19*t5*t22-4.0*t14*t19*t6*t21+t33*(-0.25E1*t34+0.15E1*t37+
1.0)-t33*(-7.0*t34+3.0*t37+0.945940263E1)+0.1E1*t1/t13*t18*t3*t49+0.5*((
-0.1268910394E3*t54+0.2537820788E3*t57-0.1268910394E3*t60)*t19-t13/t64*(
-0.225E2*t54+0.45E2*t57-0.225E2*t60))*t49/t74*t64);
  }
}


inline double phiSource ( double x) const
{
  double t1;
  double t11;
  double t13;
  double t14;
  double t18;
  double t19;
  double t21;
  double t22;
  double t3;
  double t33;
  double t34;
  double t36;
  double t37;
  double t4;
  double t45;
  double t48;
  double t49;
  double t5;
  double t54;
  double t57;
  double t6;
  double t60;
  double t64;
  double t7;
  double t74;
  double t8;
  {
    t1 = 1/delta_;
    t3 = tanh(t1*x);
    t4 = 0.5*t3;
    t5 = t4+0.5;
    t6 = t5*t5;
    t7 = t6*t6;
    t8 = t5*t7;
    t11 = t6*t5;
    t13 = 0.113810037E2-0.5075641578E2*t8+0.1268910394E3*t7-0.845940263E2*t11;
    t14 = t1*t13;
    t18 = 3.0-0.9E1*t8+0.225E2*t7-0.15E2*t11;
    t19 = 1/t18;
    t21 = 0.5-t4;
    t22 = t21*t21;
    t33 = 30.0*t7-60.0*t11+30.0*t6;
    t34 = t13*t19;
    t36 = log(t34);
    t37 = t34*t36;
    t45 = 1/t13;
    t48 = t3*t3;
    t49 = 1.0-t48;
    t54 = t7*t49*t1;
    t57 = t11*t49*t1;
    t60 = t6*t49*t1;
    t64 = t18*t18;
    t74 = t13*t13;
    return((4.0*t14*t19*t5*t22-4.0*t14*t19*t6*t21+t33*(-0.25E1*t34+0.15E1*t37+
1.0)-t33*(-7.0*t34+3.0*t37+0.945940263E1)+0.1E1*t1*t45*t18*t3*t49+0.5*((
-0.1268910394E3*t54+0.2537820788E3*t57-0.1268910394E3*t60)*t19-t13/t64*(
-0.225E2*t54+0.45E2*t57-0.225E2*t60))*t49/t74*t64)*t45*t18);
  }
}


inline double musol ( double x) const
{
  double t1;
  double t12;
  double t13;
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
  double t4;
  double t41;
  double t43;
  double t5;
  double t6;
  double t9;
  {
    t1 = 1/delta_;
    t3 = tanh(t1*x);
    t4 = 0.5*t3;
    t5 = t4+0.5;
    t6 = t5*t5;
    t9 = pow(0.5-t4,2.0);
    t12 = t6*t6;
    t13 = t12*t5;
    t14 = 6.0*t13;
    t15 = 15.0*t12;
    t16 = t6*t5;
    t17 = 10.0*t16;
    t22 = 0.113810037E2-0.5075641578E2*t13+0.1268910394E3*t12-0.845940263E2*t16
;
    t26 = -0.9E1*t13+0.225E2*t12-0.15E2*t16+3.0;
    t29 = log(t22/t26);
    t37 = t3*t3;
    t39 = pow(1.0-t37,2.0);
    t41 = t22*t22;
    t43 = t26*t26;
    return(2.0*t1*t6*t9+(t14-t15+t17)*(-1.0+0.15E1*t29)+(1.0-t14+t15-t17)*(-4.0
+3.0*t29)-0.125*t1*t39/t41*t43);
  }
}


inline double veloSource ( double x) const
{
  double t1;
  double t100;
  double t11;
  double t112;
  double t114;
  double t13;
  double t17;
  double t18;
  double t19;
  double t20;
  double t21;
  double t23;
  double t24;
  double t25;
  double t26;
  double t3;
  double t35;
  double t38;
  double t4;
  double t41;
  double t43;
  double t44;
  double t48;
  double t49;
  double t5;
  double t50;
  double t55;
  double t57;
  double t6;
  double t63;
  double t65;
  double t67;
  double t68;
  double t7;
  double t78;
  double t8;
  double t80;
  double t81;
  double t82;
  double t86;
  {
    t1 = 1/delta_;
    t3 = tanh(t1*x);
    t4 = 0.5*t3;
    t5 = t4+0.5;
    t6 = t5*t5;
    t7 = t6*t6;
    t8 = t7*t5;
    t11 = t6*t5;
    t13 = 0.113810037E2-0.5075641578E2*t8+0.1268910394E3*t7-0.845940263E2*t11;
    t17 = -0.9E1*t8+0.225E2*t7-0.15E2*t11+3.0;
    t18 = 1/t17;
    t19 = t13*t18;
    t20 = delta_*delta_;
    t21 = 1/t20;
    t23 = 0.5-t4;
    t24 = t23*t23;
    t25 = t3*t3;
    t26 = 1.0-t25;
    t35 = t7*t26*t1;
    t38 = t11*t26*t1;
    t41 = t6*t26*t1;
    t43 = 0.15E2*t35-0.3E2*t38+0.15E2*t41;
    t44 = log(t19);
    t48 = 6.0*t8;
    t49 = 15.0*t7;
    t50 = 10.0*t11;
    t55 = -0.1268910394E3*t35+0.2537820788E3*t38-0.1268910394E3*t41;
    t57 = t17*t17;
    t63 = -0.225E2*t35+0.45E2*t38-0.225E2*t41;
    t65 = t55*t18-t13/t57*t63;
    t67 = 1/t13;
    t68 = t67*t17;
    t78 = t26*t26;
    t80 = t13*t13;
    t81 = 1/t80;
    t82 = t81*t57;
    t86 = t1*t78;
    t100 = t1*t13;
    t112 = 30.0*t7-60.0*t11+30.0*t6;
    t114 = t19*t44;
    return(t19*(0.2E1*t21*t5*t24*t26-0.2E1*t21*t6*t23*t26+t43*(-1.0+0.15E1*t44)
+0.15E1*(t48-t49+t50)*t65*t68-t43*(-4.0+3.0*t44)+3.0*(1.0-t48+t49-t50)*t65*t68+
0.5*t21*t78*t82*t3+0.25*t86/t80/t13*t57*t55-0.25*t86*t81*t17*t63)-0.5*t26*t1*(
4.0*t100*t18*t5*t24-4.0*t100*t18*t6*t23+t112*(-0.25E1*t19+0.15E1*t114+1.0)-t112
*(-7.0*t19+3.0*t114+0.945940263E1)+0.1E1*t1*t67*t17*t3*t26+0.5*t65*t26*t82));
  }
}


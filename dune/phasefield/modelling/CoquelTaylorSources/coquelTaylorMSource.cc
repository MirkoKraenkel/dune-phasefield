
inline double rhsRho (double t, double x, double y ) const
{
  double t10;
  double t11;
  double t12;
  double t14;
  double t15;
  double t17;
  double t19;
  double t2;
  double t22;
  double t25;
  double t29;
  double t3;
  double t30;
  double t35;
  double t36;
  double t38;
  double t46;
  double t49;
  double t5;
  double t50;
  double t51;
  double t54;
  double t57;
  double t6;
  double t9;
  {
    t2 = 2.0*0.3141592653589793E1*t;
    t3 = cos(t2);
    t5 = 2.0*0.3141592653589793E1*x;
    t6 = cos(t5);
    t9 = 0.5*t3*t6+0.5;
    t10 = t9*t9;
    t11 = t10*t10;
    t12 = sin(t2);
    t14 = 0.3141592653589793E1*t6;
    t15 = t11*t12*t14;
    t17 = t10*t9;
    t19 = t17*t12*t14;
    t22 = t10*t12*t14;
    t25 = t11*t9;
    t29 = -0.3E1*t25+0.75E1*t11-0.5E1*t17+0.15E1;
    t30 = 1/t29;
    t35 = 0.4158883083359672E1*t25-0.1039720770839918E2*t11+
0.6931471805599453E1*t17;
    t36 = t29*t29;
    t38 = t35/t36;
    t46 = exp(t35*t30);
    t49 = sin(t5);
    t50 = t49*0.3141592653589793E1;
    t51 = t11*t3*t50;
    t54 = t17*t3*t50;
    t57 = t10*t3*t50;
    return(((-0.2079441541679836E2*t15+0.4158883083359672E2*t19
-0.2079441541679836E2*t22)*t30-t38*(0.15E2*t15-0.3E2*t19+0.15E2*t22))*t46+((
-0.2079441541679836E2*t51+0.4158883083359672E2*t54-0.2079441541679836E2*t57)*
t30-t38*(0.15E2*t51-0.3E2*t54+0.15E2*t57))*t46*t3*t49+2.0*t46*t3*t14);
  }
}


inline double rhsV1 (double t, double x, double y ) const
{
  double t10;
  double t105;
  double t106;
  double t11;
  double t117;
  double t12;
  double t15;
  double t17;
  double t2;
  double t21;
  double t22;
  double t24;
  double t25;
  double t27;
  double t28;
  double t3;
  double t32;
  double t35;
  double t38;
  double t41;
  double t42;
  double t49;
  double t5;
  double t50;
  double t51;
  double t52;
  double t53;
  double t56;
  double t59;
  double t6;
  double t68;
  double t69;
  double t7;
  double t72;
  double t73;
  double t74;
  double t9;
  {
    t2 = 2.0*0.3141592653589793E1*t;
    t3 = cos(t2);
    t5 = 2.0*0.3141592653589793E1*x;
    t6 = cos(t5);
    t7 = t3*t6;
    t9 = 0.5*t7+0.5;
    t10 = t9*t9;
    t11 = t10*t10;
    t12 = t11*t9;
    t15 = t10*t9;
    t17 = 0.4158883083359672E1*t12-0.1039720770839918E2*t11+
0.6931471805599453E1*t15;
    t21 = -0.3E1*t12+0.75E1*t11-0.5E1*t15+0.15E1;
    t22 = 1/t21;
    t24 = exp(t17*t22);
    t25 = sin(t2);
    t27 = sin(t5);
    t28 = 0.3141592653589793E1*t27;
    t32 = t11*t3*t28;
    t35 = t15*t3*t28;
    t38 = t10*t3*t28;
    t41 = (-0.2079441541679836E2*t32+0.4158883083359672E2*t35
-0.2079441541679836E2*t38)*t22;
    t42 = t21*t21;
    t49 = t17/t42*(0.15E2*t32-0.3E2*t35+0.15E2*t38);
    t50 = t41-t49;
    t51 = t50*t24;
    t52 = t3*t3;
    t53 = t27*t27;
    t56 = t3*t27;
    t59 = 0.3141592653589793E1*t6;
    t68 = -0.3E2*t32+0.6E2*t35-0.3E2*t38;
    t69 = log(t24);
    t72 = 6.0*t12;
    t73 = 15.0*t11;
    t74 = 10.0*t15;
    t105 = 30.0*t11-60.0*t15+30.0*t10;
    t106 = t24*t69;
    t117 = 0.3141592653589793E1*0.3141592653589793E1;
    return(-2.0*t24*t25*t28+t51*t52*t53-(t51*t56+2.0*t24*t3*t59)*t3*t27+t24*(
0.1*t68*t69+0.1*(t72-t73+t74)*t50-0.1*t68*(0.6931471805599453+0.15E1*t69)+0.1*(
1.0-t72+t73-t74)*(0.15E1*t41-0.15E1*t49)+0.4E1*t52*t27*t59)+0.1E1*t56*
0.3141592653589793E1*(0.2*A_*(4.0*t15-6.0*t10+0.1E1*t7+0.1E1)/delta_+0.1*t105*(-
t24+t106+0.5)-0.1*t105*(-0.8068528194400547*t24+0.15E1*t106)+0.2E1*delta_*A_*t7*
t117)+4.0*mu1Liq_*t3*t27*t117);
  }
}


inline double rhsV2 (double t, double x, double y ) const
{
  double t10;
  double t11;
  double t12;
  double t15;
  double t17;
  double t2;
  double t21;
  double t22;
  double t24;
  double t25;
  double t27;
  double t28;
  double t3;
  double t32;
  double t35;
  double t38;
  double t42;
  double t5;
  double t51;
  double t52;
  double t53;
  double t6;
  double t71;
  double t9;
  {
    t2 = 2.0*0.3141592653589793E1*t;
    t3 = cos(t2);
    t5 = 2.0*0.3141592653589793E1*x;
    t6 = cos(t5);
    t9 = 0.5*t3*t6+0.5;
    t10 = t9*t9;
    t11 = t10*t10;
    t12 = t11*t9;
    t15 = t10*t9;
    t17 = 0.4158883083359672E1*t12-0.1039720770839918E2*t11+
0.6931471805599453E1*t15;
    t21 = -0.3E1*t12+0.75E1*t11-0.5E1*t15+0.15E1;
    t22 = 1/t21;
    t24 = exp(t17*t22);
    t25 = sin(t2);
    t27 = sin(t5);
    t28 = 0.3141592653589793E1*t27;
    t32 = t11*t3*t28;
    t35 = t15*t3*t28;
    t38 = t10*t3*t28;
    t42 = t21*t21;
    t51 = ((-0.2079441541679836E2*t32+0.4158883083359672E2*t35
-0.2079441541679836E2*t38)*t22-t17/t42*(0.15E2*t32-0.3E2*t35+0.15E2*t38))*t24;
    t52 = t3*t3;
    t53 = t27*t27;
    t71 = 0.3141592653589793E1*0.3141592653589793E1;
    return(-2.0*t24*t25*t28+t51*t52*t53+4.0*t24*t52*t27*t6*0.3141592653589793E1
-(t51*t3*t27+2.0*t24*t3*t6*0.3141592653589793E1)*t3*t27+4.0*mu1Liq_*t3*t27*t71);
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
  double t34;
  double t35;
  double t46;
  double t47;
  double t48;
  double t59;
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
    t34 = 30.0*t30-60.0*t21+30.0*t20;
    t35 = t30*t19;
    t46 = exp((0.4158883083359672E1*t35-0.1039720770839918E2*t30+
0.6931471805599453E1*t21)/(-0.3E1*t35+0.75E1*t30-0.5E1*t21+0.15E1));
    t47 = log(t46);
    t48 = t46*t47;
    t59 = 0.3141592653589793E1*0.3141592653589793E1;
    return(-0.1E1*t3*0.3141592653589793E1*t7-0.1E1*t11*t13*0.3141592653589793E1
+(0.2*A_*(4.0*t21-6.0*t20+0.1E1*t17+0.1E1)/delta_+0.1*t34*(-t46+t48+0.5)-0.1*t34*
(-0.8068528194400547*t46+0.15E1*t48)+0.2E1*delta_*A_*t17*t59)/t46);
  }
}


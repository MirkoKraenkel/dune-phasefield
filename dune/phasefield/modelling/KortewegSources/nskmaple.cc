
inline double helmholtz ( double rho ) const
{
  double t4;
  {
    t4 = log(rho/(1.0-rho));
    return(rho*(-rho+0.2518518519*t4));
  }
}


inline double chemicalPotential ( double rho ,double old ) const
{
  double t1;
  double t10;
  double t12;
  double t13;
  double t14;
  double t2;
  double t20;
  double t23;
  double t24;
  double t27;
  double t28;
  double t29;
  double t3;
  double t33;
  double t4;
  double t5;
  double t7;
  double t9;
  {
    t1 = 0.5*rho;
    t2 = 0.5*old;
    t3 = t1+t2;
    t4 = 1.0-t1-t2;
    t5 = 1/t4;
    t7 = log(t3*t5);
    t9 = t4*t4;
    t10 = 1/t9;
    t12 = t5+t3*t10;
    t13 = 1/t3;
    t14 = t12*t13;
    t20 = 1/t9/t4;
    t23 = 2.0*t10+2.0*t3*t20;
    t24 = t23*t13;
    t27 = t3*t3;
    t28 = 1/t27;
    t29 = t12*t28;
    t33 = t9*t9;
    return(-t1-t2+0.2518518519*t7+t3*(-1.0+0.2518518519*t14*t4)+(0.7555555557*
t24*t4-0.7555555557*t29*t4-0.7555555557*t14+t3*(0.2518518519*(6.0*t20+6.0*t3/
t33)*t13*t4-0.5037037038*t23*t28*t4-0.5037037038*t24+0.5037037038*t12/t27/t3*t4
+0.5037037038*t29))*(rho-old)/24.0);
  }
}

inline double drhochemicalPotential ( double rho ,double old ) const
{
  double t1;
  double t12;
  double t16;
  double t17;
  double t2;
  double t22;
  double t23;
  double t25;
  double t26;
  double t29;
  double t3;
  double t30;
  double t31;
  double t32;
  double t38;
  double t39;
  double t4;
  double t40;
  double t42;
  double t43;
  double t47;
  double t48;
  double t49;
  double t51;
  double t53;
  double t57;
  double t58;
  double t59;
  double t6;
  double t64;
  double t65;
  double t66;
  double t7;
  double t8;
  double t9;
  double t92;
  double t99;
  {
    t1 = 0.5*rho;
    t2 = 0.5*old;
    t3 = 1.0-t1-t2;
    t4 = 1/t3;
    t6 = t1+t2;
    t7 = t3*t3;
    t8 = 1/t7;
    t9 = t6*t8;
    t12 = 1/t6;
    t16 = t4+t9;
    t17 = t16*t12;
    t22 = 1/t7/t3;
    t23 = t6*t22;
    t25 = 0.1E1*t8+0.1E1*t23;
    t26 = t25*t12;
    t29 = t6*t6;
    t30 = 1/t29;
    t31 = t16*t30;
    t32 = t31*t3;
    t38 = t7*t7;
    t39 = 1/t38;
    t40 = t6*t39;
    t42 = 0.3E1*t22+0.3E1*t40;
    t43 = t42*t12;
    t47 = 2.0*t8+2.0*t23;
    t48 = t47*t30;
    t49 = t48*t3;
    t51 = t47*t12;
    t53 = t25*t30;
    t57 = 1/t29/t6;
    t58 = t16*t57;
    t59 = t58*t3;
    t64 = 6.0*t22+6.0*t40;
    t65 = t64*t12;
    t66 = t65*t3;
    t92 = t29*t29;
    t99 = 0.2518518519*(0.12E2*t39+0.12E2*t6/t38/t3)*t12*t3-0.125925926*t64*t30
*t3-0.125925926*t65-0.5037037038*t42*t30*t3+0.5037037038*t47*t57*t3+
0.5037037038*t48-0.5037037038*t43+0.5037037038*t25*t57*t3-0.7555555557*t16/t92*
t3-0.7555555557*t58+0.5037037038*t53;
    return(-0.1E1+0.2518518519*(0.5*t4+0.5*t9)*t12*t3+0.125925926*t17*t3+t6*(
0.2518518519*t26*t3-0.125925926*t32-0.125925926*t17)+(0.7555555557*t43*t3
-0.6296296297*t49-0.6296296297*t51-0.7555555557*t53*t3+0.1007407408E1*t59+
0.1007407408E1*t31-0.7555555557*t26+0.125925926*t66+t6*t99)*(rho-old)/24.0+
0.3148148149E-1*t51*t3-0.3148148149E-1*t32-0.3148148149E-1*t17+t6*(0.2518518519
*t66-0.5037037038*t49-0.5037037038*t51+0.5037037038*t59+0.5037037038*t31)/24.0
);
  }
}


inline double pressure ( double rho ) const
{
  double t1;
  double t2;
  double t4;
  double t5;
  double t8;
  {
    t1 = 1.0-rho;
    t2 = 1/t1;
    t4 = log(rho*t2);
    t5 = 0.2518518519*t4;
    t8 = t1*t1;
    return(-rho*(-rho+t5)+rho*(-rho+t5+rho*(-1.0+0.2518518519*(t2+rho/t8)/rho*
t1)));
  }
}


inline double a ( double rho ) const
{
  double t1;
  double t12;
  double t8;
  {
    t1 = rho*rho;
    t8 = pow(rho-1.0,2.0);
    t12 = (0.0>-0.1E-9*(0.2E11*t1*rho-0.4E11*t1+0.2E11*rho-2518518519.0)/t8 ? 
0.0 : -0.1E-9*(0.2E11*t1*rho-0.4E11*t1+0.2E11*rho-2518518519.0)/t8);
    return(sqrt(t12));
  }
}


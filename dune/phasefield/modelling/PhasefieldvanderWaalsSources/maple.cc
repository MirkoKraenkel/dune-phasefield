
inline double helmholtz ( double rho ,double phi ) const
{
  double t1;
  double t12;
  double t13;
  double t14;
  double t2;
  double t20;
  double t21;
  double t27;
  double t4;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = t1*phi;
    t12 = 6.0*t2*phi;
    t13 = 15.0*t2;
    t14 = 10.0*t4;
    t20 = log(rho);
    t21 = rho*t20;
    t27 = rho*rho;
    return(A_*(2.0*t2-4.0*t4+2.0*t1)/delta_+(t12-t13+t14)*(-0.1351941256E-1*rho/(
rho-1.0)-0.9990832252E-2*t21-0.8058411491E-1*rho+0.2501037014E-1)+(1.0-t12+t13-
t14)*(-0.7127752726E-1*t27+0.7074554549E-1*t21+0.1028390454*rho+0.6730211801E-2
));
  }
}


inline double reactionSource ( double rho ,double phi ,double old ) const
{
  double t12;
  double t14;
  double t18;
  double t23;
  double t24;
  double t27;
  double t3;
  double t30;
  double t34;
  double t4;
  double t44;
  double t5;
  {
    t3 = 0.5*phi+0.5*old;
    t4 = t3*t3;
    t5 = t4*t3;
    t12 = 1/delta_;
    t14 = t4*t4;
    t18 = 30.0*t14-60.0*t5+30.0*t4;
    t23 = log(rho);
    t24 = rho*t23;
    t27 = -0.1351941256E-1*rho/(rho-1.0)-0.9990832252E-2*t24-0.8058411491E-1*
rho+0.2501037014E-1;
    t30 = rho*rho;
    t34 = -0.7127752726E-1*t30+0.7074554549E-1*t24+0.1028390454*rho+
0.6730211801E-2;
    t44 = 360.0*t4-0.18E3*phi-0.18E3*old+60.0;
    return(A_*(8.0*t5-12.0*t4+0.2E1*phi+0.2E1*old)*t12+t18*t27-t18*t34+(A_*(
0.24E2*phi+0.24E2*old-24.0)*t12+t44*t27-t44*t34)*(phi-old)/24.0);
  }
}


inline double dphireactionSource ( double rho ,double phi ,double old ) const
{
  double t10;
  double t17;
  double t22;
  double t23;
  double t26;
  double t29;
  double t3;
  double t33;
  double t39;
  double t4;
  double t56;
  {
    t3 = 0.5*phi+0.5*old;
    t4 = t3*t3;
    t10 = 1/delta_;
    t17 = 0.6E2*t4*t3-0.9E2*t4+0.15E2*phi+0.15E2*old;
    t22 = log(rho);
    t23 = rho*t22;
    t26 = -0.1351941256E-1*rho/(rho-1.0)-0.9990832252E-2*t23-0.8058411491E-1*
rho+0.2501037014E-1;
    t29 = rho*rho;
    t33 = -0.7127752726E-1*t29+0.7074554549E-1*t23+0.1028390454*rho+
0.6730211801E-2;
    t39 = 0.18E3*phi+0.18E3*old-0.18E3;
    t56 = 360.0*t4-0.18E3*phi-0.18E3*old+60.0;
    return(A_*(0.12E2*t4-0.6E1*phi-0.6E1*old+0.2E1)*t10+t17*t26-t17*t33+(0.24E2*
A_*t10+t39*t26-t39*t33)*(phi-old)/24.0+A_*(0.24E2*phi+0.24E2*old-24.0)*t10/24.0+
t56*t26/24.0-t56*t33/24.0);
  }
}


inline double chemicalPotential ( double rho ,double phi ,double old ) const
{
  double t1;
  double t10;
  double t11;
  double t14;
  double t15;
  double t19;
  double t2;
  double t23;
  double t32;
  double t36;
  double t37;
  double t4;
  double t5;
  double t7;
  double t8;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t7 = 10.0*t1*phi;
    t8 = t4-t5+t7;
    t9 = 0.5*rho;
    t10 = 0.5*old;
    t11 = t9+t10-1.0;
    t14 = t9+t10;
    t15 = t11*t11;
    t19 = log(t14);
    t23 = 1.0-t4+t5-t7;
    t32 = t15*t15;
    t36 = t14*t14;
    t37 = 1/t36;
    return(t8*(-0.1351941256E-1/t11+0.1351941256E-1*t14/t15-0.9990832252E-2*t19
-0.9057494716E-1)+t23*(-0.7127752725E-1*rho-0.7127752725E-1*old+0.7074554549E-1
*t19+0.1735845909)+(t8*(-0.8111647536E-1/t15/t11+0.8111647536E-1*t14/t32+
0.9990832252E-2*t37)-0.7074554549E-1*t23*t37)*(rho-old)/24.0);
  }
}

inline double drhochemicalPotential ( double rho ,double phi ,double old ) const
{
  double t1;
  double t10;
  double t11;
  double t12;
  double t15;
  double t17;
  double t2;
  double t20;
  double t24;
  double t28;
  double t29;
  double t35;
  double t37;
  double t4;
  double t5;
  double t50;
  double t7;
  double t8;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t7 = 10.0*t1*phi;
    t8 = t4-t5+t7;
    t9 = 0.5*rho;
    t10 = 0.5*old;
    t11 = t9+t10-1.0;
    t12 = t11*t11;
    t15 = t9+t10;
    t17 = 1/t12/t11;
    t20 = 1/t15;
    t24 = 1.0-t4+t5-t7;
    t28 = t12*t12;
    t29 = 1/t28;
    t35 = t15*t15;
    t37 = 1/t35/t15;
    t50 = 1/t35;
    return(t8*(0.1351941256E-1/t12-0.1351941256E-1*t15*t17-0.4995416126E-2*t20)
+t24*(-0.7127752725E-1+0.3537277274E-1*t20)+(t8*(0.1622329507*t29-0.1622329507*
t15/t28/t11-0.9990832252E-2*t37)+0.7074554549E-1*t24*t37)*(rho-old)/24.0+t8*(
-0.8111647536E-1*t17+0.8111647536E-1*t15*t29+0.9990832252E-2*t50)/24.0
-0.2947731062E-2*t24*t50);
  }
}


inline double dphichemicalPotential ( double rho ,double phi ,double old ) const
{
  double t1;
  double t10;
  double t13;
  double t14;
  double t18;
  double t2;
  double t22;
  double t31;
  double t35;
  double t36;
  double t7;
  double t8;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t7 = 30.0*t2-60.0*t1*phi+30.0*t1;
    t8 = 0.5*rho;
    t9 = 0.5*old;
    t10 = t8+t9-1.0;
    t13 = t8+t9;
    t14 = t10*t10;
    t18 = log(t13);
    t22 = -t7;
    t31 = t14*t14;
    t35 = t13*t13;
    t36 = 1/t35;
    return(t7*(-0.1351941256E-1/t10+0.1351941256E-1*t13/t14-0.9990832252E-2*t18
-0.9057494716E-1)+t22*(-0.7127752725E-1*rho-0.7127752725E-1*old+0.7074554549E-1
*t18+0.1735845909)+(t7*(-0.8111647536E-1/t14/t10+0.8111647536E-1*t13/t31+
0.9990832252E-2*t36)-0.7074554549E-1*t22*t36)*(rho-old)/24.0);
  }
}


inline double pressure ( double rho ,double phi ) const
{
  double t1;
  double t10;
  double t13;
  double t14;
  double t18;
  double t2;
  double t28;
  double t4;
  double t5;
  double t6;
  double t7;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t6 = t1*phi;
    t7 = 10.0*t6;
    t9 = rho-1.0;
    t10 = 1/t9;
    t13 = log(rho);
    t14 = rho*t13;
    t18 = t9*t9;
    t28 = rho*rho;
    return((t4-t5+t7)*(0.1351941256E-1*rho*t10+0.9990832252E-2*t14+
0.8058411491E-1*rho-0.2501037014E-1+rho*(-0.1351941256E-1*t10+0.1351941256E-1*
rho/t18-0.9990832252E-2*t13-0.9057494716E-1))+(1.0-t4+t5-t7)*(0.7127752726E-1*
t28-0.7074554549E-1*t14-0.1028390454*rho-0.6730211801E-2+rho*(-0.1425550545*rho
+0.7074554549E-1*t13+0.1735845909))-(2.0*t2-4.0*t6+2.0*t1)/delta_);
  }
}


inline double a ( double rho ,double phi ) const
{
  double t11;
  double t13;
  double t2;
  double t39;
  double t4;
  double t40;
  double t41;
  double t45;
  double t6;
  double t7;
  double t9;
  {
    t2 = rho*rho;
    t4 = t2*t2;
    t6 = phi*phi;
    t7 = t6*t6;
    t9 = t6*phi;
    t11 = t7*phi;
    t13 = t2*rho;
    t39 = 0.1773958455E11*rho-0.3199509E11*t2-7127752725.0*t4-0.605522833E11*t7
-3537277274.0+0.4036818887E11*t9+0.2422091332E11*t11+0.2492053545E11*t13+
0.200962289E12*t11*t2-0.1235409039E12*t11*rho-0.1525204624E12*t11*t13
-0.5024057225E12*t7*t2+0.3088522596E12*t7*rho+0.3813011559E12*t7*t13+
0.3349371484E12*t9*t2-0.2059015064E12*t9*rho-0.2542007706E12*t9*t13+
0.4276651635E11*t11*t4-0.1069162909E12*t7*t4+0.7127752725E11*t9*t4;
    t40 = rho-1.0;
    t41 = t40*t40;
    t45 = sqrt(t39/t41/t40);
    return(0.4472135955E-5*t45);
  }
}


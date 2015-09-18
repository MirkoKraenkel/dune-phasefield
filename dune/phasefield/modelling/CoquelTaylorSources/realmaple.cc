
inline double helmholtz ( double rho ,double phi ) const
{
  double t1;
  double t11;
  double t12;
  double t13;
  double t16;
  double t17;
  double t2;
  double t3;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = t1*phi;
    t11 = 6.0*t2*phi;
    t12 = 15.0*t2;
    t13 = 10.0*t3;
    t16 = log(rho);
    t17 = rho*t16;
    return(2.0*A_*(t2-2.0*t3+t1)/delta_+(t11-t12+t13)*(-0.2011952932E2*rho+
0.289445110522E2*t17+0.2133795654E2)+(1.0-t11+t12-t13)*(0.193942126E2*rho+
0.950205321859E1*t17+0.4540504209));
  }
}


inline double reactionSource ( double rho ,double phi ,double old ) const
{
  double t12;
  double t15;
  double t19;
  double t21;
  double t22;
  double t24;
  double t29;
  double t3;
  double t4;
  double t40;
  double t5;
  {
    t3 = 0.5*phi+0.5*old;
    t4 = t3*t3;
    t5 = t4*t3;
    t12 = 1/delta_;
    t15 = t4*t4;
    t19 = 30.0*t15-60.0*t5+30.0*t4;
    t21 = log(rho);
    t22 = rho*t21;
    t24 = -0.2011952932E2*rho+0.289445110522E2*t22+0.2133795654E2;
    t29 = 0.193942126E2*rho+0.950205321859E1*t22+0.4540504209;
    t40 = 360.0*t4-0.18E3*phi-0.18E3*old+60.0;
    return(2.0*A_*(4.0*t5-6.0*t4+0.1E1*phi+0.1E1*old)*t12+t19*t24-t19*t29+(2.0*A_
*(0.12E2*phi+0.12E2*old-12.0)*t12+t40*t24-t40*t29)*(phi-old)/24.0);
  }
}


inline double dphireactionSource ( double rho ,double phi ,double old ) const
{
  double t10;
  double t18;
  double t20;
  double t21;
  double t23;
  double t28;
  double t3;
  double t34;
  double t4;
  double t51;
  {
    t3 = 0.5*phi+0.5*old;
    t4 = t3*t3;
    t10 = 1/delta_;
    t18 = 0.6E2*t4*t3-0.9E2*t4+0.15E2*phi+0.15E2*old;
    t20 = log(rho);
    t21 = rho*t20;
    t23 = -0.2011952932E2*rho+0.289445110522E2*t21+0.2133795654E2;
    t28 = 0.193942126E2*rho+0.950205321859E1*t21+0.4540504209;
    t34 = 0.18E3*phi+0.18E3*old-0.18E3;
    t51 = 360.0*t4-0.18E3*phi-0.18E3*old+60.0;
    return(2.0*A_*(0.6E1*t4-0.3E1*phi-0.3E1*old+0.1E1)*t10+t18*t23-t18*t28+(
0.24E2*A_*t10+t34*t23-t34*t28)*(phi-old)/24.0+A_*(0.12E2*phi+0.12E2*old-12.0)*t10
/12.0+t51*t23/24.0-t51*t28/24.0);
  }
}


inline double drhoreactionSource ( double rho ,double phi ,double old ) const
{
  double t10;
  double t11;
  double t13;
  double t17;
  double t22;
  double t3;
  double t4;
  double t5;
  {
    t3 = 0.5*phi+0.5*old;
    t4 = t3*t3;
    t5 = t4*t4;
    t10 = 30.0*t5-60.0*t4*t3+30.0*t4;
    t11 = log(rho);
    t13 = 0.882498173E1+0.289445110522E2*t11;
    t17 = 0.2889626582E2+0.950205321859E1*t11;
    t22 = 360.0*t4-0.18E3*phi-0.18E3*old+60.0;
    return(t10*t13-t10*t17+(t22*t13-t22*t17)*(phi-old)/24.0);
  }
}


inline double chemicalPotential ( double rho ,double phi ,double old ) const
{
  double t1;
  double t11;
  double t12;
  double t16;
  double t2;
  double t20;
  double t21;
  double t4;
  double t5;
  double t7;
  double t8;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t7 = 10.0*t1*phi;
    t8 = t4-t5+t7;
    t11 = 0.5*rho+0.5*old;
    t12 = log(t11);
    t16 = 1.0-t4+t5-t7;
    t20 = t11*t11;
    t21 = 1/t20;
    return(t8*(0.882498173E1+0.289445110522E2*t12)+t16*(0.2889626582E2+
0.950205321859E1*t12)+(-0.2894451105E2*t8*t21-0.9502053219E1*t16*t21)*(rho-old)
/24.0);
  }
}

inline double drhochemicalPotential ( double rho ,double phi ,double old ) const
{
  double t1;
  double t11;
  double t12;
  double t15;
  double t18;
  double t2;
  double t20;
  double t29;
  double t4;
  double t5;
  double t7;
  double t8;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t7 = 10.0*t1*phi;
    t8 = t4-t5+t7;
    t11 = 0.5*rho+0.5*old;
    t12 = 1/t11;
    t15 = 1.0-t4+t5-t7;
    t18 = t11*t11;
    t20 = 1/t18/t11;
    t29 = 1/t18;
    return(0.1447225552E2*t8*t12+0.475102661E1*t15*t12+(0.2894451105E2*t8*t20+
0.9502053219E1*t15*t20)*(rho-old)/24.0-0.1206021294E1*t8*t29-0.3959188841*t15*
t29);
  }
}


inline double dphichemicalPotential ( double rho ,double phi ,double old ) const
{
  double t1;
  double t10;
  double t11;
  double t15;
  double t19;
  double t2;
  double t20;
  double t7;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t7 = 30.0*t2-60.0*t1*phi+30.0*t1;
    t10 = 0.5*rho+0.5*old;
    t11 = log(t10);
    t15 = -t7;
    t19 = t10*t10;
    t20 = 1/t19;
    return(t7*(0.882498173E1+0.289445110522E2*t11)+t15*(0.2889626582E2+
0.950205321859E1*t11)+(-0.2894451105E2*t7*t20-0.9502053219E1*t15*t20)*(rho-old)
/24.0);
  }
}


inline double pressure ( double rho ,double phi ) const
{
  double t1;
  double t10;
  double t11;
  double t2;
  double t4;
  double t5;
  double t6;
  double t7;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t6 = t1*phi;
    t7 = 10.0*t6;
    t10 = log(rho);
    t11 = rho*t10;
    return((t4-t5+t7)*(-0.2094881058E2+0.4291934897E2*rho-0.2894451105E2*t11+
rho*(-0.1397483792E2+0.289445110522E2*t10))+(1.0-t4+t5-t7)*(-0.6490446091E-1+
0.3405607049E1*rho-0.9502053219E1*t11+rho*(0.609644617E1+0.950205321859E1*t10))
-2.0*A_*(t2-2.0*t6+t1)/delta_);
  }
}


inline double a ( double rho ,double phi ) const
{
  double t1;
  double t2;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t9 = sqrt(0.116654747E11*t2*phi-0.2916368675E11*t2+0.1944245783E11*t1*phi+
950205322.0);
    return(0.1E-3*t9);
  }
}


double evalRho(double x)
{
  {
    return(0.18708E1*x+0.1292);
  }
}


inline double helmholtz ( double rho ,double phi ) const
{
  double t12;
  double t13;
  double t14;
  double t18;
  double t2;
  double t26;
  double t3;
  double t4;
  {
    t2 = phi*phi;
    t3 = t2*t2;
    t4 = t2*phi;
    t12 = 6.0*t3*phi;
    t13 = 15.0*t3;
    t14 = 10.0*t4;
    t18 = log(0.773993808E1*rho);
    t26 = log(rho/2.0);
    return(2.0*rho*A_*(t3-2.0*t4+t2)/delta_+(t12-t13+t14)*beta_*(10.0*rho*(t18-1.0
)+0.1292E1)+(1.0-t12+t13-t14)*(0.1*rho*(t26-1.0)+0.2));
  }
}


inline double reactionSource ( double rho ,double phi ) const
{
  double t1;
  double t12;
  double t16;
  double t20;
  double t23;
  double t27;
  double t31;
  double t35;
  double t4;
  double t46;
  double t5;
  double t6;
  {
    t1 = rho*A_;
    t4 = 0.5*phi+0.5*old;
    t5 = t4*t4;
    t6 = t5*t4;
    t12 = 1/delta_;
    t16 = t5*t5;
    t20 = 30.0*t16-60.0*t6+30.0*t5;
    t23 = log(0.773993808E1*rho);
    t27 = 10.0*rho*(t23-1.0)+0.1292E1;
    t31 = log(rho/2.0);
    t35 = 0.1*rho*(t31-1.0)+0.2;
    t46 = 360.0*t5-0.18E3*phi-0.18E3*old+60.0;
    return(2.0*t1*(4.0*t6-6.0*t5+0.1E1*phi+0.1E1*old)*t12+t20*beta_*t27-t20*t35+
(2.0*t1*(0.12E2*phi+0.12E2*old-12.0)*t12+t46*beta_*t27-t46*t35)*(phi-old)/24.0);
  }
}


inline double dphireactionSource ( double rho ,double phi ) const
{
  double t1;
  double t10;
  double t19;
  double t22;
  double t26;
  double t30;
  double t34;
  double t4;
  double t40;
  double t5;
  double t58;
  {
    t1 = rho*A_;
    t4 = 0.5*phi+0.5*old;
    t5 = t4*t4;
    t10 = 1/delta_;
    t19 = 0.6E2*t5*t4-0.9E2*t5+0.15E2*phi+0.15E2*old;
    t22 = log(0.773993808E1*rho);
    t26 = 10.0*rho*(t22-1.0)+0.1292E1;
    t30 = log(rho/2.0);
    t34 = 0.1*rho*(t30-1.0)+0.2;
    t40 = 0.18E3*phi+0.18E3*old-0.18E3;
    t58 = 360.0*t5-0.18E3*phi-0.18E3*old+60.0;
    return(2.0*t1*(0.6E1*t5-0.3E1*phi-0.3E1*old+0.1E1)*t10+t19*beta_*t26-t19*t34
+(0.24E2*t1*t10+t40*beta_*t26-t40*t34)*(phi-old)/24.0+t1*(0.12E2*phi+0.12E2*old
-12.0)*t10/12.0+t58*beta_*t26/24.0-t58*t34/24.0);
  }
}


inline double chemicalPotential ( double rho ,double phi ) const
{
  double t1;
  double t11;
  double t12;
  double t13;
  double t15;
  double t19;
  double t2;
  double t22;
  double t26;
  double t3;
  double t32;
  double t33;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = t1*phi;
    t11 = 6.0*t2*phi;
    t12 = 15.0*t2;
    t13 = 10.0*t3;
    t15 = (t11-t12+t13)*beta_;
    t19 = log(0.386996904E1*rho+0.386996904E1*old);
    t22 = 1.0-t11+t12-t13;
    t26 = log(0.25*rho+0.25*old);
    t32 = pow(0.5*rho+0.5*old,2.0);
    t33 = 1/t32;
    return(2.0*A_*(t2-2.0*t3+t1)/delta_+10.0*t15*t19+0.1*t22*t26+(-0.1E2*t15*t33
-0.1*t22*t33)*(rho-old)/24.0);
  }
}


inline double drhochemicalPotential ( double rho ,double phi ) const
{
  double t1;
  double t16;
  double t2;
  double t25;
  double t26;
  double t28;
  double t37;
  double t4;
  double t5;
  double t7;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t7 = 10.0*t1*phi;
    t9 = (t4-t5+t7)*beta_;
    t16 = 1.0-t4+t5-t7;
    t25 = 0.5*rho+0.5*old;
    t26 = t25*t25;
    t28 = 1/t26/t25;
    t37 = 1/t26;
    return(0.386996904E2*t9/(0.386996904E1*rho+0.386996904E1*old)+0.25E-1*t16/(
0.25*rho+0.25*old)+(0.1E2*t9*t28+0.1*t16*t28)*(rho-old)/24.0-0.4166666667*t9*
t37-0.4166666667E-2*t16*t37);
  }
}


inline double dphichemicalPotential ( double rho ,double phi ) const
{
  double t1;
  double t11;
  double t15;
  double t16;
  double t2;
  double t20;
  double t23;
  double t27;
  double t33;
  double t34;
  {
    t1 = phi*phi;
    t2 = t1*phi;
    t11 = t1*t1;
    t15 = 30.0*t11-60.0*t2+30.0*t1;
    t16 = t15*beta_;
    t20 = log(0.386996904E1*rho+0.386996904E1*old);
    t23 = -t15;
    t27 = log(0.25*rho+0.25*old);
    t33 = pow(0.5*rho+0.5*old,2.0);
    t34 = 1/t33;
    return(2.0*A_*(4.0*t2-6.0*t1+2.0*phi)/delta_+10.0*t16*t20+0.1*t23*t27+(-0.1E2
*t16*t34-0.1*t23*t34)*(rho-old)/24.0);
  }
}


inline double pressure ( double rho ,double phi ) const
{
  double t1;
  double t10;
  double t2;
  double t23;
  double t4;
  double t5;
  double t7;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t7 = 10.0*t1*phi;
    t10 = log(0.773993808E1*rho);
    t23 = log(rho/2.0);
    return(beta_*((t4-t5+t7)*(-beta_*(10.0*rho*(t10-1.0)+0.1292E1)+10.0*rho*beta_*
t10)+(1.0-t4+t5-t7)*(-0.1*rho*(t23-1.0)-0.2+0.1*rho*t23)));
  }
}


inline double a ( double rho ,double phi ) const
{
  double t1;
  double t16;
  double t2;
  double t3;
  double t8;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = t2*phi;
    t8 = t1*phi;
    t16 = sqrt(beta_*(600.0*beta_*t3-1500.0*beta_*t2+1000.0*beta_*t8+1.0-6.0*t3+
15.0*t2-10.0*t8));
    return(0.316227766*t16);
  }
}


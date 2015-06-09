double evalRho(double x)
{
  {
    return(0.4958034541*x+0.1065766549);
  }
}


inline double helmholtz ( double rho ,double phi ) const
{
  double t1;
  double t16;
  double t17;
  double t18;
  double t2;
  double t20;
  double t22;
  double t29;
  double t5;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t5 = t1*phi;
    t16 = 6.0*t2*phi;
    t17 = 15.0*t2;
    t18 = 10.0*t5;
    t20 = rho*rho;
    t22 = t20*t20;
    t29 = log(rho);
    return(2.0*A_*(t2+(2.0*alpha_-2.0)*t5+(-3.0*alpha_+1.0)*t1+alpha_)/delta_+beta_*(
(t16-t17+t18)*(0.8575855211E-1*t22*t20*rho-0.2868132559E-1*rho+0.1480890852E-1)
+(1.0-t16+t17-t18)*(0.3001504024E-1*rho*t29+0.3718535686E-1*rho+0.3198902536E-2
)));
  }
}


inline double reactionSource ( double rho ,double phi ,double old ) const
{
  double t17;
  double t20;
  double t24;
  double t25;
  double t27;
  double t3;
  double t31;
  double t34;
  double t38;
  double t4;
  double t5;
  double t52;
  {
    t3 = 0.5*phi+0.5*old;
    t4 = t3*t3;
    t5 = t4*t3;
    t17 = 1/delta_;
    t20 = t4*t4;
    t24 = 30.0*t20-60.0*t5+30.0*t4;
    t25 = rho*rho;
    t27 = t25*t25;
    t31 = 0.8575855211E-1*t27*t25*rho-0.2868132559E-1*rho+0.1480890852E-1;
    t34 = log(rho);
    t38 = 0.3001504024E-1*rho*t34+0.3718535686E-1*rho+0.3198902536E-2;
    t52 = 360.0*t4-0.18E3*phi-0.18E3*old+60.0;
    return(2.0*A_*(4.0*t5+3.0*(2.0*alpha_-2.0)*t4+2.0*(-3.0*alpha_+1.0)*t3)*t17+
beta_*(t24*t31-t24*t38)+(2.0*A_*(0.12E2*phi+0.12E2*old+12.0*alpha_-12.0)*t17+beta_*
(t52*t31-t52*t38))*(phi-old)/24.0);
  }
}


inline double dphireactionSource ( double rho ,double phi ,double old ) const
{
  double t13;
  double t21;
  double t22;
  double t24;
  double t28;
  double t3;
  double t31;
  double t35;
  double t4;
  double t43;
  double t63;
  {
    t3 = 0.5*phi+0.5*old;
    t4 = t3*t3;
    t13 = 1/delta_;
    t21 = 0.6E2*t4*t3-0.9E2*t4+0.15E2*phi+0.15E2*old;
    t22 = rho*rho;
    t24 = t22*t22;
    t28 = 0.8575855211E-1*t24*t22*rho-0.2868132559E-1*rho+0.1480890852E-1;
    t31 = log(rho);
    t35 = 0.3001504024E-1*rho*t31+0.3718535686E-1*rho+0.3198902536E-2;
    t43 = 0.18E3*phi+0.18E3*old-0.18E3;
    t63 = 360.0*t4-0.18E3*phi-0.18E3*old+60.0;
    return(2.0*A_*(0.6E1*t4+0.3E1*(2.0*alpha_-2.0)*t3-0.3E1*alpha_+0.1E1)*t13+beta_
*(t21*t28-t21*t35)+(0.24E2*A_*t13+beta_*(t43*t28-t43*t35))*(phi-old)/24.0+A_*(
0.12E2*phi+0.12E2*old+12.0*alpha_-12.0)*t13/12.0+beta_*(t63*t28-t63*t35)/24.0);
  }
}


inline double chemicalPotential ( double rho ,double phi ,double old ) const
{
  double t1;
  double t11;
  double t12;
  double t13;
  double t18;
  double t19;
  double t2;
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
    t12 = t11*t11;
    t13 = t12*t12;
    t18 = 1.0-t4+t5-t7;
    t19 = log(t11);
    return(beta_*(t8*(0.6003098648*t13*t12-0.2868132559E-1)+t18*(0.3001504024E-1
*t19+0.672003971E-1))+beta_*(0.1800929594E2*t8*t13-0.3001504024E-1*t18/t12)*(rho
-old)/24.0);
  }
}

inline double drhochemicalPotential ( double rho ,double phi ,double old ) const
{
  double t1;
  double t11;
  double t12;
  double t13;
  double t17;
  double t2;
  double t23;
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
    t12 = t11*t11;
    t13 = t12*t12;
    t17 = 1.0-t4+t5-t7;
    t23 = t12*t11;
    return(beta_*(0.1800929594E1*t8*t13*t11+0.1500752012E-1*t17/t11)+beta_*(
0.3601859188E2*t8*t23+0.3001504024E-1*t17/t23)*(rho-old)/24.0+beta_*(
0.1800929594E2*t8*t13-0.3001504024E-1*t17/t12)/24.0);
  }
}


inline double dphichemicalPotential ( double rho ,double phi ,double old ) const
{
  double t1;
  double t10;
  double t11;
  double t12;
  double t17;
  double t18;
  double t2;
  double t7;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t7 = 30.0*t2-60.0*t1*phi+30.0*t1;
    t10 = 0.5*rho+0.5*old;
    t11 = t10*t10;
    t12 = t11*t11;
    t17 = -t7;
    t18 = log(t10);
    return(beta_*(t7*(0.6003098648*t12*t11-0.2868132559E-1)+t17*(0.3001504024E-1
*t18+0.672003971E-1))+beta_*(0.1800929594E2*t7*t12-0.3001504024E-1*t17/t11)*(rho
-old)/24.0);
  }
}


inline double pressure ( double rho ,double phi ) const
{
  double t1;
  double t11;
  double t2;
  double t22;
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
    t9 = rho*rho;
    t11 = t9*t9;
    t22 = log(rho);
    return(beta_*((t4-t5+t7)*(-0.8575855211E-1*t11*t9*rho+0.2868132559E-1*rho
-0.1480890852E-1+rho*(0.6003098648*t11*t9-0.2868132559E-1))+(1.0-t4+t5-t7)*(
-0.3001504024E-1*rho*t22-0.3718535686E-1*rho-0.3198902536E-2+rho*(
0.3001504024E-1*t22+0.672003971E-1)))-2.0*A_*(t2+(2.0*alpha_-2.0)*t6+(-3.0*alpha_+
1.0)*t1+alpha_)/delta_);
  }
}


inline double a ( double rho ,double phi ) const
{
  double t1;
  double t11;
  double t19;
  double t2;
  double t3;
  double t4;
  double t5;
  double t6;
  {
    t1 = rho*rho;
    t2 = t1*t1;
    t3 = t2*t1;
    t4 = phi*phi;
    t5 = t4*t4;
    t6 = t5*phi;
    t11 = t4*phi;
    t19 = sqrt(beta_*(0.5402788784E12*t3*t6-0.1350697196E13*t3*t5+
0.9004647972E12*t3*t11+750376006.0-4502256035.0*t6+0.1125564009E11*t5
-7503760060.0*t11));
    return(0.632455532E-5*t19);
  }
}


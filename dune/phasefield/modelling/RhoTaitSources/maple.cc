double evalRho(double x)
{
  {
    return(0.4958034541*x+0.1065766549);
  }
}


inline double helmholtz ( double rho ,double phi ) const
{
  double t17;
  double t18;
  double t19;
  double t2;
  double t21;
  double t23;
  double t3;
  double t30;
  double t6;
  {
    t2 = phi*phi;
    t3 = t2*t2;
    t6 = t2*phi;
    t17 = 6.0*t3*phi;
    t18 = 15.0*t3;
    t19 = 10.0*t6;
    t21 = rho*rho;
    t23 = t21*t21;
    t30 = log(rho);
    return(2.0*rho*A_*(t3+(2.0*alpha_-2.0)*t6+(-3.0*alpha_+1.0)*t2+alpha_)/delta_+(
t17-t18+t19)*(0.8575855211E-1*t23*t21*rho-0.2868132559E-1*rho+0.1480890852E-1)+
(1.0-t17+t18-t19)*(0.3001504024E-1*rho*t30+0.3718535686E-1*rho+0.3198902536E-2)
);
  }
}


inline double reactionSource ( double rho ,double phi ) const
{
  double t18;
  double t2;
  double t22;
  double t23;
  double t25;
  double t3;
  double t32;
  {
    t2 = phi*phi;
    t3 = t2*phi;
    t18 = t2*t2;
    t22 = 30.0*t18-60.0*t3+30.0*t2;
    t23 = rho*rho;
    t25 = t23*t23;
    t32 = log(rho);
    return(2.0*rho*A_*(4.0*t3+3.0*(2.0*alpha_-2.0)*t2+2.0*(-3.0*alpha_+1.0)*phi)/
delta_+t22*(0.8575855211E-1*t25*t23*rho-0.2868132559E-1*rho+0.1480890852E-1)-t22
*(0.3001504024E-1*rho*t32+0.3718535686E-1*rho+0.3198902536E-2));
  }
}


inline double dphireactionSource ( double rho ,double phi ) const
{
  double t18;
  double t19;
  double t2;
  double t21;
  double t28;
  {
    t2 = phi*phi;
    t18 = 120.0*t2*phi-180.0*t2+60.0*phi;
    t19 = rho*rho;
    t21 = t19*t19;
    t28 = log(rho);
    return(2.0*rho*A_*(12.0*t2+6.0*(2.0*alpha_-2.0)*phi-6.0*alpha_+2.0)/delta_+t18*
(0.8575855211E-1*t21*t19*rho-0.2868132559E-1*rho+0.1480890852E-1)-t18*(
0.3001504024E-1*rho*t28+0.3718535686E-1*rho+0.3198902536E-2));
  }
}


inline double chemicalPotential ( double rho ,double phi ) const
{
  double t1;
  double t11;
  double t16;
  double t2;
  double t20;
  double t33;
  double t37;
  double t4;
  double t43;
  double t5;
  double t6;
  double t7;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = t2*phi*delta_;
    t5 = rho*rho;
    t6 = t5*t5;
    t7 = t6*t5;
    t11 = t2*delta_;
    t16 = t1*phi*delta_;
    t20 = log(rho);
    t33 = A_*t5*rho;
    t37 = A_*t5;
    t43 = 0.1800929594E12*t4*t7-0.287645168E11*t4-0.4502323986E12*t11*t7+
0.71911292E11*t11+0.3001549324E12*t16*t7-0.4794086134E11*t16+1500752012.0*delta_
*t20+3360019855.0*delta_-9004512070.0*t4*t20+0.2251128018E11*t11*t20
-0.1500752012E11*t16*t20+0.1E12*A_*t6+0.2E12*t33*alpha_-0.2E12*t33-0.3E12*t37*
alpha_+0.1E12*t37+0.1E12*A_*alpha_;
    return(0.2E-10*t43/delta_);
  }
}

inline double drhochemicalPotential ( double rho ,double phi ) const
{
  double t1;
  double t10;
  double t14;
  double t2;
  double t24;
  double t28;
  double t32;
  double t4;
  double t5;
  double t6;
  double t7;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = t2*phi*delta_;
    t5 = rho*rho;
    t6 = t5*t5;
    t7 = t6*t5;
    t10 = t2*delta_;
    t14 = t1*phi*delta_;
    t24 = A_*t5*rho;
    t28 = A_*t5;
    t32 = 0.5402788784E12*t4*t7-0.1350697196E13*t10*t7+0.9004647972E12*t14*t7+
750376006.0*delta_-4502256036.0*t4+0.1125564009E11*t10-7503760060.0*t14+0.2E12*A_
*t6+0.3E12*t24*alpha_-0.3E12*t24-0.3E12*t28*alpha_+0.1E12*t28;
    return(0.4E-10*t32/rho/delta_);
  }
}


inline double dphichemicalPotential ( double rho ,double phi ) const
{
  double t1;
  double t16;
  double t2;
  double t3;
  double t4;
  double t5;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = rho*rho;
    t4 = t3*t3;
    t5 = t4*t3;
    t9 = t1*phi;
    t16 = log(rho);
    return(0.1800929594E2*t2*t5-0.2876451681E1*t2-0.3601859189E2*t9*t5+
0.5752903361E1*t9+0.1800929594E2*t1*t5-0.2876451681E1*t1-0.9004512072*t2*t16+
0.1800902414E1*t9*t16-0.9004512072*t1*t16);
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
  double t7;
  double t9;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t7 = 10.0*t1*phi;
    t9 = rho*rho;
    t11 = t9*t9;
    t22 = log(rho);
    return((t4-t5+t7)*(-0.8575855211E-1*t11*t9*rho+0.2868132559E-1*rho
-0.1480890852E-1+rho*(0.6003098648*t11*t9-0.2868132559E-1))+(1.0-t4+t5-t7)*(
-0.3001504024E-1*rho*t22-0.3718535686E-1*rho-0.3198902536E-2+rho*(
0.3001504024E-1*t22+0.672003971E-1)));
  }
}


inline double a ( double rho ,double phi ) const
{
  double t1;
  double t11;
  double t18;
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
    t18 = sqrt(0.5402788782E13*t3*t6-0.1350697196E14*t5*t3+0.9004647972E13*t11*
t3+7503760060.0-0.4502256035E11*t6+0.1125564009E12*t5-0.750376006E11*t11);
    return(0.2E-5*t18);
  }
}



inline double helmholtz ( double rho ,double phi ) const
{
  double t1;
  double t16;
  double t17;
  double t18;
  double t2;
  double t23;
  double t5;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t5 = t1*phi;
    t16 = 6.0*t2*phi;
    t17 = 15.0*t2;
    t18 = 10.0*t5;
    t23 = log(rho);
    return(2.0*A_*(t2+(2.0*alpha_-2.0)*t5+(-3.0*alpha_+1.0)*t1+alpha_)/delta_+beta_*(
(t16-t17+t18)*((bl-cl)*rho+cl*rho*t23+dl)+(1.0-t16+t17-t18)*((bv-cv)*rho+cv*rho
*t23+dv)));
  }
}


inline double reactionSource ( double rho ,double phi ,double old ) const
{
  double t17;
  double t20;
  double t24;
  double t28;
  double t3;
  double t30;
  double t37;
  double t4;
  double t5;
  double t51;
  {
    t3 = 0.5*phi+0.5*old;
    t4 = t3*t3;
    t5 = t4*t3;
    t17 = 1/delta_;
    t20 = t4*t4;
    t24 = 30.0*t20-60.0*t5+30.0*t4;
    t28 = log(rho);
    t30 = (bl-cl)*rho+cl*rho*t28+dl;
    t37 = (bv-cv)*rho+cv*rho*t28+dv;
    t51 = 360.0*t4-0.18E3*phi-0.18E3*old+60.0;
    return(2.0*A_*(4.0*t5+3.0*(2.0*alpha_-2.0)*t4+2.0*(-3.0*alpha_+1.0)*t3)*t17+
beta_*(t24*t30-t24*t37)+(2.0*A_*(0.12E2*phi+0.12E2*old+12.0*alpha_-12.0)*t17+beta_*
(t51*t30-t51*t37))*(phi-old)/24.0);
  }
}


inline double dphireactionSource ( double rho ,double phi ,double old ) const
{
  double t13;
  double t21;
  double t25;
  double t27;
  double t3;
  double t34;
  double t4;
  double t42;
  double t62;
  {
    t3 = 0.5*phi+0.5*old;
    t4 = t3*t3;
    t13 = 1/delta_;
    t21 = 0.6E2*t4*t3-0.9E2*t4+0.15E2*phi+0.15E2*old;
    t25 = log(rho);
    t27 = (bl-cl)*rho+cl*rho*t25+dl;
    t34 = (bv-cv)*rho+cv*rho*t25+dv;
    t42 = 0.18E3*phi+0.18E3*old-0.18E3;
    t62 = 360.0*t4-0.18E3*phi-0.18E3*old+60.0;
    return(2.0*A_*(0.6E1*t4+0.3E1*(2.0*alpha_-2.0)*t3-0.3E1*alpha_+0.1E1)*t13+beta_
*(t21*t27-t21*t34)+(0.24E2*A_*t13+beta_*(t42*t27-t42*t34))*(phi-old)/24.0+A_*(
0.12E2*phi+0.12E2*old+12.0*alpha_-12.0)*t13/12.0+beta_*(t62*t27-t62*t34)/24.0);
  }
}


inline double drhoreactionSource ( double rho ,double phi ,double old ) const
{
  double t10;
  double t11;
  double t13;
  double t17;
  double t24;
  double t3;
  double t4;
  double t5;
  {
    t3 = 0.5*phi+0.5*old;
    t4 = t3*t3;
    t5 = t4*t4;
    t10 = 30.0*t5-60.0*t4*t3+30.0*t4;
    t11 = log(rho);
    t13 = bl+cl*t11;
    t17 = bv+cv*t11;
    t24 = 360.0*t4-0.18E3*phi-0.18E3*old+60.0;
    return(beta_*(t10*t13-t10*t17)+beta_*(t24*t13-t24*t17)*(phi-old)/24.0);
  }
}


inline double chemicalPotential ( double rho ,double phi ,double old ) const
{
  double t1;
  double t11;
  double t12;
  double t16;
  double t2;
  double t23;
  double t24;
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
    t23 = t11*t11;
    t24 = 1/t23;
    return(beta_*(t8*(bl+cl*t12)+t16*(bv+cv*t12))+beta_*(-t8*cl*t24-t16*cv*t24)*(
rho-old)/24.0);
  }
}

inline double drhochemicalPotential ( double rho ,double phi ,double old ) const
{
  double t1;
  double t12;
  double t13;
  double t17;
  double t2;
  double t22;
  double t24;
  double t34;
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
    t9 = (t4-t5+t7)*cl;
    t12 = 0.5*rho+0.5*old;
    t13 = 1/t12;
    t17 = (1.0-t4+t5-t7)*cv;
    t22 = t12*t12;
    t24 = 1/t22/t12;
    t34 = 1/t22;
    return(beta_*(0.5*t9*t13+0.5*t17*t13)+beta_*(0.1E1*t9*t24+0.1E1*t24*t17)*(rho
-old)/24.0+beta_*(-t9*t34-t17*t34)/24.0);
  }
}


inline double dphichemicalPotential ( double rho ,double phi ,double old ) const
{
  double t1;
  double t10;
  double t11;
  double t15;
  double t2;
  double t22;
  double t23;
  double t7;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t7 = 30.0*t2-60.0*t1*phi+30.0*t1;
    t10 = 0.5*rho+0.5*old;
    t11 = log(t10);
    t15 = -t7;
    t22 = t10*t10;
    t23 = 1/t22;
    return(beta_*(t7*(bl+cl*t11)+t15*(bv+cv*t11))+beta_*(-t7*cl*t23-t15*cv*t23)*(
rho-old)/24.0);
  }
}


inline double pressure ( double rho ,double phi ) const
{
  double t1;
  double t12;
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
    t12 = log(rho);
    return(beta_*((t4-t5+t7)*(-(bl-cl)*rho-cl*rho*t12-dl+rho*(bl+cl*t12))+(1.0-
t4+t5-t7)*(-(bv-cv)*rho-cv*rho*t12-dv+rho*(bv+cv*t12)))-2.0*A_*(t2+(2.0*alpha_
-2.0)*t6+(-3.0*alpha_+1.0)*t1+alpha_)/delta_);
  }
}


inline double a ( double rho ,double phi ) const
{
  double t1;
  double t2;
  double t3;
  double t8;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = t2*phi;
    t8 = t1*phi;
    return(sqrt(beta_*(6.0*cl*t3-15.0*cl*t2+10.0*cl*t8+cv-6.0*cv*t3+15.0*cv*t2
-10.0*cv*t8)));
  }
}


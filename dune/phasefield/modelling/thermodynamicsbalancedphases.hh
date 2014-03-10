#ifndef THERMODYNAMICSBALANCEDPHASES_HH
#define THERMODYNAMICSBALANCEDPHASES_HH

// system include
#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>

// dune-fem include
#include <dune/fem/io/parameter.hh>

using namespace Dune;
// Thermodynamics
// see MininimizerSonntag.mv

#include "thermodynamicsinterface.hh"
//thermodynamics with balanced phases, see balancedthermo.mw
class BalancedThermodynamics:
  public Thermodynamics<BalancedThermodynamics>
{
  typedef Thermodynamics<BalancedThermodynamics> BaseType;

  public:
  BalancedThermodynamics():
    BaseType(),
    delta_( Dune::Fem::Parameter::getValue<double>( "phasefield.delta" ) ),
    deltaInv_( 1./delta_ ),
    epsilon_(Dune::Fem::Parameter::getValue<double>( "phasefield.mu1" ) ),
    mu1_( epsilon_ ),
    mu2_( epsilon_ )
  {
  }

  inline void init() const 
  {
    abort();
  }

  inline double helmholtz(double& rho,double& phi) const
  {
    double t10;
    double t11;
    double t13;
    double t16;
    double t17;
    double t2;
    double t5;
    double t8;

    t2 = phi*phi;
    t5 = pow(1.0-phi,2.0);
    t8 = t2*t2;
    t10 = 6.0*t8*phi;
    t11 = 15.0*t8;
    t13 = 10.0*t2*phi;
    t16 = log(rho);
    t17 = rho*t16;
    return(2.0*deltaInv_*t2*t5+(t10-t11+t13)*(-0.25E1*rho+0.15E1*t17+1.0)+(1.0-t10+
          t11-t13)*(-7.0*rho+3.0*t17+0.945940263E1));
  }



  inline double reactionSource(double& rho,double& phi) const
  { 
    double t1;
    double t11;
    double t16;
    double t18;
    double t19;
    double t3;
    double t4;
    double t7;

    t1 = deltaInv_;
    t3 = 1.0-phi;
    t4 = t3*t3;
    t7 = phi*phi;
    t11 = t7*t7;
    t16 = 30.0*t11-60.0*t7*phi+30.0*t7;
    t18 = log(rho);
    t19 = rho*t18;
    return(4.0*t1*phi*t4-4.0*t1*t7*t3+t16*(-0.25E1*rho+0.15E1*t19+1.0)-t16*( -7.0*rho+3.0*t19+0.945940263E1));

  }

  inline double dphireactionSource( double &rho, double & phi)const
  { 
    double t1;
    double t11;
    double t16;
    double t18;
    double t19;
    double t3;
    double t4;
    double t7;

    t1 = deltaInv_;
    t3 = 1.0-phi;
    t4 = t3*t3;
    t7 = phi*phi;
    t11 = t7*t7;
    t16 = 30.0*t11-60.0*t7*phi+30.0*t7;
    t18 = log(rho);
    t19 = rho*t18;
    return(4.0*t1*phi*t4-4.0*t1*t7*t3+t16*(-0.25E1*rho+0.15E1*t19+1.0)-t16*(
          -7.0*rho+3.0*t19+0.945940263E1));

  }

  inline double chemicalPotential(double& rho,double& phi) const
  {
    double t1;
    double t11;
    double t2;
    double t3;
    double t5;
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = t2*phi;
    t5 = log(rho);
    t11 = t1*phi;
    return(18.0*t3-9.0*t3*t5-45.0*t2+0.225E2*t2*t5+30.0*t11-15.0*t11*t5-4.0+3.0
        *t5);

  }

  inline double dphichemicalPotential(double& rho,double& phi) const
  {

    double t1;
    double t2;
    double t4;
    double t7;
      t1 = phi*phi;
      t2 = t1*t1;
      t4 = log(rho);
      t7 = t1*phi;
      return(90.0*t2-45.0*t2*t4-180.0*t7+0.9E2*t7*t4+90.0*t1-45.0*t1*t4);
  }

  inline double drhochemicalPotential(double& rho,double& phi) const
  {
    double t1;
    double t2;

    t1 = phi*phi;
    t2 = t1*t1;
    return(-0.15E1*(6.0*t2*phi-15.0*t2+10.0*t1*phi-2.0)/rho);
  }

  inline double  pressure( double& rho, double& phi) const
  {
    double delta=delta_;

    double t1;
    double t12;
    double t13;
    double t2;
    double t23;
    double t4;
    double t8;
    
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = t2*phi*delta;
    t8 = t2*delta;
    t12 = t1*phi;
    t13 = t12*delta;
    t23 = 900000000.0*t4*rho-5075641578.0*t4-2250000000.0*t8*rho+
      0.1268910394E11*t8+1500000000.0*t13*rho-8459402630.0*t13-300000000.0*delta*rho+
      945940263.0*delta+200000000.0*t1-400000000.0*t12+200000000.0*t2;
    return(-0.1E-7*t23/delta);

  }



  inline double a(double rho,double phi) const
  {
    return 1.6;	
  }



  public:

  inline double delta()    const { return delta_; }
  inline double deltaInv() const { return deltaInv_; }
  inline double mu1()      const { return mu2_; }
  inline double mu2()      const { return mu1_; }

  private:
  mutable double  delta_;
  mutable double  deltaInv_;
  mutable double  epsilon_;
  mutable double  mu1_,mu2_;

};



#endif // file define

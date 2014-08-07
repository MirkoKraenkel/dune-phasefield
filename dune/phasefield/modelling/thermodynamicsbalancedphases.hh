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
    mu2_( Dune::Fem::Parameter::getValue<double>("phasefield.mu2") ),
    reaction_( Dune::Fem::Parameter::getValue<double>( "phasefield.reactionrate") )
  {
  }
#if RHOMODEL
  inline double h( double rho ) const
  {
    return rho;
  }
  
  inline double h2( double rho) const
  {
    return  1./h(rho);
  }

  inline double h2prime( double rho ) const
  {
    return -1./h(rho)*h(rho);
  }
#else
  inline double h( double rho ) const
  {
    return 1.;
  }
  
  inline double h2( double rho) const
  {
    return  1.;
  }

  inline double h2prime( double rho ) const
  {
    return 0.;
  }
#endif
  inline double reactionFactor() const
  {
    return reaction_;
  }

  inline double doubleWell(double phi) const
  {
    double    t2 = phi*phi;
    double    t5 = pow(1.0-phi,2.0);
    return t2*t5; 
  }


  inline double helmholtz(double rho,double phi) const
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
    return(2.*deltaInv_*t2*t5*h(rho)+(t10-t11+t13)*(-0.25E1*rho+0.15E1*t17+1.0)+(1.0-t10+
          t11-t13)*(-7.0*rho+3.0*t17+0.945940263E1));
  }

  inline double homSource( double rho, double phi ) const
  {
    double t11;
    double t16;
    double t18;
    double t19;
    double t7;

    t7 = phi*phi;
    t11 = t7*t7;
    t16 = 30.0*t11-60.0*t7*phi+30.0*t7;
    t18 = log(rho);
    t19 = rho*t18;
    return(t16*(-0.25E1*rho+0.15E1*t19+1.0)-t16*( -7.0*rho+3.0*t19+0.945940263E1));

  }

  inline double reactionSource(double rho,double phi) const
  { 
    double t1;
    double t3;
    double t4;
    double t7;

    t1 = deltaInv_*h(rho);
    t3 = 1.0-phi;
    t4 = t3*t3;
    t7 = phi*phi;
    return( 4.0*t1*phi*t4-4.0*t1*t7*t3 + homSource( rho, phi ));
  }

  inline double dphireactionSource( double rho, double  phi) const
  { 
    double t13;
    double t19;
    double t2;
    double t25;
    double t5;
    double t8;
    double t9;
    {
      t2 = phi*phi;
      t5 = t2*phi*delta_;
      t8 = log(rho);
      t9 = rho*t8;
      t13 = t2*delta_;
      t19 = phi*delta_;
      t25 = -20000000.0+120000000.0*phi-120000000.0*t2-2700000000.0*t5*rho+
        900000000.0*t5*t9+5075641580.0*t5+4050000000.0*t13*rho-1350000000.0*t13*t9
        -7613462365.0*t13-1350000000.0*t19*rho+450000000.0*t19*t9+2537820789.0*t19;
      return(-0.2E-6*t25/delta_);
    }





  }

  inline double chemicalPotential(double rho,double phi) const
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
    double dW;
#if RHOMODEL
    dW=2.0*deltaInv_*doubleWell(phi);
#else
    dW=0;
#endif
    return(18.0*t3-9.0*t3*t5-45.0*t2+0.225E2*t2*t5+30.0*t11-15.0*t11*t5-4.0+3.0
        *t5)+dW;

  }

  inline double dphichemicalPotential(double rho,double phi) const
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

  inline double drhochemicalPotential(double rho,double phi) const
  {
    double t1;
    double t2;

    t1 = phi*phi;
    t2 = t1*t1;
    return(-0.15E1*(6.0*t2*phi-15.0*t2+10.0*t1*phi-2.0)/rho);
  }

  inline double  pressure( double rho, double phi) const
  {
    double t1;
    double t10;
    double t11;
    double t2;
    double t4;
    double t5;
    double t7;
    double pressure;
    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t7 = 10.0*t1*phi;
    t10 = log(rho);
    t11 = rho*t10;

    pressure=((t4-t5+t7)*(0.25E1*rho-0.15E1*t11-1.0+rho*(-1.0+0.15E1*t10))+(1.0-t4
          +t5-t7)*(7.0*rho-3.0*t11-0.945940263E1+rho*(-4.0+3.0*t10)));
#if RHOMODEL
    pressure-=doubleWell(phi)*deltaInv_;
#endif
    return pressure;

  }



  inline double a(double rho,double phi) const
  {
    return 1.6;	
  }



  public:

  inline double delta()    const { return delta_; }
  inline double deltaInv() const { return deltaInv_; }
  inline double mu1()      const { return mu1_; }
  inline double mu2()      const { return mu2_; }

  private:
  mutable double  delta_;
  mutable double  deltaInv_;
  mutable double  epsilon_;
  mutable double  mu1_,mu2_;
  mutable double reaction_;
};



#endif // file define

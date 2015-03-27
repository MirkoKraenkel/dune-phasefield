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
    alpha_(Dune::Fem::Parameter::getValue<double>("phasefield.alpha",0)),
    beta_(Dune::Fem::Parameter::getValue<double>("phasefield.beta",1)),
    A_(Dune::Fem::Parameter::getValue<double>("phasefield.A")),
    mu1_(Dune::Fem::Parameter::getValue<double>( "phasefield.mu1" ) ),
    mu2_( Dune::Fem::Parameter::getValue<double>("phasefield.mu2") ),
    reaction_( Dune::Fem::Parameter::getValue<double>( "phasefield.reactionrate") )
  {
    lipschitzC_=1/(4*A_*reaction_/(delta_));
  }
#if RHO_H2CONSTMODEL
  inline double h( double rho ) const
  {
    return rho;
  }
  
  inline double h2( double rho) const
  {
    return 1.;

 }

  inline double h2prime( double rho ) const
  {
    return 0;
  }
#include "InvRhoSources/maple.cc"

#else
  inline double h( double rho ) const
  {
#if UNBALMODEL
    return rho;
#else
    return 1.;
#endif
}
  
  inline double h2( double rho) const
  {
    return  1.;
  }

  inline double h2prime( double rho ) const
  {
    return 0.;
  }
#if UNBALMODEL
#include "AltAltSources/maple.cc"
#else
#include "TaitSources/maple.cc"
#endif
#endif

  public:

  inline double reactionFactor() const { return reaction_; }
  inline double lipschitzC() const { return lipschitzC_; }
  inline double delta()    const { return delta_; }
  inline double deltaInv() const { return deltaInv_; }
  inline double mu1()      const { return mu1_; }
  inline double mu2()      const { return mu2_; }

  private:
  mutable double  delta_;
  mutable double  deltaInv_;
  mutable double  alpha_; 
  mutable double  beta_;
  mutable double  A_;
  mutable double  mu1_,mu2_;
  mutable double  reaction_;
  mutable double  lipschitzC_;
};



#endif // file define

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
    return -1./(h(rho)*h(rho));
  }
#include "InvRhoSources/maple.c"
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
#include "ConstRhoSources/maple.c"
#endif
  inline double reactionFactor() const
  {
    return reaction_;
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

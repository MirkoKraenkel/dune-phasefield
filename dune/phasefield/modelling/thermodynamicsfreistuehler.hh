#ifndef THERMODYNAMICSFREISTUEHLER_HH
#define THERMODYNAMICSFREISTUEHLER_HH

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
class ThermodynamicsFreistuehler:
  public Thermodynamics< ThermodynamicsFreistuehler >
{
  typedef Thermodynamics< ThermodynamicsFreistuehler > BaseType;

  public:
  ThermodynamicsFreistuehler():
    BaseType(),
    delta_( Dune::Fem::Parameter::getValue<double>( "phasefield.delta" ) ),
    deltaInv_( 1./delta_ ),
    epsilon_(Dune::Fem::Parameter::getValue<double>( "phasefield.mu1" ) ),
    mu1_( epsilon_ ),
    mu2_( epsilon_ )
    reaction_( Dune::Fem::Parameter::getValue<double>( "phasefield.reactionrate")),
    theta_( Dune::Fem::Parameter::getValue<double>(" phasefield.criticalTemp", 0.21))
    {
    }
    
  inline double h( double rho ) const
    {
      return rho;
    }
  
  inline double h2( double rho) const
    {
      return  h(rho);
    }

  inline double h2prime( double rho ) const
    {
      return 1;
    }

  
  #include "FreistuehlerSources/maple.cc"
  
  inline double reactionFactor() const
  {
    return reaction_;
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
  mutable double  reaction_;
  mutable double  theta_;
};



#endif // file define

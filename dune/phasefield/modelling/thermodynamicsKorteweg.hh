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
    mu1_( Dune::Fem::Parameter::getValue<double>("phasefield.mu1") ),
    mu2_( Dune::Fem::Parameter::getValue<double>("phasefield.mu2") ),
    cst_( Dune::Fem::Parameter::getValue<double>("korteweg.c") ),
    theta_( Dune::Fem::Parameter::getValue<double>( "korteweg.theta") ),
    reaction_( Dune::Fem::Parameter::getValue<double>( "phasefield.reactionrate") )
  {
  }
#include "KortewegSources/maple.c"

  public:

  inline double delta()    const { return delta_; }
  inline double deltaInv() const { return deltaInv_; }
  inline double mu1()      const { return mu1_; }
  inline double mu2()      const { return mu2_; }
  inline double theta()    const { return theta_;}
  inline double cst()      const { return cst_;}
  private:
  mutable double  delta_;
  mutable double  deltaInv_;
  mutable double  epsilon_;
  mutable double  mu1_,mu2_;
  mutable double cst_;
  mutable double theta_;
  mutable double  reaction_;

};



#endif // file define

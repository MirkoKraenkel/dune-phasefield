#ifndef THERMODYNAMICS_KORTEWEG_HH
#define THERMODYNAMICS_KORTEWEG_HH

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
class NskThermodynamics:
  public Thermodynamics<NskThermodynamics>
{
  typedef Thermodynamics<NskThermodynamics> BaseType;

  public:
  NskThermodynamics():
    BaseType(),
    delta_( Dune::Fem::Parameter::getValue<double>( "phasefield.delta" ) ),
    deltaInv_( 1./delta_ ),
    mu1Liq_(Dune::Fem::Parameter::getValue<double>( "phasefield.mu1_liq" ) ),
    mu2Liq_( Dune::Fem::Parameter::getValue<double>("phasefield.mu2_liq") ),
    mu1Vap_(Dune::Fem::Parameter::getValue<double>( "phasefield.mu1_vap" ) ),
    mu2Vap_( Dune::Fem::Parameter::getValue<double>("phasefield.mu2_vap") ),
    theta_( Dune::Fem::Parameter::getValue<double>( "korteweg.theta") )
  {
  }
#include "KortewegSources/nskmaple.cc"

  public:

  inline double delta()    const { return delta_; }
  inline double deltaInv() const { return deltaInv_; }
  inline double mu1Liq()      const { return mu1Liq_; }
  inline double mu2Liq()      const { return mu2Liq_; }
  inline double mu1Vap()      const { return mu1Vap_; }
  inline double mu2Vap()      const { return mu2Vap_; }
  inline double theta()    const { return theta_;}
  private:
  mutable double  delta_;
  mutable double  deltaInv_;
  mutable double  epsilon_;
  mutable double  mu1Liq_,mu1Vap_,mu2Liq_,mu2Vap_;
  mutable double  theta_;

};



#endif // file define

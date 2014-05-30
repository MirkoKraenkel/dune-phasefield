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
// see balancedthermoGeneric.hh
#include "thermodynamicsinterface.hh"

//thermodynamics with balanced phases, see balancedthermo.mw
class GeneralThermodynamics:
  public Thermodynamics<BalancedThermodynamics>
{
  typedef Thermodynamics<BalancedThermodynamics> BaseType;

  public:
  GeneralThermodynamics():
    BaseType(),
    delta_( Dune::Fem::Parameter::getValue<double>( "phasefield.delta" ) ),
    deltaInv_( 1./delta_ ),
    epsilon_( Dune::Fem::Parameter::getValue<double>( "phasefield.mu1" ) ),
    alpha_( Dune::Fem::Parameter::getValue<double>( "phasefield.alpha" ) ), 
    mu1_( epsilon_ ),
    mu2_( epsilon_ ),
    a1_( Dune::Fem::Parameter::getValue<double>( "phasefield.stiffenedgas.a1" ) ),
    a2_( Dune::Fem::Parameter::getValue<double>( "phasefield.stiffenedgas.a2" ) ),
    b1_( Dune::Fem::Parameter::getValue<double>( "phasefield.stiffenedgas.b1" ) ),  
    b2_( Dune::Fem::Parameter::getValue<double>( "phasefield.stiffenedgas.b2" ) ),
    p1_( Dune::Fem::Parameter::getValue<double>( "phasefield.stiffenedgas.p1" ) ),  
    p2_( Dune::Fem::Parameter::getValue<double>( "phasefield.stiffenedgas.p2" ) ),
  {
  }


  inline void init() const 
  {
    abort();
  }

  inline double h1(double rho)const
  {
    return rho;
  }

  inline double h2( double rho) const
  {
    return 1./rho;
  }

  inline double h1prime(double rho)const
  {
    return 1.;
  }

  inline double h2prime(double rho)const
  {
    return -1./(rho*rho); 
  }



  inline double doubleWell(double phi) const
  {
    double t1;
    t1 = phi*phi;
    t1*=deltaInv_;
    return(2.0*t1*(t1-2.0*phi*(1.0+alpha_)+1.0+3.0*alpha_));

  }

  inline double doubleWellprime(double phi) const
  {
    double t1;
    t1 = phi*phi;
    t1*=deltaInv_;

    return(8.0*t1*phi+6.0*(2.0*alpha-2.0)*t1+4.0*(-3.0*alpha+1.0)*phi);
  }



  inline double helmholtz(double& rho,double& phi) const
  {
    double t1;
    double t12;
    double t2;
    double t4;
    double t5;
    double t7;

    t1 = phi*phi;
    t2 = t1*t1;
    t4 = 6.0*t2*phi;
    t5 = 15.0*t2;
    t7 = 10.0*t1*phi;
    t12 = log(rho);

    double helm=((t4-t5+t7)*((b1-a1)*rho+a1*rho*t12+p1)+(1.0-t4+t5-t7)*((b2-a2)*rho+
          a2*rho*t12+p2));

    helm+=h1(rho)*doubleWell(phi);

    return helm;

  }



  inline double reactionSource(double& rho,double& phi) const
  { 
    double t1;
    double t11;
    double t2;
    double t7;

    t1 = phi*phi;
    t2 = t1*t1;
    t7 = 30.0*t2-60.0*t1*phi+30.0*t1;
    t11 = log(rho);
    double source=(t7*((b1-a1)*rho+a1*rho*t11+p1)-t7*((b2-a2)*rho+a2*rho*t11+p2));
    source+=h1(rho)*doubleWellprime(phi); 
    return source;

  }

  inline double dphireactionSource( double &rho, double & phi)const
  { 
  }

  inline double chemicalPotential(double& rho,double& phi) const
  {
  }

  inline double dphichemicalPotential(double& rho,double& phi) const
  {

  }

  inline double drhochemicalPotential(double& rho,double& phi) const
  {
  }

  inline double  pressure( double& rho, double& phi) const
  {

    double t1;
    double t15;
    double t2;
    double t3;
    double t7;

    t1 = phi*phi;
    t2 = t1*t1;
    t3 = t2*phi;
    t7 = log(rho);
    t15 = t1*phi;
    double press=(6.0*t3*b1+6.0*t3*a1*t7-15.0*t2*b1-15.0*t2*a1*t7+10.0*t15*b1+10.0*t15
        *a1*t7+b2+a2*t7-6.0*t3*b2-6.0*t3*a2*t7+15.0*t2*b2+15.0*t2*a2*t7-10.0*t15*b2
        -10.0*t15*a2*t7);

    press=(-h1(rho)+rho*h1prime(rho))*doubleWell(rho); 

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
  mutable double  a1_,a2_,b1_,b2_,p1_,p2_;
};



#endif // file define

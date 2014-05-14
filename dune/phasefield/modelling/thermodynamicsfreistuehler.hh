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
  {
  }



  inline void init() const 
  {
    abort();
  }


  inline double helmholtz(double& rho,double& phi) const
  {
    double t1;
    double t3;
    double t4;
    double t6;
    double t9;

    t1 = log(phi);
    t3 = 1.0-phi;
    t4 = log(t3);
    t6 = phi*phi;
    t9 = log(rho);
    return(phi*t1+t3*t4-t6/2.0+0.15E1*phi*rho*t9+3.0*t3*rho*t9);

  }



  inline double reactionSource(double& rho,double& phi) const
  {
    double t1;
    double t3;
    double t4;
    {
      t1 = log(phi);
      t3 = log(1.0-phi);
      t4 = log(rho);
      return(t1-t3-phi-0.15E1*rho*t4);
    }

  }

  inline double dphireactionSource( double &rho, double & phi)const
  { 
    std::cout<<"thermodynamicsfreistuehler.hh dphireactionSource not implemented!\n";
    abort();
  }

  inline double chemicalPotential(double& rho,double& phi) const
  {
    double t1;

    t1 = log(rho);
    return(-0.15E1*phi*t1-0.15E1*phi+3.0*t1+3.0);


  }

  inline double dphichemicalPotential(double& rho,double& phi) const
  {
    std::cout<<"thermodynamicsfreistuehler.hh dphichemicalPotential not implemented!\n";
    abort();
  }

  inline double drhochemicalPotential(double& rho,double& phi) const
  {
    std::cout<<"thermodynamicsfreistuehler.hh dphichemicalPotential not implemented!\n";
    abort();
  }

  inline double  pressure( double rho, double phi) const
  {
    return(-0.15E1*phi*rho+3.0*rho);
  }



  inline double a(double rho,double phi) const
  {
    return 1.6;	
  }

  private:
  inline double h(double rho)const
  {

  }
  inline double doubleWell(double phi) const
  {
    abort(); 
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

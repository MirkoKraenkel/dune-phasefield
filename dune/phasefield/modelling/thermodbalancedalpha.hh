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
// see balancedthermo2.hh
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
    epsilon_( Dune::Fem::Parameter::getValue<double>( "phasefield.mu1" ) ),
    alpha_( Dune::Fem::Parameter::getValue<double>( "phasefield.alpha" ) ), 
    mu1_( epsilon_ ),
    mu2_( epsilon_ )
  {
  }


  inline void init() const 
  {
    abort();
  }

  inline double h(double rho)const
  {
#if RHOMODEL
    return rho;
#else
    return 1.;
#endif
  }
  
  inline double h2( double rho) const
  {
#if RHOMODEL  
   return 1./rho;
#else
   return 1.;
#endif
  }

  inline double doubleWell(double phi) const
  {
    double t1;
    t1 = phi*phi;
    t1*=deltaInv_;
    return(2.0*t1*(t1-2.0*phi*(1.0+alpha_)+1.0+3.0*alpha_));
  }

  inline double helmholtz(double& rho,double& phi) const
  {

    double delta=delta_;
    double alpha=alpha_;

    double t12;
    double t14;
    double t15;
    double t17;
    double t20;
    double t21;
    double t3;

    t3 = phi*phi;
    t12 = t3*t3;
    t14 = 6.0*t12*phi;
    t15 = 15.0*t12;
    t17 = 10.0*t3*phi;
    t20 = log(rho);
    t21 = rho*t20;
    return(2.0/delta*h(rho)*t3*(t3-2.0*phi*(1.0+alpha)+1.0+3.0*alpha)+(t14-t15+t17)*
        (-0.25E1*rho+0.15E1*t21+1.0)+(1.0-t14+t15-t17)*(-7.0*rho+3.0*t21+0.945940263E1)
        );

  }

  inline double reactionSource(double& rho,double& phi) const
  { 
    double alpha=alpha_;
    double delta=delta_;
    
    double t17;
    double t2;
    double t22;
    double t24;
    double t25;
    double t3;
    
    t2 = 1/delta*h(rho);
    t3 = phi*phi;
    t17 = t3*t3;
    t22 = 30.0*t17-60.0*t3*phi+30.0*t3;
    t24 = log(rho);
    t25 = rho*t24;
    return(4.0*t2*phi*(t3-2.0*phi*(1.0+alpha)+1.0+3.0*alpha)+2.0*t2*t3*(2.0*phi
          -2.0-2.0*alpha)+t22*(-0.25E1*rho+0.15E1*t25+1.0)-t22*(-7.0*rho+3.0*t25+
            0.945940263E1));

  }

  inline double dphireactionSource( double &rho, double & phi)const
  { 
  }

  inline double chemicalPotential(double& rho,double& phi) const
  {
    double t11;
    double t2;
    double t3;
    double t5;
    double chemPot;
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = t2*phi;
    t5 = log(rho);
    t11 = t1*phi;
  
    chemPot=18.0*t3-9.0*t3*t5-45.0*t2+0.225E2*t2*t5+30.0*t11-15.0*t11*t5-4.0+3.0
        *t5;
#if RHOMODEL 
    chemPot+=doubleWell(rho)*deltaInv_;
#else
#endif
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
    double t2;
    double t3;
    double t8;
    double pressure;
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = t2*phi;
    t8 = t1*phi;

    pressure=6.0*p1*t3-15.0*p1*t2+10.0*p1*t8+p0-6.0*p0*t3+15.0*p0*t2-10.0*p0*t8);

#if RHOMODEL
#else
    pressure-=doubleWell(phi)*deltaInv;
#endif

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

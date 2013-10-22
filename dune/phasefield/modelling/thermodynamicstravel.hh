#ifndef THERMODYNAMICSBALANCEDPHASES_HH
#define THERMODYNAMICSBALANCEDPHASES_HH

// system include
#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>

// dune-fem include
#include <dune/fem/io/parameter.hh>
#include "sourceterms2.hh"
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
   mu2_( epsilon_ ),
   gamma_(Dune::Fem::Parameter::getValue<double>( "gamma" ))
  {
  }

  inline void init() const 
  {
  abort();
  }
  
  inline double helmholtz(double& rho,double& phi) const
  {
		double delta=delta_;
		return 0.;
	}


	

	inline double reactionSource(double& rho,double& phi) const
	{ 
    double delta=delta_;
    double t1;
    double t11;
    double t16;
    double t18;
    double t22;
    double t3;
    double t4;
    double t7;
  
    t1 = 1/delta;
    t3 = 1.0-phi;
    t4 = t3*t3;
    t7 = phi*phi;
    t11 = t7*t7;
    t16 = 30.0*t11-60.0*t7*phi+30.0*t7;
    t18 = pow(-1.0+rho,2.0);
    t22 = pow(-2.0+rho,2.0);
    return(4.0*t1*phi*t4-4.0*t1*t7*t3+t16*t18-2.0*t16*t22);
	}
  


	inline double chemicalPotential(double& rho,double& phi) const
	{
	  double t1;
    double t10;
    double t2;
    double t3;
    
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = t2*phi;
    t10 = t1*phi;
    return(36.0*t3-12.0*t3*rho-90.0*t2+30.0*t2*rho+60.0*t10-20.0*t10*rho-8.0+4.0*rho);
	}
	
	inline double  pressure( double& rho, double& phi) const
	{
    double t1;
    double t11;
    double t2;
    double t3;
    double t5;

    t1 = phi*phi;
    t2 = t1*t1;
    t3 = t2*phi;
    t5 = rho*rho;
    t11 = t1*phi;
    return(42.0*t3-6.0*t3*t5-105.0*t2+15.0*t2*t5
          +70.0*t11-10.0*t11*t5-8.0+2.0*t5);
    
  }

	inline double a(double rho,double phi) const
	{
		return 1.6;	
	}
  
 

public:

	inline double delta()    const { return delta_;}
	inline double deltaInv() const { return deltaInv_;}
	inline double mu1()      const { return mu2_; }
  inline double mu2()      const { return mu1_; }
  inline double velo()     const { return gamma_;}

private:
	mutable double  delta_;
	mutable double  deltaInv_;
	mutable  double epsilon_;
  mutable double mu1_,mu2_;
  mutable double gamma_;
  };



#endif // file define

#ifndef THERMODYNAMICSBALANCEDPHASES_HH
#define THERMODYNAMICSBALANCEDPHASES_HH

// system include
#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>

// dune-fem include
#include <dune/fem/io/parameter.hh>
#include "sourcetermsonepressure.hh"
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
	  double t2;
    double t5;
    double t8;
    {
      t2 = phi*phi;
      t5 = pow(1.0-phi,2.0);
      t8 = log(rho);
      return(2.0/delta*t2*t5+rho*t8+1.0);
    }
  
  }


	

	inline double reactionSource(double& rho,double& phi) const
	{ 

    double t1;
    double t3;
    double t4;
    double t6;
    {
      t1 = 1/delta_;
      t3 = 1.0-phi;
      t4 = t3*t3;
      t6 = phi*phi;

      return(4.0*t1*phi*t4-4.0*t1*t6*t3);
    }
  
  }
  


	inline double chemicalPotential(double& rho,double& phi) const
	{
	  double t1;
   
    t1 = log(rho);
    return(t1+1.0);
   
  }
	
	inline double  pressure( double& rho, double& phi) const
	{
    return rho; 
  }

	inline double a(double rho,double phi) const
	{
		return velo();	
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

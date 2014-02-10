#ifndef THERMODYNAMICSPERFECTGAS_HH
#define THERMODYNAMICSPERFECTGAS_HH

// system include
#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>

// dune-fem include
#include <dune/fem/io/parameter.hh>

using namespace Dune;
// Thermodynamics
// see simplehelmholtz.mv
#include "sourcetermsSinusSolution.hh"
#include "thermodynamicsinterface.hh"


class TestThermodynamics:
	public Thermodynamics<TestThermodynamics>
{
	typedef Thermodynamics<TestThermodynamics> BaseType;

public:
  TestThermodynamics():
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
    return deltaInv_*2*phi*phi*(1-phi)*(1-phi);
  }
  
	

	inline double reactionSource(double& rho,double& phi) const
	{ 
    //    return phi*phi;
    return deltaInv_*4*(2*phi*phi*phi-3*phi*phi+phi);
   // return 0;
  }
  


	inline double chemicalPotential(double& rho,double& phi) const
	{
    return 0;
  }

	
	inline double  pressure( double& rho, double& phi) const
	{
    //return -helmholtz(rho,phi);
    return -helmholtz(rho,phi);
}
		

	inline double a(double rho,double phi) const
	{
    return 1.6;	
	}
  
 

  public:

  	inline double delta() const { return delta_; }
		inline double deltaInv() const { return deltaInv_; }
    inline double mu1() const { return mu1_; }
	 	inline double mu2() const { return mu2_; }
    inline double velo() const {return gamma_;} 
	private:
		mutable double  delta_;
    mutable double  deltaInv_;
		mutable  double epsilon_;
    mutable double mu1_, mu2_;
    mutable double gamma_;
		};



#endif // file define

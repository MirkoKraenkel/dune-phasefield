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

#include "thermodynamicsinterface.hh"

class PGThermodynamics:
	public Thermodynamics<PGThermodynamics>
{
	typedef Thermodynamics<PGThermodynamics> BaseType;

public:
  PGThermodynamics():
    BaseType()
  {
    init();
  }
  inline void init() const 
  {
    delta_=Dune::Fem::Parameter::getValue<double>("phasefield.delta");
    delta_inv_=1./delta_;
    epsilon_=Dune::Fem::Parameter::getValue<double>("phasefield.delta");
    std::cout<<"init epsilon="<<epsilon_<<"\n"; 
  }
  
  inline double helmholtz(double& rho,double& phi) const
  {
    return delta_inv_*2*phi*phi*(1-phi)*(1-phi);
  }
  
	

	inline double reactionSource(double& rho,double& phi) const
	{ 

    return delta_inv_*4*(2*phi*phi*phi-3*phi*phi+phi);
  }
  


	inline double chemicalPotential(double& rho,double& phi) const
	{
    return 0;
  }

	
	inline double  pressure( double& rho, double& phi) const
	{
    return -helmholtz(rho,phi);
  }
		

	inline double a(double rho,double phi) const
	{
    return 1.6;	
	}
  
 

  public:

  	inline double delta()const {return delta_;}
		inline double delta_inv()const {return delta_inv_;}
    inline double mu1() const{ return epsilon_;}
	 	inline double mu2()const {  return epsilon_;}
   
	private:
		mutable double  delta_;
    mutable double  delta_inv_;
		mutable  double epsilon_;

		};



#endif // file define

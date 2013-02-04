#ifndef THERMODYNAMICSPERFECTGAS2_HH
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
// see MininimizerSonntag.mv

#include "thermodynamicsinterface.hh"

class PG2Thermodynamics:
	public Thermodynamics<PG2Thermodynamics>
{
	typedef Thermodynamics<PG2Thermodynamics> BaseType;

public:
  PG2Thermodynamics():
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
		double t10;
		double t11;
		double t13;
		double t16;
		double t17;
		double t2;
		double t5;
		double t8;
		
		t2 = phi*phi;
		t5 = pow(1.0-phi,2.0);
		t8 = t2*t2;
		t10 = 6.0*t8*phi;
		t11 = 15.0*t8;
		t13 = 10.0*t2*phi;
		t16 = log(rho);
		t17 = rho*t16;
		return(2.0/delta_*t2*t5+(t10-t11+t13)*(-0.5*rho+0.15E1*t17)+(1.0-t10+t11-t13)*(0.16E1*t17+0.5-0.15E1*rho));

	}
  
	

	inline double reactionSource(double& rho,double& phi) const
	{ 
		double t1;
		double t11;
		double t16;
		double t18;
		double t19;
		double t3;
		double t4;
		double t7;
				
		t1 = 1/delta_;
		t3 = 1.0-phi;
		t4 = t3*t3;
		t7 = phi*phi;
		t11 = t7*t7;
		t16 = 30.0*t11-60.0*t7*phi+30.0*t7;
		t18 = log(rho);
		t19 = rho*t18;
		return(4.0*t1*phi*t4-4.0*t1*t7*t3+t16*(-0.5*rho+0.15E1*t19)-t16*(0.16E1*t19+0.5-0.15E1*rho));

	}
  


	inline double chemicalPotential(double& rho,double& phi) const
	{
		double t1;
		double t11;
		double t2;
		double t3;
		double t5;
						
		t1 = phi*phi;
		t2 = t1*t1;
		t3 = t2*phi;
		t5 = log(rho);
		t11 = t1*phi;
		return(0.54E1*t3-0.6*t3*t5-0.135E2*t2+0.15E1*t2*t5+9.0*t11-1.0*t11*t5+
					 0.16E1*t5+0.1);
	}

	
	inline double  pressure( double& rho, double& phi) const
	{
		double delta=delta_;
double t1;
  double t10;
  double t2;
  double t26;
  double t3;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = t2*phi;
    t10 = t1*phi;
    t26 = 6.0*t3*rho*delta-15.0*t2*rho*delta+10.0*t10*rho*delta+5.0*delta-16.0*
delta*rho-30.0*delta*t3+75.0*delta*t2-50.0*delta*t10+20.0*t1-40.0*t10+20.0*t2;
    return(-0.1*t26/delta);
  }

   

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

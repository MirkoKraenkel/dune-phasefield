#ifndef SIMPLETHERMODYNAMICS_HH
#define SIMPLETHERMODYNAMICS_HH

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

class SimpleThermodynamics:
public Thermodynamics<SimpleThermodynamics>
{
 typedef Thermodynamics<SimpleThermodynamics> BaseType;

public:
  SimpleThermodynamics():
    BaseType(),
    delta_ (Fem::Parameter::getValue< double >( "phasefield.delta" )),
    epsilon_ (Fem::Parameter::getValue< double >( "epsilon",1. )),
    b_(Fem::Parameter::getValue<double>("pressure")),
    gamma_(Fem::Parameter::getValue<double>("gamma"))
  {
    delta_inv_     = 1. / delta_;
	}

  inline void init() const {std::cout<<"init me!\n";abort();}
  inline double helmholtz(double& rho,double& phi) const
  {
		return 0.;
	}



	

  inline double reactionSource(double& rho,double& phi) const
  { 
		
		double t1;
		double t10;
		double t2;
		double t4;
		double t9;
		
		t1 = phi-1.0;
		t2 = t1*t1;
		t4 = pow(rho-c2_,2);
		t9 = pow(rho-c1_,2);
		t10 = phi*phi;
		return(2.0*phi*(t2+t4)+2.0*(t9+t10)*t1);
	}
  


  inline double chemicalPotential(double& rho,double& phi) const
  {
		double t1;
		double t10;
		double t4;
		double t5;
		double t6;
		double t9;
		
		t1 = rho-c1_;
		t4 = pow(phi-1.0,2.0);
		t5 = rho-c2_;
		t6 = t5*t5;
		t9 = t1*t1;
		t10 = phi*phi;
		return(2.0*t1*(t4+t6)+2.0*(t9+t10)*t5);
		
	}

	inline double  pressure( double& rho, double& phi) const
  {
		double t1;
		double t10;
		double t11;
		double t4;
		double t5;
		double t6;
		double t7;
		double t9;
		
    t1 = rho-c1_;
    t4 = pow(phi-1.0,2.0);
    t5 = rho-c2_;
    t6 = t5*t5;
    t7 = t4+t6;
    t9 = t1*t1;
    t10 = phi*phi;
    t11 = t9+t10;
    return(rho*(2.0*t1*t7+2.0*t11*t5)-t11*t7);
  
	}


 

public:

  inline double delta()const {return delta_;}
  inline double delta_inv(){return delta_inv_;}
  inline double c1(){return c1_;}
  inline double c2(){return c2_;}
private:
  double c1_; // Parameter for rho minima of double well
  double c2_;
  const double delta_; double epsilon_;
 
  const double b_;
  const double gamma_;
  double delta_inv_;

};



#endif // file define

#ifndef THERMODYNAMICS_HH
#define THERMODYNAMICS_HH

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
    double t2 = phi*phi;
    double t4 = 1.0-phi;
    double t5 = t4*t4;
    double t9 = pow(rho-1.0,2.0);
    double t12 = pow(rho-2.0,2.0);
    double t13 = t12*t12;
  
    return(2.0/delta_*t2*t5+phi*t9+t4*t13);
  }
  
	

	inline double reactionSource(double& rho,double& phi) const
	{ 
    double t1 = 1/delta_;
    double t3 = 1.0-phi;
    double t4 = t3*t3;
    double t7 = phi*phi;
    double t12 = pow(rho-1.0,2.0);
    double t14 = pow(rho-2.0,2.0);
    double t15 = t14*t14;
   
    return(4.0*t1*phi*t4-4.0*t1*t7*t3+t12-t15);
  }
  


	inline double chemicalPotential(double& rho,double& phi) const
	{
   double t4 = rho*rho;
   double t5 = t4*rho;
   
   return(30.0*phi-46.0*phi*rho+4.0*t5-24.0*t4+48.0*rho-32.0-4.0*t5*phi+24.0*phi*t4);
  }

	
  inline double  pressure( double& rho, double& phi) const
	{
    double t1 = phi*delta;
    double t2 = rho*rho;
    double t6 = t2*t2;
    double t9 = t2*rho;
    double t19 = phi*phi;
    double t23 = t19*t19;
    double t25 = 23.0*t1*t2-15.0*t1-3.0*delta*t6+16.0*delta*t9-24.0*delta*t2+16.0*delta+3.0*t1*t6-16.0*t1*t9+2.0*t19-4.0*t19*phi+2.0*t23;
    return(-t25/delta);

  }
		

	inline double a(double rho,double phi) const
	{
    return 1.6;	
	}
  
  inline double acResiduum(double x, double t) const
  { 
    double t10;
    double t2;
    double t3;
    double t5;
    double t6;
    
    t2 = 0.6283185308E1*x+t;
    t3 = sin(t2);
    t5 = cos(t2);
    t6 = t5*t5;
    t10 = t6*t6;
  
     return(0.2973920881E1*t3+0.11875E1-0.125*t6-10.0*t3*t6-0.625E-1*t10);

  }

  inline double momResiduum(double x, double t) const
  { 
    
    double t13;
    double t17;
    double t22;
    double t28;
    double t3;
    double t31;
    double t35;
    double t36;
    double t38;
    double t4;
    double t5;
    double t7;
    double t8;
    t3 = 2.0*0.3141592653589793E1*x+t;
    t4 = sin(t3);
    t5 = 0.5*t4;
    t7 = cos(t3);
    t8 = t7*0.3141592653589793E1;
    t13 = t4/2.0+1.0/2.0;
    t17 = t4*t4;
    t22 = 1.0/2.0-t4/2.0;
    t28 = t22*t22;
    t31 = t13*t13;
    t35 = pow(t5+1.0,2.0);
    t36 = t17*t17;
    t38 = 0.3141592653589793E1*0.3141592653589793E1;
    return((2.0+t5)*(t8*(2.0+0.1E1*t4)+0.2E1*t13*t7*0.3141592653589793E1-0.5*t8
           *t17*t4+0.3E1*t22*t17*t8)-t8*(0.4E2*t13*t28-0.4E2*t31*t22+t35-0.625E-1*t36+0.2* t4*t38));
  
  }
  public:

  	inline double delta() const { return delta_; }
		inline double deltaInv() const { return deltaInv_; }
    inline double mu1() const { return mu1_; }
	 	inline double mu2() const { return mu2_; }
   
	private:
		mutable double  delta_;
    mutable double  deltaInv_;
		mutable  double epsilon_;
    mutable double mu1_, mu2_;
		};



#endif // file define

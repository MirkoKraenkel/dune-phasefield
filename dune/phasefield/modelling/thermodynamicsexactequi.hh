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
#include "thermodynamicsinterface.hh"
//thermodynamics with balanced phases, see balancedthermo_correctionterm.mw
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
   mu2_( epsilon_ )
    {
    }

  inline void init() const 
  {
    std::cout<<"thermo.init()\n";
  abort();
  }
  
  inline double helmholtz(double& rho,double& phi) const
  {
		double delta=delta_;
    double t10;
    double t11;
		double t12;
		double t13;
		double t14;
		double t2;
		double t25;
		double t27;
		double t28;
		double t32;
		double t38;
		double t39;
		double t5;
		double t8;
		double t9;
		{
			t2 = phi*phi;
			t5 = pow(1.0-phi,2.0);
			t8 = t2*t2;
			t9 = t8*phi;
			t10 = 6.0*t9;
			t11 = 15.0*t8;
			t12 = t2*phi;
			t13 = 10.0*t12;
			t14 = t10-t11+t13;
			t25 = exp((4.0-18.0*t9+45.0*t8-30.0*t12)/(-0.9E1*t9+0.225E2*t8-0.15E2*t12+3.0));
			t27 = log(t25);
			t28 = t25*t27;
			t32 = 1.0-t10+t11-t13;
			t38 = log(rho);
			t39 = rho*t38;
			return(2.0/delta*t2*t5-t14*(-0.25E1*t25+0.15E1*t28+1.0)-t32*(-7.0*t25+3.0*t28+0.945940263E1)+t14*(-0.25E1*rho+0.15E1*t39+1.0)
          +t32*(-7.0*rho+3.0*t39+0.945940263E1));
		}
	}


	

	inline double reactionSource(double& rho,double& phi) const
	{ 
    double delta=delta_;
  
    double t1;
		double t11;
		double t13;
		double t16;
		double t17;
		double t21;
		double t25;
		double t26;
		double t28;
		double t3;
		double t30;
		double t31;
		double t35;
		double t36;
		double t37;
		double t4;
		double t44;
		double t53;
		double t55;
		double t59;
		double t7;
		double t70;
		double t71;
		{
			t1 = 1/delta;
			t3 = 1.0-phi;
			t4 = t3*t3;
			t7 = phi*phi;
			t11 = t7*t7;
			t13 = t7*phi;
			t16 = 30.0*t11-60.0*t13+30.0*t7;
			t17 = t11*phi;
			t21 = 4.0-18.0*t17+45.0*t11-30.0*t13;
			t25 = -0.9E1*t17+0.225E2*t11-0.15E2*t13+3.0;
			t26 = 1/t25;
			t28 = exp(t21*t26);
			t30 = log(t28);
			t31 = t28*t30;
			t35 = 6.0*t17;
			t36 = 15.0*t11;
			t37 = 10.0*t13;
			t44 = t25*t25;
			t53 = ((-90.0*t11+180.0*t13-90.0*t7)*t26-t21/t44*(-0.45E2*t11+0.9E2*t13-0.45E2*t7))*t28;
			t55 = t53*t30;
			t59 = -t16;
			t70 = log(rho);
			t71 = rho*t70;
			return(4.0*t1*phi*t4-4.0*t1*t7*t3-t16*(-0.25E1*t28+0.15E1*t31+1.0)
          -(t35-t36+t37)*(-0.1E1*t53+0.15E1*t55)-t59*(0.945940263E1-7.0*t28+3.0*t31)
          -(1.0-t35+t36-t37)*(-4.0*t53+3.0*t55)+t16*(1.0-0.25E1*rho+0.15E1*t71)+t59*(0.945940263E1-7.0*rho+3.0*t71));
		}
 

	  
	}
  


	inline double chemicalPotential(double& rho,double& phi) const
	{
    double t1;
    double t11;
    double t2;
    double t3;
    double t5;
    {
      t1 = phi*phi;
      t2 = t1*t1;
      t3 = t2*phi;
      t5 = log(rho);
      t11 = t1*phi;
      return(18.0*t3-9.0*t3*t5-45.0*t2+0.225E2*t2*t5+30.0*t11-15.0*t11*t5-4.0+3.0
						 *t5);
    }  
  
	}

	
	inline double  pressure( double& rho, double& phi) const
	{
    double delta=delta_;

		double t1;
		double t12;
		double t13;
		double t2;
		double t23;
		double t4;
		double t8;
		{
			t1 = phi*phi;
			t2 = t1*t1;
			t4 = t2*phi*delta;
			t8 = t2*delta;
			t12 = t1*phi;
			t13 = t12*delta;
			t23 = 900000000.0*t4*rho-5075641578.0*t4-2250000000.0*t8*rho+
				0.1268910394E11*t8+1500000000.0*t13*rho-8459402630.0*t13-300000000.0*delta*rho+
				945940263.0*delta+200000000.0*t1-400000000.0*t12+200000000.0*t2;
			return(-0.1E-7*t23/delta);
		}

  
  
  
  
  
  
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


private:
	mutable double  delta_;
	mutable double  deltaInv_;
	mutable  double epsilon_;
  mutable double mu1_,mu2_;
};



#endif // file define

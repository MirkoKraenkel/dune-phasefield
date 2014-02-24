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
// see MininimizerSonntag.mv

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
   alpha_( Dune::Fem::Parameter::getValue<double>( "phasefield.alpha" ) )
  {
  }

  inline void init() const 
  {
    abort();
    deltaInv_=1./delta_;
    epsilon_=Dune::Fem::Parameter::getValue<double>("phasefield.epsilon");
  }
  
  inline double helmholtz(double& rho,double& phi) const
  {
    double delta=delta_;
    double alpha=alpha_;
  
    double t15;
    double t16;
    double t17;
    double t2;
    double t20;
    double t21;
    double t3;
    double t6;
    {
      t2 = phi*phi;
      t3 = t2*t2;
      t6 = t2*phi;
      t15 = 6.0*t3*phi;
      t16 = 15.0*t3;
      t17 = 10.0*t6;
      t20 = log(rho);
      t21 = rho*t20;
      return(1/delta*(2.0*t3+2.0*(2.0*alpha-2.0)*t6+2.0*(-3.0*alpha+1.0)*t2+2.0*alpha)
      +(t15-t16+t17)*(-0.25E1*rho+0.15E1*t21+1.0)+(1.0-t15+t16-t17)*(-7.0*rho+3.0*t21+0.945940263E1));
    }	

  }
  
	

  inline double reactionSource(double& rho,double& phi) const
    { 
      double delta=delta_;
      double alpha=alpha_;
      double t15;
      double t19;
      double t2;
      double t21;
      double t22;
      double t3;
      {
        t2 = phi*phi;
        t3 = t2*phi;
        t15 = t2*t2;
        t19 = 30.0*t15-60.0*t3+30.0*t2;
        t21 = log(rho);
        t22 = rho*t21;

        return(1/delta*(8.0*t3+6.0*(2.0*alpha-2.0)*t2+4.0*(-3.0*alpha+1.0)*phi)+t19
          *(-0.25E1*rho+0.15E1*t22+1.0)-t19*(-7.0*rho+3.0*t22+0.945940263E1));
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
      double alpha=alpha_;
      double t1;
  double t10;
  double t14;
  double t2;
  double t29;
  double t3;
  {
    t1 = phi*phi;
    t2 = t1*t1;
    t3 = t2*phi;
    t10 = t1*phi;
    t14 = p0*delta;
    t29 = -6.0*t3*p1*delta+15.0*t2*p1*delta-10.0*t10*p1*delta-t14+6.0*t14*t3
-15.0*t14*t2+10.0*t14*t10+2.0*t2+4.0*t10*alpha-4.0*t10-6.0*t1*alpha+2.0*t1+2.0*
alpha;
    return(-t29/delta);
  }


 
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
  mutable double  alpha_;

};



#endif // file define

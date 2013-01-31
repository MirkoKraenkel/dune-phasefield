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
   return

     2.0/delta_*phi*phi*pow(1.0-phi,2.0)+(6.0*phi*phi*phi*phi*phi-15.0*phi
*phi*phi*phi+10.0*phi*phi*phi)*(-0.4*rho+0.14E1*rho*log(rho)+1.0)+(1.0-6.0*phi*
phi*phi*phi*phi+15.0*phi*phi*phi*phi-10.0*phi*phi*phi)*(0.16E1*rho*log(rho)+
0.15E1-0.15E1*rho);

 
  }

	

	inline double reactionSource(double& rho,double& phi) const
	{ 
		double s=delta_inv_;
		double t11;
		double t16;
		double t17;
		double t18;
		double t2;
		double t3;
		double t6;
		{
			t2 = 1.0-phi;
			t3 = t2*t2;
			t6 = phi*phi;
			t11 = t6*t6;
			t16 = 30.0*t11-60.0*t6*phi+30.0*t6;
			t17 = log(rho);
			t18 = rho*t17;
			return(rho*s*(2.0*phi*t3-2.0*t6*t2-0.1E-1)+t16*(0.15E1*t18+0.1693147181E1*
																											rho)-t16*(1.0*t18-rho+0.5));
		}

	}
  


	inline double chemicalPotential(double& rho,double& phi) const
	{
	  double s=delta_inv_;
		double t1;
		double t10;
		double t11;
		double t13;
		double t15;
		double t3;
		double t8;
		{
			t1 = phi*phi;
			t3 = pow(1.0-phi,2.0);
			t8 = t1*t1;
			t10 = 6.0*t8*phi;
			t11 = 15.0*t8;
			t13 = 10.0*t1*phi;
			t15 = log(rho);
			return(s*(t1*t3-0.1E-1*phi+0.1E-1)+(t10-t11+t13)*(0.15E1*t15+0.3193147181E1
																												)+1.0*(1.0-t10+t11-t13)*t15);

		}
	}
	inline double  pressure( double& rho, double& phi) const
	{
		double t1;
		double t2;
		double t3;
		double t8;
		{
			t1 = phi*phi;
			t2 = t1*t1;
			t3 = t2*phi;
			t8 = t1*phi;
			return(3.0*t3*rho-0.75E1*t2*rho+5.0*t8*rho-0.5+rho+3.0*t3-0.75E1*t2+5.0*t8);
			
		}
	}
	inline double a(double rho,double phi) const
	{
    return 1.5;	
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

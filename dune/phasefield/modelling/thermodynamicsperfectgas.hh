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
    double s=delta_inv_; 
    double t1;
    double t12;
    double t13;
    double t14;
    double t18;
    double t2;
    double t22;
    double t4;
    double t7;
    {
      t1 = rho*s;
      t2 = phi*phi;
      t4 = t2*phi;
      t7 = t2*t2;
      t12 = t7*phi;
      t13 = t12*rho;
      t14 = log(rho);
      t18 = t7*rho;
      t22 = t4*rho;
    
      return(t1*t2-2.0*t1*t4+t1*t7-0.1E-1*t1*phi+0.1E-1*t1+3.0*t13*t14+
0.1615888309E2*t13-0.75E1*t18*t14-0.4039720772E2*t18+5.0*t22*t14+0.2693147181E2
*t22+rho*t14-1.0*rho+0.5-3.0*t12+0.75E1*t7-5.0*t4);
	
 	}

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

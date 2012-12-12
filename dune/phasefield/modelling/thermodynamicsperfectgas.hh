
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
		double s=delta_inv_;
	  double t1;
		double t11;
		double t12;
		double t14;
		double t15;
		double t17;
		double t19;
		double t20;
		double t3;
		double t5;
		{
			t1 = 0.9834594901*phi;
			t3 = pow(t1+0.6055746688E-1,2.0);
			t5 = pow(0.9394425331-t1,2.0);
			t11 = phi*phi;
			t12 = t11*t11;
			t14 = 6.0*t12*phi;
			t15 = 15.0*t12;
			t17 = 10.0*t11*phi;
			t19 = log(rho);
			t20 = rho*t19;
			return(s*(t3*t5-0.9834594901E-1*phi+0.9394425331E-1)*rho+(t14-t15+t17)
						 *(0.15E1*t20+0.6931471806*rho+1.0)+(1.0-t14+t15-t17)*(t20+1.0));

		}

	}

	

	inline double reactionSource(double& rho,double& phi) const
	{ 
		double s=delta_inv_;
		double t1;
		double t13;
		double t14;
		double t19;
		double t2;
		double t20;
		double t21;
		double t3;
		double t4;
		double t7;
		{
			t1 = 0.9834594901*phi;
			t2 = t1+0.6055746688E-1;
			t3 = 0.9394425331-t1;
			t4 = t3*t3;
			t7 = t2*t2;
			t13 = phi*phi;
			t14 = t13*t13;
			t19 = 30.0*t14-60.0*t13*phi+30.0*t13;
			t20 = log(rho);
			t21 = rho*t20;
			return(s*(0.196691898E1*t2*t4-0.196691898E1*t7*t3-0.9834594901E-1)*rho+t19*
						 (0.15E1*t21+0.6931471806*rho+1.0)-t19*(t21+1.0));
		}
	
	}
  


	inline double chemicalPotential(double& rho,double& phi) const
	{
    double s=delta_inv_;
		double t1;
		double t10;
		double t11;
		double t13;
		double t14;
		double t16;
		double t18;
		double t3;
		double t5;
		{
			t1 = 0.9834594901*phi;
			t3 = pow(t1+0.6055746688E-1,2.0);
			t5 = pow(0.9394425331-t1,2.0);
			t10 = phi*phi;
			t11 = t10*t10;
			t13 = 6.0*t11*phi;
			t14 = 15.0*t11;
			t16 = 10.0*t10*phi;
			t18 = log(rho);
			return(s*(t3*t5-0.9834594901E-1*phi+0.9394425331E-1)+(t13-t14+t16)
						 *(0.15E1*t18+0.2193147181E1)+(1.0-t13+t14-t16)*(t18+1.0));
		}
	
	}

	inline double  pressure( double& rho, double& phi) const
	{
		double s=delta_inv_;
		double t1;
		double t11;
		double t12;
		double t14;
		double t15;
		double t17;
		double t18;
		double t19;
		double t20;
		double t25;
		double t3;
		double t5;
		double t9;
		{
			t1 = 0.9834594901*phi;
			t3 = pow(t1+0.6055746688E-1,2.0);
			t5 = pow(0.9394425331-t1,2.0);
			t9 = s*(t3*t5-0.9834594901E-1*phi+0.9394425331E-1);
			t11 = phi*phi;
			t12 = t11*t11;
			t14 = 6.0*t12*phi;
			t15 = 15.0*t12;
			t17 = 10.0*t11*phi;
			t18 = t14-t15+t17;
			t19 = log(rho);
			t20 = rho*t19;
			t25 = 1.0-t14+t15-t17;
			return(-t9*rho-t18*(0.15E1*t20+0.6931471806*rho+1.0)-t25*(t20+1.0)
						 -rho*(t9+t18*(0.15E1*t19+0.2193147181E1)+t25*(t19+1.0)));
		}
	
  
	}

 inline double a(double rho,double phi) const
  {
    return 1.5; 
  }
  
 

public:

	inline double delta()const {return delta_;}
	inline double delta_inv()const {return delta_inv_;}
	inline double mu1() const{return epsilon_;}
	inline double mu2()const {return epsilon_;}
private:
	double c1_; // Parameter for rho minima of double well
	double c2_;
	const double delta_; double epsilon_;
 
	const double b_;
	const double gamma_;
	double delta_inv_;

};



#endif // file define

#ifndef DUNE_THERMODYNAMICS_HH
#define DUNE_THERMODYNAMICS_HH

// system include
#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>

// dune-fem include
#include <dune/fem/io/parameter.hh>

using namespace Dune;
// Thermodynamics
// --------------

/** \class Thermodynamics
 *  \brief deals with physics used for the atmosphere, in particular
 *         physical constants, equation of state etc.
 *
 *  \tparam dimDomain dimension of the domain
 */
template< int dimDomain >
class Thermodynamics
{
  enum{ energyId = dimDomain+1 };

public:
  Thermodynamics():
    delta_ (Parameter::getValue< double >( "femhowto.delta" )),
    epsilon_ (Parameter::getValue< double >( "epsilon",1. )),
    b_(Parameter::getValue<double>("pressure"))
  {
    delta_inv_     = 1. / delta_;
  }

  inline void init();
  /** \brief calculate the pressure and the temperature assuming the energy form
   *         for conservative variables: \f$[\rho,\rho\boldsymbol{v},\rho e]\f$
   *
   *  Calculate the pressure and temperature from the conservative variables
   *  \f$[\rho,\rho\boldsymbol{v},\rho e]\f$, where \f$ e \f$ is the sum of 
   *  the internal and the kinetic energy.
   *
   *  \param[in] cons Conervative variables
   *  \param[out] p Pressure W-Wrho
   *  \param[out] Wphi Source Part of Phasefield eq.  
   *
   *  \tparam RangeType Type of the range value
   */
  

  inline double Sopt(double& rho,double& phi) const
  { 
    double a=delta_inv_;
    double b=b_;
    double res;
    double t1;
    double t13;
    double t14;
    double t19;
    double t2;
    double t22;
    double t27;
    double t3;
    double t4;
    double t7;
    
    t1 = 0.9834594901*phi;
    t2 = t1+0.6055746688E-1;
    t3 = 0.9394425331-t1;
    t4 = t3*t3;
    t7 = t2*t2;
    t13 = phi*phi;
    t14 = t13*t13;
    t19 = 30.0*t14-60.0*t13*phi+30.0*t13;
    t22 = log(0.6666666667*rho);
    t27 = log(rho);
    res=(a*(0.196691898E1*t2*t4-0.196691898E1*t7*t3-0.9834594901E-1)*rho+b
	 *(2.0*t19*rho*t22-0.12E1*t19*rho*t27));
    
    //    res*=delta_inv_;
    return res; 
    
  }
  


  inline double Popt(double& rho,double& phi)
    const
  {
    double a=delta_inv_;
    double b=b_;
    double res;
    double t1;
    double t11;
    double t12;
    double t13;
    double t14;
    double t15;
    double t16;
    double t17;
    double t18;
    double t21;
    double t24;
    double t26;
    double t3;
    double t5;
    double t9;
  
    t1 = 0.9834594901*phi;
    t3 = pow(t1+0.6055746688E-1,2.0);
    t5 = pow(0.9394425331-t1,2.0);
    t9 = a*(t3*t5-0.9834594901E-1*phi+0.9394425331E-1);
    t11 = phi*phi;
    t12 = t11*t11;
    t13 = t12*phi;
    t14 = 6.0*t13;
    t15 = 15.0*t12;
    t16 = t11*phi;
    t17 = 10.0*t16;
    t18 = t14-t15+t17;
    t21 = log(0.6666666667*rho);
    t24 = 1.0-t14+t15-t17;
    t26 = log(rho);
    
    res=(-t9*rho-b*(2.0*t18*rho*t21+0.12E1*t24*rho*t26)+rho*
	 (t9+b*(2.0*t18*t21+0.48E1*t13-0.12E2*t12+0.8E1*t16+0.12E1+0.12E1*t24*t26)));
    return res;
  
  }




  template< class RangeType >
  inline void pressAndTempEnergyForm( const RangeType& cons, 
                                      double& p, double& Wphi ) const
  {
    assert( cons[0] > 1e-20 );
    ///double delta=delta_;
    double rho=cons[0];
    double phi=cons[dimDomain+1];
    phi/=rho;
    
    p=Popt(rho,phi);
    
    Wphi=0.;
    
    Wphi=Sopt(rho,phi); 
    Wphi*=-1.;




  
   //  assert( p > 1e-20 );
    // assert( T > 1e-20 );
  }

  /** \brief calculate the pressure and the temperature assuming the theta form
   *         for conservative variables: \f$[\rho,\rho\boldsymbol{v},\rho\theta]\f$
   *
   *  \param[in] cons Conervative variables
   *  \param[out] p Pressure
   *  \param[out] T temperature
   *
   *  \tparam RangeType Type of the range value
   */
  template< class RangeType >
  inline void pressAndTempThetaForm( const RangeType& cons, 
                                     double& p, double& T ) const
  {
    // cons = [rho, rho*v, rho*theta]  
    assert( cons[0] > 1e-20 );
    assert( cons[energyId] > 1e-20 ); 

    
   //  std::cout<<"pressure="<<p<<"\n";
    assert( p > 1e-20 );
    assert( T > 1e-20 );
  }

  template< class RangeType >
  void conservativeToPrimitiveThetaForm( const RangeType& cons, RangeType& prim ) const;

  template< class RangeType >
  void conservativeToPrimitiveEnergyForm( const RangeType& cons, RangeType& prim ) const;

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
  double delta_inv_;

};



template< int dimDomain >
inline void Thermodynamics< dimDomain >
:: init()
{
  c1_      = Parameter::getValue< double >( "c1",0.1 );
  c2_      = Parameter::getValue< double >( "c2",0.9 );
  delta_   = Parameter::getValue< double >( "femhowto.delta" );
  epsilon_ = Parameter::getValue< double >( "epsilon",1. );
  b_= Parameter::getValue< double >( "pessure" );
  delta_inv_     = 1. / delta_;
  
}



/** \brief converts conservative variables in the energy form to primitive ones
 *
 *  Converts conservative variables \f$[\rho\boldsymbol{v},p,\rho e\f$
 *  to primitive ones \f$[\boldsymbol{v},p,\theta]\f$, where \f$e\f$ is the sum
 *  of internal and kinetic energy, and \f$\theta\f$ potential temperature.
 *
 *  \param[in] cons Conservative variables
 *  \param[out] prim Primitive variables
 *
 *  \tparam dimDomain dimension of the domain
 *  \tparam RangeType type of the range value
 */
template< int dimDomain >
template< class RangeType >
void Thermodynamics< dimDomain >
:: conservativeToPrimitiveEnergyForm( const RangeType& cons, RangeType& prim ) const
{
  std::cerr <<"conservativeToPrimitiveEnergyForm not implemented" <<std::endl;
  abort();

  //const double rho_inv = 1./ cons[0];

  //double p, T;
  //pressAndTempEnergyForm( cons, p, T );

  //prim[energyId-1] = p;
  // this is not pot. temp !!!!!!!
  //prim[energyId] = cons[energyId]/cons[0];
}


/** \brief converts conservative variables in the theta form to primitive ones
 *
 *  Converts conservative variables \f$[\rho\boldsymbol{v},p,\rho\theta]\f$
 *  to primitive ones \f$[\boldsymbol{v},p,\theta]\f$, where \f$\theta\f$ is
 *  potential temperature
 *
 *  \param[in] cons Conservative variables
 *  \param[out] prim Primitive variables
 *
 *  \tparam dimDomain dimension of the domain
 *  \tparam RangeType type of the range value
 */
template< int dimDomain >
template< class RangeType >
void Thermodynamics< dimDomain >
:: conservativeToPrimitiveThetaForm( const RangeType& cons, RangeType& prim ) const
{
  
  
   assert( cons[0] > 0. );
  
  double p, T;
  pressAndTempEnergyForm( cons, p, T );

  for( int i = 0; i < dimDomain; ++i )
    prim[i] = cons[i+1]/cons[0];

  prim[energyId-1] = p;
  prim[energyId] = cons[energyId] / cons[0];
}

#endif // file define

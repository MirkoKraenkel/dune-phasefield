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
    epsilon_ (Parameter::getValue< double >( "epsilon",1. ))
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

  inline double Sopt( double& rho, double& phi) const
  {
    double t1;
    double t10;
    double t2;
    double t4;
    double t9;
    {
      t1 = phi-1.0;
      t2 = t1*t1;
      t4 = pow(rho-c2_,2);
      t9 = pow(rho-c1_,2);
      t10 = phi*phi;
      return(2.0*phi*(t2+t4)+2.0*(t9+t10)*t1);
    }
  }


  inline double  Popt( double& rho, double& phi) const
  {
  double t1;
  double t10;
  double t11;
  double t4;
  double t5;
  double t6;
  double t7;
  double t9;
  {
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
}

  template< class RangeType >
  inline void pressAndTempEnergyForm( const RangeType& cons, 
                                      double& p, double& Wphi ) const
  {
    assert( cons[0] > 1e-20 );
    double delta=delta_;
    double rho=cons[0];
    double phi=cons[dimDomain+1];
    phi/=rho;


  //   p=4.0*c1_*rho*rho*c2_-phi*phi*phi*phi+2.0*phi*phi*phi-phi*phi+rho*rho-c1_
//       *c1_+2.0*rho*rho*phi*phi-2.0*rho*rho*phi-4.0*rho*rho*rho*c2_+rho*rho*c2_*c2_-4.0*c1_
//       *rho*rho*rho+c1_*c1_*rho*rho-c1_*c1_*phi*phi+2.0*c1_*c1_*phi-c1_*c1_*c2_*c2_-phi*phi*c2_*
//       c2_+3.0*rho*rho*rho*rho;
  
    p=Popt(rho,phi);
    p*=delta_inv_;

      

    Wphi=0.;
    
   //  Wphi=2.0*phi*(pow(phi-1.0,2.0)+pow(rho-c2_,2.0))+(pow(rho-c1_,2.0)+phi*phi)
//       *(2.0*phi-2.0);
    Wphi=Sopt(rho,phi);
         Wphi*=delta_inv_;

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
  const double delta_;
  double delta_inv_;
  double epsilon_;
};



template< int dimDomain >
inline void Thermodynamics< dimDomain >
:: init()
{
  c1_      = Parameter::getValue< double >( "c1",0.1 );
  c2_      = Parameter::getValue< double >( "c2",0.9 );
  delta_   = Parameter::getValue< double >( "femhowto.delta" );
  epsilon_ = Parameter::getValue< double >( "epsilon",1. );

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
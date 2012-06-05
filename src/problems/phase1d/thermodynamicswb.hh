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
// see simplehelmholtz.mv

template< int dimDomain >
class Thermodynamics
{
  enum{ phaseId = dimDomain+1 };

public:
  Thermodynamics():
    delta_ (Parameter::getValue< double >( "femhowto.delta" )),
    epsilon_ (Parameter::getValue< double >( "epsilon",1. )),
    b_(Parameter::getValue<double>("pressure")),
    gamma_(Parameter::getValue<double>("gamma"))
  {
    delta_inv_     = 1. / delta_;
    
  }

  inline void init();
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
  


  inline double chemicalPotential(double& rho,double& phi)
    const
  {
	 double t1;
  double t10;
  double t4;
  double t5;
  double t6;
  double t9;
  {
    t1 = rho-c1_;
    t4 = pow(phi-1.0,2.0);
    t5 = rho-c2_;
    t6 = t5*t5;
    t9 = t1*t1;
    t10 = phi*phi;
    return(2.0*t1*(t4+t6)+2.0*(t9+t10)*t5);
  }
	
			
	
	}




  template< class RangeType >
  inline void chemPotAndReaction( const RangeType& cons, 
																	double& mu,double& reaction ) const
  {
    assert( cons[0] > 1e-20 );
  
		double rho=cons[0];
    double phi=cons[dimDomain+1];
    
		phi/=rho;
		
    mu=chemicalPotential(rho,phi);
		reaction=reactionSource(rho,phi); 
    reaction*=-1.;
  
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
    assert( cons[phaseId] > 1e-20 ); 

    
   //  std::cout<<"pressure="<<p<<"\n";
    assert( p > 1e-20 );
    assert( T > 1e-20 );
  }

  template< class RangeType >
  void conservativeToPrimitiveThetaForm( const RangeType& cons, RangeType& prim ) const;
  
  template< class RangeType,class GradRangeType >
  void totalEnergy( const RangeType& cons, const GradRangeType& grad, double& res ) const;
  
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
  const double gamma_;
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

  //prim[phaseId-1] = p;
  // this is not pot. temp !!!!!!!
  //prim[phaseId] = cons[phaseId]/cons[0];
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
 //  pressAndTempEnergyForm( cons, p, T );

  for( int i = 0; i < dimDomain; ++i )
    prim[i] = cons[i+1]/cons[0];

  prim[phaseId-1] = p;
  prim[phaseId] = cons[phaseId] / cons[0];
}


template< int dimDomain >
template< class RangeType,class GradRangeType >
void Thermodynamics< dimDomain >
:: totalEnergy( const RangeType& cons,const GradRangeType& grad,double& res ) const
{
	abort();
  assert( cons[0] > 0. );
  double rho = cons[0];
  double rho_inv = 1. /rho;
  double phi=cons[dimDomain+1]*rho_inv;
  
  const double v = cons[1]*rho_inv;
  double kin= cons[1]*v*0.5;
  const double dxrho     = grad[0]; //drho/dx

   const double dxrhophi  = grad[2]; //d(rho*phi)/dx
   
  const double rhodxphi = (dxrhophi - phi*dxrho);
  //double tens=gamma_*rhodxphi*rhodxphi*rho_inv;
  double tens=gamma_*delta_*rhodxphi*rhodxphi*rho_inv*rho_inv;
  
  //double tens=gamma_*rhodxphi*rhodxphi*rho_inv*rho_inv*rho_inv;
  double helm=helmholtz( rho,phi);
  res=helm+tens+kin;

}


#endif // file define

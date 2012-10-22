
#ifndef DUNEPHASEFIELD_WBPHYSICS_HH
#define DUNEPHASEFIELD_WBPHYSICS_HH

// system include
#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>

// dune-fem include
#include <dune/fem/io/parameter.hh>

#include "thermodynamics.hh"

using namespace Dune;

template<int dimDomain>
class PhasfieldPhysics
{
public:
	enum{ phaseId = dimDomain+1} ;
private:
	Thermodynamics& thermoDynamics_;

public:
	PhasefieldPhysics(const Thermodynamics& thermodyn):
		thermoDynamics_(thermodyn)
	{
	}
	
	 
  template< class RangeType >
  void conservativeToPrimitive( const RangeType& cons, RangeType& prim ) const;
  
  template< class RangeType,class GradRangeType >
  void totalEnergy( const RangeType& cons, const GradRangeType& grad, double& res ) const;
	




	
	template< class RangeType >
	inline void chemPotAndReaction( const RangeType& cons, 
																	double& mu,
																	double& reaction ) const
	{
		assert( cons[0] > 1e-20 );
		
		double rho=cons[0];
		double phi=cons[phaseId];
		phi/=rho;
		
		mu=thermoDynamics_.chemicalPotential(rho,phi);
		reaction=thermoDynamics_.reactionSource(rho,phi); 
		reaction*=-1.;
	}
	
	template< class RangeType >
	inline void pressureAndReaction( const RangeType& cons, 
																	 double& p,
																	 double& reaction ) const
	{
		assert( cons[0] > 1e-20 );
	  
		double rho=cons[0];
		double phi=cons[dimDomain+1];
		phi/=rho;
		
		p=thermodynamics_.pressure(rho,phi);
		reaction=thermodynamics_.reactionSource(rho,phi); 
		reaction*=-1.;
			
		assert( p > 1e-20 );
	}



	inline void nonConProduct(const RangeType & uL, 
														const RangeType & uR,
														const ThetaRangeType& thetaL,
														const ThetaRangeType& thetaR,
														RangeType& ret) const;


	
	template <class JacobianRangeImp>
	inline void diffusion( const RangeType& u,
												 const JacobianRangeImp& du,
												 JacobianRangeType& f ) const;
	
	
	inline void allenCahn( const RangeType& u,
												 const GradientRangeType& du,
												 ThetaJacobianRangeType& f ) const;
	
  
	template< class JacobianRangeImp >  
	inline void allenCahn(const RangeType& u,
												const JacobianRangeImp& du,
												ThetaJacobianRangeType& f ) const;

	

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
	void conservativeToPrimitive( const RangeType& cons, RangeType& prim ) const;
  
	template< class RangeType,class GradRangeType >
	void totalEnergy( const RangeType& cons, const GradRangeType& grad, double& res ) const;
  

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
	c1_      = Fem::Parameter::getValue< double >( "c1",0.1 );
	c2_      = Fem::Parameter::getValue< double >( "c2",0.9 );
	delta_   = Fem::Parameter::getValue< double >( "femhowto.delta" );
	epsilon_ = Fem::Parameter::getValue< double >( "epsilon",1. );
	b_= Fem::Parameter::getValue< double >( "pessure" );
	delta_inv_     = 1. / delta_;
  
}

template< int dimDomain >
template< class RangeType >
void Thermodynamics< dimDomain >
:: conservativeToPrimitive( const RangeType& cons, RangeType& prim ) const
{
  
  
	assert( cons[0] > 0. );
  
	double p,rho,phi;
	rho=cons[0];
	phi=cons[phaseId]/rho;
	p=pressure(rho,phi);
	

	for( int i = 0; i < dimDomain; ++i )
		prim[i] = cons[i+1]/cons[0];

	prim[phaseId-1] = p;
	prim[phaseId] = phi;
}


template< int dimDomain >
template< class RangeType,class GradRangeType >
void Thermodynamics< dimDomain >
:: totalEnergy( const RangeType& cons,const GradRangeType& grad,double& res ) const
{
	assert( cons[0] > 0. );
	double rho = cons[0];
	double rho_inv = 1. /rho;
	double phi=cons[dimDomain+1];
  
	const double v = cons[1]*rho_inv;
	double kin= cons[1]*v*0.5;
	const double dxrho     = grad[0]; //drho/dx
	const double dxphi     = grad[2]; //d(rho*phi)/dx
	 
	double tens=gamma_*delta_*dxphi*dxphi;
	double helm=helmholtz( rho,phi);
	res=helm+tens+kin;

}


#endif // DUNEPHASEFIELD_WBPHYSICS_HH

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

namespace Dune{

// 
template<int dimDomain>
class PhasefieldPhysics
{
  public:
	  enum{ phaseId = dimDomain+1} ;
    enum { dimRange = dimDomain + 2 };
    enum { dimThetaRange =  2 };
    enum { dimThetaGradRange = dimThetaRange*dimDomain };
    enum { dimGradRange = dimRange * dimDomain };
    enum { dimGradient = dimDomain + 1 };
    
    typedef double RangeFieldType;

    typedef FieldVector< double, dimDomain >                  DomainType;
    typedef FieldVector< double, dimDomain - 1 >              FaceDomainType;
    typedef FieldVector< double, dimRange >                   RangeType;
    typedef FieldVector< double, dimThetaRange >              ThetaRangeType;
    typedef FieldVector< double, dimGradRange >               GradientType;
    typedef FieldVector< double, dimThetaGradRange >          ThetaGradientRangeType;
    typedef FieldMatrix< double, dimRange, dimDomain >        JacobianRangeType;                          
    typedef FieldMatrix< double, dimRange, dimDomain >        FluxRangeType;
    typedef FieldVector< double, dimGradRange >               GradientRangeType;

  //  typedef FieldMatrix< double, dimRange, 
   typedef FieldMatrix< double, dimThetaRange, dimDomain >    ThetaJacobianRangeType;
    typedef FieldMatrix< double, dimGradRange, dimDomain >    JacobianFluxRangeType;

  private:
    const Thermodynamics& thermoDynamics_;
  public:
	PhasefieldPhysics(const Thermodynamics& thermodyn):
    thermoDynamics_(thermodyn),
    delta_(Dune::Fem::Parameter::getValue<double>("phasefield.delta")),
    delta_inv_(1./delta_)
  {
	}
	
  inline void analyticalFlux( const RangeType& u, JacobianRangeType& f ) const;
  
  inline void jacobian( const RangeType& u, JacobianFluxRangeType& a) const;

  inline double maxSpeed( const DomainType& n, const RangeType& u ) const;
  
  inline void conservativeToPrimitive( const RangeType& cons, RangeType& prim ) const;
  
  inline void totalEnergy( const RangeType& cons, 
                           const GradientRangeType& grad,
                          double& res ) const;

  inline void chemPotAndReaction( const RangeType& cons, 
																	double& mu,
																	double& reaction ) const;

	inline void pressureAndReaction( const RangeType& cons, 
																	 double& p,
																	 double& reaction ) const;

  inline void nonConProduct(const RangeType & uL, 
														const RangeType & uR,
														const ThetaRangeType& thetaL,
														const ThetaRangeType& thetaR,
														RangeType& ret) const;

  inline double stiffSource(const RangeType& u,
								            const GradientRangeType& du,
								            const ThetaRangeType& theta,
								            const ThetaJacobianRangeType& dtheta,
								            RangeType& f) const;
 
  template< class JacobianRangeImp >
	inline void diffusion( const RangeType& u,
												 const JacobianRangeImp& du,
												 JacobianRangeType& f ) const;
  
  template< class JacobianRangeImp >
	inline void allenCahn( const RangeType& u,
												 const JacobianRangeImp& du,
												 ThetaJacobianRangeType& f ) const;
	
  


public:

	inline double delta()const {return delta_;}
	inline double delta_inv(){return delta_inv_;}
private:
	const double delta_; 
	double delta_inv_;

};
#if 0
//================================================
//
//for all dim
//
//================================================
  template<int dimDomain >
  inline void PhasefieldPhysics< dimDomain >
  :: conservativeToPrimitive( const RangeType& cons, RangeType& prim ) const
  {
  	assert( cons[0] > 0. );
  
  	double p,rho,phi;
  	rho=cons[0];
  	phi=cons[phaseId]/rho;
  	for( int i = 0; i < dimDomain; ++i )
	  	prim[i] = cons[i+1]/cons[0];

  	prim[phaseId-1] = p;
  	prim[phaseId] = phi;
  
  }

  template< int dimDomain >
  inline void PhasefieldPhysics< dimDomain >
  :: totalEnergy( const RangeType& cons, 
                  const GradientRangeType& grad , 
                  double& res ) const
  {
  	assert( cons[0] > 0. );
	  double rho = cons[0];
	  double rho_inv = 1. /rho;
	  double phi = cons[dimDomain+1];
    FieldVector<RangeFieldType, dimDomain> v;
    double kineticEnergy,surfaceEnergy;
    for(int i = 0; i < dimDomain; i++)
    { 
      kineticEnergy+=cons[i]*cons[i];
      surfaceEnergy+=grad[phaseId][i]*grad[phaseId];
    }
  
    kineticEnergy*=0.5*rho_inv;
    surfaceEnergy*=delta_*0.5*rho_inv;


	  double freeEnergy = thermoDynamics_.helmholtz( rho, phi );

	  res = freeEnergy + surfaceEnergy + kineticEnergy;

  }

  template< int dimDomain >
  inline void PhasefieldPhysics< dimDomain >
  ::chemPotAndReaction( const RangeType& cons, 
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

  template< int dimDomain >
	inline void PhasefieldPhysics< dimDomain >
  ::pressureAndReaction( const RangeType& cons, 
                         double& p,
												 double& reaction ) const
	{
		assert( cons[0] > 1e-20 );
	  
		double rho=cons[0];
		double phi=cons[phaseId];
		phi/=rho;
		
		p=thermoDynamics_.pressure(rho,phi);
		reaction=thermoDynamics_.reactionSource(rho,phi); 
		reaction*=-1.;
			
		assert( p > 1e-20 );
	}



//================================================
//
//for dim=1
//
//================================================

  template< >
  inline void  PhasefieldPhysics<1>
  ::analyticalFlux( const RangeType& u, JacobianRangeType& f ) const
  {
		assert( u[0] > 1e-10 );
		double rho_inv = 1. / u[0];
		const double v = u[1]*rho_inv;
		double p;
		double W;
   
		JacobianRangeType phiT(0.);
		f[0][0] = u[1];
		f[1][0] = v*u[1];
		f[e][0] = u[phaseId]*v;

  }



	template<>
 	inline void   PhasePhysics<1>
	::nonConProduct(const RangeType & uL, 
									const RangeType & uR,
									const ThetaRangeType& thetaL,
									const ThetaRangeType& thetaR,
									RangeType& ret) const
	{
		ret[1]=0.5*(uL[0]+uR[0])*(thetaL[0]-thetaR[0]);
 	}
	
	template<>
	inline double PhasePhysics<1>
	::stiffSource(const RangeType& u,
								const GradientRangeType& du,
								const ThetaRangeType& theta,
								const ThetaJacobianRangeType& dtheta,
								RangeType& f) const
	{
			f[0]=0;
			f[1]=dtheta[0]*u[0]+du[2]*theta[1];
			f[2]=theta[1];

	
	}
#endif
}	//end namspace DUNE

#include "physics_inline.hh"









#endif // DUNEPHASEFIELD_WBPHYSICS_HH

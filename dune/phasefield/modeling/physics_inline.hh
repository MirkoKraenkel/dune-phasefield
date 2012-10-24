#ifndef PHYSICS_INLINE_HH
#define PHYSICS_INLINE_HH
namespace Dune{

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
		f[phaseId][0] = u[phaseId]*v;
  }

template<>
  inline void PhasefieldPhysics<1>
  ::jacobian( const RangeType& u, JacobianFluxRangeType& a) const
  {
    assert(u[0] > 1e-10);

    a[0][0] = u[0]; //rho
    a[1][0] = u[1]/u[0];//(rho v)/rho
    a[2][0] = u[2]/u[0];//(rho phi)/rho
  }

	template<>
 	inline void   PhasefieldPhysics<1>
	::nonConProduct(const RangeType & uL, 
									const RangeType & uR,
									const ThetaRangeType& thetaL,
									const ThetaRangeType& thetaR,
									RangeType& ret) const
	{
		ret[1]=0.5*(uL[0]+uR[0])*(thetaL[0]-thetaR[0]);
 	}
	
	template<>
	inline double PhasefieldPhysics<1>
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

  template<>
  template< class JacobianRangeImp >
  inline void PhasefieldPhysics<1>
  ::diffusion( const RangeType& u,
               const JacobianRangeImp& du,
               JacobianRangeType& diff) const
  {
    assert( u[0] > 1e-10 );
		double rho_inv = 1. / u[0];
		double p;
		double T;
		const double muLoc = 1;
  	const double dxv   = du[1][0]; //dv/dx
	
		diff[0][0]=0.;
		diff[2][0]=0.;
		diff[1][0]=muLoc*dxv;
   
  }
  template<>
  template< class JacobianRangeImp >
  inline void PhasefieldPhysics<1>
  ::allenCahn( const RangeType& u,
               const JacobianRangeImp& du,
               ThetaJacobianRangeType& diff ) const
  {
   assert( u[0] > 1e-10 );
	 
	 diff[0][0]=0.;
	 diff[1][0]=delta_*du[2][0];

  }

 template<>
 inline double PhasefieldPhysics<1>
 ::maxSpeed( const DomainType& n, const RangeType& u) const
 {
   abort();
  return 0.;
 } 



}//end namespace Dune
#endif// PHYSICS_INLINE_HH


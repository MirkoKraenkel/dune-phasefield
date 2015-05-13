#ifndef PHYSICSWB_INLINE2D_HH
#define PHYSICSWB_INLINE2D_HH
namespace Dune{
//================================================
//
//
//================================================


template<class Thermodynamics>
class PhasefieldPhysics<2,Thermodynamics>
{
   typedef Thermodynamics ThermodynamicsType;
 
  public:
    enum { dimDomain = 2 };
    enum { phaseId = dimDomain+1} ;
    enum { dimRange = dimDomain + 2 };
    enum { dimThetaRange =  2 };
    enum { dimThetaGradRange = dimThetaRange*dimDomain };
    enum { dimGradRange = dimRange * dimDomain };
    
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

    typedef FieldMatrix< double, dimThetaRange, dimDomain >    ThetaJacobianRangeType;
    typedef FieldMatrix< double, dimGradRange, dimDomain >    JacobianFluxRangeType;
  protected:
    const ThermodynamicsType thermoDynamics_;
  public:
  PhasefieldPhysics(const ThermodynamicsType& thermodyn):
    thermoDynamics_(thermodyn)
  {
 
  }
   
  inline void conservativeToPrimitive( const RangeType& cons, RangeType& prim ) const;
 
  template< class JacobianRangeImp >
  inline void totalEnergy( const RangeType& cons, 
                           const JacobianRangeImp& grad,
                           double& kin,
                           double& therm,
                           double& surf,
                           double& total ) const;

  inline void chemPotAndReaction( const RangeType& cons, 
																	const JacobianRangeType du,
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

 
  inline void analyticalFlux( const RangeType& u, JacobianRangeType& f ) const;
  
  inline void jacobian( const RangeType& u, JacobianFluxRangeType& a) const;

  inline double maxSpeed( const DomainType& n, const RangeType& u ) const;
  
  inline double stiffSource(const DomainType& x,
                            const double time,
                            const RangeType& u,
								            const GradientRangeType& du,
								            const ThetaRangeType& theta,
								            const ThetaJacobianRangeType& dtheta,
	 							            const JacobianRangeType& uJac,
                            RangeType& f) const;
 
  template< class JacobianRangeImp >
	inline void diffusion( const RangeType& u,
												 const JacobianRangeImp& du,
												 JacobianRangeType& f ) const;
  
  template< class JacobianRangeImp >
	inline void allenCahn( const RangeType& u,
												 const JacobianRangeImp& du,
												 ThetaJacobianRangeType& f ) const;
	
  template< class JacobianRangeImp >
	inline void boundarydiffusion( const RangeType& u,
												 const JacobianRangeImp& du,
												 JacobianRangeType& f ) const;
  
  template< class JacobianRangeImp >
	inline void boundaryallenCahn( const RangeType& u,
												 const JacobianRangeImp& du,
												 ThetaJacobianRangeType& f ) const;

  inline double nstkSource(const DomainType xglobal, const double time) const ;  
  inline double acSource(const DomainType xglobal, const double time) const ;
 
  
public:

	inline double delta()const  { return thermoDynamics_.delta(); }
	inline double deltaInv()const{ return thermoDynamics_.deltaInv(); }
  inline double mu1() const { return thermoDynamics_.mu1();}
 	inline double mu2() const { return thermoDynamics_.mu2();}


};



template< class Thermodynamics >
  inline void PhasefieldPhysics<2,Thermodynamics>
::conservativeToPrimitive( const RangeType& cons, RangeType& prim ) const
 {
  	assert( cons[0] > 0. );
  	double rho,rho_inv,phi;
  	rho=cons[0];
    rho_inv=1./rho;
    phi=cons[phaseId];
    
    //velocity 
    for(int i=0;int<dimDomain;++i)
      prim[i] = cons[i+1]*rho_inv;
    //pressure  
  	prim[phaseId-1] = thermoDynamics_.pressure(rho,phi);
    //phasefield
    prim[phaseId] = phi;
  
  }

 template< class Thermodynamics > 
 template<class JacobianRangeImp>   
  inline void PhasefieldPhysics< 2, Thermodynamics >
  :: totalEnergy( const RangeType& cons, 
                  const JacobianRangeImp& grad , 
                  double& kin, 
                  double& therm,
                  double& surf,
                  double& total ) const
  {
	  double rho = cons[0];
	  double rho_inv = 1./rho;
	  double phi = cons[phaseId];
    FieldVector<RangeFieldType, dimDomain> v;
    double kineticEnergy,surfaceEnergy;
    
    for(int i=0;i<dimDomain;++i)
      kineticEnergy+=cons[i]*cons[i];

   
    for(int i=0;i<dimDomain;++i)
      surfaceEnergy+=grad[phaseId][i]*grad[phaseId][i];
    
    kineticEnergy*=0.5*rho_inv;
    surfaceEnergy*=delta()*0.5*thermoDynamics_.h2(rho);
   
	  therm = thermoDynamics_.helmholtz( rho, phi );
	  therm +=surfaceEnergy;
    kin  = kineticEnergy;
    total = therm+kineticEnergy; 
    surf = surfaceEnergy;
  }

  template< class Thermodynamics >
  inline void PhasefieldPhysics< 2,Thermodynamics >
  ::chemPotAndReaction( const RangeType& cons, 
												const JacobianRangeType& du,
                        double& mu,
												double& reaction ) const
	{
		assert( cons[0] > 1e-20 );

		double rho=cons[0];
		double phi=cons[phaseId];
		phi/=rho;
    double dxphi=du[3][0];
    double dyphi=du[3][1];

  	mu=thermoDynamics_.chemicalPotential(rho,phi);
		mu+=delta()*thermoDynamics_.h2prime(rho)*0.5*(dxphi*dxphi+dyphi*dyphi);
    reaction=thermoDynamics_.reactionSource(rho,phi); 
  }

  template< class Thermodynamics >
	inline void PhasefieldPhysics< 2,Thermodynamics >
  ::pressureAndReaction( const RangeType& cons, 
                         double& p,
												 double& reaction ) const
	{
    std::cout<<"Don't call this with WB\n";
	}



  template< typename Thermodynamics >
  inline void  PhasefieldPhysics< 2,Thermodynamics >
  ::analyticalFlux( const RangeType& u, JacobianRangeType& f ) const
  {
		assert( u[0] > 1e-10 );
		const double rho_inv = 1. / u[0];
		const double vx = u[1]*rho_inv;
		const double vy = u[2]*rho_inv;
    const double phi= u[3];
		
		f[0][0] = u[1];  // \rho*vx
	  f[0][1] = u[2];  // \rho*vy
    f[1][0] = vx*u[1];  // \rho*vx*vx
	  f[1][1] = vx*u[2];  // \rho*vx*vy
    f[2][0] = vy*u[1];  // \rho*vy*vx
    f[2][1] = vy*u[2];  // \rho*vy*vy
    f[phaseId][0]=0.;
    f[phaseId][1]=0.;
  }

  template< class Thermodynamics > 
  inline void PhasefieldPhysics< 2,Thermodynamics>
  ::jacobian( const RangeType & u, JacobianFluxRangeType& a) const
  {
    assert(u[0] > 1e-10);
    a = 0;  
    const double rho_inv=1./u[0];
    const double vx = u[1]*rho_inv;
		const double vy = u[2]*rho_inv;
    const double phi=u[phaseId];
 
    a[0][0] = u[0];
    a[1][1] = u[0];
    a[2][0] = vx;
    a[3][1] = vx;
    a[4][0] = vy;
    a[5][1] = vy;
    a[6][0] = phi;
    a[7][1] = phi;
  }

  

	template< class Thermodynamics >
	inline double PhasefieldPhysics< 2, Thermodynamics  >
	::stiffSource(const DomainType& x,
                const double time,
                const RangeType& u,
								const GradientRangeType& du,
					  		const ThetaRangeType& theta,
								const ThetaJacobianRangeType& dtheta,
								const JacobianRangeType& jacU,
                RangeType& f) const
	{
    double rho_inv=1./u[0];
    double phi=u[phaseId];
    double dphix=jacU[phaseId][0];
    double dphiy=jacU[phaseId][1];
    double vx=u[1]*rho_inv;
    double vy=u[2]*rho_inv;
    
    double reactionFac=thermoDynamics_.reactionFactor();

    f[0]=0;
    //-(\rho\nabla\mu-\tau\nabla\phi) 
   	f[1]=-dtheta[0][0]*u[0]+dphix*theta[1];
		f[2]=-dtheta[0][1]*u[0]+dphiy*theta[1];
    f[phaseId]=theta[1]
    f[phaseId]*=-reactionFac;
    f[phaseId]*=rho_inv;
    //nonconservative Discretization of transport term
    f[phaseId]-=vx*dphix+vy*dphiy;
    return 0.4*reactionFac*deltaInv(); 
  }
  
  template< class Thermodynamics >
  template< class JacobianRangeImp >
  inline void PhasefieldPhysics< 2, Thermodynamics >
  ::diffusion( const RangeType& u,
               const JacobianRangeImp& du,
               JacobianRangeType& diff) const
  {
    // du isd(u) which is 4x2 matrix (for 2d case)
    assert( u[0] > 1e-10 );
		const double muLoc = mu1();
    const double lambdaLoc =mu2();
  
    // get dx_u, dz_u, dx_w, dz_w, dx_T, dz_T (in 2d case) for du
    const double du10=du[1][0];//dx U_1
    const double du11=du[1][1];//dy U_1
    const double du20=du[2][0];//dx U_2
    const double du21=du[2][1];//dy U_2
  
    const double dxu = du10;
    const double dzu = du11;
    const double dxw = du20; 
    const double dzw = du21; 
   
    const double tau00 = (2.*muLoc+lambdaLoc)*dxu + lambdaLoc*dzw;
    const double tau01 = muLoc*(dxw + dzu);
    const double tau10 = tau01;
    const double tau11 = lambdaLoc*dxu + (2.*muLoc+lambdaLoc)*dzw;

    // 1st row
    diff[0][0] = 0.;                   diff[0][1] = 0.;

    // 2nd row
    diff[1][0] = tau00;                diff[1][1] = tau01;

    // 3rd row
    diff[2][0] = tau10;                diff[2][1] = tau11;

    // 4th row
    diff[3][0] =0.;                diff[3][1] = 0.;
    
  }
  
  template<class Thermodynamics>
  template< class JacobianRangeImp >
  inline void PhasefieldPhysics< 2, Thermodynamics >
  ::boundarydiffusion( const RangeType& u,
               const JacobianRangeImp& du,
               JacobianRangeType& diff) const
  {
    diff = 0.;                

    
  } 

  template< class Thermodynamics >
  template< class JacobianRangeImp >
  inline void PhasefieldPhysics< 2, Thermodynamics >
  ::allenCahn( const RangeType& u,
               const JacobianRangeImp& du,
               ThetaJacobianRangeType& diff ) const
  {
   
	 
	  diff[0][0]=0.;
    diff[0][1]=0.;
    diff[1][0]=-delta()*thermoDynamics_.h2( u[0] )*du[3][0];
    diff[1][1]=-delta()*thermoDynamics_.h2( u[0] )*du[3][1];
    
  }
  
  template<class Thermodynamics>
  template< class JacobianRangeImp >
  inline void PhasefieldPhysics< 2, Thermodynamics >
  ::boundaryallenCahn( const RangeType& u,
               const JacobianRangeImp& du,
               ThetaJacobianRangeType& diff ) const
  {
	
	  diff[0][0]=0.;
    diff[0][1]=0.;
	  diff[1][0]=0.;
    diff[1][1]=0.;
  }


 template< class Thermodynamics >
 inline double PhasefieldPhysics< 2, Thermodynamics >
 ::maxSpeed( const DomainType& n, const RangeType& u) const
 {
  assert( u[0] > 1e-10 );
  RangeFieldType u_normal = (u[1]*n[0]+u[2]*n[1]) / u[0];
  RangeFieldType c = thermoDynamics_.a(u[0],u[3])* n.two_norm2();
  return std::abs(u_normal) + std::sqrt(c);
 } 

}//end namespace Dune
#endif// PHYSICS_INLINE_HH


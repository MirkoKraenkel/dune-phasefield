#ifndef PHYSICSWB_INLINE1D_HH
#define PHYSICSWB_INLINE1D_HH
namespace Dune{
//================================================
//
//
//================================================


template<class Thermodynamics>
class PhasefieldPhysics< 1,Thermodynamics>
{
  typedef Thermodynamics ThermodynamicsType;
 
  public:
    enum{dimDomain=1};  
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
                           double& total) const;

  inline void chemPotAndReaction( const RangeType& cons, 
																  const JacobianRangeType& du,	
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
inline void PhasefieldPhysics< 1, Thermodynamics >
::conservativeToPrimitive( const RangeType& cons, RangeType& prim ) const
 {
  	assert( cons[0] > 0. );
  	double rho,rho_inv,phi;
  	rho=cons[0];
    rho_inv=1./rho;
    phi=cons[phaseId];
    
    //velocity 
    prim[0] = cons[1]*rho_inv;
    //pressure  
  	prim[phaseId-1] = thermoDynamics_.pressure(rho,phi);
    //phasefield
    prim[phaseId] = phi;
  
  }

 template< class Thermodynamics > 
 template<class JacobianRangeImp>   
 inline void PhasefieldPhysics< 1,Thermodynamics >
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
    
    double kineticEnergy,surfaceEnergy;
    kineticEnergy=cons[1]*cons[1];
    double gradphi=grad[2][0];
    surfaceEnergy=thermoDynamics_.h2(rho)*gradphi*gradphi;
    
    kineticEnergy*=0.5*rho_inv;
    surfaceEnergy*=delta()*0.5;
   
	  therm = thermoDynamics_.helmholtz( rho, phi );
	  therm +=surfaceEnergy;
    kin  = kineticEnergy;
    total = therm+kineticEnergy; 
    surf = surfaceEnergy;
  }

  template< class Thermodynamics >
  inline void PhasefieldPhysics< 1,Thermodynamics >
  ::chemPotAndReaction( const RangeType& cons, 
												const JacobianRangeType& du,
                        double& mu,
												double& reaction ) const
	{
		assert( cons[0] > 1e-20 );

		double rho=cons[0];
		double phi=cons[phaseId];
   
    double dxphi=du[2][0];
  	mu=thermoDynamics_.chemicalPotential(rho,phi);
		mu+=delta()*thermoDynamics_.h2prime(rho)*0.5*dxphi*dxphi;
    reaction=thermoDynamics_.reactionSource(rho,phi); 
  }

  template< class Thermodynamics >
	inline void PhasefieldPhysics< 1,Thermodynamics >
  ::pressureAndReaction( const RangeType& cons, 
                         double& p,
												 double& reaction ) const
	{
		assert( cons[0] > 1e-20 );
	  
		double rho=cons[0];
		double phi=cons[phaseId];
		
		p=thermoDynamics_.pressure(rho,phi);
		reaction=thermoDynamics_.reactionSource(rho,phi); 
			
		assert( p > 1e-20 );
	}

  template< typename Thermodynamics >
  inline void  PhasefieldPhysics< 1, Thermodynamics>
  ::analyticalFlux( const RangeType& u, JacobianRangeType& f ) const
  {
		assert( u[0] > 1e-10 );
		double rho_inv = 1. / u[0];
		const double v = u[1]*rho_inv;
    
    f[0][0] = u[1];
    f[1][0] = v*u[1];
    f[2][0] = 0.;
  }

  template< class Thermodynamics > 
  inline void PhasefieldPhysics< 1, Thermodynamics>
  ::jacobian( const RangeType & u, JacobianFluxRangeType& a) const
  {
    //assert(u[0] > 1e-10);
    double rho=u[0];
    a[0][0] = rho;     //rho
    a[1][0] = u[1]/rho;//(rho v)/rho
    a[2][0] = u[2];    // phi
  }

	template< class Thermodynamics >
	inline double PhasefieldPhysics< 1, Thermodynamics  >
	::stiffSource(const DomainType& xglobal, //model already gives globla coordinate
                const double time,
                const RangeType& u,
								const GradientRangeType& du,
					  		const ThetaRangeType& theta,
								const ThetaJacobianRangeType& dtheta,
								const JacobianRangeType& jacU,
                RangeType& f) const
	{
    
    const double rho_inv=1./u[0];
    double dphi=jacU[phaseId][0];
    double v=u[1]*rho_inv;
  
    double reactionFac=thermoDynamics_.reactionFactor();

    f[0]=0;
    //-(\rho\nabla\mu-\tau\nabla\phi) 
    f[1]=-dtheta[0]*u[0]+dphi*theta[1];
    f[phaseId]=theta[1];
    f[phaseId]*=-reactionFac;
    f[phaseId]*=rho_inv;
    //nonconservative Discretization of transport term
    f[phaseId]-=v*dphi;

    return 0.4*reactionFac*deltaInv(); 
  }
  
  template< class Thermodynamics >
  template< class JacobianRangeImp >
  inline void PhasefieldPhysics< 1 ,Thermodynamics >
  ::diffusion( const RangeType& u,
               const JacobianRangeImp& du,
               JacobianRangeType& diff) const
  {
    assert( u[0] > 1e-10 );
		const double muLoc = mu1();
    const double dxv   = du[1][0]; //dv/dx
		diff[0][0]=0.;
		diff[1][0]=muLoc*dxv;
  	diff[2][0]=0.;
  }
  
  template<class Thermodynamics>
  template< class JacobianRangeImp >
  inline void PhasefieldPhysics< 1, Thermodynamics>
  ::allenCahn( const RangeType& u,
               const JacobianRangeImp& du,
               ThetaJacobianRangeType& diff ) const
  {
	
    diff[0][0]=0.;
    diff[1][0]=-delta()*thermoDynamics_.h2( u[0] )*du[2][0];
  }

  template< class Thermodynamics >
  template< class JacobianRangeImp >
  inline void PhasefieldPhysics< 1 ,Thermodynamics >
  ::boundarydiffusion( const RangeType& u,
               const JacobianRangeImp& du,
               JacobianRangeType& diff) const
  {
     diffusion(u,du,diff); 
  }
  
  template<class Thermodynamics>
  template< class JacobianRangeImp >
  inline void PhasefieldPhysics< 1, Thermodynamics>
  ::boundaryallenCahn( const RangeType& u,
               const JacobianRangeImp& du,
               ThetaJacobianRangeType& diff ) const
  {
	  diff[0][0]=0.;
	  diff[1][0]=0.;
  }


 template< class Thermodynamics >
 inline double PhasefieldPhysics< 1, Thermodynamics>
 ::maxSpeed( const DomainType& n, const RangeType& u) const
 {
  assert(u[0] > 1e-20);
  double u_normal=(u[1]*n[0])/u[0];
  double c=thermoDynamics_.a( u[0] , u[2] );
  return std::abs(u_normal)+sqrt(c);
 } 

}//end namespace Dune
#endif// PHYSICS_INLINE_HH


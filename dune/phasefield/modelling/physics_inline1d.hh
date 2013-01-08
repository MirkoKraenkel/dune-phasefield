#ifndef PHYSICS_INLINE1D_HH
#define PHYSICS_INLINE1D_HH
namespace Dune{
//================================================
//
//for dim=1
//
//================================================


template<class Thermodynamics>
class PhasefieldPhysics<1,Thermodynamics>
{
  typedef Thermodynamics ThermodynamicsType;
 
  public:
    enum{dimDomain=1};  
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

   typedef FieldMatrix< double, dimThetaRange, dimDomain >    ThetaJacobianRangeType;
   typedef FieldMatrix< double, dimGradRange, dimDomain >    JacobianFluxRangeType;
  protected:
    const ThermodynamicsType thermoDynamics_;
  public:
  PhasefieldPhysics(const ThermodynamicsType& thermodyn):
    thermoDynamics_(thermodyn),
    delta_(Dune::Fem::Parameter::getValue<double>("phasefield.delta")),
    delta_inv_(1./delta_)
  {
   std::cout<<"PHYSICS";
   std::cout<<"================"<<thermoDynamics_.delta()<<"==========\n";
  }
 
  inline void conservativeToPrimitive( const RangeType& cons, RangeType& prim ) const;
 
  template < class JacobianRangeImp >
  inline void totalEnergy( const RangeType& cons, 
                           const JacobianRangeImp& grad,
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

 
  inline void analyticalFlux( const RangeType& u, JacobianRangeType& f ) const;
  
  inline void jacobian( const RangeType& u, JacobianFluxRangeType& a) const;

  inline double maxSpeed( const DomainType& n, const RangeType& u ) const;
  
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
	
  template< class JacobianRangeImp >
  inline void tension( const RangeType& u,
                       const JacobianRangeImp& du,
                       GradientRangeType& tens) const;



public:

	inline double delta()const {return delta_;}
	inline double delta_inv(){return delta_inv_;}
	inline double mu1()const {return thermoDynamics_.mu1();}
	inline double mu2(){return 1;}



protected:
	const double delta_; 
	double delta_inv_;
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
    prim[0] = cons[1]/rho;
    //pressure  
  	prim[phaseId-1] = thermoDynamics_.pressure(rho,phi);
    //phasefield
    prim[phaseId]   = phi;
  
  }

  template< class Thermodynamics > 
  template< class JacobianRangeImp >
  inline void PhasefieldPhysics<1,Thermodynamics >
  :: totalEnergy( const RangeType& cons, 
                  const JacobianRangeImp& grad , 
                  double& res ) const
  {
    assert( cons[0] > 0. );
	  double rho = cons[0];
	  double rho_inv = 1. /rho;
	  double phi = cons[phaseId];
    phi*=rho_inv;
    
    double kineticEnergy,surfaceEnergy;
    kineticEnergy=cons[1]*cons[1];
    surfaceEnergy=grad[phaseId][0]*grad[phaseId][0];
    
  
    kineticEnergy*=0.5*rho_inv;
    surfaceEnergy*=delta_*0.5*rho_inv;

 
	  double freeEnergy = thermoDynamics_.helmholtz( rho, phi );

	  res = freeEnergy + surfaceEnergy + kineticEnergy;

  }

  template< class Thermodynamics >
  inline void PhasefieldPhysics<1,Thermodynamics >
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

  template< class Thermodynamics >
	inline void PhasefieldPhysics<1,Thermodynamics >
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



  template< typename Thermodynamics >
  inline void  PhasefieldPhysics< 1, Thermodynamics>
  ::analyticalFlux( const RangeType& u, JacobianRangeType& f ) const
  {
		assert( u[0] > 1e-10 );
		double rho_inv = 1. / u[0];
		const double v = u[1]*rho_inv;
		double p;
  	double rho=u[0];
    double phi=u[phaseId];
    phi*=rho_inv;
    p=thermoDynamics_.pressure(rho,phi);
 
  
 		f[0][0] = u[1];
		f[1][0] = v*u[1]+p;
		f[2][0] = u[2]*v;
  }

  template< class Thermodynamics > 
  inline void PhasefieldPhysics< 1, Thermodynamics>
  ::jacobian( const RangeType & u, JacobianFluxRangeType& a) const
  {
    assert(u[0] > 1e-10);

    a[0][0] = u[0]; //rho
    a[1][0] = u[1];//(rho v) 
    a[2][0] = u[2];//(rho phi)
  }

	template< class Thermodynamics >
 	inline void   PhasefieldPhysics< 1, Thermodynamics >
	::nonConProduct(const RangeType & uL, 
									const RangeType & uR,
									const ThetaRangeType& thetaL,
									const ThetaRangeType& thetaR,
									RangeType& ret) const
	{
	abort();
 	}
	
	template< class Thermodynamics >
	inline double PhasefieldPhysics< 1, Thermodynamics  >
	::stiffSource(const RangeType& u,
								const GradientRangeType& du,
								const ThetaRangeType& theta,
								const ThetaJacobianRangeType& dtheta,
								RangeType& f) const
	{
			f[0]=0;
			f[1]=dtheta[0]*u[0]+du[2]*theta[1];
			f[2]=theta[1];
	  return 1.;
  }

  template< class Thermodynamics >
  template< class JacobianRangeImp >
  inline void PhasefieldPhysics< 1 ,Thermodynamics >
  ::diffusion( const RangeType& u,
               const JacobianRangeImp& du,
               JacobianRangeType& diff) const
  {
    assert( u[0] > 1e-10 );
    double rho=u[0];
    double rho_inv = 1. / rho;
 
    const double muLoc = mu1 ();

    const double v   =  u[1]*rho_inv;
   double phi =  u[2]*rho_inv;
    const double dxrho     = du[0][0]; //drho/dx
    const double dxrhou    = du[1][0]; //d(rho*v)/dx
    const double dxrhophi  = du[2][0]; //d(rho*phi)/dx
    thermoDynamics_.pressure( rho,phi);
  
    const double dxv   = rho_inv*(dxrhou - v*dxrho);
    const double rhodxphi = (dxrhophi - phi*dxrho);
    double delta_inv=1./delta_;

    diff[0][0]=0.;
    diff[1][0]=muLoc*dxv;
    diff[2][0]=rhodxphi*delta_inv;
    
  }
  
  template< class Thermodynamics>
  template< class JacobianRangeImp >
  inline void PhasefieldPhysics< 1, Thermodynamics>
  ::allenCahn( const RangeType& u,
               const JacobianRangeImp& du,
               ThetaJacobianRangeType& diff ) const
  {
   abort(); 
  }

 template< class Thermodynamics >
 inline double PhasefieldPhysics< 1, Thermodynamics>
 ::maxSpeed( const DomainType& n, const RangeType& u) const
 {
  assert(u[0] > MIN_RHO);
  double u_normal=u[1]*n[0]/u[0];
  double c=thermoDynamics_.a(u[0],u[3]);   
  return std::abs(u_normal)+sqrt(c);
 
 } 


#if 1 
// template< class Thermodynamics >
// template< class JacobianRangeImp>
 template< class Thermodynamics>
 template< class JacobianRangeImp >
 inline void PhasefieldPhysics< 1, Thermodynamics>
 ::tension( const RangeType& u,
            const JacobianRangeImp& du,
            GradientRangeType& diff ) const
  { 
    assert( u[0] > 1e-10 );
    double rho_inv = 1. / u[0];
      
    const double phi =  u[2]*rho_inv;
    const double dxrho     = du[0][0]; //drho/dx
    const double dxrhophi  = du[2][0]; //d(rho*phi)/dx
           
 //   const double rhodxphi = (dxrhophi - phi*dxrho);
    const double dxphi = rho_inv*(dxrhophi - phi*dxrho);
    double tension =delta_*(dxphi*dxphi*rho_inv-0.5*dxphi*dxphi*rho_inv);
                  
    diff[0]=0.;
    diff[1]=tension;
    diff[2]=0;
                          
  }
#endif

}//end namespace Dune
#endif// PHYSICS_INLINE_HH


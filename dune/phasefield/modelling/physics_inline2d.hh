#ifndef PHYSICS_INLINE_HH
#define PHYSICS_INLINE_HH
namespace Dune{
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
	}
 
  inline void conservativeToPrimitive( const RangeType& cons, RangeType& prim ) const;
 
  template<class JacobianRangeImp>
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
	
  template< class JacobianRangeImp>
  inline void tension( const RangeType& u,
                       const JacobianRangeImp& du,
                       GradientRangeType& tens) const;



public:

	inline double delta()const {return delta_;}
	inline double delta_inv()const {return delta_inv_;}
	inline double mu1()const {return thermoDynamics_.mu1();}
	inline double mu2()const {return thermoDynamics_.mu2();}


protected:
	const double delta_; 
	double delta_inv_;
 };


//================================================
//
//for all dim
//
//================================================
  template< class Thermodynamics >
  inline void PhasefieldPhysics<2,Thermodynamics>
  :: conservativeToPrimitive( const RangeType& cons, RangeType& prim ) const
  {
  	assert( cons[0] > 0. );
  
  	double rho,rho_inv,phi;
  	rho=cons[0];
    rho_inv=1./rho; 
    phi=cons[phaseId];
    phi*=rho_inv;
    
    //velocity
    prim[0] = cons[1]*rho_inv;
    prim[1] = cons[2]*rho_inv;
    //pressure
  	prim[phaseId-1] = thermoDynamics_.pressure(rho,phi);
  	//phasefield
    prim[phaseId]   = phi;
  
  }

  template< class Thermodynamics > 
  template< class JacobianRangeImp>
  inline void PhasefieldPhysics< 2, Thermodynamics >
  :: totalEnergy( const RangeType& cons, 
                  const JacobianRangeImp& grad , 
                  double& res ) const
  {
    assert( cons[0] > 0. );
	  double rho = cons[0];
	  double rho_inv = 1. /rho;
	  double phi = cons[dimDomain+1];
    FieldVector<RangeFieldType, dimDomain> v;
    double kineticEnergy,surfaceEnergy;
    
    kineticEnergy=cons[1]*cons[1]+cons[2]*cons[2];
    surfaceEnergy=grad[phaseId][0]*grad[phaseId][0]+grad[phaseId][1]*grad[phaseId][1];
  
    kineticEnergy*=0.5*rho_inv;
    surfaceEnergy*=delta_*0.5*rho_inv;

 
	  double freeEnergy = thermoDynamics_.helmholtz( rho, phi );

	  res = freeEnergy + surfaceEnergy + kineticEnergy;

  }

  template< class Thermodynamics >
  inline void PhasefieldPhysics< 2,Thermodynamics >
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

  template<class Thermodynamics>
	inline void PhasefieldPhysics< 2,Thermodynamics >
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
			
//		assert( p > 1e-20 );
	}



  template< typename Thermodynamics >
  inline void  PhasefieldPhysics< 2,Thermodynamics >
  ::analyticalFlux( const RangeType& u, JacobianRangeType& f ) const
  {
		assert( u[0] > 1e-10 );
		const double rho_inv = 1. / u[0];
		const double vx = u[1]*rho_inv;
		const double vy = u[2]*rho_inv;

   
		f[0][0] = u[1];
	  f[0][1] = u[2]; 
    f[1][0] = vx*u[1];
	  f[1][1] = vx*u[2];
    f[2][0] = vy*u[1];
    f[2][1] = vy*u[2];
    f[3][0] = u[phaseId];
  }

  template<class Thermodynamics> 
  inline void PhasefieldPhysics< 2,Thermodynamics>
  ::jacobian( const RangeType & u, JacobianFluxRangeType& a) const
  {
    assert(u[0] > 1e-10);
    
    const double rho_inv=1./u[0];
    const double vx = u[1]*rho_inv;
		const double vy = u[2]*rho_inv;
    const double phi=u[phaseId]*rho_inv;








   a[0][0] = u[0];       a[0][1] = 0.;
   a[1][0] = 0.;         a[1][1] = u[0];
   a[2][0] = u[1];       a[2][1] = 0.;
   a[3][0] = 0.;         a[3][1] = u[1];
   a[4][0] = u[2];       a[4][1] = 0.;
   a[5][0] = 0.;         a[5][1] = u[2];
   a[6][0] = u[3];       a[6][1] = 0.;
   a[7][0] = 0.;         a[7][1] = u[3];
#if 0
    a[0][0] = u[0];
    a[1][1] = u[0];
    a[2][0] = vx;
    a[2][1] = vx;
    a[3][0] = vy;
    a[3][0] = vy;
    a[4][0] = phi;
    a[4][0] = phi;
#endif
  }

	template< class Thermodynamics >
 	inline void   PhasefieldPhysics< 2, Thermodynamics >
	::nonConProduct(const RangeType & uL, 
									const RangeType & uR,
									const ThetaRangeType& thetaL,
									const ThetaRangeType& thetaR,
									RangeType& ret) const
	{
    std::cout<<"Checkme physics_inline2d.hh nonConProduct\n";
	 ret[1]=0.5*(uL[0]+uR[0])*(thetaL[0]-thetaR[0]);
 	 ret[2]=0.5*(uL[0]+uR[0])*(thetaL[0]-thetaR[0]);
 

  }
	
	template<class Thermodynamics>
	inline double PhasefieldPhysics< 2, Thermodynamics  >
	::stiffSource(const RangeType& u,
								const GradientRangeType& du,
								const ThetaRangeType& theta,
								const ThetaJacobianRangeType& dtheta,
								RangeType& f) const
	{
		    std::cout<<"Checkme physics_inline2d.hh stiffSourc\n";
		
      f[0]=0;
			f[1]=dtheta[0]*u[0]+du[2]*theta[1];
			f[2]=0;
	    f[3]=theta[1];
    return 1.;
  }
  
  template<class Thermodynamics>
  template< class JacobianRangeImp >
  inline void PhasefieldPhysics< 2, Thermodynamics >
  ::diffusion( const RangeType& u,
               const JacobianRangeImp& du,
               JacobianRangeType& diff) const
  {
    // du is grad(u) which is 4x2 matrix (for 2d case)
    assert( u[0] > 1e-10 );
    double rho_inv = 1. / u[0];
    const double v[2] = { u[1]*rho_inv, u[2]*rho_inv };
  
    const double muLoc = 1.;
    const double lambdaLoc =1.;
  
    // get dx_u, dz_u, dx_w, dz_w, dx_T, dz_T (in 2d case) for du
    const double du00=du[0][0];
    const double du01=du[0][1];
    const double du10=du[1][0];//dx rho*U_1
    const double du11=du[1][1];//dy rho*U_1
    const double du20=du[2][0];//dx rho*U_2
    const double du21=du[2][1];//dy rho*U_2
  
    const double dxu = rho_inv*(du10 - v[0]*du00);//=1/rho(dx(rho*v1)-v1*dx(rho))=dx(v1);
    const double dzu = rho_inv*(du11 - v[0]*du01);
    const double dxw = rho_inv*(du20 - v[1]*du00);
    const double dzw = rho_inv*(du21 - v[1]*du01);
  
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
  

  template< class Thermodynamics >
  template< class JacobianRangeImp >
  inline void PhasefieldPhysics< 2, Thermodynamics >
  ::allenCahn( const RangeType& u,
               const JacobianRangeImp& du,
               ThetaJacobianRangeType& diff ) const
  {
   assert( u[0] > 1e-10 );
	 
	 diff[0][0]=0.;
	 diff[1][0]=delta_*du[2][0];

  }

  template< class Thermodynamics >
  template< class JacobianRangeImp >
  inline void ::PhasefieldPhysics< 2,Thermodynamics >
  ::tension( const RangeType& u,
	           const JacobianRangeImp& du,
	           GradientRangeType& diff ) const 
  {
    assert( u[0] > 1e-10 );
    double rho_inv = 1. / u[0];
  
    const double phi =  u[2]*rho_inv;
    const double dxrho     = du[0][0]; //drho/dx
    const double dyrho     = du[0][1]; //drho/dy
    const double dxrhophi  = du[3][0]; //d(rho*phi)/dx
    const double dyrhophi  = du[3][1]; //d(rho*phi)/dy
  
    const double rhodxphi = (dxrhophi - phi*dxrho);
    const double rhodyphi = (dyrhophi - phi*dyrho);
  
    const double dxphi = rho_inv*(dxrhophi - phi*dxrho);
    const double dyphi = rho_inv*(dyrhophi - phi*dyrho);

  // double tension = delta_*rhodxphi*dxphi;
  
    diff[0]=0; 
    diff[1]=0;
    diff[2]=rhodxphi*dxphi;
    diff[3]=rhodyphi*dxphi;
    diff[4]=rhodyphi*dxphi;
    diff[5]=rhodyphi*dyphi;
    diff[6]=0;
    diff[7]=0;
    diff*=delta_;
}


 template< class Thermodynamics >
 inline double PhasefieldPhysics< 2, Thermodynamics >
 ::maxSpeed( const DomainType& n, const RangeType& u) const
 {
//   abort();
  return 0.;
 } 



}//end namespace Dune
#endif// PHYSICS_INLINE_HH


#ifndef PHYSICS_INLINE2D_HH
#define PHYSICS_INLINE2D_HH
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
 
  template<class JacobianRangeImp>
  inline void totalEnergy( const RangeType& cons, 
                           const JacobianRangeImp& grad,
                           double& kin,
                           double& total ) const;

  inline void chemPotAndReaction( const RangeType& cons, 
																	double& mu,
																	double& reaction ) const;

	inline void pressureAndReaction( const RangeType& cons, 
																	 double& p,
																	 double& reaction ) const;

 
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
	inline void boundarydiffusion( const RangeType& u,
												 const JacobianRangeImp& du,
												 JacobianRangeType& f ) const;
  
  template< class JacobianRangeImp >
	inline void baoundaryallenCahn( const RangeType& u,
												 const JacobianRangeImp& du,
												 ThetaJacobianRangeType& f ) const;
	
  template< class JacobianRangeImp>
  inline void tension( const RangeType& u,
                       const JacobianRangeImp& du,
                       GradientRangeType& tens) const;



public:

	inline double delta() const    {return thermoDynamics_.delta();}
	inline double deltaInv() const { return thermoDynamics_.deltaInv();}
	inline double mu1() const      { return thermoDynamics_.mu1();}
	inline double mu2() const      { return thermoDynamics_.mu2();}


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
                  double& kin,
                  double& total ) const
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
    surfaceEnergy*=delta()*0.5*rho_inv;

 
	  double freeEnergy = thermoDynamics_.helmholtz( rho, phi );
    kin = kineticEnergy;
	  total = freeEnergy + surfaceEnergy + kineticEnergy;

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
			
//		assert( p > 1e-20 );
	}



  template< typename Thermodynamics >
  inline void  PhasefieldPhysics< 2,Thermodynamics >
  ::analyticalFlux( const RangeType& u, JacobianRangeType& f ) const
  {
    f=0.;
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
		double rho=u[0];
    double phi=u[phaseId];
    phi/=rho;
    double reaction=thermoDynamics_.reaction(rho,phi);
  
   f[0]=0;
   f[1]=0;
	 f[2]=0;
   f[3]=0;-reaction;
    return 0.4*deltaInv*deltaInv();
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
  
    const double phi = u[3]*rho_inv;
    const double muLoc = mu1();
    const double lambdaLoc = mu2();
  
    // get dx_u, dz_u, dx_w, dz_w, dx_T, dz_T (in 2d case) for du
    const double du00=du[0][0];
    const double du01=du[0][1];
    const double du10=du[1][0];//dx rho*U_1
    const double du11=du[1][1];//dy rho*U_1
    const double du20=du[2][0];//dx rho*U_2
    const double du21=du[2][1];//dy rho*U_2
    const double du30=du[3][0];//dx rho*phi
    const double du31=du[3][1];;//dz rho*phi
 
    const double dxu = rho_inv*(du10 - v[0]*du00);//=1/rho(dx(rho*v1)-v1*dx(rho))=dx(v1);
    const double dzu = rho_inv*(du11 - v[0]*du01);
    const double dxw = rho_inv*(du20 - v[1]*du00);
    const double dzw = rho_inv*(du21 - v[1]*du01);
    const double dxphi = rho_inv*( du30 - phi*du00); //dx(rho*phi)-phi*dx(rho)=rho*dx(phi);
    const double dzphi = rho_inv*( du31 - phi*du01);
 
    const double tau00 = (2.*muLoc+lambdaLoc)*dxu + lambdaLoc*dzw;
    const double tau01 = muLoc*(dxw + dzu);
    const double tau10 = tau01;
    const double tau11 = lambdaLoc*dxu + (2.*muLoc+lambdaLoc)*dzw;

    // 1st row
    diff[0][0] = 0.;                   diff[0][1] = 0.;
    // 2nd row
    diff[1][0] =0;                diff[1][1] = 0.; 
   // 3rd row
    diff[2][0] = 0.;                diff[2][1] =0.;
  // 4th row
    diff[3][0] = delta()*dxphi;                diff[3][1] =delta()*dzphi;
    
   }
   template<class Thermodynamics>
  template< class JacobianRangeImp >
  inline void PhasefieldPhysics< 2, Thermodynamics >
  ::boundarydiffusion( const RangeType& u,
               const JacobianRangeImp& du,
               JacobianRangeType& diff) const
  {
    // du is grad(u) which is 4x2 matrix (for 2d case)
    assert( u[0] > 1e-10 );
    double rho_inv = 1. / u[0];
    const double v[2] = { u[1]*rho_inv, u[2]*rho_inv };

    const double phi=u[3]*rho_inv;

    const double muLoc = mu1();
    const double lambdaLoc = mu2();
  
    // get dx_u, dz_u, dx_w, dz_w, dx_T, dz_T (in 2d case) for du
    const double du00=du[0][0];
    const double du01=du[0][1];
    const double du10=du[1][0];//dx rho*U_1
    const double du11=du[1][1];//dy rho*U_1
    const double du20=du[2][0];//dx rho*U_2
    const double du21=du[2][1];//dy rho*U_2
    const double du30=du[3][0];//dx rho*phi
    const double du31=du[3][1];;//dz rho*phi
 
    const double dxu = rho_inv*(du10 - v[0]*du00);//=1/rho(dx(rho*v1)-v1*dx(rho))=dx(v1);
    const double dzu = rho_inv*(du11 - v[0]*du01);
    const double dxw = rho_inv*(du20 - v[1]*du00);
    const double dzw = rho_inv*(du21 - v[1]*du01);
    const double dxphi = rho_inv*( du30 - phi*du00); //dx(rho*phi)-phi*dx(rho)=rho*dx(phi);
    const double dzphi = rho_inv*( du31 - phi*du01);
 
    const double tau00 = (2.*muLoc+lambdaLoc)*dxu + lambdaLoc*dzw;
    const double tau01 = muLoc*(dxw + dzu);
    const double tau10 = tau01;
    const double tau11 = lambdaLoc*dxu + (2.*muLoc+lambdaLoc)*dzw;

    // 1st row
    diff[0][0] = 0.;                   diff[0][1] = 0.;
    // 2nd row
    // diff[1][0] = tau00;                diff[1][1] = tau01;
    diff[1][0] = 0.;                diff[1][1] = 0.; 
      // 3rd row
  //  diff[2][0] = tau10;                diff[2][1] = tau11;
      diff[2][0] = 0.;                diff[2][1] = 0.;
  // 4th row
    diff[3][0] = 0;                diff[3][1] = 0;
    
   }
 

  template< class Thermodynamics >
  template< class JacobianRangeImp >
  inline void PhasefieldPhysics< 2, Thermodynamics >
  ::allenCahn( const RangeType& u,
               const JacobianRangeImp& du,
               ThetaJacobianRangeType& diff ) const
  {
	   abort(); 
  }

  template< class Thermodynamics >
  template< class JacobianRangeImp >
  inline void ::PhasefieldPhysics< 2,Thermodynamics >
  ::tension( const RangeType& u,
	           const JacobianRangeImp& du,
	           GradientRangeType& diff ) const 
  {
    assert( u[0] > 1e-10 );
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


#ifndef PHYSICS_INLINE1D_HH
#define PHYSICS_INLINE1D_HH
namespace Dune{
//================================================
//
//for dim=1
//
//================================================


template<class Thermodynamics>
class PhasefieldPhysics< 1,Thermodynamics>
{
   typedef Thermodynamics ThermodynamicsType;
 
  public:
    enum { dimDomain = 1 };  
    enum { phaseId = dimDomain + 1 };
    enum { dimRange = dimDomain + 2 };
    enum { dimThetaRange =  2 };
    enum { dimThetaGradRange = dimThetaRange*dimDomain };
    enum { dimGradRange = dimRange * dimDomain };
    
    typedef double RangeFieldType;

    using DomainType             = FieldVector< double, dimDomain >;
    using FaceDomainType         = FieldVector< double, dimDomain - 1 >;
    using RangeType              = FieldVector< double, dimRange >;
    using ThetaRangeType         = FieldVector< double, dimThetaRange >;
    using GradientType           = FieldVector< double, dimGradRange >;
    using ThetaGradientRangeType = FieldVector< double, dimThetaGradRange >;
    using GradientRangeType           = FieldVector< double, dimGradRange >;

    using JacobianRangeType      = FieldMatrix< double, dimRange, dimDomain >;
    using FluxRangeType          = FieldMatrix< double, dimRange, dimDomain >;
    using ThetaJacobianRangeType = FieldMatrix< double, dimThetaRange, dimDomain >;
    using JacobianFluxRangeType  = FieldMatrix< double, dimGradRange, dimDomain >;

  protected:
    const ThermodynamicsType thermoDynamics_;
  public:
    PhasefieldPhysics(const ThermodynamicsType& thermodyn):
      thermoDynamics_(thermodyn)
      {}
 
    inline void conservativeToPrimitive ( const RangeType& cons,
                                          RangeType& prim ) const;
 
    template< class JacobianRangeImp >
    inline void totalEnergy ( const RangeType& cons,
                              const JacobianRangeImp& grad,
                              double& kin,
                              double& therm,
                              double& surf,
                              double& total) const;

	  inline void pressureAndReaction ( const RangeType& cons,
                                      double& p,
                                      double& reaction ) const;
  
    inline void analyticalFlux ( const RangeType& u, JacobianRangeType& f ) const;
  
    inline void jacobian ( const RangeType& u, JacobianFluxRangeType& a) const;

 
    inline double stiffSource ( const DomainType& xglobal,
                                const double time,
                                const RangeType& u,
                                const GradientRangeType& du,
                                RangeType& f) const;
 
    template< class JacobianRangeImp >
	  inline void diffusion ( const RangeType& u,
                            const JacobianRangeImp& du,
                            JacobianRangeType& f ) const;

    template< class JacobianRangeImp >
	  inline void boundarydiffusion ( const RangeType& u,
                                    const JacobianRangeImp& du,
                                    JacobianRangeType& f ) const;

    template< class JacobianRangeImp >
    inline void tension( const RangeType& u,
                       const JacobianRangeImp& du,
                       GradientRangeType& tens) const;

    inline double maxSpeed( const DomainType& n, const RangeType& u ) const;

	  inline double delta() const { return thermoDynamics_.delta(); }
	  inline double deltaInv() const { return thermoDynamics_.deltaInv(); }
	  inline double mu1() const { return thermoDynamics_.mu1Liq(); }
	  inline double mu2() const { return thermoDynamics_.mu2Liq(); }
 
 };



  template< class Thermodynamics >
  inline void PhasefieldPhysics< 1, Thermodynamics >
  ::conservativeToPrimitive ( const RangeType& cons, RangeType& prim ) const
  {
    assert( cons[0] > 0. );
 
  	double rho,rho_inv,phi;
  	rho=cons[0];
    rho_inv=1./rho;
  	phi=cons[phaseId];
    phi*=rho_inv;
    //velocity 
    for( int i = 0 ; i < dimDomain ; ++i )
      prim[i] = cons[i+1]*rho_inv;

    //pressure  
  	prim[phaseId-1] = thermoDynamics_.pressure(rho,phi);
    //phasefield
    prim[phaseId]   = phi;
  
  }

  template< class Thermodynamics > 
  template< class JacobianRangeImp >
  inline void PhasefieldPhysics< 1,Thermodynamics >
  ::totalEnergy ( const RangeType& cons, 
                  const JacobianRangeImp& grad , 
                  double& kin,
                  double& therm,
                  double& surf,
                  double& total) const
  {
	  double rho = cons[0];
	  double rho_inv = 1. /rho;
	  double phi = cons[phaseId];
    phi*=rho_inv;
    
    double kineticEnergy,surfaceEnergy;
    
    for(int i=0;i<dimDomain;++i)
      kineticEnergy+=cons[i]*cons[i];
    
    kineticEnergy*=0.5*rho_inv;
    //recontsruction of gradphi
    double gradphi=grad[2][0];
    gradphi-=phi*grad[0][0];
    gradphi*=rho_inv;
   
    surfaceEnergy=gradphi*gradphi;
    surfaceEnergy*=delta()*0.5*thermoDynamics_.h2(rho);
	  therm = thermoDynamics_.helmholtz( rho, phi );
    therm+=surfaceEnergy;
    kin = kineticEnergy;
	  total = therm+kineticEnergy;
    surf = surfaceEnergy;
  }

  template< class Thermodynamics >
  inline void PhasefieldPhysics< 1,Thermodynamics >
  ::pressureAndReaction ( const RangeType& cons, 
                          double& p,
                          double& reaction ) const
  {
    assert( cons[0] > 1e-20 );
	  
		double rho=cons[0];
		double phi=cons[phaseId];
		phi/=rho;
		
		p=thermoDynamics_.pressure(rho,phi);
		reaction=thermoDynamics_.reactionSource(rho,phi); 
	}

  template< typename Thermodynamics >
  inline void  PhasefieldPhysics< 1, Thermodynamics>
  ::analyticalFlux ( const RangeType& u, JacobianRangeType& f ) const
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
  ::jacobian ( const RangeType & u, JacobianFluxRangeType& a) const
  {
    assert(u[0] > 1e-10);

    a[0][0] = u[0]; //rho
    a[1][0] = u[1];//(rho v) 
    a[2][0] = u[2];//(rho phi)
  }

	
  template< class Thermodynamics >
  inline double PhasefieldPhysics< 1, Thermodynamics  >
  ::stiffSource ( const DomainType& xglobal,
                  const double time,
                  const RangeType& u,
							    const GradientRangeType& du,
							    RangeType& f) const
	{

    double rho=u[0];
    double phi=u[phaseId];
    phi/=rho;
    double reaction=thermoDynamics_.reactionSource(rho,phi);
    for( int ii = 0 ; ii<=dimDomain ; ++ii)
      f[ii]=0;
    
    f[phaseId]=-reaction*thermoDynamics_.reactionFactor();
    
    return 4*thermoDynamics_.reactionFactor()*deltaInv();
  }

  template< class Thermodynamics >
  template< class JacobianRangeImp >
  inline void PhasefieldPhysics< 1 ,Thermodynamics >
  ::diffusion ( const RangeType& u,
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
  
    const double dxv   = rho_inv*(dxrhou - v*dxrho);
    const double dxphi = rho_inv*(dxrhophi - phi*dxrho);

    const double reactionfactor=thermoDynamics_.reactionFactor();
    diff[0][0]=0.;
    diff[1][0]=muLoc*dxv;
    diff[2][0]=delta()*thermoDynamics_.h2(rho)*dxphi;
    diff[2][0]*=reactionfactor;
  }

  template< class Thermodynamics >
  template< class JacobianRangeImp >
  inline void PhasefieldPhysics< 1 ,Thermodynamics >
  ::boundarydiffusion ( const RangeType& u,
                        const JacobianRangeImp& du,
                        JacobianRangeType& diff) const
  {
    assert( u[0] > 1e-10 );
    double rho=u[0];
    double rho_inv = 1. / rho;
    
    const double muLoc = mu1 ();
    const double v   =  u[1]*rho_inv;
    const double dxrho     = du[0][0]; //drho/dx
    const double dxrhou    = du[1][0]; //d(rho*v)/dx
  
    const double dxv   = rho_inv*(dxrhou - v*dxrho);

    diff[0][0]=0.;
    diff[1][0]=muLoc*dxv;
    diff[2][0]=0;
 
  }


  template< class Thermodynamics>
  template< class JacobianRangeImp>
  inline void PhasefieldPhysics< 1, Thermodynamics>
  ::tension ( const RangeType& u,
              const JacobianRangeImp& du,
              GradientRangeType& diff ) const
  { 
    assert( u[0] > 1e-10 );
    double rho_inv = 1. / u[0];
    double rho = u[0];  
    const double phi =  u[dimDomain+1]*rho_inv;
    const double dxrho     = du[0][0]; //drho/dx
    const double dxrhophi  = du[dimDomain+1][0]; //d(rho*phi)/dx
           
    const double dxphi = rho_inv*(dxrhophi - phi*dxrho);
    //1/2 ph2(rho)+h2(rho)=1/2 h2(rho)+rho h2prime(rho)
    double tension =delta()*0.5*(rho*thermoDynamics_.h2prime(rho)+thermoDynamics_.h2(rho))*dxphi*dxphi;
                  
    diff[0]=0.;
    diff[1]=tension;
    diff[2]=0;

}

  template< class Thermodynamics >
  inline double PhasefieldPhysics< 1, Thermodynamics>
  ::maxSpeed ( const DomainType& n, const RangeType& u) const
  {
    double u_normal=u[1]*n[0]/u[0];
    double phi=u[phaseId]/u[0];
    double c=thermoDynamics_.a(u[0],phi)*n.two_norm2();

    return std::abs(u_normal)+sqrt(c);
  } 

}//end namespace Dune
#endif// PHYSICS_INLINE_HH


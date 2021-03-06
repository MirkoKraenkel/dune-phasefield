#ifndef PHYSICS_INLINE2D_HH
#define PHYSICS_INLINE2D_HH
namespace Dune{
//================================================
//
//for dim=2
//
//================================================


template<class Thermodynamics>
class PhasefieldPhysics<2,Thermodynamics>
{
   typedef Thermodynamics ThermodynamicsType;
 
  public:
    enum { dimDomain = 2 };
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
  inline void PhasefieldPhysics< 2 , Thermodynamics >
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
  inline void PhasefieldPhysics< 2, Thermodynamics >
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
   
    for(int i=1;i<dimDomain+1;++i)
      kineticEnergy+=cons[i]*cons[i];
  
    kineticEnergy*=0.5*rho_inv;
    //recontsruction of gradphi 
    double dxphi=grad[phaseId][0];
    dxphi-=phi*grad[0][0];
    double dyphi=grad[phaseId][1];
    dyphi-=phi*grad[0][1];
 
    surfaceEnergy=dxphi*dxphi+dyphi*dyphi;
    surfaceEnergy*=rho_inv;
    surfaceEnergy*=rho_inv;
    surfaceEnergy*=delta()*0.5*thermoDynamics_.h2(rho);
	  therm = thermoDynamics_.helmholtz( rho, phi );
	  therm+= surfaceEnergy;
    kin = kineticEnergy;
    total = therm+kineticEnergy; 
    surf = surfaceEnergy;
  }

  template< class Thermodynamics >
  inline void PhasefieldPhysics< 2,Thermodynamics >
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
  inline void  PhasefieldPhysics< 2,Thermodynamics >
  ::analyticalFlux ( const RangeType& u, JacobianRangeType& f ) const
  {
		assert( u[0] > 1e-10 );
		const double rho_inv = 1. / u[0];
		const double vx = u[1]*rho_inv;
		const double vy = u[2]*rho_inv;
		double p;
  	double rho=u[0];
    double phi=u[phaseId];
    phi*=rho_inv;
    p=thermoDynamics_.pressure(rho,phi);
 
  
    f[0][0] = u[1]; //rho*vx
	  f[0][1] = u[2]; //rho*vy
    f[1][0] = vx*u[1]+p; //rho*vx^2+p
	  f[1][1] = vx*u[2];   //rho*vx*vy
    f[2][0] = f[1][1]; 
    f[2][1] = vy*u[2]+p;
    f[3][0] = u[phaseId]*vx;
    f[3][1] = u[phaseId]*vy;
 
		
  }
  
  template< class Thermodynamics > 
  inline void PhasefieldPhysics< 2,Thermodynamics>
  ::jacobian ( const RangeType & u, JacobianFluxRangeType& a) const
  {
    assert(u[0] > 1e-10);


	
    a[0][0] = u[0];       a[0][1] = 0.;
    a[1][0] = 0.;         a[1][1] = u[0];
    a[2][0] = u[1];       a[2][1] = 0.;
    a[3][0] = 0.;         a[3][1] = u[1];
    a[4][0] = u[2];       a[4][1] = 0.;
    a[5][0] = 0.;         a[5][1] = u[2];
    a[6][0] = u[3];       a[6][1] = 0.;
    a[7][0] = 0.;         a[7][1] = u[3];
  }

  template< class Thermodynamics >
  inline double PhasefieldPhysics< 2, Thermodynamics  >
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
    
    f[ phaseId ]=-reaction*thermoDynamics_.reactionFactor();
    
    return 4*thermoDynamics_.reactionFactor()*deltaInv();
  }

  template< class Thermodynamics >
  template< class JacobianRangeImp >
  inline void PhasefieldPhysics< 2, Thermodynamics >
  ::diffusion ( const RangeType& u,
                const JacobianRangeImp& du,
                JacobianRangeType& diff) const
  {
    // du is grad(u) which is 4x2 matrix (for 2d case)
    assert( u[0] > 1e-10 );
    double rho=u[0];
    double rho_inv = 1. / rho;

    const double v[2] = { u[1]*rho_inv, u[2]*rho_inv };
  
    const double phi = u[3]*rho_inv;
    double reactionFactor=thermoDynamics_.reactionFactor();
 
    const double muLoc = mu1 ();
    const double lambdaLoc =mu2();
  
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
#if LAPLACE
    const double tau00 = muLoc*dxu;
    const double tau01 = 0;
    const double tau10 = 0;
    const double tau11 = muLoc*dzw;
#else
    const double tau00 = (2.*muLoc+lambdaLoc)*dxu + lambdaLoc*dzw;
    const double tau01 = muLoc*(dxw + dzu);
    const double tau10 = tau01;
    const double tau11 = lambdaLoc*dxu + (2.*muLoc+lambdaLoc)*dzw;
#endif
    // 1st row
    diff[0][0] = 0.;                   diff[0][1] = 0.;

    // 2nd row
    diff[1][0] = tau00;                diff[1][1] = tau01;

    // 3rd row
    diff[2][0] = tau10;                diff[2][1] = tau11;

    // 4th row
    diff[3][0] = reactionFactor*thermoDynamics_.h2(rho)*delta()*dxphi;        
    diff[3][1] = reactionFactor*thermoDynamics_.h2(rho)*delta()*dzphi;

  }

  template< class Thermodynamics >
  template< class JacobianRangeImp >
  inline void PhasefieldPhysics< 2, Thermodynamics >
  ::boundarydiffusion ( const RangeType& u,
                        const JacobianRangeImp& du,
                        JacobianRangeType& diff) const
  {
    // du is grad(u) which is 4x2 matrix (for 2d case)
    assert( u[0] > 1e-10 );
    double rho_inv = 1. / u[0];
    const double v[2] = { u[1]*rho_inv, u[2]*rho_inv };

    
    const double muLoc = mu1 ();
    const double lambdaLoc = mu2();
  
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
#if LAPLACE   
    const double tau00 = muLoc*dxu;
    const double tau01 = 0;
    const double tau10 = 0;
    const double tau11 = muLoc*dzw;
#else
    const double tau00 = (2.*muLoc+lambdaLoc)*dxu + lambdaLoc*dzw;
    const double tau01 = muLoc*(dxw + dzu);
    const double tau10 = tau01;
    const double tau11 = lambdaLoc*dxu + (2.*muLoc+lambdaLoc)*dzw;
#endif
    // 1st row
    diff[0][0] = 0.;                   diff[0][1] = 0.;
    // 2nd row
    diff[1][0] = tau00;                diff[1][1] = tau01;
    // 3rd row
    diff[2][0] = tau10;                diff[2][1] = tau11;
    // 4th row
    diff[3][0] = 0.;                   diff[3][1] = 0.;
 
  }


  template< class Thermodynamics>
  template< class JacobianRangeImp>
  inline void PhasefieldPhysics< 2, Thermodynamics >
  ::tension ( const RangeType& u,
              const JacobianRangeImp& du,
              GradientRangeType& diff ) const
  { 
    assert( u[0] > 1e-10 );
    double rho_inv = 1. / u[0];
    double rho = u[0];  
    const double phi =  u[dimDomain+1]*rho_inv;
    const double dxrho     = du[0][0]; //drho/dx
    const double dyrho     = du[0][1]; //drho/dy
   
    const double dxrhophi  = du[3][0]; //d(rho*phi)/dx
    const double dyrhophi  = du[3][1]; //d(rho*phi)/dy
  
           
    const double dxphi = rho_inv*(dxrhophi - phi*dxrho);
    const double dyphi = rho_inv*(dyrhophi - phi*dyrho);
                  
    //(-h_2+rho*h2prime)*|\nabla\phi|^2*0.5
    double gradphisquarehalf=0.5*(dxphi*dxphi+dyphi*dyphi);
    gradphisquarehalf*=rho*thermoDynamics_.h2prime(rho)-thermoDynamics_.h2(rho);
    diff[0]=0; 
    diff[1]=0;
    diff[2]=dxphi*dxphi;
    diff[3]=dyphi*dxphi;
    diff[4]=dyphi*dxphi;
    diff[5]=dyphi*dyphi;
    diff[6]=0;
    diff[7]=0;
    diff*=thermoDynamics_.h2(rho);
    diff[2]+=gradphisquarehalf;
    diff[5]+=gradphisquarehalf;

    diff*=delta();
}

  template< class Thermodynamics >
  inline double PhasefieldPhysics< 2, Thermodynamics >
  ::maxSpeed ( const DomainType& n, const RangeType& u) const
  {
    double u_normal=(u[1]*n[0]+u[2]*n[1])/u[0];
    double phi=u[phaseId]/u[0];
    double c=thermoDynamics_.a(u[0],phi)*n.two_norm2();

    return std::abs(u_normal)+sqrt(c);
  } 

}//end namespace Dune
#endif// PHYSICS_INLINE_HH


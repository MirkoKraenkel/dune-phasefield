#ifndef PHYSICSWB_INLINE1D_HH
#define PHYSICSWB_INLINE1D_HH
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

    inline void chemPotAndReaction ( const RangeType& cons,
                                     const JacobianRangeType& du,
                                     double& mu,
                                      double& reaction ) const;
  
    inline void analyticalFlux ( const RangeType& u, JacobianRangeType& f ) const;
  
    inline void jacobian ( const RangeType& u, JacobianFluxRangeType& a) const;

 
    inline double stiffSource ( const DomainType& xglobal,
                                const double time,
                                const RangeType& u,
                                const GradientRangeType& du,
                                const ThetaRangeType& theta,
                                const ThetaJacobianRangeType& dtheta,
                                const JacobianRangeType& uJac,
                                RangeType& f) const;
 
    template< class JacobianRangeImp >
	  inline void diffusion ( const RangeType& u,
                            const JacobianRangeImp& du,
                            JacobianRangeType& f ) const;

    template< class JacobianRangeImp >
	  inline void boundarydiffusion ( const RangeType& u,
                                    const JacobianRangeImp& du,
                                    JacobianRangeType& f ) const;

    //f|phi-div(f|nabla phi)
    template< class JacobianRangeImp >
	  inline void allenCahn ( const RangeType& u,
                       const JacobianRangeImp& du,
                            ThetaJacobianRangeType& f ) const;

    //f|phi-div(f|nabla phi)
    template< class JacobianRangeImp >
	  inline void boundaryallenCahn ( const RangeType& u,
                                    const JacobianRangeImp& du,
                                    ThetaJacobianRangeType& f ) const;

    inline double maxSpeed( const DomainType& n, const RangeType& u ) const;

  inline double nstkSource(const DomainType xglobal, const double time) const ;  
  inline double acSource(const DomainType xglobal, const double time) const ;
 
  
public:

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
    prim[0] = cons[1]*rho_inv;
    

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
    
    
    kineticEnergy=cons[1]*cons[1];
    double gradphi=grad[2][0];
    surfaceEnergy=gradphi*gradphi;
    
    kineticEnergy*=0.5*rho_inv;
    surfaceEnergy*=delta()*0.5*thermoDynamics_.h2(rho);
	  therm = thermoDynamics_.helmholtz( rho, phi );
    therm+=surfaceEnergy;
    kin = kineticEnergy;
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
		phi/=rho;
		
    double dxphi=du[2][0];
    mu = thermoDynamics_.chemicalPotential(rho,phi);
		mu+= delta()*thermoDynamics_.h2prime(rho)*0.5*dxphi*dxphi;
    
		reaction=thermoDynamics_.reactionSource(rho,phi); 
	}

  template< typename Thermodynamics >
  inline void  PhasefieldPhysics< 1, Thermodynamics>
  ::analyticalFlux ( const RangeType& u, JacobianRangeType& f ) const
  {
		assert( u[0] > 1e-10 );
		double rho_inv = 1. / u[0];
		const double v = u[1]*rho_inv;
 
 		f[0][0] = u[1];
    f[1][0] = v*u[1];
		f[2][0] = u[2]*v;
  }
  
  template< class Thermodynamics > 
  inline void PhasefieldPhysics< 1, Thermodynamics>
  ::jacobian ( const RangeType & u, JacobianFluxRangeType& a) const
  {
    //assert(u[0] > 1e-10);

    a[0][0] = u[0]; //rho
    a[1][0] = u[1]/u[0];//(rho v)/rho
    a[2][0] = u[2]/u[0];//(rho phi)/rho

  }
	
  template< class Thermodynamics >
  inline double PhasefieldPhysics< 1, Thermodynamics  >
  ::stiffSource ( const DomainType& xglobal, //model already gives global coordinate
                  const double time,
                  const RangeType& u,
							    const GradientRangeType& du,
					  	    const ThetaRangeType& theta,
							    const ThetaJacobianRangeType& dtheta,
							    const JacobianRangeType& jacU,
							    RangeType& f) const
	{
    double rho_inv=1./u[0];
    double phi=u[2]*rho_inv;
    //dphi=1/rho*(d(rho*phi)-phi*drho)   
    double dphi=jacU[2][0]-phi*jacU[0][0];
    dphi*=rho_inv;
    double reactionFac=thermoDynamics_.reactionFactor();

    f[0]=0;
    //-(\rho\nabla\mu-\tau\nabla\phi) 
    f[1]=-dtheta[0]*u[0]+dphi*theta[1];
    
    
    f[phaseId]=theta[1];
    f[phaseId]*=-reactionFac;

    return 4*reactionFac*deltaInv(); 
  }

  template< class Thermodynamics >
  template< class JacobianRangeImp >
  inline void PhasefieldPhysics< 1 ,Thermodynamics >
  ::diffusion ( const RangeType& u,
                const JacobianRangeImp& du,
                JacobianRangeType& diff) const
  {
    assert( u[0] > 1e-10 );
    const double muLoc = mu1 ();
    const double dxv   = du[1][0]; //dv/dx
    diff[0][0]=0.;
    diff[1][0]=muLoc*dxv;
  	diff[2][0]=0.;
  }

  template< class Thermodynamics >
  template< class JacobianRangeImp >
  inline void PhasefieldPhysics< 1 ,Thermodynamics >
  ::boundarydiffusion ( const RangeType& u,
                        const JacobianRangeImp& du,
                        JacobianRangeType& diff) const
  {
    diffusion(u,du,diff); 
  }

  template<class Thermodynamics>
  template< class JacobianRangeImp >
  inline void PhasefieldPhysics< 1, Thermodynamics>
  ::allenCahn ( const RangeType& u,
                const JacobianRangeImp& du,
                ThetaJacobianRangeType& diff ) const
  {
    diff[0][0]=0.;
    diff[1][0]=-delta()*thermoDynamics_.h2( u[0] )*du[phaseId][0];
  }


  template< class Thermodynamics>
  template< class JacobianRangeImp>
  inline void PhasefieldPhysics< 1, Thermodynamics>
  ::boundaryallenCahn ( const RangeType& u,
                        const JacobianRangeImp& du,
                        ThetaJacobianRangeType& diff ) const
  { 
    diff[0][0]=0.;
	  diff[1][0]=0.;
  }

  template< class Thermodynamics >
  inline double PhasefieldPhysics< 1, Thermodynamics>
  ::maxSpeed ( const DomainType& n, const RangeType& u) const
  {
    assert(u[0] > 1e-20);
    double u_normal=(u[1]*n[0])/u[0];
    double c = thermoDynamics_.a( u[0] , u[2]/u[0] )*n.two_norm2();
    return std::abs(u_normal) + std::sqrt(c);
  } 

}//end namespace Dune
#endif// PHYSICS_INLINE_HH


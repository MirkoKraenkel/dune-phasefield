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
    enum { dimGradRange = dimRange * dimDomain };//=8
    
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
   typedef FieldMatrix< double, dimGradRange, dimDomain >    JacobianFluxRangeType;//8*2

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
                           double& total ) const;

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
  
  inline double stiffSource(const DomainType& x,
                            const double time,
                            const RangeType& u,
								            const GradientRangeType& du,
								            const ThetaRangeType& theta,
								            const ThetaJacobianRangeType& dtheta,
								            const JacobianRangeType& jacU,
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
	


public:
//	inline double delta()const { return 0;}
		inline double delta()const { return thermoDynamics_.delta();}
	inline double deltaInv()  const  { return thermoDynamics_.deltaInv();}
  inline double mu1()const   { return thermoDynamics_.mu1();}
 	inline double mu2()const   { return thermoDynamics_.mu2();}


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
  
  	double rho,rho_inv,phi;
  	rho=cons[0];
    rho_inv=1.; 
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
  template< class JacobianRangeImp >
  inline void PhasefieldPhysics< 2, Thermodynamics >
  :: totalEnergy( const RangeType& cons, 
                  const JacobianRangeImp& grad , 
                  double& kin,
                  double& total ) const
  {
	  double rho = cons[0];
	  double rho_inv = 1.;
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
		
		double rho=cons[0];
		double phi=cons[phaseId];
//		phi/=rho;
	  rho=1;	
		mu=thermoDynamics_.chemicalPotential(rho,phi);
		reaction=thermoDynamics_.reactionSource(rho,phi); 
	}

  template<class Thermodynamics>
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
    f=0.;
  }
  template<class Thermodynamics> 
  inline void PhasefieldPhysics< 2,Thermodynamics>
  ::jacobian( const RangeType & u, JacobianFluxRangeType& a) const
  {
   // std::cout<<"RHHHHO "<<u[0]<<"\n";

   a=0;
#if 0
   for(int r=0;r<dimRange;r++)
     for(int d=0;d<dimDomain;d++)
       a[dimDomain*r+d][d]=u[r];
    
#endif    
    
    
    
    
#if 1    
    a[6][0] = u[phaseId];
    a[7][1] = u[phaseId];
#endif  
  }



	
	template<class Thermodynamics>
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
      double mu, reaction;
      chemPotAndReaction(u,mu,reaction);

      f[0]=0;
			f[1]=-du[6]*theta[1];
			f[2]=-du[7]*theta[1];
	    f[3]=-theta[1];

      return deltaInv();  
  }

  template<class Thermodynamics>
  template< class JacobianRangeImp >
  inline void PhasefieldPhysics< 2, Thermodynamics >
  ::diffusion( const RangeType& u,
               const JacobianRangeImp& du,
               JacobianRangeType& diff) const
  {
  diff=0;
#if 0 

    // 1st row
    diff[0][0] = 0.;                   diff[0][1] = 0.;

    // 2nd row
    diff[1][0] = 0.;                   diff[1][1] = 0.;

    // 3rd row
   diff[2][0] = 0.;                    diff[2][1] = 0.;

    // 4th row
    diff[3][0] =delta()*du[3][0];      diff[3][1] =delta()*du[3][1];
#endif
  }

  template<class Thermodynamics>
  template< class JacobianRangeImp >
  inline void PhasefieldPhysics< 2, Thermodynamics >
  ::boundarydiffusion( const RangeType& u,
               const JacobianRangeImp& du,
               JacobianRangeType& diff) const
  {
    diffusion(u,du,diff);
//    diff=0;
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
    diff[1][0]=-delta()*du[3][0];
    diff[1][1]=-delta()*du[3][1];
  }
  template< class Thermodynamics >
  template< class JacobianRangeImp >
  inline void PhasefieldPhysics< 2, Thermodynamics >
  ::boundaryallenCahn( const RangeType& u,
               const JacobianRangeImp& du,
               ThetaJacobianRangeType& diff ) const
  {
   allenCahn(u,du,diff); 
  }


 template< class Thermodynamics >
 inline double PhasefieldPhysics< 2, Thermodynamics >
 ::maxSpeed( const DomainType& n, const RangeType& u) const
 {
  return 0.;
 } 



}//end namespace Dune
#endif// PHYSICS_INLINE_HH


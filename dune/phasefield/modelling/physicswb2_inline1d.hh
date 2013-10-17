#ifndef PHYSICSWB_INLINE1D_HH
#define PHYSICSWB_INLINE1D_HH
namespace Dune{
//================================================
//
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
 	inline double mu2() const { abort();return 0.;}


protected:
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
 inline void PhasefieldPhysics<1,Thermodynamics >
  :: totalEnergy( const RangeType& cons, 
                  const JacobianRangeImp& grad , 
                  double& kin, 
                  double& therm,
                  double& total ) const
  {
    assert( cons[0] > 0. );
	  double rho = cons[0];
	  double rho_inv = 1. /rho;
	  double phi = cons[phaseId];
    
    double kineticEnergy,surfaceEnergy;
    kineticEnergy=cons[1]*cons[1];
    double gradphi=grad[2][0];
    surfaceEnergy=gradphi*gradphi;
 
    kineticEnergy*=0.5*rho_inv;
    surfaceEnergy*=delta()*0.5;
   
	  therm = thermoDynamics_.helmholtz( rho, phi );
	  therm +=surfaceEnergy;
    kin  = kineticEnergy;
    total = therm+kineticEnergy; 
    
  }

  template< class Thermodynamics >
  inline void PhasefieldPhysics<1,Thermodynamics >
  ::chemPotAndReaction( const RangeType& cons, 
												double& mu,
												double& reaction ) const
	{
		assert( cons[0] > 1e-8 );

		double rho=cons[0];
		double phi=cons[phaseId];
    
	  assert( phi > -1e-8);
    assert( phi < 1.+(1e-8));
  	mu=thermoDynamics_.chemicalPotential(rho,phi);
		reaction=thermoDynamics_.reactionSource(rho,phi); 
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
		f[2][0] = 0;
  }

  template< class Thermodynamics > 
  inline void PhasefieldPhysics< 1, Thermodynamics>
  ::jacobian( const RangeType & u, JacobianFluxRangeType& a) const
  {
    //assert(u[0] > 1e-10);

    a[0][0] = u[0]; //rho
    a[1][0] = u[1]/u[0];//(rho v)/rho
    a[2][0] = u[2];//(rho phi)/rho
  }

	template< class Thermodynamics >
 	inline void   PhasefieldPhysics< 1, Thermodynamics >
	::nonConProduct(const RangeType & uL, 
									const RangeType & uR,
									const ThetaRangeType& thetaL,
									const ThetaRangeType& thetaR,
									RangeType& ret) const
	{ abort();
  }
	
	template< class Thermodynamics >
	inline double PhasefieldPhysics< 1, Thermodynamics  >
	::stiffSource(const DomainType& xglobal,
                const double time,
                const RangeType& u,
								const GradientRangeType& du,
					  		const ThetaRangeType& theta,
								const ThetaJacobianRangeType& dtheta,
								const JacobianRangeType& jacU,
                RangeType& f) const
	{

    RangeType nstksource,acsource;
    SourceTerms::nstkSource(xglobal,
                            time,
                            thermoDynamics_.delta(),
                            thermoDynamics_.velo(),
                            nstksource);
    
      SourceTerms::acSource(xglobal,
                            time,
                            thermoDynamics_.delta(),
                            thermoDynamics_.velo(),
                            acsource);

#if USEJACOBIAN
   double rho_inv=1./u[0];
//   double phi=u[2];
   double dphi=jacU[2][0];
   dphi*=-1;
#else
   double dphi=du[2];
#endif  
    double v=u[1]*rho_inv;
  	f[0]=0;
  //  f[1]=-dtheta[0]*u[0]+dphi*theta[1]+nstkSource(xglobal,time);
    f[1]=-dtheta[0]*u[0]+dphi*theta[1];
   //nonconservative Discretization of transport term
 //  f[2]=-theta[1]*deltaInv()/*rho_inv*/+v*dphi+acSource(xglobal,time);
    f[2]=-theta[1]*deltaInv()+v*dphi;
    f+=nstksource;
    f+=acsource;
 //f[2]=v*dphi;	
     
    return 0.4*deltaInv()*deltaInv(); 
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
	 //diff[1][0]=0.;
   diff[1][0]=-delta()*du[2][0];



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
  inline double PhasefieldPhysics< 1, Thermodynamics >
  ::nstkSource(const DomainType xglobal,const double t) const
  {
		double x=xglobal;
 		double c=thermoDynamics_.velo();
		double delta=thermoDynamics_.delta();
	
    double t1;
    double t10;
    double t11;
    double t21;
    double t4;
    double t5;
    double t6;
    t1 = 1/delta;
    t4 = (x+c*t)*t1;
    t5 = sinh(t4);
    t6 = cosh(t4);
    t10 = t6*t6;
    t11 = t10*t10;
    t21 = t11*t11;
    return(-0.125E-2*t1*(t5+15.0*t6)*
          (8.0*t5*t11+4.0*t5*t10+18.0*t5-24.0*t11*t6-225.0*t6)/t21);
                              
  }




	
  template< class Thermodynamics >
  inline double PhasefieldPhysics< 1, Thermodynamics >
  ::acSource(const DomainType xglobal,const double t) const
  {
    double x=xglobal;
		double delta=thermoDynamics_.delta();
		double c=thermoDynamics_.velo();
    double t12;
	  double t4;
		double t5;
		double t6;
    double t8;
		
    t4 = (x+c*t)/delta;
    t5 = cosh(t4);
    t6 = t5*t5;
    t8 = sinh(t4);
    t12 = t6*t6;
    return(0.1875E-1*(-26.0*t6+1.0+30.0*t8*t5)/t12/t6);
                          
//return  0.;
  }

 template< class Thermodynamics >
 inline double PhasefieldPhysics< 1, Thermodynamics>
 ::maxSpeed( const DomainType& n, const RangeType& u) const
 {
  assert(u[0] > 1e-8);
  double u_normal=(u[1])/u[0];
  double c=thermoDynamics_.a(u[0],u[2]);
//  std::cout<<"physics maxSpeed"<< std::abs(u_normal) <<std::endl;
//  return std::abs(u_normal)+sqrt(c);
    return std::abs(thermoDynamics_.velo());

 } 



}//end namespace Dune
#endif// PHYSICS_INLINE_HH


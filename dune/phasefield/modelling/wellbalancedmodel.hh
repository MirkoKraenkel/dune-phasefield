#ifndef WB_MODEL_HH
#define WB_MODEL_HH



// DUNE includes
#include <dune/common/version.hh>
#include <dune/fem/misc/fmatrixconverter.hh>

// local includes
#include <dune/fem-dg/models/defaultmodel.hh>
#include "physics.hh"

namespace Dune {

	// forward declaration of PhaseFlux
	template< int >
	struct PhaseFlux;

	//////////////////////////////////////////////////////
	//
	//                 NAVIER-STOKES EQNS 
	//
	//////////////////////////////////////////////////////
	template< class GridPart > 
	class PhaseModelTraits 
	{
  public:
		typedef GridPart GridPartType;
		typedef typename GridPart :: GridType                     GridType;
		enum { dimDomain = GridType :: dimensionworld };
		enum { dimRange = dimDomain + 2 };
		enum { dimThetaRange =  2 };
		enum { dimThetaGradRange = dimThetaRange*dimDomain };
		enum { dimGradRange = dimRange * dimDomain };
		enum { dimGradient = dimDomain + 1 };
  
		typedef FieldVector< double, dimDomain >                  DomainType;
		typedef FieldVector< double, dimDomain - 1 >              FaceDomainType;
		typedef FieldVector< double, dimRange >                   RangeType;
		typedef FieldVector< double, dimThetaRange >              ThetaRangeType;
		typedef FieldVector< double, dimGradRange >               GradientType;
		typedef FieldVector< double, dimThetaGradRange >          ThetaGradientRangeType;
		
		typedef FieldMatrix< double, dimRange, dimDomain >        FluxRangeType;
		typedef FieldVector< double, dimGradRange >               GradientRangeType;
		typedef FieldMatrix< double, dimRange, dimDomain >        JacobianRangeType;
		
		//For the Ac-equation 1dim reaction diffusion
		typedef FieldMatrix< double, dimThetaRange, dimDomain >    ThetaJacobianRangeType;
		

		typedef FieldMatrix< double, dimGradRange, dimDomain >    JacobianFluxRangeType;

		typedef typename GridPart :: IntersectionIteratorType     IntersectionIterator;
		typedef typename IntersectionIterator::Intersection       Intersection;
		typedef Intersection       IntersectionType;
		typedef typename GridType :: template Codim<0> :: Entity  EntityType;
		typedef typename EntityType :: EntityPointer              EntityPointerType;

	  
  };


	template< class GridPartType , class Thermodynamics>
	class PhaseModel : public DefaultModel < PhaseModelTraits< GridPartType > >
	{
  public:
		typedef Thermodynamics          ThermodynamicsType;
    typedef typename GridPartType::GridType                   GridType;
  	typedef PhaseModelTraits< GridPartType >                   Traits;
    
  
 		enum { dimDomain = Traits :: dimDomain };
		enum { dimRange  = Traits :: dimRange  };
		enum { dimGradRange = Traits::dimGradRange } ;
    
    typedef PhasefieldPhysics< dimDomain, ThermodynamicsType > PhysicsType;


		typedef typename Traits :: EntityType                     EntityType;
		typedef typename Traits :: EntityPointerType              EntityPointerType;
		typedef typename Traits :: IntersectionIterator           IntersectionIterator;
		typedef typename Traits :: Intersection                   IntersectionType;
		typedef typename Traits :: FaceDomainType                 FaceDomainType;

		typedef typename Traits :: DomainType                     DomainType;

		typedef typename Traits :: RangeType                      RangeType;
		typedef typename Traits :: GradientRangeType              GradientRangeType;
		typedef typename Traits :: JacobianRangeType              JacobianRangeType;
	
		typedef typename Traits :: ThetaRangeType                 ThetaRangeType;
		typedef typename Traits :: ThetaGradientRangeType         ThetaGradientRangeType;
		typedef typename Traits :: ThetaJacobianRangeType         ThetaJacobianRangeType;
	
		

		typedef typename Traits :: JacobianFluxRangeType          JacobianFluxRangeType;

	public:
		PhaseModel( const ThermodynamicsType& thermo ):
      phasefieldPhysics_( thermo )
		{
    }
		inline bool hasStiffSource() const { return true; }
		inline bool hasNonStiffSource() const { return false; }
		inline bool hasFlux() const { return true ; }

	 
		inline void nonConProduct(const RangeType & uL, 
															const RangeType & uR,
															const ThetaRangeType& thetaL,
															const ThetaRangeType& thetaR,
															RangeType& ret) const
		{
			phasefieldPhysics_.nonConProduct(uL, uR,thetaL,thetaR,ret);
		}


		inline double stiffSource( const EntityType& en,
															 const double time,
															 const DomainType& x,
															 const RangeType& u,
															 const GradientRangeType& du,
															 const ThetaRangeType& theta,
															 const ThetaJacobianRangeType& dtheta,
															 const JacobianRangeType& jacU,
                               RangeType & s) const
		{	
     return phasefieldPhysics_.stiffSource(u, du,theta,dtheta,jacU,s);
 		}


		inline double thetaSource( const EntityType& en,
															 const double time,
															 const DomainType& x,
															 const RangeType& u,
															 const GradientRangeType& du,
															 ThetaRangeType & s) const
		{
			return thetaSource( en, time, x, u, s );
		}


		inline double thetaSource( const EntityType& en
															 , const double time
															 , const DomainType& x
															 , const RangeType& u
															 , ThetaRangeType& s ) const 
		{
			double mu,reaction;
      
      phasefieldPhysics_.chemPotAndReaction(u,mu,reaction);
      //stheta[0]=dF/drho stheta[1]=dF/dphi
      s[0]=mu;
			s[1]=reaction;
      double deltaInv=phasefieldPhysics_.deltaInv();
			return deltaInv*deltaInv*0.4;
		}

    
		inline void advection( const EntityType& en,
													 const double time,
													 const DomainType& x,
													 const RangeType& u,
													 JacobianRangeType& f ) const 
		{
			phasefieldPhysics_.analyticalFlux(u, f);
		}


		inline double diffusionTimeStep( const IntersectionType& it,
																		 const double enVolume,
																		 const double circumEstimate,
																		 const double time,
																		 const FaceDomainType& x,
																		 const RangeType& u ) const
		{
      const double mu = phasefieldPhysics_.mu1();
   
      return mu * circumEstimate * 1 / (0.25 * u[0] * enVolume);
  	}


		inline void jacobian( const EntityType& en
													, const double time
													, const DomainType& x
													, const RangeType& u 
													, JacobianFluxRangeType& a ) const 
		{
			phasefieldPhysics_.jacobian( u, a );
		}
  

		inline bool hasBoundaryValue( const IntersectionType& it
																	, const double time
																	, const FaceDomainType& x ) const 
		{ 
			return true;
		}
 

		inline double boundaryFlux( const IntersectionType& it
																, const double time
																, const FaceDomainType& x
																, const RangeType& uLeft
																, const GradientRangeType& duLeft
																, RangeType& gLeft ) const; 
  

		inline double diffusionBoundaryFlux( const IntersectionType& it,
																				 const double time,
																				 const FaceDomainType& x,
																				 const RangeType& uLeft,
																				 const GradientRangeType& gradLeft,
																				 RangeType& gLeft ) const  
		{
			Fem::FieldMatrixConverter< GradientRangeType, JacobianRangeType> jacLeft( gradLeft );
			return diffusionBoundaryFlux( it, time, x, uLeft, jacLeft, gLeft );
		}


		template <class JacobianRangeImp>
		inline double diffusionBoundaryFlux( const IntersectionType& it,
																				 const double time,
																				 const FaceDomainType& x,
																				 const RangeType& uLeft,
																				 const JacobianRangeImp& jacLeft,
																				 RangeType& gLeft ) const  
		{
			return 1.;
		}


		inline void boundaryValue( const IntersectionType& it,
															 const double time,
															 const FaceDomainType& x,
															 const RangeType& uLeft,
															 RangeType& uRight ) const; 
 
  
		// here x is in global coordinates
		inline void maxSpeed( const DomainType& normal,
													const double time,
													const DomainType& x,
													const RangeType& u,
													double& advspeed,
													double& totalspeed ) const 
		{
			advspeed = phasefieldPhysics_.maxSpeed( normal , u );
			totalspeed=advspeed;
		}


     inline void diffusion( const EntityType& en,
						  	      			const double time,
							  		      	const DomainType& x,
									        	const RangeType& u,
								        		const GradientRangeType& v,
										        JacobianRangeType& diff ) const
		{
			Fem::FieldMatrixConverter< GradientRangeType, JacobianRangeType> jac( v );
			diffusion( en, time, x, u, jac,diff );
		}


    inline double visc() const 
    {
      std::cout<<"REVISE ME\n!"; 
      return 1.;
    }
  

		template <class JacobianRangeImp>
    inline void diffusion( const EntityType& en,
										const double time,
										const DomainType& x,
										const RangeType& u,
										const JacobianRangeImp& jac,
										JacobianRangeType& diff ) const
		{
			phasefieldPhysics_.diffusion( u, jac, diff );
		}
	  inline void boundarydiffusion( const EntityType& en,
						  	      			const double time,
							  		      	const DomainType& x,
									        	const RangeType& u,
								        		const GradientRangeType& v,
										        JacobianRangeType& diff ) const
		{
			Fem::FieldMatrixConverter< GradientRangeType, JacobianRangeType> jac( v );
			boundarydiffusion( en, time, x, u, jac,diff );
		}



		template <class JacobianRangeImp>
    inline void boundarydiffusion( const EntityType& en,
										const double time,
										const DomainType& x,
										const RangeType& u,
										const JacobianRangeImp& jac,
										JacobianRangeType& diff ) const
		{
			phasefieldPhysics_.boundarydiffusion( u, jac, diff );
		}
	



		
  	inline void allenCahnDiffusion(const RangeType& u,const GradientRangeType& du,ThetaJacobianRangeType& dv ) const
		{

      Fem::FieldMatrixConverter< GradientRangeType, JacobianRangeType> jac( du );
			allenCahnDiffusion(u,jac,dv);
		}	

    
    template <class JacobianRangeImp>
	  inline	void allenCahnDiffusion(const RangeType& u,const JacobianRangeImp& du,ThetaJacobianRangeType& dv ) const
		{
			phasefieldPhysics_.allenCahn(u,du,dv);
		}



  	inline void boundaryallenCahnDiffusion(const RangeType& u,const GradientRangeType& du,ThetaJacobianRangeType& dv ) const
		{
			Fem::FieldMatrixConverter< GradientRangeType, JacobianRangeType> jac( du );
			boundaryallenCahnDiffusion(u,jac,dv);
		}	

    
    template <class JacobianRangeImp>
	  inline	void boundaryallenCahnDiffusion(const RangeType& u,const JacobianRangeImp& du,ThetaJacobianRangeType& dv ) const
		{
			phasefieldPhysics_.boundaryallenCahn(u,du,dv);
		}



		inline double boundaryFlux( const IntersectionType& it
																, const double time
																, const FaceDomainType& x
																, const RangeType& uLeft
																, RangeType& gLeft ) const  
		{
			gLeft = 0.;
			return 0.;
		}
		

		inline void conservativeToPrimitive( const DomainType& xgl,
																				 const RangeType& cons, 
																				 RangeType& prim ) const
		{
			phasefieldPhysics_.conservativeToPrimitive( cons, prim );
		}
		
		
    inline void totalEnergy( const DomainType& xgl,
														 const RangeType& cons, 
														 const GradientRangeType& grad, 
                             FieldVector<double,1>& kin, 
														 FieldVector<double,1>& total ) const
		{
			Fem::FieldMatrixConverter< GradientRangeType, JacobianRangeType> jac( grad);
		  totalEnergy(cons,jac, kin,total);
		}

	  template <class JacobianRangeImp>
		inline void totalEnergy( const RangeType& cons, 
														 const JacobianRangeImp& grad, 
                             FieldVector<double,1>& kin, 
														 FieldVector<double,1>& total ) const
		{
		  phasefieldPhysics_.totalEnergy(cons,grad,kin[0],total[0] );
		}

 inline double delta()const
  {
    return phasefieldPhysics_.delta();
  }


		
	protected:
		//const ThermoynamicsType thermo_;
    const PhysicsType phasefieldPhysics_;
};




	template< class GridPartType, class ProblemImp >
	inline double PhaseModel< GridPartType, ProblemImp >
	:: boundaryFlux( const IntersectionType& it
									 , const double time
									 , const FaceDomainType& x
									 , const RangeType& uLeft
									 , const GradientRangeType& duLeft
									 , RangeType& gLeft ) const  
	{
		abort();
		DomainType xgl=it.intersectionGlobal().global(x);
		const typename Traits :: DomainType normal = it.integrationOuterNormal(x); 
		double p;
		double T;
		pressAndTemp( uLeft, p, T );
		gLeft = 0;

		// bnd. cond. from euler part
		for (int i=0;i<dimDomain; ++i)
			gLeft[i+1] = 0;


		return 0.;
	}



	template< class GridPartType, class ProblemImp >
	inline void PhaseModel< GridPartType, ProblemImp >
	:: boundaryValue( const IntersectionType& it
										, const double time
										, const FaceDomainType& x
										, const RangeType& uLeft
										, RangeType& uRight ) const 
	{
  
		DomainType xgl = it.geometry().global( x );
		// uRight[0]=uLeft[0];     
  
		//v=0 
		for(int i=1;i<dimDomain+1;i++)
			uRight[i]=0.;      
		
		//Neumann Boundary for \phi and \rho
		uRight[0]=uLeft[0];
		uRight[dimDomain+1]=uLeft[dimDomain+1];
	
 
	}

} // end namespace Dune

#endif // file definition

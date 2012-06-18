#ifndef WB_MODEL_HH
#define WB_MODEL_HH



// DUNE includes
#include <dune/common/version.hh>
#include <dune/fem/misc/fmatrixconverter.hh>

// local includes
#include <dune/fem-dg/models/defaultmodel.hh>
#include "wellbalancedphysics.hh"

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

		typedef Thermodynamics< dimDomain >                       ThermodynamicsType;
	};



	template< class GridPartType , class ProblemImp >
	class PhaseModel : public DefaultModel < PhaseModelTraits< GridPartType > >
	{
  public:
		typedef ProblemImp                                        ProblemType;
		typedef typename GridPartType::GridType                   GridType;
		typedef PhaseModelTraits< GridPartType >                   Traits;

		enum { dimDomain = Traits :: dimDomain };
		enum { dimRange  = Traits :: dimRange  };
		enum { dimGradRange = Traits::dimGradRange } ;

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
		typedef typename Traits :: ThermodynamicsType             ThermodynamicsType;

	public:
		PhaseModel( const ProblemType& problem ) 
			: thermodynamics_( problem.thermodynamics() )
			, problem_( problem )
			, nsFlux_( problem )
			, visc_(Parameter::getValue<double>("visc"))
			, alpha1_(Parameter::getValue<double>("alpha1",1.))
		
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
			//{rho}[theta1]
			nsFlux_.nonConProduct(uL, uR,thetaL,thetaR,ret);
		}
		
		
		inline double stiffSource( const EntityType& en,
															 const double time,
															 const DomainType& x,
															 const RangeType& u,
															 const GradientRangeType& du,
															 const ThetaRangeType& theta,
															 const ThetaJacobianRangeType& dtheta,
															 RangeType & s) const
		{	
			nsFlux_.stiffSource(u,du,theta,dtheta,s);
			return 1;
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
			const DomainType& xgl = en.geometry().global(x);
			double mu,reaction;
			thermodynamics_.chemPotAndReaction(u,mu,reaction);
			s[0]=mu;
			s[1]=reaction;
			return 1.;
		}

// 		inline double nonStiffSource( const EntityType& en,
// 																	const double time,
// 																	const DomainType& x,
// 																	const RangeType& u,
// 																	const GradientRangeType& du,
// 																	RangeType & s) const
// 		{
// 			FieldMatrixConverter< GradientRangeType, JacobianRangeType > jac( du );
// 			return nonStiffSource( en, time, x, u, jac, s );
// 		}

// 		template< class JacobianRangeTypeImp >
// 		inline double nonStiffSource( const EntityType& en,
// 																	const double time,
// 																	const DomainType& x,
// 																	const RangeType& u,
// 																	const JacobianRangeTypeImp& jac,
// 																	RangeType& s) const
// 		{
// 			abort();
// 			return 0.;
// 		}

  

		inline void advection( const EntityType& en,
													 const double time,
													 const DomainType& x,
													 const RangeType& u,
													 JacobianRangeType& f ) const 
		{
  
			nsFlux_.analyticalFlux(u, f);
    
		}

		inline double diffusionTimeStep( const IntersectionType& it,
																		 const double enVolume,
																		 const double circumEstimate,
																		 const double time,
																		 const FaceDomainType& x,
																		 const RangeType& u ) const
		{
			// look at Ch. Merkle Diplom thesis, pg. 38

			// get value of mu at temperature T
			double p;
			double T;
// 			pressAndTemp( u, p, T );
 			// get mu 
			const double mu = problem_.mu(  );

			// ksi = 0.25 
    
			return mu * circumEstimate * 1 / (0.25 * u[0] * enVolume);
		}


		inline void jacobian( const EntityType& en
													, const double time
													, const DomainType& x
													, const RangeType& u 
													, JacobianFluxRangeType& a ) const 
		{
    
			nsFlux_.jacobian( u, a );
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
			FieldMatrixConverter< GradientRangeType, JacobianRangeType> jacLeft( gradLeft );
			return diffusionBoundaryFlux( it, time, x, uLeft, jacLeft, gLeft );
		}

		/** \brief boundary flux for the diffusion part
		 */
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
			advspeed = nsFlux_.maxSpeed( normal , u );
			totalspeed=advspeed;
		}


		void diffusion( const EntityType& en,
										const double time,
										const DomainType& x,
										const RangeType& u,
										const GradientRangeType& v,
										JacobianRangeType& diff ) const
		{
			FieldMatrixConverter< GradientRangeType, JacobianRangeType> jac( v );
			diffusion( en, time, x, u, jac,diff );
		}


		template <class JacobianRangeImp>
		void diffusion( const EntityType& en,
										const double time,
										const DomainType& x,
										const RangeType& u,
										const JacobianRangeImp& jac,
										JacobianRangeType& diff ) const
		{
			nsFlux_.diffusion( u, jac, diff );
		}
	

		
		void allenCahnDiffusion(const RangeType& u,const GradientRangeType& du,ThetaJacobianRangeType& dv ) const
		{
			FieldMatrixConverter< GradientRangeType, JacobianRangeType> jac( du );
			allenCahnDiffusion(u,jac,dv);
		}	
		template <class JacobianRangeImp>
		void allenCahnDiffusion(const RangeType& u,const JacobianRangeImp& du,ThetaJacobianRangeType& dv ) const
		{
			nsFlux_.allenCahn(u,du,dv);
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
			thermodynamics_.conservativeToPrimitive( cons, prim );
		}
		
	
		inline void totalEnergy( const DomainType& xgl,
														 const RangeType& cons, 
														 const GradientRangeType& grad,
														 FieldVector<double,1>& energy ) const
		{
    
			thermodynamics_.totalEnergy(cons,grad,energy );
		}


		inline const ProblemType& problem() const { return problem_; }
		inline const double mu( const double T ) const { return problem_.mu(T); }
		inline const double delta( ) const {return problem_.delta();}
		inline const double visc( ) const {return visc_;}
		inline const double alpha1() const {return alpha1_;}
		
	protected:
		const ThermodynamicsType& thermodynamics_;
		const ProblemType& problem_;
		const PhasePhysics< dimDomain > nsFlux_;
		const double visc_;
		const double alpha1_;
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
		for(int i=0;i<dimDomain+1;i++)
			uRight[i]=0.;      
		
		//Neumann Boundary for \phi and \rho
		uRight[0]=uLeft[0];
		uRight[dimDomain+1]=uLeft[dimDomain+1];
	
 
	}

} // end namespace Dune

#endif // file definition

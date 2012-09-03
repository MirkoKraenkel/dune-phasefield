#ifndef DUNE_PHASEFIELD_WELLBALANCEDDISCRETEMODEL_HH
#define DUNE_PHASEFIELD_WELLBALANCEDDISCRETEMODEL_HH

// Dune includes
#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>

// Dune-Fem includes
#include <dune/fem/space/dgspace.hh>
#include <dune/fem/pass/dgdiscretemodel.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/misc/boundaryidentifier.hh>
#include <dune/fem/misc/fmatrixconverter.hh>

// local includes
#include "projectiondiscretemodelcommon.hh"
#include <dune/fem/fluxes/ldgfluxwellbalanced.hh>

//*************************************************************
namespace Dune {

  // GradientModel
  //--------------
  
  template <class Model, class NumFlux, int polOrd, int passUId> /*@LST1S@*/
  class GradientModel;


  // GradientTraits
  //---------------  
  template <class Model, class NumFlux, int polOrd, int passUId >
  struct GradientTraits
  {
    typedef typename Model :: Traits                                 ModelTraits;
    typedef typename ModelTraits :: GridType                         GridType;

    enum { dimRange  = ModelTraits::dimGradRange };
    enum { dimDomain = ModelTraits::dimDomain };

    typedef PassTraits< Model, dimRange, polOrd >                    Traits;
    typedef typename Traits :: FunctionSpaceType                     FunctionSpaceType;

    typedef typename Traits :: VolumeQuadratureType                  VolumeQuadratureType;
    typedef typename Traits :: FaceQuadratureType                    FaceQuadratureType;
    typedef typename Traits :: GridPartType                          GridPartType;

    typedef typename Traits :: DiscreteFunctionSpaceType             DiscreteFunctionSpaceType;
    typedef typename Traits :: DestinationType                       DestinationType;
    typedef DestinationType                                          DiscreteFunctionType;

    typedef typename ModelTraits :: DomainType                       DomainType;
    typedef typename DestinationType :: RangeType                    RangeType;
		typedef typename DestinationType :: JacobianRangeType            JacobianRangeType;

    typedef GradientModel< Model, NumFlux, polOrd, passUId >         DGDiscreteModelType;
  };


  // GradientModel
  //--------------

  template <class Model, class NumFlux, int polOrd, int passUId> 
  class GradientModel :
    public Fem::DGDiscreteModelDefaultWithInsideOutside
		<GradientTraits<Model, NumFlux, polOrd, passUId>, passUId >
  {
    typedef Fem::DGDiscreteModelDefaultWithInsideOutside
		<GradientTraits< Model, NumFlux, polOrd, passUId >,passUId > BaseType;

    using BaseType :: inside;
    using BaseType :: outside;

    integral_constant< int, passUId > uVar;

  public:
    typedef GradientTraits< Model, NumFlux, polOrd, passUId >           Traits;

    typedef FieldVector< double, Traits :: dimDomain >                  DomainType;
    typedef FieldVector< double, Traits :: dimDomain-1 >                FaceDomainType;
    typedef typename Traits :: RangeType                                RangeType;
    typedef typename Traits :: GridType                                 GridType;
    typedef typename Traits :: JacobianRangeType                        JacobianRangeType;
    typedef typename Traits :: GridPartType:: IntersectionIteratorType  IntersectionIterator;
    typedef typename IntersectionIterator :: Intersection               Intersection;
    typedef typename GridType :: template Codim< 0 > :: Entity          EntityType;

    enum { evaluateJacobian = NumFlux :: evaluateJacobian };

    // necessary for TESTOPERATOR
    // not sure how it works for dual operators!
    static const bool ApplyInverseMassOperator = true;

  public:
    /**
     * @brief constructor
     */
    GradientModel(const Model& mod,        /*@LST1S@*/
                  const NumFlux& numf) :
      model_( mod ),
      gradientFlux_( numf ),
      cflDiffinv_( 2.0 * ( polOrd + 1) )
    {
    }

    bool hasSource() const { return false; }
    bool hasFlux() const { return true; }  /*@LST1E@*/

    template< class ArgumentTuple, class JacobianTuple > 
    inline double source( const EntityType& en,
                          const double time, 
                          const DomainType& x,
                          const ArgumentTuple& u, 
                          const JacobianTuple& jac, 
                          RangeType& s ) const
    {
      return 0.;
    }

    template <class QuadratureImp, class ArgumentTupleVector > 
    void initializeIntersection(const Intersection& it,
                                const double time,
                                const QuadratureImp& quadInner, 
                                const QuadratureImp& quadOuter,
                                const ArgumentTupleVector& uLeftVec,
                                const ArgumentTupleVector& uRightVec) 
    {
    }

    template <class QuadratureImp, class ArgumentTupleVector > 
    void initializeBoundary(const Intersection& it,
                            const double time,
                            const QuadratureImp& quadInner, 
                            const ArgumentTupleVector& uLeftVec)
    {
    }

    //! dummy method 
    void switchUpwind() const {}

    /**
     * @brief flux function on interfaces between cells for advectionand diffusion
     *
     * @param[in] it intersection
     * @param[in] time current time given by TimeProvider
     * @param[in] x coordinate of required evaluation local to \c it
     * @param[in] uLeft DOF evaluation on this side of \c it
     * @param[in] uRight DOF evaluation on the other side of \c it
     * @param[out] gLeft num. flux projected on normal on this side
     *             of \c it for multiplication with \f$ \phi \f$
     * @param[out] gRight advection flux projected on normal for the other side 
     *             of \c it for multiplication with \f$ \phi \f$
     * @param[out] gDiffLeft num. flux projected on normal on this side
     *             of \c it for multiplication with \f$ \nabla\phi \f$
     * @param[out] gDiffRight advection flux projected on normal for the other side 
     *             of \c it for multiplication with \f$ \nabla\phi \f$
     *
     * @note For dual operators we have \c gDiffLeft = 0 and \c gDiffRight = 0.
     *
     * @return wave speed estimate (multiplied with the integration element of the intersection),
     *              to estimate the time step |T|/wave.
		 */
    template <class QuadratureImp,
              class ArgumentTuple, 
              class JacobianTuple >          /*@LST0S@*/
    double numericalFlux(const Intersection& it,
                         const double time,
                         const QuadratureImp& faceQuadInner,
                         const QuadratureImp& faceQuadOuter,
                         const int quadPoint, 
                         const ArgumentTuple& uLeft,
                         const ArgumentTuple& uRight,
                         const JacobianTuple& jacLeft,
                         const JacobianTuple& jacRight,
                         RangeType& gLeft,
                         RangeType& gRight,
                         JacobianRangeType& gDiffLeft,
                         JacobianRangeType& gDiffRight ) const
    {
      return gradientFlux_.gradientNumericalFlux( it, inside(), outside(), time,
																									faceQuadInner, faceQuadOuter, quadPoint,
																									uLeft[ uVar ], uRight[ uVar ], 
																									gLeft, gRight, gDiffLeft, gDiffRight);
    }

    /**
     * @brief same as numericalFlux() but for the boundary
     */
    template <class QuadratureImp, 
              class ArgumentTuple, class JacobianTuple>
    double boundaryFlux(const Intersection& it,
                        double time, 
                        const QuadratureImp& faceQuadInner,
                        const int quadPoint,
                        const ArgumentTuple& uLeft,
                        const JacobianTuple& jacLeft,
                        RangeType& gLeft,
                        JacobianRangeType& gDiffLeft ) const   /*@LST0E@*/
    {
      const FaceDomainType& x = faceQuadInner.localPoint( quadPoint );

      typedef typename ArgumentTuple::template Get<passUId>::Type UType;
      UType uRight;
   
			if( model_.hasBoundaryValue(it,time,x) )
				{
	
					model_.boundaryValue(it, time, x, uLeft[ uVar ], uRight);
				}
			else 
				{
					uRight = uLeft[ uVar ];
				}
      
      return gradientFlux_.gradientBoundaryFlux(it, inside(),
																								time, faceQuadInner, quadPoint,
																								uLeft[ uVar ],
																								uRight, 
																								gLeft,
																								gDiffLeft );
    }

    /**
     * @brief method required by LocalDGPass
     */
    template <class ArgumentTuple, class JacobianTuple>    /*@LST1S@*/
    void analyticalFlux(const EntityType& en,
                        const double time, const DomainType& x,
                        const ArgumentTuple& u, 
                        const JacobianTuple& jac,
                        JacobianRangeType& f)
    {
      model_.jacobian(en, time, x, u[ uVar ], f);
    }                                /*@LST1E@*/

  private:
    const Model& model_;
    const NumFlux& gradientFlux_;
    const double cflDiffinv_;
  };







	// AdvectionDiffusionLDGModel
	//---------------------------

	template< class Model, 
						class NumFlux, 
						int polOrd, int passUId,int passProjId,int passGradId,
						bool advectionPartExists, bool diffusionPartExists >
	class AdvectionDiffusionLDGModel;


	// AdvectionDiffusionLDGTraits
	//----------------------------
  
	template <class Model, class NumFlux,
						int polOrd, int passUId, int passProjId, int passGradId, 
						bool advectionPartExists, bool diffusionPartExists >
	struct AdvectionDiffusionLDGTraits 
		: public AdvectionProjTraits< Model, NumFlux, polOrd, passUId,passProjId, passGradId, advectionPartExists>
          
	{
		typedef AdvectionDiffusionLDGModel
		< Model, NumFlux, polOrd, passUId,passProjId, passGradId, 
			advectionPartExists, diffusionPartExists >                   DGDiscreteModelType;
	};


	// AdvectionDiffusionLDGModel
	//---------------------------
	template< class Model, 
						class NumFlux, 
						int polOrd, int passUId, int passProjId,int passGradId,
						bool advectionPartExists, bool diffusionPartExists >
	class AdvectionDiffusionLDGModel :
		public AdvectionProjModel< Model, NumFlux, polOrd, passUId, passProjId, passGradId, advectionPartExists > 
	{
	public:
		typedef AdvectionDiffusionLDGTraits< Model, 
																				 NumFlux, 
																				 polOrd,
																				 passUId, 
																				 passProjId,
																				 passGradId,
																				 advectionPartExists, 
																				 diffusionPartExists >  Traits; 
		
		typedef AdvectionProjModel< Model, 
																NumFlux, 
																polOrd, 
																passUId,
																passProjId,
																passGradId, 
																advectionPartExists>    BaseType;
		
		using BaseType :: inside ;
		using BaseType :: outside ;
		using BaseType :: model_;
		using BaseType :: uBnd_;

		using BaseType :: uVar;
    
		// defined in AdvectionModel 
		using BaseType :: maxAdvTimeStep_ ;
		using BaseType :: maxDiffTimeStep_ ;


		integral_constant< int, passGradId> sigmaVar;

		integral_constant< int, passProjId> thetaVar;
    
	public:
		enum { dimDomain = Traits :: dimDomain };
		enum { dimRange  = Traits :: dimRange };

		enum { advection = advectionPartExists };
		enum { diffusion = diffusionPartExists };

		typedef FieldVector< double, dimDomain >               DomainType;
		typedef FieldVector< double, dimDomain-1 >             FaceDomainType;

#if defined TESTOPERATOR
		static const bool ApplyInverseMassOperator = false;
#else
		static const bool ApplyInverseMassOperator = true;
#endif

		typedef typename Traits :: GridPartType                            GridPartType;
		typedef typename Traits :: GridType                                GridType;
		typedef typename GridPartType :: IntersectionIteratorType          IntersectionIterator;
		typedef typename IntersectionIterator :: Intersection              Intersection;
		typedef typename GridType :: template Codim< 0 > :: Entity         EntityType;
		typedef typename GridType :: template Codim< 0 > :: EntityPointer  EntityPointerType;
		typedef typename Traits :: RangeFieldType                          RangeFieldType;
		typedef typename Traits :: DomainFieldType                         DomainFieldType;
		typedef typename Traits :: RangeType                               RangeType;
		typedef typename Traits :: JacobianRangeType                       JacobianRangeType;
// 		typedef typename Traits :: RangeType                               RangeType;
// 		typedef typename Traits :: JacobianRangeType                       JacobianRangeType;

	
		typedef typename Traits :: DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

		typedef LDGDiffusionFlux< DiscreteFunctionSpaceType, Model> DiffusionFluxType;
		enum { evaluateJacobian = false };
	public:
		/**
		 * @brief constructor
		 */
		AdvectionDiffusionLDGModel(const Model& mod,
															 const NumFlux& numf,
															 DiffusionFluxType& diffflux)
			: BaseType( mod, numf ),
				diffFlux_( diffflux ),
				penalty_( 1.0 ),
				cflDiffinv_( 8.0 * ( polOrd + 1) )
		{}

		bool hasSource() const
		{                
			return ( model_.hasNonStiffSource() || model_.hasStiffSource() );
		} 

		bool hasFlux() const { return advection || diffusion; };      
		/**
		 * @brief analytical flux function$
		 */
		template <class ArgumentTuple, class JacobianTuple >
		double source( const EntityType& en,
									 const double time, 
									 const DomainType& x,
									 const ArgumentTuple& u, 
									 const JacobianTuple& jac, 
									 RangeType& s ) const
		{
			s = 0;
			//we need \sigma_phi= u[sigmaVar][dimDomain+1]
			//        \theta_1=u[thetaVar][1]
			//        \nablan\theta_2 jac[thetatVar]
			double dtEst = std::numeric_limits< double > :: max();

			const double dtStiff = model_.stiffSource( en, time, x, u[uVar],u[sigmaVar],u[thetaVar],jac[thetaVar], s );


			dtEst = ( dtStiff > 0 ) ? dtStiff : dtEst;
			maxDiffTimeStep_ = std::max( dtStiff, maxDiffTimeStep_ );
	
			
			// return the fastest wave from source terms
			return dtEst;
		}

		void switchUpwind() const 
		{ 
			BaseType :: switchUpwind(); 
			diffFlux_.switchUpwind(); 
		}

	public:
		template< class QuadratureImp,
							class ArgumentTuple, 
							class JacobianTuple >          /*@LST0S@*/
		double numericalFlux(const Intersection& it,
												 const double time,
												 const QuadratureImp& faceQuadInner,
												 const QuadratureImp& faceQuadOuter,
												 const int quadPoint, 
												 const ArgumentTuple& uLeft,
												 const ArgumentTuple& uRight,
												 const JacobianTuple& jacLeft,
												 const JacobianTuple& jacRight,
												 RangeType& gLeft,
												 RangeType& gRight,
												 JacobianRangeType& gDiffLeft,
												 JacobianRangeType& gDiffRight )
		{
			// advection
								
			const double wave = BaseType :: numericalFlux( it, time, faceQuadInner, faceQuadOuter,
																										 quadPoint, uLeft, uRight, jacLeft, jacRight, 
																										 gLeft, gRight, gDiffLeft, gDiffRight );
      
			// diffusion

			double diffTimeStep = 0.0;
			if( diffusion ) 
				{
					RangeType dLeft, dRight;
					diffTimeStep = diffFlux_.numericalFlux(it, *this,
																								 time, faceQuadInner, faceQuadOuter, quadPoint,
																								 uLeft[ uVar ], uRight[ uVar ], 
																								 uLeft[ sigmaVar ], uRight[ sigmaVar ], 
																								 dLeft, dRight,
																								 gDiffLeft, gDiffRight);

					gLeft  += dLeft;
					gRight += dRight;
				}

			RangeType nonCons;
			
			model_.nonConProduct(uLeft[ uVar ], uRight[ uVar ], uLeft[ thetaVar ], uRight[ thetaVar ],nonCons); 
	
			gLeft+=nonCons;
			gRight-=nonCons;
			
			gDiffLeft  = 0;
			gDiffRight = 0;

			maxAdvTimeStep_  = std::max( wave, maxAdvTimeStep_ );
			maxDiffTimeStep_ = std::max( diffTimeStep, maxDiffTimeStep_ );

			return std::max( wave, diffTimeStep );
		}


		/**
		 * @brief same as numericalFlux() but for fluxes over boundary interfaces
		 */
		template <class QuadratureImp, 
							class ArgumentTuple, class JacobianTuple>
		double boundaryFlux(const Intersection& it,
												const double time, 
												const QuadratureImp& faceQuadInner,
												const int quadPoint,
												const ArgumentTuple& uLeft,
												const JacobianTuple& jacLeft,
												RangeType& gLeft,
												JacobianRangeType& gDiffLeft ) const   /*@LST0E@*/
		{

#if 0
			// advection

			const double wave = BaseType :: 
				boundaryFlux( it, time, faceQuadInner, quadPoint,
											uLeft, jacLeft, gLeft, gDiffLeft );
                                  
			// diffusion
      
			double diffTimeStep = 0.0;

			bool hasBoundaryValue = 
				model_.hasBoundaryValue( it, time, faceQuadInner.localPoint(0) );

			if( diffusion && hasBoundaryValue ) 
				{
					// diffusion boundary flux for Dirichlet boundaries 
					RangeType dLeft;
					diffTimeStep = diffFlux_.boundaryFlux(it, 
																								*this, 
																								time, faceQuadInner, quadPoint,
																								uLeft[ uVar ], uBnd_, // is set during call of  BaseType::boundaryFlux
																								uLeft[ sigmaVar ], 
																								uLeft[ thetaVar ],
																								dLeft,
																								gDiffLeft);
					gLeft += dLeft;
				}
			else if ( diffusion )
				{
					RangeType diffBndFlux;
					model_.diffusionBoundaryFlux( it, time, faceQuadInner.localPoint(quadPoint),
																				uLeft[uVar], jacLeft[uVar], diffBndFlux );
					gLeft += diffBndFlux;
				}
			else
				gDiffLeft = 0;

			maxAdvTimeStep_  = std::max( wave, maxAdvTimeStep_ );
			maxDiffTimeStep_ = std::max( diffTimeStep, maxDiffTimeStep_ );

			return std::max( wave, diffTimeStep );
#endif
			return 0;
		}
		/*@LST0S@*/
		/**
		 * @brief analytical flux function$
		 */
		template <class ArgumentTuple, class JacobianTuple >
		void analyticalFlux( const EntityType& en,
												 const double time, 
												 const DomainType& x,
												 const ArgumentTuple& u, 
												 const JacobianTuple& jac, 
												 JacobianRangeType& f ) const
		{
			// advection
    
			BaseType :: analyticalFlux( en, time, x, u, jac, f );

			// diffusion
   
      
			if( diffusion ) 
				{
					JacobianRangeType diffmatrix;
					JacobianRangeType tensionmatrix;
					model_.diffusion(en, time, x, u[ uVar ],u[sigmaVar], diffmatrix);
					// ldg case 
					f += diffmatrix;
					f += tensionmatrix;
				}
 
		}
	protected:
		mutable DiffusionFluxType& diffFlux_;
		const double penalty_;
		const double cflDiffinv_;
	};                                              /*@LST0E@*/

}
#include "thetadiscretemodel.hh"
#endif

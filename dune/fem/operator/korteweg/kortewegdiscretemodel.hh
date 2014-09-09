#ifndef DUNE_PHASEFIELD_WELLBALANCEDDISCRETEMODEL_HH
#define DUNE_PHASEFIELD_WELLBALANCEDDISCRETEMODEL_HH

// Dune includes
#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>

// Dune-Fem includes
#include <dune/fem/space/discontinuousgalerkin.hh>
#include <dune/fem/pass/localdg/discretemodel.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/misc/boundaryidentifier.hh>
#include <dune/fem/misc/fmatrixconverter.hh>

// local includes
#include "projectiondiscretemodelcommon.hh"


//*************************************************************
namespace Dune {

  // GradientModel
   // Discrete Model for firt approximating \sigma 
  template <class Model, class NumFlux, int polOrd, int passUId> 
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

    typedef MyPassTraits< Model, dimRange, polOrd >                    Traits;
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

    static const bool ApplyInverseMassOperator = true;

  public:
    GradientModel(const Model& mod,        
                  const NumFlux& numf) :
      model_( mod ),
      gradientFlux_( numf ),
      cflDiffinv_( 2.0 * ( polOrd + 1) )
    {
    }
 
    void setTime (double time) {}

    bool hasSource() const { return false; }
    bool hasFlux() const { return true; } 

    template< class ArgumentTuple, class JacobianTuple > 
    inline double source ( const EntityType& en,
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

    template <class QuadratureImp,
              class ArgumentTuple, 
              class JacobianTuple >         
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
       return gradientFlux_.gradientNumericalFlux(it, inside(), outside(), time,
																									faceQuadInner, faceQuadOuter, quadPoint,
																									uLeft[ uVar ], uRight[ uVar ], 
																									gLeft, gRight, gDiffLeft, gDiffRight);
    }

    template <class ArgumentTuple, class JacobianTuple> 
    void analyticalFlux(const EntityType& en,
                        const double time, const DomainType& x,
                        const ArgumentTuple& u, 
                        const JacobianTuple& jac,
                        JacobianRangeType& f)
    {
      model_.jacobian(en, time, x, u[ uVar ], f);
    }    

    template <class QuadratureImp, 
              class ArgumentTuple, class JacobianTuple>
    double boundaryFlux(const Intersection& it,
                        double time, 
                        const QuadratureImp& faceQuadInner,
                        const int quadPoint,
                        const ArgumentTuple& uLeft,
                        const JacobianTuple& jacLeft,
                        RangeType& gLeft,
                        JacobianRangeType& gDiffLeft ) const   
    {
      return boundaryFluxImpl( it, time, faceQuadInner, quadPoint,
                               uLeft[ uVar ], gLeft, gDiffLeft );
    }

  protected:
    template<class QuadratureImp,
             class UType>
    double boundaryFluxImpl( const Intersection& it,
                             double time,
                             const QuadratureImp& faceQuadInner,
                             const int quadPoint,
                             const UType& uLeft,
                             RangeType& gLeft,
                             JacobianRangeType& gDiffLeft ) const
    {
      const FaceDomainType& x=faceQuadInner.localPoint( quadPoint );
      UType uRight(0.);
      
      if( model_.hasBoundaryValue( it , time , x) )
      {
         model_.boundaryValue(it, time, x, uLeft, uRight);
      }  
      else
      {
        uRight=uLeft;
      }
    
      
      return gradientFlux_.gradientBoundaryFlux(it, inside(),
																								time, faceQuadInner, quadPoint,
																								uLeft,
																								uRight, 
																								gLeft,
																								gDiffLeft );
    }

  private:
    const Model& model_;
    const NumFlux& gradientFlux_;
    const double cflDiffinv_;
  };

	// AdvectionDiffusionLDGModel
	//---------------------------
  // for the 3rd Pass
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

		static const bool ApplyInverseMassOperator = true;

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

	
		typedef typename Traits :: DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

		typedef LDGDiffusionFlux< DiscreteFunctionSpaceType, Model> DiffusionFluxType;
		enum { evaluateJacobian = true };

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
        switch_(Fem::Parameter::getValue<double>("phasefield.thetaswitch")),
				cflDiffinv_( 8.0 * ( polOrd + 1) )
		{}
	
    void setTime (double setTime){}
    
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

			const double dtStiff = model_.stiffSource( en, time, x, u[uVar],u[sigmaVar],u[thetaVar],jac[thetaVar],jac[uVar], s );
      
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
		  const FaceDomainType& x = faceQuadInner.localPoint( quadPoint );
     
      const DomainType normal = it.integrationOuterNormal( x );
  
      double wave = BaseType :: numericalFlux( it, time, faceQuadInner, faceQuadOuter,
																										 quadPoint, uLeft, uRight, jacLeft, jacRight, 
																										 gLeft, gRight, gDiffLeft, gDiffRight );
      
			// diffusion
			double diffTimeStep = 0.0;
			if( diffusion) 
				{ 
					RangeType dLeft(0.), dRight(0.);
					diffTimeStep = diffFlux_.numericalFlux(it, *this,
																								 time, faceQuadInner, faceQuadOuter, quadPoint,
																								 uLeft[ uVar ], uRight[ uVar ], 
																								 uLeft[ sigmaVar ], uRight[ sigmaVar ], 
																								 dLeft, dRight,
																								 gDiffLeft, gDiffRight);

					gLeft  += dLeft;
					gRight += dRight;
				}

 
    
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

			// advection
			const double wave = BaseType ::boundaryFlux( it, 
                                                   time, 
                                                   faceQuadInner, 
                                                   quadPoint,
											                             uLeft, 
                                                   jacLeft, 
                                                   gLeft, 
                                                   gDiffLeft );
                                 
		  // diffusion
      
			double diffTimeStep = 0.0;

			bool hasBoundaryValue = model_.hasBoundaryValue( it, 
                                                       time, 
                                                       faceQuadInner.localPoint(0) );
	
      if( diffusion && hasBoundaryValue )
      {
        // diffusion boundary flux for Dirichlet boundaries 
				RangeType dLeft;
				diffTimeStep = diffFlux_.boundaryFlux(it, 
																							*this, 
																							time, faceQuadInner, quadPoint,
																							uLeft[ uVar ], uBnd_, // is set during call of  BaseType::boundaryFlux
																							uLeft[ sigmaVar ], 
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
      f=0;
      // advection
			BaseType :: analyticalFlux( en, time, x, u, jac, f );
			// diffusion
			if( diffusion ) 
				{
					JacobianRangeType diffmatrix(0.);
					model_.diffusion(en, time, x, u[ uVar ],u[sigmaVar], diffmatrix);
					// ldg case 
					f += diffmatrix;
				}
		}
	protected:
	 DiffusionFluxType& diffFlux_;
		const double penalty_;
    const double switch_; 
    const double cflDiffinv_;
	};                                              /*@LST0E@*/

}
#include "thetadiscretemodel.hh"
#endif

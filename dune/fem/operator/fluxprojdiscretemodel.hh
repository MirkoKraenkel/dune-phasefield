#ifndef DUNE_FEM_DG_FLUXDISCRETEMODEL_HH
#define DUNE_FEM_DG_FLUXDISCRETEMODEL_HH

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
#include <dune/fem/fluxes/ldgfluxtheta.hh>

//*************************************************************
namespace Dune {

  // GradientModel
  // Discrete Model for firt approximating \sigma 
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

    typedef MyPassTraits< Model, dimRange, polOrd >                   Traits;
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

  template <class Model, class NumFlux, int polOrd, int passUId> /*@LST1S@*/
  class GradientModel :
    public Dune::Fem::DGDiscreteModelDefaultWithInsideOutside
      <GradientTraits<Model, NumFlux, polOrd, passUId>, passUId >
  {
    typedef Dune::Fem::DGDiscreteModelDefaultWithInsideOutside
              < GradientTraits< Model, NumFlux, polOrd, passUId >,passUId > BaseType;

    using BaseType :: inside;
    using BaseType :: outside;

    // This type definition allows a convenient access to arguments of passes.
    integral_constant< int, passUId > uVar;

  public:
    typedef GradientTraits< Model, NumFlux, polOrd, passUId >        Traits;

    typedef FieldVector< double, Traits :: dimDomain >               DomainType;
    typedef FieldVector< double, Traits :: dimDomain-1 >             FaceDomainType;
    typedef typename Traits :: RangeType                             RangeType;
    typedef typename Traits :: GridType                              GridType;
    typedef typename Traits :: JacobianRangeType                     JacobianRangeType;
    typedef typename Traits :: GridPartType
              :: IntersectionIteratorType                            IntersectionIterator;
    typedef typename IntersectionIterator :: Intersection            Intersection;
    typedef typename GridType :: template Codim< 0 > :: Entity       EntityType;

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
    
       if( false /*model_.hasBoundaryValue(it,time,x)*/ )
       {
        model_.boundaryValue(it, time, x, uLeft[ uVar ], uRight);
       }
       else 
       {
        uRight = uLeft[ uVar ];
  	   } 

 double tmp= gradientFlux_.gradientBoundaryFlux(it, inside(),
						time, faceQuadInner, quadPoint,
						uLeft[ uVar ],
						uRight, 
						gLeft,
						gDiffLeft );
    

    return tmp;
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
  private:
    const Model& model_;
    const NumFlux& gradientFlux_;
    const double cflDiffinv_;
  };




  //////////////////////////////////////////////////////
  //
  // ProjectionModel
  //
  //////////////////////////////////////////////////////

	    
  template< class Model, 
	    int polOrd, int passUId, int passGradId> 
  class ProjectionModel;



 template <class Model,int polOrd, int passUId, int passGradId>
  struct ProjectionTraits
  {
    typedef typename Model :: Traits                                 ModelTraits;
    typedef typename ModelTraits :: GridType                         GridType;

    enum { dimRange = ModelTraits::dimGradRange };
    enum { dimDomain = ModelTraits::dimDomain };

    typedef MyPassTraits< Model, dimRange, polOrd >                    Traits;
    typedef typename Traits :: FunctionSpaceType                     FunctionSpaceType;

    typedef typename Traits :: VolumeQuadratureType                  VolumeQuadratureType;
    typedef typename Traits :: FaceQuadratureType                    FaceQuadratureType;
    typedef typename Traits :: GridPartType                          GridPartType;
    typedef typename Traits :: DiscreteFunctionSpaceType             DiscreteFunctionSpaceType;
    typedef typename Traits :: DestinationType                       DestinationType;
    typedef DestinationType                                          DiscreteFunctionType;
    typedef typename Traits :: IndicatorType                         IndicatorType;

    typedef typename DestinationType :: DomainType                   DomainType;
    typedef typename DestinationType :: RangeType                    RangeType;
    typedef typename DestinationType :: RangeFieldType               RangeFieldType;
    typedef typename DestinationType :: DomainFieldType              DomainFieldType;
    typedef typename DestinationType :: JacobianRangeType            JacobianRangeType;

    typedef typename Traits :: AdaptationHandlerType  AdaptationHandlerType ;

    typedef ProjectionModel< Model, polOrd, passUId, passGradId >       DGDiscreteModelType;
  };



  template<class Model,int polOrd, int passUId,int passGradId>
  class ProjectionModel:
    public Dune::Fem::DGDiscreteModelDefault
      <ProjectionTraits<Model,polOrd,passUId,passGradId>,passUId,passGradId>
  {
    typedef Dune::Fem::DGDiscreteModelDefault< ProjectionTraits< Model,polOrd,passUId,passGradId> , passUId,passGradId > BaseType ;  
    integral_constant< int, passUId > uVar;
    integral_constant< int, passGradId > sigmaVar;

  public:
    typedef ProjectionTraits< Model, polOrd, passUId,passGradId >        Traits;
    typedef FieldVector< double, Traits :: dimDomain >               DomainType;
    typedef typename Traits :: RangeType                             RangeType;
    typedef typename Traits :: JacobianRangeType                     JacobianRangeType;
    
    typedef typename Traits :: GridType                              GridType; 
    typedef typename GridType :: template Codim< 0 > :: Entity       EntityType;
 
    typedef typename Traits :: GridPartType
              :: IntersectionIteratorType                            IntersectionIterator;
    typedef typename IntersectionIterator :: Intersection            Intersection;

    enum { evaluateJacobian = false };
    static const bool ApplyInverseMassOperator = true;

 public:
    
    ProjectionModel(const Model& mod):
      model_(mod){}
      
    bool hasSource() const { return true; } 
    template <class ArgumentTuple, class JacobianTuple>
    double
    source( const EntityType& en, 
            const double time,
            const DomainType& x,
            const ArgumentTuple& u,
            const JacobianTuple& jac, 
            RangeType& s )
    { 
      model_.tension(en,time,x,u[uVar],u[sigmaVar],s);
      return 0.0;
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

    template <class QuadratureImp,class ArgumentTuple, class JacobianTuple>
    double boundaryFlux(const Intersection& it,
			double time, 
			const QuadratureImp& faceQuadInner,
                        const int quadPoint,
			const ArgumentTuple& uLeft,
			const JacobianTuple& jacLeft,
			RangeType& gLeft,
			JacobianRangeType& gDiffLeft ) const   /*@LST0E@*/
    {
      return 0;
    }

    template <class ArgumentTuple, class JacobianTuple>    /*@LST1S@*/
    void analyticalFlux(const EntityType& en,
                        const double time, const DomainType& x,
                        const ArgumentTuple& u, 
                        const JacobianTuple& jac,
                        JacobianRangeType& f)
    {
    
    }                                
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

      return 0.;
    }





  private:
    const Model& model_;
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
    typedef AdvectionDiffusionLDGTraits
    < Model, NumFlux, polOrd, passUId, passProjId,passGradId, 
          advectionPartExists, diffusionPartExists >  Traits; /*@LST0E@*/

    typedef AdvectionProjModel< Model, NumFlux, polOrd, 
                            passUId,passProjId, passGradId, advectionPartExists>    BaseType;

    using BaseType :: inside ;
    using BaseType :: outside ;
    using BaseType :: model_;
    using BaseType :: uBnd_;

    // These type definitions allow a convenient access to arguments of pass.
    using BaseType :: uVar;
    
    // defined in AdvectionModel 
    using BaseType :: maxAdvTimeStep_ ;
    using BaseType :: maxDiffTimeStep_ ;

#if DUNE_VERSION_NEWER_REV(DUNE_COMMON,2,1,0)
    integral_constant< int, passGradId> sigmaVar;
#else
    Int2Type<passGradId> sigmaVar;
#endif

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

    typedef typename Traits :: DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

    typedef ThetaDiffusionFlux< DiscreteFunctionSpaceType, Model> DiffusionFluxType;
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
    } /*@\label{dm:hasSource}@*/

    bool hasFlux() const { return advection || diffusion; };       /*@LST0S@*/

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

      double dtEst = std::numeric_limits< double > :: max();

      if (diffusion)
      {
        const double dtStiff = 
          model_.stiffSource( en, time, x, u[uVar], u[sigmaVar], s );
        dtEst = ( dtStiff > 0 ) ? dtStiff : dtEst;
        maxDiffTimeStep_ = std::max( dtStiff, maxDiffTimeStep_ );
      }

      if (model_.hasNonStiffSource())
      {

	      RangeType sNonStiff(0);
        const double dtNon = 
          model_.nonStiffSource( en, time, x, u[uVar], u[sigmaVar], sNonStiff );

        s += sNonStiff;

        dtEst = ( dtNon > 0 ) ? std::min( dtEst, dtNon ) : dtEst;
        maxAdvTimeStep_  = std::max( dtNon, maxAdvTimeStep_ );
      }

      // return the fastest wave from source terms
      return dtEst;
    }

    void switchUpwind() const 
    { 
      BaseType :: switchUpwind(); 
      diffFlux_.switchUpwind(); 
    }

  public:
    /**
     * @brief flux function on interfaces between cells for advection and diffusion
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
    
      const double wave = BaseType :: 
				numericalFlux( it, time, faceQuadInner, faceQuadOuter,
                       quadPoint, uLeft, uRight, jacLeft, jacRight, 
                       gLeft, gRight, gDiffLeft, gDiffRight );
      
      // diffusion

      double diffTimeStep = 0.0;
      if( diffusion ) 
      {
        RangeType dLeft, dRight;
        diffTimeStep = 
          diffFlux_.numericalFlux(it, *this,
                                  time, faceQuadInner, faceQuadOuter, quadPoint,
                                  uLeft[ uVar ], uRight[ uVar ], 
																	uLeft[ sigmaVar ], uRight[ sigmaVar ], 
																	uLeft[ thetaVar],uRight[thetaVar],
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
      {  abort();
        RangeType diffBndFlux;
        model_.diffusionBoundaryFlux( it, time, faceQuadInner.localPoint(quadPoint),
                                      uLeft[uVar], jacLeft[uVar], diffBndFlux );
        gLeft += diffBndFlux;
      }
      else
      {abort();
        gDiffLeft = 0;
      }
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
      // advection
    
      BaseType :: analyticalFlux( en, time, x, u, jac, f );

      // diffusion
   
      
      if( diffusion ) 
      {
        JacobianRangeType diffmatrix;
	JacobianRangeType tensionmatrix;
	model_.diffusion(en, time, x, u[ uVar ],u[sigmaVar], diffmatrix);
 	model_.tension(en,time,x,u[uVar],u[thetaVar],tensionmatrix);

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
#endif

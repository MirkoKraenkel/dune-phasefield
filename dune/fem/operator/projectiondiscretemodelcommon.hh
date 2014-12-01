#ifndef DUNE_THETA_DISCRETEMODELCOMMON_HH
#define DUNE_THETA_DISCRETEMODELCOMMON_HH

#include<limits>

// Dune includes
#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>

// Dune-Fem includes
#include <dune/fem/space/finitevolume.hh>
#include <dune/fem/space/discontinuousgalerkin.hh>
//#include <dune/fem/pass/localdg/discretemodel.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/misc/boundaryidentifier.hh>
#include <dune/fem/misc/fmatrixconverter.hh>

// local includes
//#include <dune/fem-dg/operator/limiter/limiter.hh>
#if WELLBALANCED || NSK
#include <dune/fem/fluxes/meanfluxwellbalanced.hh>
#else
#include <dune/fem/fluxes/ldgfluxtheta.hh>
#endif
#include <dune/fem-dg/operator/adaptation/adaptation.hh>

namespace Dune {

  //PassTraits
  //----------

  template <class Model,int dimRange,int polOrd>
  class MyPassTraits
  {
  public:
    typedef typename Model :: Traits                                 ModelTraits;
    typedef typename ModelTraits :: GridPartType                     GridPartType;
    typedef typename GridPartType :: GridType                        GridType;
    typedef typename GridType :: ctype                               ctype;
    static const int dimDomain = Model :: Traits :: dimDomain;

    typedef Fem::CachingQuadrature< GridPartType, 0 >                     VolumeQuadratureType;
    typedef Fem::CachingQuadrature< GridPartType, 1 >                     FaceQuadratureType;
 
    // Allow generalization to systems
    typedef Fem::FunctionSpace< ctype, double, dimDomain, dimRange >      FunctionSpaceType;
    typedef Fem::FunctionSpace< ctype, double, dimDomain, 1 >             ScalarFunctionSpaceType;
    
    typedef Fem::DiscontinuousGalerkinSpace< FunctionSpaceType,GridPartType, polOrd,Fem::CachingStorage >       DiscreteFunctionSpaceType;
    typedef Fem::DiscontinuousGalerkinSpace< FunctionSpaceType,GridPartType, polOrd,Fem::CachingStorage >       ProjectionDiscreteFunctionSpaceType;
		
    typedef Fem::DiscontinuousGalerkinSpace< ScalarFunctionSpaceType,GridPartType, polOrd,Fem::CachingStorage >   ScalarDiscreteSpaceType;
		typedef ScalarDiscreteSpaceType ScalarDiscreteFunctionSpaceType;
		typedef Fem::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType >         DestinationType;
		typedef Fem::AdaptiveDiscreteFunction< ScalarDiscreteFunctionSpaceType >   DiscreteScalarType;
		typedef DiscreteScalarType ScalarDFType;
		// Indicator for Limiter
		typedef Fem::FunctionSpace< ctype, double, dimDomain, 3> FVFunctionSpaceType;
		typedef Fem::FiniteVolumeSpace<FVFunctionSpaceType,GridPartType, 0, Fem::SimpleStorage> IndicatorSpaceType;
		typedef Fem::AdaptiveDiscreteFunction<IndicatorSpaceType> IndicatorType;
		typedef AdaptationHandler< GridType, FunctionSpaceType >  AdaptationHandlerType ;
};

// AdvectionProjModel
//---------------

template< class Model, 
	  class NumFlux, 
	  int polOrd, int passUId,int passProjId, int passGradId,
	  bool returnAdvectionPart> 
class AdvectionProjModel;


// AdvectionTraits
//----------------

template <class Model, class NumFlux,int polOrd, int passUId, int passProjId,int passGradId, bool returnAdvectionPart>
struct AdvectionProjTraits
{
  typedef typename Model :: Traits                                 ModelTraits;
  typedef typename ModelTraits :: GridType                         GridType;

  enum { dimRange = ModelTraits::dimRange };
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

  typedef AdvectionProjModel
  < Model, NumFlux, polOrd, passUId,passProjId, passGradId, returnAdvectionPart >       DGDiscreteModelType;
};


// AdvectionProjModel
//---------------
  
/*  \class AdvectionProjModel
 *
 *  \tparam Model Mathematical model
 *  \tparam NumFlux Numerical flux
 *  \tparam polOrd Polynomial degree
 *  \tparam passUId The id of a pass whose value is used here
 *  \tparam passGradId The id of a pass whose value is used here
 *  \tparam returnAdvectionPart Switch on/off the advection
 */
template< class Model, 
          class NumFlux, 
          int polOrd, int passUId, int passProjId,int passGradId,
          bool returnAdvectionPart> 
class AdvectionProjModel :
		public Fem::DGDiscreteModelDefaultWithInsideOutside
  <AdvectionProjTraits<Model, NumFlux, polOrd, passUId, passProjId, passGradId, returnAdvectionPart>,
   passUId, passProjId, passGradId>
{
public:
  typedef AdvectionProjTraits <Model, NumFlux, polOrd, passUId,  passProjId,passGradId, returnAdvectionPart> Traits;

  typedef Model   ModelType ;
  typedef NumFlux NumFluxType ;

  typedef Fem::DGDiscreteModelDefaultWithInsideOutside
  < Traits, passUId, passProjId, passGradId >                          BaseType;

  // These type definitions allow a convenient access to arguments of paesss.
  integral_constant< int, passUId >    uVar;
  integral_constant< int, passGradId > sigmaVar;    
  integral_constant< int, passProjId > thetaVar;

public:
  enum { dimDomain = Traits :: dimDomain };
  enum { dimRange  = Traits :: dimRange };

  enum { advection = returnAdvectionPart  };
  enum { evaluateJacobian = false };

  typedef FieldVector< double, dimDomain >               DomainType;
  typedef FieldVector< double, dimDomain-1 >             FaceDomainType;


  typedef typename Traits :: GridPartType                            GridPartType;
  typedef typename Traits :: GridType                                GridType;
  typedef typename GridPartType :: IntersectionIteratorType          IntersectionIterator;
  typedef typename IntersectionIterator :: Intersection              Intersection;
  typedef typename BaseType :: EntityType                            EntityType;
  typedef typename EntityType :: EntityPointer                       EntityPointerType;
  typedef typename Traits :: RangeFieldType                          RangeFieldType;
  typedef typename Traits :: DomainFieldType                         DomainFieldType;
  typedef typename Traits :: RangeType                               RangeType;
  typedef typename Traits :: JacobianRangeType                       JacobianRangeType;

  // discrete function storing the adaptation indicator information 
  typedef typename Traits :: IndicatorType          IndicatorTpye;

  // discrete function storing the adaptation indicator information 
  typedef typename Traits :: AdaptationHandlerType  AdaptationHandlerType ;

public:
  /**
   * @brief constructor
   */
  AdvectionProjModel(const Model& mod,
		     const NumFlux& numf)
    : model_(mod),
      numflux_( numf ),
      maxAdvTimeStep_( 0.0 ),
      maxDiffTimeStep_( 0.0 )
  {
  }

  //! copy constructor (for thread parallel progs mainly)
  AdvectionProjModel( const AdvectionProjModel& other )
    : BaseType( other ),
      model_( other.model_ ),
      numflux_( other.numflux_ ),
      maxAdvTimeStep_( other.maxAdvTimeStep_ ),
      maxDiffTimeStep_( other.maxDiffTimeStep_ )
  {
  }

  void setTime( double time) {}
  //! dummy method 
  void switchUpwind() const 
  {
    maxAdvTimeStep_ = 0;
    maxDiffTimeStep_ = 0;
  } 

  // cummy methods doing nothing 
  void setAdaptationHandler( AdaptationHandlerType& adaptation, double weight ) 
  {
  }

  //! remove pointer to adaptation indicator 
  void removeAdaptationHandler() 
  {
  }

  inline bool hasSource() const 
  { 
    return model_.hasNonStiffSource(); 
  }  

  inline bool hasFlux() const { return advection; }

  //! return true if diffusion time step is defining the time step size (default is false)
  double maxAdvectionTimeStep() const { return maxAdvTimeStep_; }
  double maxDiffusionTimeStep() const { return maxDiffTimeStep_; }

  //! this method is needed for thread pass 
  void setMaxTimeSteps( const double advStep, const double diffStep ) 
  {
    maxAdvTimeStep_  = advStep ;
    maxDiffTimeStep_ = diffStep ;
  }

  /**
   * @brief Stiff source associated with advection
   */
  template <class ArgumentTuple, class JacobianTuple >
  inline double source( const EntityType& en,
			const double time, 
			const DomainType& x,
			const ArgumentTuple& u, 
			const JacobianTuple& jac, 
			RangeType& s ) const
  {
    abort();
    return model_.nonStiffSource( en, time, x, u[ uVar ], s );
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

public:
  template < class QuadratureImp,
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

    gDiffLeft  = 0;
    gDiffRight = 0;

    if( advection ) 
      {
#if WELLBALANCED || NSK
				double ldt = numflux_.numericalFlux(it, 
                                            this->inside(), 
                                            this->outside(),
																						time, 
                                            faceQuadInner, 
                                            faceQuadOuter, 
                                            quadPoint, 
																						uLeft[ uVar ], 
                                            uRight[ uVar ],
                                            uLeft[ thetaVar ],
                                            uRight[ thetaVar ],
                                            gLeft, 
                                            gRight);
#else
				double ldt = numflux_.numericalFlux(it, 
                                            this->inside(), 
                                            this->outside(),
																						time, 
                                            faceQuadInner, 
                                            faceQuadOuter, 
                                            quadPoint, 
																						uLeft[ uVar ], 
                                            uRight[ uVar ],
                                            gLeft, 
                                            gRight);
#endif
				return ldt ;
      }
    else 
      {
        gLeft = 0;
        gRight = 0;
        return 0.0;
      }
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
    const FaceDomainType& x = faceQuadInner.localPoint( quadPoint );

    //see if there is a boundaryValue, if yes calculate uBnd_ by calling the model    
    const bool hasBndValue = boundaryValue( it, 
                                            time, 
                                            faceQuadInner, 
                                            quadPoint, 
                                            uLeft );

    // make sure user sets specific boundary implementation
    gLeft = std::numeric_limits< double >::quiet_NaN();
    gDiffLeft = 0;

    if (advection)
    {
      if( hasBndValue )
	    { 
        RangeType gRight;
        
#if WELLBALANCED || NSK
        double fluxret=numflux_.numericalFlux(it, this->inside(), this->inside(),
																		          time, faceQuadInner, faceQuadInner, quadPoint, 
																			        uLeft[ uVar ], uBnd_,uLeft[thetaVar],uLeft[thetaVar],gLeft, gRight);
        return fluxret;
#else
        double fluxret=numflux_.numericalFlux(it, this->inside(), this->inside(),
					                                    time, faceQuadInner, faceQuadInner, quadPoint, 
                                              uLeft[ uVar ], uBnd_, gLeft, gRight);
        return fluxret;
#endif
      }
      else 
	    {
	      return model_.boundaryFlux( it, time, x, uLeft[uVar], gLeft );
	    }
    }
    else
    {
      gLeft=0;
      return 0.;
    } 
    return 0.;
  }
  
  
  /**
   * @brief analytical flux function for advection only
   */
  template <class ArgumentTuple, class JacobianTuple >
  void analyticalFlux( const EntityType& en,
            		       const double time, 
            		       const DomainType& x,
            		       const ArgumentTuple& u, 
            		       const JacobianTuple& jac, 
            		       JacobianRangeType& f ) const
  {
    if( advection ) 
    {
      model_.advection(en, time, x, u[ uVar ], f);
    }
    else 
    {
      f = 0;
    }
  }


protected:
  template <class QuadratureImp, 
	    class ArgumentTuple>
  bool boundaryValue(const Intersection& it,
		     const double time, 
		     const QuadratureImp& faceQuadInner,
		     const int quadPoint,
		     const ArgumentTuple& uLeft ) const 
  {
    const FaceDomainType& x = faceQuadInner.localPoint( quadPoint );
    const bool hasBndValue = model_.hasBoundaryValue(it, time, x);
    if( hasBndValue ) 
      { 
        model_.boundaryValue(it, time, x, uLeft[ uVar ], uBnd_ );
      }
    else 
      // do something bad to uBnd_ as it shouldn't be used
      uBnd_ = std::numeric_limits< double >::quiet_NaN();

    return hasBndValue;
  }

  const Model&   model_;
  const NumFlux& numflux_;
  mutable RangeType uBnd_;
  mutable double maxAdvTimeStep_;
  mutable double maxDiffTimeStep_;
};                                              /*@LST0E@*/

///END ADVPROMOD//////////////////////////////////////////////////////////


  // AdvectionModel
  //---------------

template< class Model, 
	  class NumFlux, 
	  int polOrd, int passUId, int passGradId,
	  bool returnAdvectionPart> 
class AdvectionModel;


// AdvectionTraits
//----------------

template <class Model, class NumFlux,
	  int polOrd, int passUId, int passGradId, bool returnAdvectionPart>
struct AdvectionTraits
{
  typedef typename Model :: Traits                                 ModelTraits;
  typedef typename ModelTraits :: GridType                         GridType;

  enum { dimRange = ModelTraits::dimRange };
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

  typedef AdvectionModel
  < Model, NumFlux, polOrd, passUId, passGradId, returnAdvectionPart >       DGDiscreteModelType;
};


// AdvectionModel
//---------------
  
/*  \class AdvectionModel
 *
 *  \tparam Model Mathematical model
 *  \tparam NumFlux Numerical flux
 *  \tparam polOrd Polynomial degree
 *  \tparam passUId The id of a pass whose value is used here
 *  \tparam passGradId The id of a pass whose value is used here
 *  \tparam returnAdvectionPart Switch on/off the advection
 */
template< class Model, 
	  class NumFlux, 
	  int polOrd, int passUId, int passGradId,
	  bool returnAdvectionPart> 
class AdvectionModel :
  public Fem::DGDiscreteModelDefaultWithInsideOutside
  <AdvectionTraits<Model, NumFlux, polOrd, passUId, passGradId, returnAdvectionPart>,
   passUId, passGradId>
{
public:
  typedef AdvectionTraits 
  <Model, NumFlux, polOrd, passUId, passGradId, returnAdvectionPart> Traits;

  typedef Model   ModelType ;
  typedef NumFlux NumFluxType ;

  typedef Fem::DGDiscreteModelDefaultWithInsideOutside
  < Traits, passUId, passGradId >                          BaseType;

  // These type definitions allow a convenient access to arguments of paesss.
  integral_constant< int, passUId > uVar;
  integral_constant< int, passGradId > sigmaVar;    
public:
  enum { dimDomain = Traits :: dimDomain };
  enum { dimRange  = Traits :: dimRange };

  enum { advection = returnAdvectionPart  };
  enum { evaluateJacobian = false };

  typedef FieldVector< double, dimDomain >               DomainType;
  typedef FieldVector< double, dimDomain-1 >             FaceDomainType;


  typedef typename Traits :: GridPartType                            GridPartType;
  typedef typename Traits :: GridType                                GridType;
  typedef typename GridPartType :: IntersectionIteratorType          IntersectionIterator;
  typedef typename IntersectionIterator :: Intersection              Intersection;
  typedef typename BaseType :: EntityType                            EntityType;
  typedef typename EntityType :: EntityPointer                       EntityPointerType;
  typedef typename Traits :: RangeFieldType                          RangeFieldType;
  typedef typename Traits :: DomainFieldType                         DomainFieldType;
  typedef typename Traits :: RangeType                               RangeType;
  typedef typename Traits :: JacobianRangeType                       JacobianRangeType;

  // discrete function storing the adaptation indicator information 
  typedef typename Traits :: IndicatorType          IndicatorTpye;

  // discrete function storing the adaptation indicator information 
  typedef typename Traits :: AdaptationHandlerType  AdaptationHandlerType ;

public:
  /**
   * @brief constructor
   */
  AdvectionModel(const Model& mod,
		 const NumFlux& numf)
    : model_(mod),
      numflux_( numf ),
      maxAdvTimeStep_( 0.0 ),
      maxDiffTimeStep_( 0.0 )
  {
  }

  //! copy constructor (for thread parallel progs mainly)
  AdvectionModel( const AdvectionModel& other )
    : BaseType( other ),
      model_( other.model_ ),
      numflux_( other.numflux_ ),
      maxAdvTimeStep_( other.maxAdvTimeStep_ ),
      maxDiffTimeStep_( other.maxDiffTimeStep_ )
  {
  }

  //! dummy method 
  void switchUpwind() const {
    maxAdvTimeStep_ = 0;
    maxDiffTimeStep_ = 0;
  } 

  // dummy methods doing nothing 
  void setAdaptationHandler( AdaptationHandlerType& adaptation, double weight ) 
  {
  }

  //! remove pointer to adaptation indicator 
  void removeAdaptationHandler() 
  {
  }

  inline bool hasSource() const 
  { 
    return model_.hasNonStiffSource(); 
  }  /*@\label{dm:hasSource}@*/

  inline bool hasFlux() const { return advection; }

  //! return true if diffusion time step is defining the time step size (default is false)
  double maxAdvectionTimeStep() const { return maxAdvTimeStep_; }
  double maxDiffusionTimeStep() const { return maxDiffTimeStep_; }

  //! this method is needed for thread pass 
  void setMaxTimeSteps( const double advStep, const double diffStep ) 
  {
    maxAdvTimeStep_  = advStep ;
    maxDiffTimeStep_ = diffStep ;
  }

  /**
   * @brief Stiff source associated with advection
   */
  template <class ArgumentTuple, class JacobianTuple >
  inline double source( const EntityType& en,
			const double time, 
			const DomainType& x,
			const ArgumentTuple& u, 
			const JacobianTuple& jac, 
			RangeType& s ) const
  {
   
    return model_.nonStiffSource( en, time, x, u[ uVar ], s );
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
    gDiffLeft = 0;
    gDiffRight = 0;

    if( advection ) 
      {
    	  double ldt = numflux_.numericalFlux(it, 
                                            this->inside(), 
                                            this->outside(),
                                            time, 
                                            faceQuadInner, 
                                            faceQuadOuter, 
                                            quadPoint, 
                                     				uLeft[ uVar ], 
                                            uRight[ uVar ], 
                                            uLeft[sigmaVar],
                                            uRight[sigmaVar],
                                            gLeft, 
                                            gRight);

	      return ldt ;
      }
    else 
      {
        gLeft  = 0;
        gRight = 0;
        return 0.0;
      }
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
		                  JacobianRangeType& gDiffLeft ) const   
  {
     
    const FaceDomainType& x = faceQuadInner.localPoint( quadPoint );

    const bool hasBndValue = boundaryValue( it, 
                                            time,
                                            faceQuadInner, 
                                            quadPoint,
					                                  uLeft );
      
    // make sure user sets specific boundary implementation
    gLeft = std::numeric_limits< double >::quiet_NaN();
   
    gDiffLeft = 0;
    
    if (advection)
    {
      if( hasBndValue )
      {
        RangeType gRight;

        return numflux_.numericalFlux(it, 
                                      this->inside(), 
                                      this->inside(),
                                      time, 
                                      faceQuadInner, 
                                      faceQuadInner, 
                                      quadPoint, 
					                            uLeft[ uVar ], 
                                      uBnd_,
                                      uLeft[sigmaVar],
                                      uLeft[sigmaVar], 
                                      gLeft, 
                                      gRight);
      }
      else 
      {
        return model_.boundaryFlux( it, time, x, uLeft[uVar], gLeft );
      }
    }
    else
    {
      gLeft = 0.;
      return 0.;
    }
    
    return 0.;
  }




  template <class ArgumentTuple, class JacobianTuple >
  void analyticalFlux( const EntityType& en,
											 const double time, 
            		       const DomainType& x,
            		       const ArgumentTuple& u, 
            		       const JacobianTuple& jac, 
            		       JacobianRangeType& f ) const
  {
      
    if( advection ) 
      {
				model_.advection(en, time, x, u[ uVar ],u[sigmaVar], f);
      }
    else 
      {
				f = 0;
      }
  }


protected:
  template <class QuadratureImp,class ArgumentTuple>
  bool boundaryValue(const Intersection& it,
		     const double time, 
		     const QuadratureImp& faceQuadInner,
		     const int quadPoint,
		     const ArgumentTuple& uLeft ) const 
  {
    const FaceDomainType& x = faceQuadInner.localPoint( quadPoint );
    const bool hasBndValue = model_.hasBoundaryValue(it, time, x);
    if( hasBndValue ) 
      { 
        model_.boundaryValue(it, time, x, uLeft[ uVar ], uBnd_ );
      }
    else 
      // do something bad to uBnd_ as it shouldn't be used
      uBnd_ = std::numeric_limits< double >::quiet_NaN();

    return hasBndValue;
  }

  const Model&   model_;
  const NumFlux& numflux_;
  mutable RangeType uBnd_;
  mutable double maxAdvTimeStep_;
  mutable double maxDiffTimeStep_;
};                                             

} // end namespace Dune

#endif

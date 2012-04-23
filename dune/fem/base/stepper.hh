#ifndef NAVIER_STOKES_STEPPER_HH
#define NAVIER_STOKES_STEPPER_HH
#warning "Orginal Stepper"
#include "stepperbase.hh"

template <class GridImp,
          class ProblemTraits, 
          int polynomialOrder>             
struct Stepper 
  : public StepperBase< GridImp, ProblemTraits, polynomialOrder >
{
  typedef StepperBase< GridImp, ProblemTraits, polynomialOrder > BaseType ;

  // type of Grid
  typedef typename BaseType :: GridType                 GridType;

  // Choose a suitable GridView
  typedef typename BaseType :: GridPartType             GridPartType;

  // initial data type 
  typedef typename BaseType :: InitialDataType          InitialDataType;

  // An analytical version of our model
  typedef typename BaseType :: ModelType                 ModelType;

  // The flux for the discretization of advection terms
  typedef typename BaseType :: FluxType                  FluxType;

  // The DG space operator
  // The first operator is sum of the other two
  // The other two are needed for semi-implicit time discretization
  typedef typename BaseType :: DgType                    DgType;
  typedef typename BaseType :: DgAdvectionType           DgAdvectionType;
  typedef typename BaseType :: DgDiffusionType           DgDiffusionType;

  // The discrete function for the unknown solution is defined in the DgOperator
  typedef typename BaseType :: DiscreteFunctionType      DiscreteFunctionType;

  // The indicator function in case of limiting 
  //typedef typename BaseType :: IndicatorType             IndicatorType;

  // ... as well as the Space type
  typedef typename BaseType :: DiscreteSpaceType         DiscreteSpaceType;

  // The ODE Solvers
  typedef typename BaseType :: OdeSolverType     OdeSolverType;

  typedef typename BaseType :: TimeProviderType       TimeProviderType;
  typedef typename BaseType :: AdaptationManagerType  AdaptationManagerType;
  typedef typename BaseType :: AdaptationHandlerType  AdaptationHandlerType; 

  static const Dune::DGDiffusionFluxIdentifier DiffusionFluxId =
    BaseType::Traits::DiffusionFluxId ;

  // advection = true , diffusion = true
  typedef Dune :: DGAdaptationIndicatorOperator< ModelType, FluxType,
            DiffusionFluxId, polynomialOrder, true, true >  DGIndicatorType;

  //typedef Estimator< DiscreteFunctionType, InitialDataType >        EstimatorType;

  using BaseType :: grid_;
  using BaseType :: space;
  using BaseType :: convectionFlux_ ;
  using BaseType :: problem;
  using BaseType :: solution_;
  using BaseType :: adaptationHandler_ ;
  using BaseType :: adaptationParameters_;
  using BaseType :: adaptive_ ;

  Stepper(GridType& grid) :
    BaseType( grid ),
    dgOperator_(grid_, convectionFlux_),
    dgAdvectionOperator_(grid_, convectionFlux_),
    dgDiffusionOperator_(grid_, convectionFlux_),
    dgIndicator_( grid_, convectionFlux_ )
  {
  }                                                                        /*@LST1E@*/

  virtual size_t gridSize() const 
  {
    // one of them is not zero, 
    // use int because the unintialized size_t is the largest 
    size_t advSize  = dgAdvectionOperator_.numberOfElements();
    size_t diffSize = dgDiffusionOperator_.numberOfElements();
    size_t grSize   = std::max( advSize, diffSize );
    grSize = std::max( dgOperator_.numberOfElements(), grSize );
    return grid_.comm().sum( grSize );
  }

  virtual OdeSolverType* createOdeSolver(TimeProviderType& tp) 
  {
    if( adaptive_ )
    {
      if( ! adaptationHandler_ && adaptationParameters_.aposterioriIndicator() )
      {
        adaptationHandler_ = new AdaptationHandlerType( grid_, tp );
        dgIndicator_.setAdaptationHandler( *adaptationHandler_ );
      }
    }

    typedef SmartOdeSolver< DgType, DgAdvectionType, DgDiffusionType > OdeSolverImpl;
    return new OdeSolverImpl( tp, dgOperator_, 
                              dgAdvectionOperator_,
                              dgDiffusionOperator_ );
  }

  //! estimate and mark solution 
  virtual void estimateMarkAdapt( AdaptationManagerType& am ) 
  {
    if( adaptationHandler_ )
    {
      DiscreteFunctionType tmp( solution_ );
      // call operator once to calculate indicator 
      dgIndicator_( solution_, tmp );
      adaptationHandler_->adapt( am );
    }
    else 
    {
      //EstimatorType estimator( solution_, problem() );
      //estimator.estimateAndMark();
      am.adapt();
    }
  }

protected:
  DgType                  dgOperator_;
  DgAdvectionType         dgAdvectionOperator_;
  DgDiffusionType         dgDiffusionOperator_;
  DGIndicatorType         dgIndicator_;
};
#endif // FEMHOWTO_STEPPER_HH

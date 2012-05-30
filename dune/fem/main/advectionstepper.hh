#ifndef NAVIER_STOKES_STEPPER_HH
#define NAVIER_STOKES_STEPPER_HH

#include <dune/fem-dg/operator/adaptation/estimator.hh>
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
  typedef typename BaseType :: DgAdvectionType           DgType;
  typedef DgType  DgAdvectionType;
  typedef DgType  DgDiffusionType;

  // The discrete function for the unknown solution is defined in the DgOperator
  typedef typename BaseType :: DiscreteFunctionType      DiscreteFunctionType;

  // ... as well as the Space type
  typedef typename BaseType :: DiscreteSpaceType         DiscreteSpaceType;

  // The ODE Solvers
  typedef typename BaseType :: OdeSolverType     OdeSolverType;

  typedef typename BaseType :: TimeProviderType        TimeProviderType;
  typedef typename BaseType :: AdaptationManagerType   AdaptationManagerType;
  typedef typename BaseType :: AdaptationHandlerType   AdaptationHandlerType;

  static const Dune::DGDiffusionFluxIdentifier DiffusionFluxId =
    BaseType::Traits::DiffusionFluxId ;
  // advection = true , diffusion = false 
  typedef Dune :: DGAdaptationIndicatorOperator< ModelType, FluxType,
            DiffusionFluxId, polynomialOrder, true, false >  DGIndicatorType;

  //typedef Estimator< DiscreteFunctionType, InitialDataType >        EstimatorType;

  using BaseType :: grid_;
  using BaseType :: convectionFlux_ ;
  using BaseType :: problem;
  using BaseType :: adaptationHandler_ ;
  using BaseType :: solution_ ;
  using BaseType :: adaptive_ ;
  using BaseType :: adaptationParameters_;

  Stepper(GridType& grid) :
    BaseType( grid ),
    dgAdvectionOperator_(grid_, convectionFlux_),
    dgIndicator_( grid_, convectionFlux_ )
  {
  }                                                                        /*@LST1E@*/

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

    typedef SmartOdeSolver< DgAdvectionType, DgAdvectionType, DgAdvectionType > OdeSolverImpl;
    return new OdeSolverImpl( tp, dgAdvectionOperator_, 
                              dgAdvectionOperator_,
                              dgAdvectionOperator_ );
  }

  //! call limiter (only if dgAdvectionOperator_ is DGLimitedAdvectionOperator)
  void limitSolution() 
  { 
    dgAdvectionOperator_.limit( solution_ );
  }

  //! estimate and mark solution 
  virtual void estimateMarkAdapt( AdaptationManagerType& am )
  {
    Estimator< DiscreteFunctionType, InitialDataType > gradientIndicator( solution_, problem() );
    doEstimateMarkAdapt( dgIndicator_, gradientIndicator, am );
  }

protected:
  DgAdvectionType         dgAdvectionOperator_;
  DGIndicatorType         dgIndicator_;
};
#endif // FEMHOWTO_STEPPER_HH

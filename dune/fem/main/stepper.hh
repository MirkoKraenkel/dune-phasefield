#ifndef NAVIER_STOKES_STEPPER_HH
#define NAVIER_STOKES_STEPPER_HH

#define NS_USE_ADAPTATION 0


#include "stepperbase.hh"
#include <dune/fem/adaptation/estimator2.hh>
//#include "estimator1.hh"




struct MyDataOutputParameters:
  public LocalParameter< DataOutputParameters, MyDataOutputParameters >
{
  virtual std::string path() const
  {
    return Parameter::getValue< std::string >( "theta.io.path", "./" );
  }


};







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
  typedef typename BaseType :: DiscreteGradientType      DiscreteGradientType;
  typedef typename BaseType :: ScalarDFType      ScalarDFType;
  
  // The indicator function in case of limiting 
  //typedef typename BaseType :: IndicatorType             IndicatorType;

  // ... as well as the Space type
  typedef typename BaseType :: DiscreteSpaceType         DiscreteSpaceType;
  typedef typename BaseType :: ScalarDiscreteSpaceType         ScalarDiscreteSpaceType;
  
  // The ODE Solvers
  typedef typename BaseType :: OdeSolverType     OdeSolverType;

  typedef typename BaseType :: TimeProviderType       TimeProviderType;
  typedef typename BaseType :: AdaptationManagerType  AdaptationManagerType;
  typedef typename BaseType :: AdaptationHandlerType  AdaptationHandlerType; 

  static const Dune::DGDiffusionFluxIdentifier DiffusionFluxId =    BaseType::Traits::DiffusionFluxId ;
  
  
  typedef Estimator1< DiscreteFunctionType >        EstimatorType;

  using BaseType :: grid_;
  using BaseType :: space;
  using BaseType :: convectionFlux_ ;
  using BaseType :: problem;
  using BaseType :: model;
  using BaseType :: solution_;
	using BaseType :: theta_;
	using BaseType :: energy_;
  using BaseType :: adaptationHandler_ ;
  using BaseType :: overallTimer_;
  using BaseType :: adaptationParameters_;
  using BaseType :: odeSolver_;
  using BaseType :: odeSolverMonitor_;
  using BaseType :: adaptive_ ;


  Stepper(GridType& grid) :
    BaseType( grid ),
    dgOperator_(grid_, convectionFlux_),
    dgAdvectionOperator_(grid_, convectionFlux_),
    dgDiffusionOperator_(grid_, convectionFlux_),
    tolerance_(Dune::Parameter::getValue<double>("tolerance",1e-2) ) 
  {
  }                                                                        /*@LST1E@*/
  
  virtual DiscreteFunctionType& solution() { return solution_; }
	
	virtual DiscreteGradientType& theta() { return theta_; }
	
	virtual ScalarDFType& energy() { return energy_; }
  
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
#if NS_USE_ADAPTATION
    if( adaptive_ )
      {
				if( ! adaptationHandler_ && adaptationParameters_.aposterioriIndicator() )
					{
						adaptationHandler_ = new AdaptationHandlerType( grid_, tp );
						dgIndicator_.setAdaptationHandler( *adaptationHandler_ );
					}
      }
#endif
		
    typedef SmartOdeSolver< DgType, DgAdvectionType, DgDiffusionType > OdeSolverImpl;
    return new OdeSolverImpl( tp, dgOperator_, 
															dgAdvectionOperator_,
															dgDiffusionOperator_ );
  }

  
	void step(TimeProviderType& tp, //DiscreteFunctionType& U,
            int& newton_iterations, 
            int& ils_iterations,
            int& max_newton_iterations,
            int& max_ils_iterations) 
  {
    DiscreteFunctionType& U = solution();
		DiscreteGradientType& theta1=theta();
		ScalarDFType& energy1=energy();

    // reset overall timer
    overallTimer_.reset();

    assert(odeSolver_);
    odeSolver_->solve( U, odeSolverMonitor_ );
    
		dgOperator_.gradient(U,theta1);
		dgOperator_.energy(U,theta1, energy1);
   

#if 0
    typedef Dune::tuple< DiscreteFunctionType* > MyIOTupleType;
   
    MyIOTupleType iotup(  &theta1 );
    
    typedef Dune::DataOutput< GridType, MyIOTupleType >    MyDataWriterType;
    MyDataWriterType mywrite(grid_, iotup,MyDataOutputParameters());
  
    mywrite.write(tp);   
		//limitSolution ();
#endif
    newton_iterations     = odeSolverMonitor_.newtonIterations_;
    ils_iterations        = odeSolverMonitor_.linearSolverIterations_;
    max_newton_iterations = odeSolverMonitor_.maxNewtonIterations_ ;
    max_ils_iterations    = odeSolverMonitor_.maxLinearSolverIterations_;

  }


  //! estimate and mark solution 
  virtual void estimateMarkAdapt( AdaptationManagerType& am ) 
  {

    EstimatorType estimator( solution_,grid_ );
    estimator.estimateAndMark(tolerance_);
    am.adapt();


  }

protected:
  DgType                  dgOperator_;
  DgAdvectionType         dgAdvectionOperator_;
  DgDiffusionType         dgDiffusionOperator_;
#if NS_USE_ADAPTATION
  DGIndicatorType         dgIndicator_;
#endif
  double tolerance_;

};
#endif // FEMHOWTO_STEPPER_HH

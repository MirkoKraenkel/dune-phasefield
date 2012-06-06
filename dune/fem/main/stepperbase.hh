#ifndef FEMHOWTO_DGSTEPPERBASE_HH
#define FEMHOWTO_DGSTEPPERBASE_HH
#include <config.h>

#ifdef LOCALDEBUG
static double sum_ = 0.;
static double sum2_ = 0.;
static double localMaxRatio_ = 0.;
static double localMinRatio_ = 1e+100;
static double maxRatioOfSums = 0.;
static double minRatioOfSums = 1e+100;
#endif


#ifndef NDEBUG 
// enable fvector and fmatrix checking
#define DUNE_ISTL_WITH_CHECKING
#endif

// include std libs
#include <iostream>
#include <string>

// Dune includes
#include <dune/fem/misc/l2error.hh>
#include <dune/fem/operator/projection/l2projection.hh>
#include <dune/fem/gridpart/common/gridpart.hh>

#include <dune/fem-dg/solver/smartodesolver.hh>

#ifdef BASEFUNCTIONSET_CODEGEN_GENERATE
#include <dune/fem/space/basefunctions/codegen.hh>
#endif

#if NONCON
#include <dune/fem/operator/fluxoperator.hh>
#elif WELLBALANCED
#include <dune/fem/operator/wellbalancedoperator.hh>
#else
#include <dune/fem/operator/fluxprojoperator.hh>
#endif
// include local header files
#include <dune/fem/base/baseevolution.hh>

#include "steppertraits.hh"

#include <dune/fem-dg/misc/runfile.hh>
#include <dune/fem-dg/operator/adaptation/estimatorbase.hh>

using namespace Dune;                                        


template <class GridImp,
          class ProblemTraits, 
          int polynomialOrder>             
struct StepperBase
  : public AlgorithmBase< StepperTraits< GridImp, ProblemTraits, polynomialOrder> > 
{
  // my traits class 
  typedef StepperTraits< GridImp, ProblemTraits, polynomialOrder> Traits ;
  
	typedef AlgorithmBase< Traits > BaseType;


  // type of Grid
  typedef typename Traits :: GridType                 GridType;

  // Choose a suitable GridView
  typedef typename Traits :: GridPartType             GridPartType;

  // initial data type 
  typedef typename Traits :: InitialDataType          InitialDataType;

  // An analytical version of our model
  typedef typename Traits :: ModelType                 ModelType;

  // The flux for the discretization of advection terms
  typedef typename Traits :: FluxType                  FluxType;

  // The DG space operator
  // The first operator is sum of the other two
  // The other two are needed for semi-implicit time discretization
  typedef typename Traits :: DgType                    DgType;
  typedef typename Traits :: DgAdvectionType           DgAdvectionType;
  typedef typename Traits :: DgDiffusionType           DgDiffusionType;

  // The discrete function for the unknown solution is defined in the DgOperator
  typedef typename Traits :: DiscreteFunctionType      DiscreteFunctionType;
  typedef typename Traits :: DiscreteGradientType      DiscreteGradientType;
  typedef typename Traits :: ScalarDFType      ScalarDFType;
  
  // ... as well as the Space type
 
  typedef typename Traits :: DiscreteSpaceType         DiscreteSpaceType;
  typedef typename Traits :: ScalarDiscreteSpaceType    ScalarDiscreteSpaceType;
  

  // The ODE Solvers
  typedef typename Traits :: OdeSolverType       OdeSolverType;
  typedef typename OdeSolverType :: MonitorType  OdeSolverMonitorType ;

  typedef typename BaseType :: TimeProviderType       TimeProviderType;
  typedef typename BaseType :: AdaptationManagerType  AdaptationManagerType;

  typedef AdaptationHandler< GridType, 
           typename DiscreteSpaceType::FunctionSpaceType >  AdaptationHandlerType;

  typedef RunFile< GridType >  RunFileType;

  using BaseType :: grid_;
  using BaseType :: gridPart_;
  using BaseType :: solution;
	using BaseType :: theta;
	using BaseType :: energy;
  using BaseType :: space;
	using BaseType :: energyspace;
  using BaseType :: limitSolution ;

  // constructor taking grid 
  StepperBase(GridType& grid) :
    BaseType( grid ),
    solution_( "solution", space() ),
		theta_( "theta", space() ),
		energy_("energy",energyspace()),
    additionalVariables_( Parameter :: getValue< bool >("femhowto.additionalvariables", false) ? 
													new DiscreteFunctionType("additional", space() ) : 0 ),
    problem_( ProblemTraits::problem() ),
    model_( new ModelType( problem() ) ),
    convectionFlux_( *model_ ),
    adaptationHandler_( 0 ),
    runfile_( grid.comm(), true ),
    overallTimer_(),
    eocId_( FemEoc::addEntry(std::string("$L^2$-error")) ),
    odeSolver_( 0 ),
    adaptive_( Dune::AdaptationMethod< GridType >( grid_ ).adaptive() ),
    adaptationParameters_( )
  {
  }                                                                      

  //! destructor 
  virtual ~StepperBase()
  {
    delete odeSolver_;
    odeSolver_ = 0;
    delete problem_ ;
    problem_ = 0;
    delete adaptationHandler_ ;
    adaptationHandler_ = 0;
    delete additionalVariables_; 
    additionalVariables_ = 0;
  }

  // return reference to discrete function holding solution 
  virtual DiscreteFunctionType& solution() { return solution_; }

	virtual DiscreteGradientType& theta() { return theta_; }
  virtual ScalarDFType& energy(){return energy_;}

  // return pointer to additional variables, can be zero 
  virtual DiscreteFunctionType* additionalVariables() { return additionalVariables_; }

  // function creating the ode solvers 
  virtual OdeSolverType* createOdeSolver( TimeProviderType& ) = 0;

  virtual void writeCheckPoint(TimeProviderType& tp,
                               AdaptationManagerType& am ) const 
  {
    assert( odeSolver_ );

    const double ldt = tp.deltaT();
    // write times to run file 
    #if 0  
    runfile_.write( tp.time() + ldt, ldt, 
                    odeSolverMonitor_.numberOfElements_,
                    odeSolverMonitor_.operatorTime_, 
                    odeSolverMonitor_.odeSolveTime_, 
                    am.adaptationTime(),
                    am.loadBalanceTime(), 
                    overallTimer_.elapsed());

#endif
  }

  //! returns data prefix for EOC loops ( default is loop )
  virtual std::string dataPrefix() const 
  {
    return problem_->dataPrefix();
  }

  // gather information from the space operator, the time integratior
  // and the problem to output before each table in tex file
  std::string description() const
  {
    std::string latexInfo = odeSolver_->description();
    latexInfo +=  problem_->description() + "\n\n";
    return latexInfo;
  }

  // before first step, do data initialization 
  void initializeStep(TimeProviderType& tp) 
    //, DiscreteFunctionType& U)
  {
    DiscreteFunctionType& U = solution();

    if( odeSolver_ == 0 ) odeSolver_ = this->createOdeSolver( tp );

#if defined TESTOPERATOR
//    analyseOperator( dgDiffusionOperator_ );
#endif
    assert( odeSolver_ );

    L2Projection< double, double,                         /*@\label{dg:l2pro0}@*/
                  InitialDataType, DiscreteFunctionType > l2pro;
    l2pro(problem(), U );                                   /*@\label{dg:l2pro1}@*/

    // ode.initialize applies the DG Operator once to get an initial
    // estimate on the time step. This does not change the initial data u.
    odeSolver_->initialize( U );                               /*@\label{dg:odeinit}@*/
  }
;

    void step(TimeProviderType& tp, 
            int& newton_iterations, 
            int& ils_iterations,
            int& max_newton_iterations,
	      int& max_ils_iterations) =0;
#if 0

  void step(TimeProviderType& tp, //DiscreteFunctionType& U,
            int& newton_iterations, 
            int& ils_iterations,
            int& max_newton_iterations,
            int& max_ils_iterations) 
  {
    DiscreteFunctionType& U = solution();

    // reset overall timer
    overallTimer_.reset();

    assert(odeSolver_);
    odeSolver_->solve( U, odeSolverMonitor_ );
						
    
    
    limitSolution ();

    newton_iterations     = odeSolverMonitor_.newtonIterations_;
    ils_iterations        = odeSolverMonitor_.linearSolverIterations_;
    max_newton_iterations = odeSolverMonitor_.maxNewtonIterations_ ;
    max_ils_iterations    = odeSolverMonitor_.maxLinearSolverIterations_;

  }
#endif
  
  inline double error(TimeProviderType& tp, DiscreteFunctionType& u)
  {
    typedef typename DiscreteFunctionType :: RangeType RangeType;
    L2Error<DiscreteFunctionType> L2err;
    // Compute L2 error of discretized solution ...
    RangeType error = L2err.norm(problem(), u, tp.time());
    return error.two_norm();
  }

  void finalizeStep(TimeProviderType& tp)
  {
    DiscreteFunctionType& u = solution();
    // write run file (in writeatonce mode)
    runfile_.flush(); 

    bool doFemEoc = problem().calculateEOC( tp, u, eocId_ );

    // ... and print the statistics out to a file
    if( doFemEoc )
      FemEoc::setErrors(eocId_, error(tp, u ));

#ifdef LOCALDEBUG
    std::cout <<"maxRatioOfSums: " <<maxRatioOfSums <<std::endl;
    std::cout <<"minRatioOfSums: " <<minRatioOfSums <<std::endl;
#endif    

    // delete ode solver
    delete odeSolver_;
    odeSolver_ = 0;
    delete adaptationHandler_;
    adaptationHandler_ = 0;
  }                                                       /*@LST1E@*/

  const InitialDataType& problem() const 
  { 
    assert( problem_ ); 
    return *problem_; 
  }


  const ModelType& model() const 
  { 
    assert( model_ ); 
    return *model_;
  }

protected:
#if 0
  template <class IndicatorOperator, class GradientEstimator>
  void doEstimateMarkAdapt( const IndicatorOperator& dgIndicator,
                            GradientEstimator& gradientEstimator,
                            AdaptationManagerType& am ) 
  {
    if( adaptive_ )
    {
      // get grid sequence before adaptation 
      const int sequence = solution_.space().sequence();

      if( adaptationHandler_ )
      {
        // call operator once to calculate indicator 
        dgIndicator.evaluateOnly( solution_ );

        // do marking and adaptation 
        adaptationHandler_->adapt( am );
      }
      else if( adaptationParameters_.gradientBasedIndicator() )
      {
        gradientEstimator.estimateAndMark();
        am.adapt();
      }
      else if( adaptationParameters_.shockIndicator() )
      {
        // marking has been done by limiter 
        am.adapt();
      }

      // if grid has changed then limit solution again
      if( sequence != solution_.space().sequence() )
      {
        limitSolution();
      }
    }
  }

#endif
  // the solution 
  DiscreteFunctionType   solution_;
	DiscreteGradientType   theta_;
	ScalarDFType   energy_;
  
  DiscreteFunctionType*  additionalVariables_;

  // InitialDataType is a Dune::Operator that evaluates to $u_0$ and also has a
  // method that gives you the exact solution.
  const InitialDataType*  problem_;
  ModelType*              model_;
  // Initial flux for advection discretization (UpwindFlux)
  FluxType                convectionFlux_;
  AdaptationHandlerType*  adaptationHandler_;
  mutable RunFileType     runfile_;
  Timer                   overallTimer_;
  double                  odeSolve_;
  const unsigned int      eocId_;
  OdeSolverType*          odeSolver_;
  OdeSolverMonitorType    odeSolverMonitor_;
  int                     odeSolverType_;

  const bool              adaptive_;
  AdaptationParameters    adaptationParameters_;
};
#endif // FEMHOWTO_STEPPER_HH

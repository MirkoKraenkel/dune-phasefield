#ifndef PHASEFIELD_ALGORITHM_HH
#define PHASEFIELD_ALGORITHM_HH
// system includes
#include <sstream>
// dune-fem includes
#include <dune/fem/io/file/datawriter.hh>
#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/base/base.hh>

//Note Problen shoul be independent of Operator/Scheme







template <class GridImp,
					class AlgorithmTraits,
					class ProblemCreatorImp,
					int polynomilaOrder> 
class PhasefieldAlgorithm
{
	//type of Grid
	typedef GridImp GridType;
	

	
	
	






	//MemberVaribles
	
private:
	GridType&               grid_;
	GridPartType            gridPart_;
	DiscreteSpaceType       space_;
	ThetaSpaceType*         thetaSpace_;
	ScalarDiscreteSpaceType* energySpace_;
	Dune::IOInterface*       eocLoopData_;
	IOTupleType             eocDataTup_; 
	unsigned int            timeStepTimer_; 
	unsigned int            loop_ ; 
	// use fixed time step if fixedTimeStep>0
	double                  fixedTimeStep_;        
	double                  fixedTimeStepEocLoopFactor_;       
	DiscreteFunctionType    solution_;
	DiscreteThetaType*      theta_;
	ScalarDFType*           energy_;
	DiscreteFunctionType*   additionalVariables_;
	const InitialDataType*  problem_;
  ModelType*              model_;
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
	DiscreteOperatorType    dgOperator_;
#if PF_USE_ADAPTATION
  DGIndicatorType         dgIndicator_;
#endif
  double tolerance_;

	//Constructor
	PhasefieldAlgorithm(GridType& grid):
		grid_(grid)	
		gridPart_( grid_ ),
		space_( gridPart_ ),
		thetaSpace_(Parameter :: getValue< bool >("phasefield.theta", false) ? new ThetaSpaceType(gridPart_ ) : 0 ),
		energySpace_(Parameter :: getValue< bool >("phasefield.energy", false) ? new ScalarDiscreteSpaceType(gridPart_ ) : 0 ),
 		eocLoopData_( 0 ),
 		eocDataTup_(),
 		timeStepTimer_( Dune::FemTimer::addTo("max time/timestep") ),
 		loop_( 0 ),
 		fixedTimeStep_( Dune::Parameter::getValue<double>("fixedTimeStep",0) ),
 		fixedTimeStepEocLoopFactor_( Dune::Parameter::getValue<double>("fixedTimeStepEocLoopFactor",1.) ), // algorithmbase
		solution_( "solution", space() ),
		theta_(Parameter :: getValue< bool >("phasefield.theta", false) ? new DiscreteThetaType("theta",thetaSpace_ ) : 0 ),
		energy_(Parameter :: getValue< bool >("phasefield.energy", false) ? new ScalarDFType("energy",energySpace_):0),
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
    adaptationParameters_( ),
		dgOperator_(grid,convectionFlux_)
	{}

//! destructor 
  virtual ~PhasefieldAlgorithm()
  {
		delete thetaSpace_;
		thetaSpace_=0;
		delete energySapce_;
		energySpace_=0;
		delete theta_;
		theta_=0;
		delete energy_;
		energy_=0;
		delete odeSolver_;
		odeSolver_ = 0;
		delete problem_ ;
		problem_ = 0;
		delete adaptationHandler_ ;
		adaptationHandler_ = 0;
		delete additionalVariables_; 
		additionalVariables_ = 0;
  }


	virtual void initializeStep(TimeProviderType& tp)
	{
		DiscreteFunctionType& U = solution();
		//Create OdeSolver if necessary
    if( odeSolver_ == 0 ) odeSolver_ = createOdeSolver( tp );    
		assert( odeSolver_ );
		L2Projection< double, double,InitialDataType, DiscreteFunctionType > l2pro;
    l2pro(problem(), U );          
    odeSolver_->initialize( U );  
	}
	
	
 
	virtual void step(TimeProviderType& tp)
	{
		DiscreteFunctionType& U = solution();
		// reset overall timer
    overallTimer_.reset(); 
		assert(odeSolver_);
    odeSolver_->solve( U, odeSolverMonitor_ );
	}
		
	virtual void operator()()
	{	
		//some setup stuff
		const bool verbose = Dune::Parameter :: verbose ();
		int printCount = Dune::Parameter::getValue<int>("femhowto.printCount", -1);
		
		double maxTimeStep =Dune::Parameter::getValue("femhowto.maxTimeStep", std::numeric_limits<double>::max());
		const double startTime = Dune::Parameter::getValue<double>("femhowto.startTime", 0.0);
		const double endTime   = Dune::Parameter::getValue<double>("femhowto.endTime");	
	
		TimeProviderType tp(startTime, grid_ ); 

		tp.provideTimeStepEstimate(maxTimeStep);
		if ( fixedTimeStep_ > 1e-20 )
				tp.init( fixedTimeStep_ );
		else
			tp.init();
		
		for( ; tp.time() < endTime; )   
			{
			}
		writeData( eocDataOutput, tp, true );
		averagedt /= double(tp.timeStep());
		if(verbose)
			{
				std::cout << "Minimum dt: " << mindt
									<< "\nMaximum dt: " << maxdt
									<< "\nAverage dt: " << averagedt << std::endl;
			}
		finalizeStep( tp );                                  
    
		// increase loop counter
		++loop_;
		
		// prepare the fixed time step for the next eoc loop
		fixedTimeStep_ /= fixedTimeStepEocLoopFactor_; 
	}
	
	virtual void finalizeStep(TimeProviderType& tp)
	{ 
		DiscreteFunctionType& u = solution();
		
		bool doFemEoc = problem().calculateEOC( tp, u, eocId_ );

    // ... and print the statistics out to a file
    if( doFemEoc )
      FemEoc::setErrors(eocId_, error(tp, u ));
		delete odeSolver_;
    odeSolver_ = 0;
    delete adaptationHandler_;
    adaptationHandler_ = 0;
	}






	//! restore all data from check point (overload to do something)
	virtual void restoreFromCheckPoint(TimeProviderType& tp) {} 

	//! write a check point (overload to do something)
	virtual void writeCheckPoint(TimeProviderType& tp,
															 AdaptationManagerType& am ) const {}






	virtual DiscreteFunctionType& solution() { return solution_; }
	virtual ThetaType& theta() { return theta_; }


	virtual OdeSolverType* createOdeSolver(TimeProviderType& tp) 
  {

#if PF_USE_ADAPTATION
    if( adaptive_ )
      {
				if( ! adaptationHandler_ && adaptationParameters_.aposterioriIndicator() )
					{
						adaptationHandler_ = new AdaptationHandlerType( grid_, tp );
						dgIndicator_.setAdaptationHandler( *adaptationHandler_ );
					}
      }
#endif
		
    typedef PhaseFieldOdeSolver< DgType > OdeSolverImpl;
    return new OdeSolverImpl( tp, dgOperator_ );
  }

	virtual void finalize( const int eocloop ) 
		{
			DiscreteFunctionType& U = solution(); 
			DiscreteGradientType& theta1 = theta(); 
			
			//    theta1.print(cout);    
			DiscreteFunctionType* addVars = additionalVariables();

			if( eocLoopData_ == 0 ) 
				{
					eocDataTup_ = IOTupleType( &U,addVars, indicator() ); 
					eocLoopData_ = new DataWriterType( grid_, eocDataTup_ );
				}

			if( addVars ) 
				{
					problem().finalizeSimulation( *addVars, eocloop );
				}
			else 
				{
					problem().finalizeSimulation( U,  eocloop );
				}
			eocLoopData_->writeData( eocloop );
		}




};












#endif

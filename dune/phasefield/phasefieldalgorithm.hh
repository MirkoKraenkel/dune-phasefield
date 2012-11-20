#ifndef PHASEFIELD_ALGORITHM_HH
#define PHASEFIELD_ALGORITHM_HH
// system includes
#include <sstream>
// dune-fem includes
#include <dune/fem/io/file/datawriter.hh>
#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/base/base.hh>


//Note Problen should be independent of Operator/Scheme
#include<dune/phasefield/algorithmtraits.hh>


// include std libs
#include <iostream>
#include <string>
// fem includes

#include <dune/fem/misc/l2error.hh>
#include <dune/fem/operator/projection/l2projection.hh>
#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/util/phasefieldodesolver.hh>
#if WELLBALANCED
#include <dune/fem/operator/wellbalancedoperator.hh>
#else
#include <dune/fem/operator/fluxprojoperator.hh>
#endif
#include <dune/fem-dg/misc/runfile.hh>
#include <dune/fem-dg/operator/adaptation/estimatorbase.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/space/dgspace/localrestrictprolong.hh>
#include <dune/fem/solver/odesolver.hh>
#if WELLBALANCED
#include <dune/phasefield/util/wb_energyconverter.hh>
#else
#include <dune/phasefield/util/energyconverter.hh>
#endif
#include <dune/fem-dg/operator/adaptation/adaptation.hh>
#include <dune/fem/adaptation/estimator1.hh>
#include <dune/phasefield/util/cons2prim.hh>
/////////////////////////////////////////////////////////////////////////////
//
//  EOC output parameter class 
//
/////////////////////////////////////////////////////////////////////////////
namespace Dune{

struct EocDataOutputParameters :   /*@LST1S@*/
  public Dune::Fem::LocalParameter<Dune::Fem::DataWriterParameters,EocDataOutputParameters>
{
  std::string loop_;
  EocDataOutputParameters(int loop, const std::string& name) {
    std::stringstream ss;
    ss << name << loop;
    loop_ = ss.str();
  }
  EocDataOutputParameters(const EocDataOutputParameters& other)
    : loop_(other.loop_) {}

  std::string path() const {
    return loop_;
  }
};

template <class GridImp,
					class AlgorithmTraits,
					int polynomilaOrder> 
class PhasefieldAlgorithm
{
public:
	//type of Grid
	typedef GridImp GridType;
	//traits class gathers types depending on the problem and the operator which is specified there 
	typedef AlgorithmTraits Traits;
	//problem depending types
	typedef typename Traits::ProblemGeneratorType ProblemGeneratorType;
	typedef typename Traits::GridPartType GridPartType;
	typedef typename Traits::DiscreteOperatorType DiscreteOperatorType;
	//discrete spaces
	typedef typename Traits::DiscreteSpaceType       DiscreteSpaceType;
	typedef typename Traits::SigmaDiscreteSpaceType  SigmaDiscreteSpaceType;
	typedef typename Traits::ThetaDiscreteSpaceType  ThetaDiscreteSpaceType;
	typedef typename Traits::ScalarDiscreteSpaceType ScalarDiscreteSpaceType;
	//discrete functions
	typedef typename Traits::DiscreteFunctionType DiscreteFunctionType;
	typedef typename Traits::DiscreteSigmaType DiscreteSigmaType;
	typedef typename Traits::DiscreteThetaType DiscreteThetaType;
	typedef typename Traits::DiscreteScalarType DiscreteScalarType;
	//additional types
	typedef typename Traits :: RestrictionProlongationType RestrictionProlongationType;
	typedef typename Traits :: InitialDataType InitialDataType;
	typedef typename Traits :: ModelType ModelType;
	typedef typename Traits :: FluxType FluxType;


	typedef typename Traits :: OdeSolverType OdeSolverType;
  typedef typename OdeSolverType :: MonitorType  OdeSolverMonitorType ;	

	// type of adaptation manager 
	typedef Dune::Fem::AdaptationManager< GridType, RestrictionProlongationType > AdaptationManagerType;

	// type of IOTuple 
	typedef Dune::tuple< DiscreteFunctionType*,DiscreteFunctionType*,DiscreteScalarType* > IOTupleType; 
	// type of data writer 
	typedef Dune::Fem::DataWriter< GridType, IOTupleType >    DataWriterType;
	
	// type of time provider organizing time for time loops 
	typedef Dune::Fem::GridTimeProvider< GridType >                 TimeProviderType;
	

  typedef AdaptationHandler< GridType,typename DiscreteSpaceType::FunctionSpaceType >  AdaptationHandlerType;

  typedef Dune::RunFile< GridType >  RunFileType;
	
	//MemberVaribles
	
private:
	GridType&               grid_;
	GridPartType            gridPart_;
	DiscreteSpaceType       space_;
	ThetaDiscreteSpaceType*  thetaSpace_;
	SigmaDiscreteSpaceType*  sigmaSpace_;
	ScalarDiscreteSpaceType* energySpace_;
	Dune::Fem::IOInterface*       eocLoopData_;
	IOTupleType             eocDataTup_; 
	unsigned int            timeStepTimer_; 
	unsigned int            loop_ ; 
	// use fixed time step if fixedTimeStep>0
	double                  fixedTimeStep_;        
	double                  fixedTimeStepEocLoopFactor_;       
	DiscreteFunctionType    solution_;
	DiscreteSigmaType*       sigma_;
	DiscreteThetaType*      theta_;
	DiscreteScalarType*      energy_;
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
#if PF_USE_ADAPTATION
  AdaptationParameters    adaptationParameters_;
#endif
	DiscreteOperatorType    dgOperator_;
#if PF_USE_ADAPTATION
  DGIndicatorType         dgIndicator_;
#endif
  double tolerance_;
public:
	//Constructor
	PhasefieldAlgorithm(GridType& grid):
		grid_(grid)	,
		gridPart_( grid_ ),
		space_( gridPart_ ),
		sigmaSpace_( Fem::Parameter :: getValue< bool >("phasefield.energy", false) ? new SigmaDiscreteSpaceType(gridPart_ ) : 0 ),
		thetaSpace_( Fem::Parameter :: getValue< bool >("phasefield.theta", false)  ? new ThetaDiscreteSpaceType(gridPart_ ) : 0 ),
		energySpace_(Fem::Parameter :: getValue< bool >("phasefield.energy", false) ? new ScalarDiscreteSpaceType(gridPart_ ) : 0 ),
 		eocLoopData_( 0 ),
 		eocDataTup_(),
 		timeStepTimer_( Dune::FemTimer::addTo("max time/timestep") ),
 		loop_( 0 ),
 		fixedTimeStep_( Dune::Fem::Parameter::getValue<double>("fixedTimeStep",0) ),
 		fixedTimeStepEocLoopFactor_( Dune::Fem::Parameter::getValue<double>("fixedTimeStepEocLoopFactor",1.) ), // algorithmbase
		solution_( "solution", space() ),
		sigma_(Fem::Parameter :: getValue< bool >("phasefield.energy", false) ? new DiscreteSigmaType("sigma",sigmaspace()) : 0),
		theta_(Fem::Parameter :: getValue< bool >("phasefield.theta", false) ? new DiscreteThetaType("theta",thetaspace() ) : 0 ),
		energy_(Fem::Parameter :: getValue< bool >("phasefield.energy", false) ? new DiscreteScalarType("energy",energyspace()) : 0),
		additionalVariables_( Fem::Parameter :: getValue< bool >("phasefield.additionalvariables", false) ? 
													new DiscreteFunctionType("additional", space() ) : 0 ),
		problem_( ProblemGeneratorType ::problem() ),
    model_( new ModelType( problem() ) ),
    convectionFlux_( *model_ ),
    adaptationHandler_( 0 ),
    runfile_( grid.comm(), true ),
    overallTimer_(),
    eocId_( Fem::FemEoc::addEntry(std::string("$L^2$-error")) ),
    odeSolver_( 0 ),
    adaptive_( Dune::Fem::AdaptationMethod< GridType >( grid_ ).adaptive() ),
		//     adaptationParameters_( ),
		dgOperator_(grid,convectionFlux_)
	{}

	//! destructor 
  virtual ~PhasefieldAlgorithm()
  {
		delete sigma_;
		sigma_=0;
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

	//some acces methods
	//spacs
	DiscreteSpaceType& space()
	{
		return space_;
	}
	SigmaDiscreteSpaceType& sigmaspace()
	{
		assert(sigmaSpace_);
		return *sigmaSpace_;
	}

	
	ThetaDiscreteSpaceType& thetaspace()
	{
		assert(thetaSpace_);
		return *thetaSpace_;
	}

	ScalarDiscreteSpaceType& energyspace()
	{
		assert(energySpace_);
		return *energySpace_;
	}
	size_t gridSize() const
	{
		size_t grSize=grid_.size(0);
		return grid_.comm().sum(grSize);
	}

	DiscreteSigmaType*  sigma() {return sigma_;}
	DiscreteThetaType*  theta() { return theta_; }
	DiscreteScalarType* energy(){ return energy_;}

  // return pointer to additional variables, can be zero 
  virtual DiscreteFunctionType* additionalVariables() { return additionalVariables_; }	
 
	std::string dataPrefix()
	{
		assert( problem_ );  
		return problem_->dataPrefix(); 
	}
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
#if 1
  void step(TimeProviderType& tp,
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
		newton_iterations     = odeSolverMonitor_.newtonIterations_;
    ils_iterations        = odeSolverMonitor_.linearSolverIterations_;
    max_newton_iterations = odeSolverMonitor_.maxNewtonIterations_ ;
    max_ils_iterations    = odeSolverMonitor_.maxLinearSolverIterations_;

	}
#endif
	void step(TimeProviderType& tp)
	{
		DiscreteFunctionType& U = solution();
		// reset overall timer
    overallTimer_.reset(); 
		assert(odeSolver_);
    odeSolver_->solve( U, odeSolverMonitor_ );
	}

	//! estimate and mark solution 
  virtual void estimateMarkAdapt( AdaptationManagerType& am )
  {
#if 0
    EstimatorType estimator( solution_,grid_ );
    estimator.estimateAndMark(tolerance_);
    am.adapt();
#endif

  }
	//! write data, if pointer to additionalVariables is true, they are calculated first 
	void writeData( DataWriterType& eocDataOutput,
									TimeProviderType& tp,
									const bool reallyWrite )
	{


		if( reallyWrite )
			{
				 
				DiscreteFunctionType* addVariables = additionalVariables();
				// calculate DG-projection of additional variables
				if ( addVariables )
					{
						// calculate additional variables from the current num. solution
					  setupAdditionalVariables( solution(), model(), *addVariables );
					}

				//	energyconverter(solution(),sigma(),model(),energy());
			}

		// write the data 
		eocDataOutput.write( tp );
	}
	
	
	
	virtual void operator()()
	{	

		int counter;
		//some setup stuff
		const bool verbose = Dune::Fem::Parameter :: verbose ();
		int printCount = Dune::Fem::Parameter::getValue<int>("femhowto.printCount", -1);
		// if adaptCount is 0 then no dynamics grid adaptation
		int adaptCount = 0;
		int maxAdaptationLevel = 0;
		Dune::Fem::AdaptationMethod< GridType > am( grid_ );
		if( am.adaptive() )
			{
				adaptCount = Dune::Fem::Parameter::getValue<int>("fem.adaptation.adaptcount");
				maxAdaptationLevel = Dune::Fem::Parameter::getValue<int>("fem.adaptation.finestLevel");
			}

		double maxTimeStep =Dune::Fem::Parameter::getValue("femhowto.maxTimeStep", std::numeric_limits<double>::max());
		const double startTime = Dune::Fem::Parameter::getValue<double>("femhowto.startTime", 0.0);
		const double endTime   = Dune::Fem::Parameter::getValue<double>("femhowto.endTime");	
		const int maximalTimeSteps =
			Dune::Fem::Parameter::getValue("femhowto.maximaltimesteps", std::numeric_limits<int>::max());

		TimeProviderType tp(startTime, grid_ ); 

		DiscreteFunctionType& U = solution();
    RestrictionProlongationType rp( U );
		
    
		rp.setFatherChildWeight( Dune::DGFGridInfo<GridType> :: refineWeight() );
 
		// create adaptation manager 
		AdaptationManagerType adaptManager( grid_ , rp );

		// restoreData if checkpointing is enabled (default is disabled)
		//restoreFromCheckPoint( tp );
		
		// tuple with additionalVariables 
		IOTupleType dataTuple( &U,   this->additionalVariables() );
		
		// type of the data writer
		DataWriterType eocDataOutput( grid_, dataTuple, tp, EocDataOutputParameters( loop_ , dataPrefix() ) );
		
		// set initial data (and create ode solver)
		initializeStep( tp );
		
		// start first time step with prescribed fixed time step 
		// if it is not 0 otherwise use the internal estimate
		
		tp.provideTimeStepEstimate(maxTimeStep);
		if ( fixedTimeStep_ > 1e-20 )
			tp.init( fixedTimeStep_ );
		else
			tp.init();
		


		// adapt the grid to the initial data
		int startCount = 0;
		if( adaptCount > 0 )
			while( startCount < maxAdaptationLevel )
				{
					estimateMarkAdapt( adaptManager );
					
					initializeStep( tp );

					if( verbose )
						std::cout << "start: " << startCount << " grid size: " << grid_.size(0)
											<< std::endl;
					++startCount;
				}
		
		

		writeData( eocDataOutput, tp, eocDataOutput.willWrite( tp ) );
		


		
		for( ; tp.time() < endTime; )   
			{ 
				tp.provideTimeStepEstimate(maxTimeStep);                                         
			
				const double tnow  = tp.time();
				const double ldt   = tp.deltaT();
				
				int newton_iterations;
				int ils_iterations;
				counter  = tp.timeStep();
					
				Dune::FemTimer::start(timeStepTimer_);
				// grid adaptation (including marking of elements)
				if( (adaptCount > 0) && (counter % adaptCount) == 0 )
					estimateMarkAdapt( adaptManager );
		
				//this is where the magic happens
				step( tp);
				
				// Check that no NAN have been generated
				if (! U.dofsValid()) 
					{
						std::cout << "Loop(" << loop_ << "): Invalid DOFs" << std::endl;
						eocDataOutput.write(tp);
						abort();
					}
					
				if( (printCount > 0) && (counter % printCount == 0))
					{
						size_t grSize = gridSize();
						if( grid_.comm().rank() == 0 )
							std::cout << "step: " << counter << "  time = " << tnow << ", dt = " << ldt
												<<" grid size: " << grSize << std::endl;
					}


				writeData( eocDataOutput, tp, eocDataOutput.willWrite( tp ) );
				writeCheckPoint( tp, adaptManager );
							 
				// next time step is prescribed by fixedTimeStep
				// it fixedTimeStep is not 0
				if ( fixedTimeStep_ > 1e-20 )
					tp.next( fixedTimeStep_ );
				else
					tp.next();

			}
		/*end of timeloop*/
		// write last time step  
		writeData( eocDataOutput, tp, true );


	
	
		// 		writeData( eocDataOutput, tp, true );
	
		finalizeStep( tp );                                  
    
	
		++loop_;
	
		// prepare the fixed time step for the next eoc loop
		fixedTimeStep_ /= fixedTimeStepEocLoopFactor_; 
	}
	inline double error(TimeProviderType& tp, DiscreteFunctionType& u)
	{
		typedef typename DiscreteFunctionType :: RangeType RangeType;
		Fem::L2Error<DiscreteFunctionType> L2err;
		// Compute L2 error of discretized solution ...
		RangeType error = L2err.norm(problem(), u, tp.time());
		return error.two_norm();
	}
	
	virtual void finalizeStep(TimeProviderType& tp)
	{ 
		DiscreteFunctionType& u = solution();
		
		bool doFemEoc = problem().calculateEOC( tp, u, eocId_ );

		// ... and print the statistics out to a file
		if( doFemEoc )
			Fem::FemEoc::setErrors(eocId_, error(tp, u));
	
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
	

	virtual OdeSolverType* createOdeSolver(TimeProviderType& tp) 
	{

#if 0
		if( adaptive_ )
			{
				if( ! adaptationHandler_ && adaptationParameters_.aposterioriIndicator() )
					{
						adaptationHandler_ = new AdaptationHandlerType( grid_, tp );
						dgIndicator_.setAdaptationHandler( *adaptationHandler_ );
					}
			}
#endif
		
		typedef PhaseFieldOdeSolver< DiscreteOperatorType > OdeSolverImpl;
		return new OdeSolverImpl( tp, dgOperator_ );
	}

	virtual void finalize( const int eocloop ) 
	{
		DiscreteFunctionType& U = solution(); 
		DiscreteThetaType* theta1 = theta(); 
			
  
		DiscreteFunctionType* addVars = additionalVariables();

		if( eocLoopData_ == 0 ) 
			{
				eocDataTup_ = IOTupleType( &U,addVars ); 
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



}//end namespace Dune








#endif
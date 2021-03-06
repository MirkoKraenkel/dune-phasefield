#ifndef PHASEFIELD_ALGORITHM_HH
#define PHASEFIELD_ALGORITHM_HH
// system includes
#include <sstream>
#include <dune/common/nullptr.hh>

// dune-fem includes
#include <dune/fem/io/file/datawriter.hh>
#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/base/base.hh>

//Note Problen should be independent of Operator/Scheme
#include<dune/phasefield/allencahntraits.hh>


// include std libs
#include <iostream>
#include <string>
// fem includes

#include <dune/fem/operator/projection/l2projection.hh>
#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/util/phasefieldodesolver.hh>

#if WELLBALANCED
//#include <dune/fem/operator/wellbalancedoperator.hh>
//#include <dune/fem/operator/wbsplitoperator.hh>
#else
//#include <dune/fem/operator/fluxprojoperator.hh>
#endif
//#include <dune/fem-dg/misc/runfile.hh>
#include <dune/fem-dg/operator/adaptation/estimatorbase.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/space/discontinuousgalerkin/localrestrictprolong.hh>

#if WELLBALANCED
#include <dune/phasefield/util/wb_energyconverter.hh>
#else
#include <dune/phasefield/util/energyconverter.hh>
#endif
#include <dune/fem-dg/operator/adaptation/adaptation.hh>


#include <dune/fem/adaptation/estimator2.hh>
#include <dune/phasefield/util/cons2prim.hh>
namespace Dune{

/////////////////////////////////////////////////////////////////////////////
//
//  EOC output parameter class 
//
/////////////////////////////////////////////////////////////////////////////
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


/////////////////////////////////////////////////////////////////////////////
//
//  EOC output parameter class 
//
/////////////////////////////////////////////////////////////////////////////


template <class GridImp,
					class AlgorithmTraits,
					int polynomilaOrder> 
class AllenCahnAlgorithm
{
public:
  //type of Grid
	typedef GridImp GridType;

  //traits class gathers types depending on the problem and the operator which is specified there 
	typedef AlgorithmTraits Traits;
	//problem depending types
	
  typedef typename Traits::ProblemGeneratorType ProblemGeneratorType;
  typedef typename ProblemGeneratorType::VelocityType VelocityType;
  typedef typename Traits::GridPartType GridPartType;
	typedef typename Traits::DiscreteOperatorType DiscreteOperatorType;
 
  //for interpolation of initial Data 
	typedef typename Traits::LagrangeGridPartType LagrangeGridPartType;


  //discrete spaces
	typedef typename Traits::DiscreteSpaceType       DiscreteSpaceType;
  typedef typename Traits::DiscreteVelocitySpaceType DiscreteVelocitySpaceType;	
  //discrete functions
	typedef typename Traits::DiscreteFunctionType DiscreteFunctionType;
  typedef typename Traits::DiscreteVelocityType DiscreteVelocityType;  
  
  //additional types
	typedef typename Traits :: RestrictionProlongationType RestrictionProlongationType;
	typedef typename Traits :: InitialDataType             InitialDataType;
	typedef typename Traits :: ModelType                   ModelType;
	typedef typename Traits :: FluxType                    FluxType;

  typedef typename DiscreteSpaceType::FunctionSpaceType FunctionSpaceType;

	typedef typename Traits :: OdeSolverType       OdeSolverType;
  typedef typename OdeSolverType :: MonitorType  OdeSolverMonitorType ;	

	// type of adaptation manager 
	typedef Dune::Fem::AdaptationManager< GridType, RestrictionProlongationType > AdaptationManagerType;

	// type of IOTuple 
	//typedef Dune::tuple< DiscreteFunctionType*, DiscreteSigmaType*,DiscreteThetaType* > IOTupleType; 
  typedef Dune::tuple< DiscreteFunctionType* > IOTupleType; 
	
  // type of data 
  // writer 
  typedef Dune::Fem::DataWriter< GridType, IOTupleType >    DataWriterType;

	// type of ime provider organizing time for time loops 
	typedef Dune::Fem::GridTimeProvider< GridType >                 TimeProviderType;
	

  typedef AdaptationHandler< GridType,typename DiscreteSpaceType::FunctionSpaceType >  AdaptationHandlerType;
  typedef   Estimator1<DiscreteFunctionType> EstimatorType;
	typedef typename Dune::Fem::LagrangeDiscreteFunctionSpace< FunctionSpaceType, LagrangeGridPartType, 2> InterpolationSpaceType;
  typedef typename Dune::Fem::AdaptiveDiscreteFunction< InterpolationSpaceType > InterpolationFunctionType;


	//MemberVaribles
	
private:
	GridType&               grid_;
	GridPartType            gridPart_;
	DiscreteSpaceType       space_;
	DiscreteVelocitySpaceType    veloSpace_;
  Dune::Fem::IOInterface* eocLoopData_;
	IOTupleType             eocDataTup_; 
	unsigned int            timeStepTimer_; 
	unsigned int            loop_ ; 
	double                  fixedTimeStep_;        
	double                  fixedTimeStepEocLoopFactor_;       
  std::string             energyFilename_;
  DiscreteFunctionType    solution_;
	DiscreteVelocityType    velocity_;
  DiscreteFunctionType*   oldsolution_; 
	const InitialDataType*  problem_;
  ModelType*              model_;
  FluxType                convectionFlux_;
  AdaptationHandlerType*  adaptationHandler_;
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
  bool   interpolateData_;
  bool   calcresidual_;
  double timeStepTolerance_;
public:
	//Constructor
	AllenCahnAlgorithm(GridType& grid):
		grid_(grid)	,
		gridPart_( grid_ ),
		space_( gridPart_ ),
    veloSpace_(gridPart_),
 		eocLoopData_( 0 ),
 		eocDataTup_(),
 		timeStepTimer_( Dune::FemTimer::addTo("max time/timestep") ),
 		loop_( 0 ),
 		fixedTimeStep_( Dune::Fem::Parameter::getValue<double>("fixedTimeStep",0) ),
 		fixedTimeStepEocLoopFactor_( Dune::Fem::Parameter::getValue<double>("fixedTimeStepEocLoopFactor",1.) ), // algorithmbase
    energyFilename_(Dune::Fem::Parameter::getValue< std::string >("phasefield.energyfile","./energy.gnu")),
    solution_( "solution", space() ),
		velocity_("velocity" , velospace()),
    oldsolution_( Fem::Parameter :: getValue< bool >("phasefield.storelaststep", false) ? 
													new DiscreteFunctionType("oldsolution", space() ) : nullptr ),
		problem_( ProblemGeneratorType::problem() ),
    model_( new ModelType( problem() ) ),
    convectionFlux_( *model_ ),
    adaptationHandler_( 0 ),
    overallTimer_(),
    eocId_( Fem::FemEoc::addEntry(std::string("L2error")) ),
    odeSolver_( 0 ),
    adaptive_( Dune::Fem::AdaptationMethod< GridType >( grid_ ).adaptive() ),
		//     adaptationParameters_( ),
		dgOperator_(grid,convectionFlux_),
    tolerance_(Fem::Parameter :: getValue< double >( "phasefield.adaptTol", 100) ),
    interpolateData_(Fem::Parameter::getValue<bool>( "phasefield.interpolinitial", false) ),
    calcresidual_( Fem :: Parameter :: getValue< bool >("phasefield.calcresidual"  ) ),
    timeStepTolerance_( Fem :: Parameter :: getValue< double >( "phasefield.timesteptolerance",-1. ) )
    {
    }

  //! destructor 
  virtual ~AllenCahnAlgorithm()
  {
		delete oldsolution_;
    oldsolution_= nullptr;
		delete odeSolver_;
		odeSolver_ = 0;
		delete problem_ ;
		problem_ = 0;
		delete model_;
    model_=0;
    delete adaptationHandler_ ;
		adaptationHandler_ = 0;
  }

	//some acces methods
	//spacs
	DiscreteSpaceType& space()
	{
		return space_;
	}

  DiscreteVelocitySpaceType& velospace()
  {
    return veloSpace_;
  } 
 
  size_t gridSize() const
	{
		size_t grSize=grid_.size(0);
		return grid_.comm().sum(grSize);
	}

  DiscreteFunctionType* oldsolution() {return oldsolution_;}

 
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
   
   
    if( interpolateData_ )
      {
        LagrangeGridPartType lagGridPart(grid_); 
        InterpolationSpaceType interpolSpace(lagGridPart); 
        InterpolationFunctionType interpolSol("uinterpol",interpolSpace); 
  
        Dune::Fem::LagrangeInterpolation<InitialDataType,InterpolationFunctionType>::interpolateFunction( problem(),interpolSol);
        Dune::Fem::DGL2ProjectionImpl::project(interpolSol, U);
      }
    else 
    {   
       Dune::Fem::DGL2ProjectionImpl::project(problem().fixedTimeFunction(tp.time()), U);
     }
       odeSolver_->initialize( U );  

 
  }

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
	
  void step(TimeProviderType& p)
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
    EstimatorType estimator( solution_,grid_ );
    estimator.estimateAndMark(tolerance_);
    am.adapt();
  }
  template<class Stream>
  void writeEnergy( TimeProviderType& tp,
                    Stream& str)
  {
#if 0
    DiscreteSigmaType* gradient = sigma();
    DiscreteScalarType* totalenergy = energy();

    if(gradient != nullptr && gradient != nullptr)
    { 
      gradient->clear();
     
      dgOperator_.gradient(solution(),*gradient);
   
      double kineticEnergy;
      
      double chemicalEnergy; 
      
      double energyIntegral =energyconverter(solution(),*gradient,model(),*totalenergy,kineticEnergy,chemicalEnergy);
      str<<std::setprecision(10)<<tp.time()<<"\t"<<energyIntegral<<"\t"<<chemicalEnergy<<"\t"<<kineticEnergy<<"\n";
    
    }
#endif
  }
  
  //! write data, if pointer to additionalVariables is true, they are calculated first 
  void writeData( DataWriterType& eocDataOutput,
									TimeProviderType& tp,
                  const bool reallyWrite )
	{
    if( reallyWrite )
		{
				 
         // calculate DG-projection of additional variables
		
            
    }

		// write the data 
		eocDataOutput.write( tp );
	}
	
	
	
	virtual void operator()(int loopNumber,double& averagedt, double& mindt, double& maxdt,
                          size_t& counter, int& total_newton_iterations, int& total_ils_iterations,
                          int& max_newton_iterations, int& max_ils_iterations)
	{	

		double timeStepError=std::numeric_limits<double>::max();
 
     //some setup stuff
		const bool verbose = Dune::Fem::Parameter :: verbose ();
		int printCount = Dune::Fem::Parameter::getValue<int>("phasefield.printCount", -1);

    // if adaptCount is 0 then no dynamics grid adaptation
		int adaptCount = 0;
		int maxAdaptationLevel = 0;
		
    Dune::Fem::AdaptationMethod< GridType > am( grid_ );
		if( am.adaptive() )
			{
				adaptCount = Dune::Fem::Parameter::getValue<int>("fem.adaptation.adaptcount");
				maxAdaptationLevel = Dune::Fem::Parameter::getValue<int>("fem.adaptation.finestLevel");
			}

		double maxTimeStep =Dune::Fem::Parameter::getValue("phasefield.maxTimeStep", std::numeric_limits<double>::max());
		const double startTime = Dune::Fem::Parameter::getValue<double>("phasefield.startTime", 0.0);
		const double endTime   = Dune::Fem::Parameter::getValue<double>("phasefield.endTime",1.);	
    const int maximalTimeSteps =Dune::Fem::Parameter::getValue("phasefield.maximaltimesteps", std::numeric_limits<int>::max());

		//statistics
    maxdt = 0.;
    mindt = std::numeric_limits<double>::max();
    averagedt=0.;
    total_newton_iterations = 0;
    total_ils_iterations = 0;
    max_newton_iterations = 0;
    max_ils_iterations = 0;max_ils_iterations = 0;  
    
    //Initialize Tp
    TimeProviderType tp(startTime, grid_ ); 

		DiscreteFunctionType& U = solution();
    DiscreteFunctionType* Uold = oldsolution(); 

    VelocityType velo;
    Dune::Fem::DGL2ProjectionImpl::project(velo, velocity_);

 //   typedef Fem::GridFunctionAdapter<Velocity,GridPartType> GridVeloType;
   // GridVeloType  gridvelo("grid velocity", velo,gridPart_,solution().space().order()+1);

    dgOperator_.setVelocity(velocity_);
   
    RestrictionProlongationType rp( U );
    
		rp.setFatherChildWeight( Dune::DGFGridInfo<GridType> :: refineWeight() );
 
		// create adaptation manager 
    AdaptationManagerType adaptManager(grid_,rp);
		
    // restoreData if checkpointing is enabled (default is disabled)
		//restoreFromCheckPoint( tp );
		
		// tuple with additionalVariables 
	  IOTupleType dataTuple( &U );
	
   // IOTupleType dataTuple( &U, this->sigma(),this->theta() );
    std::ofstream energyfile;
    std::ostringstream convert;
    convert<<loopNumber;
    std::string filename=energyFilename_;
    filename.append(convert.str()); 
    filename.append(".gnu");
    energyfile.open(filename.c_str());

	
		// type of the data writer
		DataWriterType eocDataOutput( grid_, dataTuple, tp, EocDataOutputParameters( loop_ , dataPrefix() ) );
		
		// set initial data (and create ode solver)
		initializeStep( tp );
	   if(Uold!=nullptr)
      { 
        Uold->assign(U);
      }
 
 		// start first time step with prescribed fixed time step 
		// if it is not 0 otherwise use the internal estimate
		
	
		// adapt the grid to the initial data
		int startCount = 0;
		if( adaptCount > 0 )
    {
      while( startCount < maxAdaptationLevel )
				{
					estimateMarkAdapt( adaptManager );
					
					initializeStep( tp );

					if( verbose )
						std::cout << "start: " << startCount << " grid size: " << grid_.size(0)<<std::endl;
          ++startCount;
				
        }
    }
    tp.provideTimeStepEstimate(maxTimeStep);
		if ( fixedTimeStep_ > 1e-20 )
			tp.init( fixedTimeStep_ );
		else
			tp.init();
		


    tp.provideTimeStepEstimate(maxTimeStep);                                         
		
    std::cout<<"deltaT "<<tp.deltaT()<<" estimate "<<dgOperator_.timeStepEstimate()<<std::endl;
					


		writeData( eocDataOutput, tp, eocDataOutput.willWrite( tp ) );

    if( calcresidual_)
        {
          std::cout<<"Residual\n";
          Uold->clear();
      dgOperator_(U,*Uold);
      U.assign(*Uold);
 			writeData( eocDataOutput , tp , eocDataOutput.willWrite( tp ) );
	
    }
    else
    {
    for( ; tp.time() < endTime && tp.timeStep() < maximalTimeSteps /* && timeStepError > timeStepTolerance_*/;  )   
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
		
			step( tp, newton_iterations, ils_iterations, 
            max_newton_iterations, max_ils_iterations);
			// Check that no NAN have been generated
			if (! U.dofsValid()) 
			{
  			std::cout << "Loop(" << loop_ << "): Invalid DOFs" << std::endl;
	   		eocDataOutput.write(tp);
        energyfile.close();
        abort();
			}
		
      if(Uold!=nullptr)
      { 
        timeStepError = stepError(U,*Uold);
        timeStepError/=ldt;
        Uold->assign(U);
      }
      double timeStepEstimate=dgOperator_.timeStepEstimate();	
//      double diffTimeStep=dgOperator_.maxDiffusionTimeStep();
  //    double advTimeStep=dgOperator_.maxAdvectionTimeStep();
     if( (printCount > 0) && (counter % printCount == 0))
			{
	
        if( grid_.comm().rank() == 0 )
        {
          std::cout <<"step: " << counter << "  time = " << tnow << ", dt = " << ldt<<" ,timeStepEstimate " <<timeStepEstimate;
          if(Uold!=nullptr)
           std::cout<< " ,Error between timesteps="<< timeStepError;
           std::cout<<std::endl;
               
        }
         writeEnergy( tp , energyfile);
        
     }


			writeData( eocDataOutput , tp , eocDataOutput.willWrite( tp ) );
			writeCheckPoint( tp, adaptManager );
      //statistics
      mindt = (ldt<mindt) ? ldt : mindt;
      maxdt = (ldt>maxdt) ? ldt : maxdt;
      averagedt += ldt;
      total_newton_iterations+=newton_iterations;
      total_ils_iterations+=ils_iterations;
  
			// next time step is prescribed by fixedTimeStep
			// it fixedTimeStep is not 0
			if ( fixedTimeStep_ > 1e-20 )
				tp.next( fixedTimeStep_ );
			else
				tp.next();

		}		/*end of timeloop*/
		
    }    // write last time step  
		writeData( eocDataOutput, tp, true );

    energyfile.close();
	  averagedt /= double(tp.timeStep());
	
	// 		writeData( eocDataOutput, tp, true );
	
		finalizeStep( tp );                                  
    
	
		++loop_;
	
		// prepare the fixed time step for the next eoc loop
		fixedTimeStep_ /= fixedTimeStepEocLoopFactor_; 
	}
	
  inline double error(TimeProviderType& tp, DiscreteFunctionType& u)
	{
		typedef typename DiscreteFunctionType :: RangeType RangeType;
		Fem::L2Norm< GridPartType > l2norm(gridPart_);
    
    //typedef Fem::GridFunctionAdapter<InitialDataType,GridPartType> GridFunctionType;
    //GridFunctionType exactsolution("exact solution",problem(),gridPart_,space().order()+1);   
		// Compute L2 error of discretized solution ...
		//RangeType error = L2err.norm(problem(), u, tp.time());
    //double error = l2norm.distance(exactsolution,u);  
   // double error = l2norm.distance(problem(),u);  
   double error = l2norm.distance(problem().fixedTimeFunction(tp.time()),u);

    
    return error;
	}
  //compute Error between old and NewTimeStep
  inline double stepError(DiscreteFunctionType uOld, DiscreteFunctionType& uNew)
	{
		typedef typename DiscreteFunctionType :: RangeType RangeType;
		Fem::L2Norm<GridPartType> l2norm(gridPart_);
	  return l2norm.distance(uOld, uNew);

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
	
    if( eocLoopData_ == 0 ) 
			{
			
        eocDataTup_ = IOTupleType( &U ); 
        //eocDataTup_ = IOTupleType( &U,sig ); 
        eocLoopData_ = new DataWriterType( grid_, eocDataTup_ );
			
      }

				problem().finalizeSimulation( U,  eocloop );
		eocLoopData_->writeData( eocloop );
	}
};



}//end namespace Dune








#endif

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
#include<dune/phasefield/algorithmtraits.hh>


// include std libs
#include <iostream>
#include <string>
// fem includes

#include <dune/fem/operator/projection/l2projection.hh>
#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
//solver
#include <dune/fem/util/phasefieldodesolver.hh>



//adaptation
#include <dune/fem-dg/operator/adaptation/estimatorbase.hh>
#include <dune/fem/space/discontinuousgalerkin/localrestrictprolong.hh>
#include <dune/fem-dg/operator/adaptation/adaptation.hh>
#include <dune/fem/adaptation/estimator2.hh>
//post processing
#if WELLBALANCED
#include <dune/phasefield/util/wb_energyconverter.hh>
#else
#warning "CONSERVATIVE!"
#include <dune/phasefield/util/energyconverter.hh>
#endif
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
//  Algorithm class 
//
/////////////////////////////////////////////////////////////////////////////


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
  //for interpolation of initial Data 
	typedef typename Traits::LagrangeGridPartType LagrangeGridPartType;
	
  //Problemdependent Types
	typedef typename Traits :: InitialDataType             InitialDataType;
	typedef typename Traits :: ModelType                   ModelType;
	typedef typename Traits :: FluxType                    FluxType;
  
  //discrete spaces
  typedef typename Traits::DiscreteSpaceType       DiscreteSpaceType;
	typedef typename Traits::SigmaDiscreteSpaceType  SigmaDiscreteSpaceType;
	typedef typename Traits::ThetaDiscreteSpaceType  ThetaDiscreteSpaceType;
	typedef typename Traits::ScalarDiscreteSpaceType ScalarDiscreteSpaceType;
  
  //discrete functions
	typedef typename Traits::DiscreteFunctionType DiscreteFunctionType;
  typedef typename Traits::DiscreteSigmaType    DiscreteSigmaType;
	typedef typename Traits::DiscreteThetaType    DiscreteThetaType;
	typedef typename Traits::DiscreteScalarType   DiscreteScalarType;
 
  
  
  typedef typename DiscreteSpaceType::FunctionSpaceType FunctionSpaceType;

  //Operator
  typedef typename Traits :: DiscreteOperatorType DiscreteOperatorType;
	
  
 
  //Solver
  typedef typename Traits :: OdeSolverType       OdeSolverType;
  typedef typename OdeSolverType :: MonitorType  OdeSolverMonitorType ;	
  //Adaptation
  typedef typename Traits :: RestrictionProlongationType RestrictionProlongationType;
 	// type of adaptation manager 
	typedef Dune::Fem::AdaptationManager< GridType, RestrictionProlongationType > AdaptationManagerType;
  typedef AdaptationHandler< GridType,typename DiscreteSpaceType::FunctionSpaceType >  AdaptationHandlerType;
  typedef  Estimator1<DiscreteFunctionType> EstimatorType;
	
  
  typedef typename Dune::Fem::LagrangeDiscreteFunctionSpace< FunctionSpaceType, LagrangeGridPartType, 2> InterpolationSpaceType;
  typedef typename Dune::Fem::AdaptiveDiscreteFunction< InterpolationSpaceType > InterpolationFunctionType;


	//type of IOTuple 
	//typedef Dune::tuple< DiscreteFunctionType*, DiscreteSigmaType*,DiscreteThetaType* > IOTupleType; 
  //Pointers for (rho,rho v,rho phi),(v,p,phi),(totalenergy),(\theta)
  typedef Dune::tuple< DiscreteFunctionType*,DiscreteFunctionType*,DiscreteScalarType*,DiscreteThetaType* > IOTupleType; 

  // type of data 
  // writer 
  typedef Dune::Fem::DataOutput< GridType, IOTupleType >    DataWriterType;
  typedef Dune::Fem::CheckPointer< GridType >   CheckPointerType;

	// type of ime provider organizing time for time loops 
	typedef Dune::Fem::GridTimeProvider< GridType >                 TimeProviderType;

	//MemberVaribles
private:
	GridType&               grid_;
	GridPartType            gridPart_;
	DiscreteSpaceType       space_;
	ScalarDiscreteSpaceType energySpace_;
	Dune::Fem::IOInterface* eocLoopData_;
	IOTupleType             eocDataTup_; 
	unsigned int            timeStepTimer_; 
	unsigned int            loop_ ; 
	double                  fixedTimeStep_;        
	double                  fixedTimeStepEocLoopFactor_;       
  std::string             energyFilename_;
  DiscreteFunctionType    solution_;
	DiscreteFunctionType*   oldsolution_; 
  DiscreteScalarType*     energy_;
  const InitialDataType*  problem_;
  ModelType*              model_;
  AdaptationHandlerType*  adaptationHandler_;
  Timer                   overallTimer_;
  const unsigned int      eocId_;
  double tolerance_;
  bool interpolateInitialData_;
  bool calcresidual_;
  double timeStepTolerance_;
  /////////////////////////////////////////////// 
  FluxType                convectionFlux_;
  DiscreteOperatorType    dgOperator_;
  SigmaDiscreteSpaceType  sigmaSpace_;
  ThetaDiscreteSpaceType  thetaSpace_;
  DiscreteSigmaType*      sigma_;
	DiscreteThetaType*      theta_;
  DiscreteFunctionType*   additionalVariables_;
  OdeSolverType*          odeSolver_;
  OdeSolverMonitorType    odeSolverMonitor_;
  int                     odeSolverType_;
public:
	//Constructor
	PhasefieldAlgorithm(GridType& grid):
		grid_(grid)	,
		gridPart_( grid_ ),
		space_( gridPart_ ),
		energySpace_(gridPart_ ),
 		eocLoopData_( 0 ),
 		eocDataTup_(),
 		timeStepTimer_( Dune::FemTimer::addTo("max time/timestep") ),
 		loop_( 0 ),
 		fixedTimeStep_( Dune::Fem::Parameter::getValue<double>("fixedTimeStep",0) ),
 		fixedTimeStepEocLoopFactor_( Dune::Fem::Parameter::getValue<double>("fixedTimeStepEocLoopFactor",1.) ), // algorithmbase
    energyFilename_(Dune::Fem::Parameter::getValue< std::string >("phasefield.energyfile","./energy.gnu")),
    solution_( "solution", space() ),
		oldsolution_( Fem::Parameter :: getValue< bool >("phasefield.storelaststep", false) ?	new DiscreteFunctionType("oldsolution", space() ) : nullptr ),
  	energy_( Fem :: Parameter :: getValue< bool >("phasefield.energy", false) ? new DiscreteScalarType("energy",energyspace()) : nullptr),
		problem_( ProblemGeneratorType::problem() ),
    model_( new ModelType( problem() ) ),
    adaptationHandler_( 0 ),
    overallTimer_(),
    eocId_( Fem::FemEoc::addEntry(std::string("L2error")) ),
    tolerance_(Fem::Parameter :: getValue< double >("phasefield.adaptTol", 100)),
    interpolateInitialData_( Fem :: Parameter :: getValue< bool >("phasefield.interpolinitial" , false ) ),
    calcresidual_( Fem :: Parameter :: getValue< bool >("phasefield.calcresidual" , false ) ),
    timeStepTolerance_( Fem :: Parameter :: getValue< double >( "phasefield.timesteptolerance",-1. ) ),
    convectionFlux_( *model_ ),
 	  dgOperator_(grid,convectionFlux_),
   	thetaSpace_( gridPart_ ),
    sigmaSpace_(gridPart_ ),
	  sigma_( Fem :: Parameter :: getValue< bool >("phasefield.sigma", false) ? new DiscreteSigmaType("sigma",sigmaspace()) : nullptr),
		theta_( Fem :: Parameter :: getValue< bool >("phasefield.theta", false) ? new DiscreteThetaType("theta",thetaspace() ) : nullptr),
		additionalVariables_( Fem::Parameter :: getValue< bool >("phasefield.additionalvariables", false) ? 
													new DiscreteFunctionType("additional", space() ) : 0 ),
		odeSolver_( 0 )
    {
    }

  //! destructor 
  virtual ~PhasefieldAlgorithm()
  {
		delete energy_;
		energy_=0;
		delete problem_ ;
		problem_ = 0;
		delete model_;
    model_=0;
    delete adaptationHandler_ ;
		adaptationHandler_ = 0;
		delete oldsolution_;
    oldsolution_= nullptr;
    delete sigma_;
		sigma_=0;
    delete theta_;
		theta_=0;
  	delete additionalVariables_; 
		additionalVariables_ = 0;
  	delete odeSolver_;
		odeSolver_ = 0;
  }

	//some acces methods
	//space
	DiscreteSpaceType& space () { return space_; }	
	

	ScalarDiscreteSpaceType& energyspace () { return energySpace_; }
	
  size_t gridSize() const
	{
		size_t grSize=grid_.size(0);
		return grid_.comm().sum(grSize);
	}


  SigmaDiscreteSpaceType& sigmaspace()
	{
		return sigmaSpace_;
	}
	
	ThetaDiscreteSpaceType& thetaspace()
	{
		return thetaSpace_;
	}


  DiscreteFunctionType* oldsolution() { return oldsolution_; }
  DiscreteScalarType* energy()        { return energy_; }
	DiscreteSigmaType*  sigma() {return sigma_;}
	DiscreteThetaType*  theta() { return theta_; }
  virtual DiscreteFunctionType* additionalVariables() { return additionalVariables_; }	
 
	std::string dataPrefix()
	{
		assert( problem_ );  
		return problem_->dataPrefix(); 
	}
	
  const InitialDataType& problem () const 
  { 
    assert( problem_ ); 
    return *problem_; 
  }

  const ModelType& model () const 
  { 
    assert( model_ ); 
    return *model_;
  }

	virtual void initializeStep(TimeProviderType& timeProvider)
	{
		DiscreteFunctionType& U = solution();
    //Create OdeSolver if necessary
    if( odeSolver_ == 0 ) odeSolver_ = createOdeSolver( timeProvider );    
		assert( odeSolver_ );
   
    
    if( interpolateInitialData_ )
      {
        LagrangeGridPartType lagGridPart(grid_); 
        InterpolationSpaceType interpolSpace(lagGridPart); 
        InterpolationFunctionType interpolSol("uinterpol",interpolSpace); 
  
        Dune::Fem::LagrangeInterpolation<InitialDataType,InterpolationFunctionType>::interpolateFunction( problem(),interpolSol);
        Dune::Fem::DGL2ProjectionImpl::project(interpolSol, U);
      }
     else 
      {   
       Dune::Fem::DGL2ProjectionImpl::project(problem().fixedTimeFunction( timeProvider.time()), U);
      }
      //MOVE IN SEPERATE METHOD
      odeSolver_->initialize( U );  


   }

  void step ( TimeProviderType& timeProvider,
						  int& newton_iterations,
						  int& ils_iterations,
						  int& max_newton_iterations,
						  int& max_ils_iterations)
	{
    DiscreteFunctionType& U=solution();
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
  void writeEnergy( TimeProviderType& timeProvider,
                    Stream& str)
  {
    DiscreteSigmaType* gradient = sigma();
    DiscreteScalarType* totalenergy = energy();

    if(gradient != nullptr  && totalenergy != nullptr)
    { 
      gradient->clear();
      
      dgOperator_.gradient(solution(),*gradient);
  
      totalenergy->clear();

      double kineticEnergy;
      
#if WELLBALANCED    
      double chemicalEnergy; 
      double energyIntegral =energyconverter(solution(),*gradient,model(),*totalenergy,kineticEnergy,chemicalEnergy);
      str<<std::setprecision(20)<< timeProvider.time()<<"\t"<<energyIntegral<<"\t"<<chemicalEnergy<<"\t"<<kineticEnergy<<"\n";
#else
      double energyIntegral =energyconverter(solution(),*gradient,model(),*totalenergy,kineticEnergy);
      str<<std::setprecision(20)<<timeProvider.time()<<"\t"<<energyIntegral<<"\t"<<kineticEnergy<<"\n";
#endif

    }
  
  }
  
  //! write data, if pointer to additionalVariables is true, they are calculated first 
  void writeData( DataWriterType& eocDataOutput,
									TimeProviderType& timeProvider,
                  const bool reallyWrite )
	{
    if( reallyWrite )
		{
				 
      DiscreteFunctionType* addVariables = additionalVariables();
      DiscreteSigmaType* gradient=sigma();
      DiscreteThetaType* theta1=theta();
  
      // calculate DG-projection of additional variables
			if ( addVariables && gradient)
			{
        gradient->clear();
  
        dgOperator_.gradient(solution(),*gradient);

        if(theta1)
        { 
        
          theta1->clear();
          dgOperator_.theta(solution(),*theta1);

        }

        // calculate additional variables from the current num. solution
        setupAdditionalVariables( solution(), *gradient,model(), *addVariables );
			}

            
    }
	// write the data 
		eocDataOutput.write( timeProvider );
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
    TimeProviderType timeProvider(startTime, grid_ ); 

    DiscreteFunctionType& U = solution();
    DiscreteFunctionType* Uold = oldsolution(); 


    RestrictionProlongationType rp( U );
    
    rp.setFatherChildWeight( Dune::DGFGridInfo<GridType> :: refineWeight() );
 
	// create adaptation manager 
    AdaptationManagerType adaptManager(grid_,rp);
		

    // restoreData if checkpointing is enabled (default is disabled)
		//restoreFromCheckPoint(timeProvider );
		
	// tuple with additionalVariables 
	
  
    IOTupleType dataTuple(&U, this->additionalVariables
        (),this->energy(),this->theta() );
	
   // IOTupleType dataTuple( &U, this->sigma(),this->theta() );
    std::ofstream energyfile;
    std::ostringstream convert;
    convert<<loopNumber;
    std::string filename=energyFilename_;
    filename.append(convert.str()); 
    filename.append(".gnu");
    energyfile.open(filename.c_str());

	
		// type of the data writer
		DataWriterType eocDataOutput( grid_, 
                                  dataTuple, 
                                  timeProvider, 
                                  EocDataOutputParameters( loop_ , dataPrefix() ) );
	  bool restart = Dune::Fem::Parameter::getValue<bool>("phasefield.restart",false);
		
    if(!restart)
    {
		// set initial data (and create ode solver)
		initializeStep(timeProvider );
	 		writeData( eocDataOutput, timeProvider, eocDataOutput.willWrite( timeProvider ) );
	
		// adapt the grid to the initial data
		int startCount = 0;
		if( adaptCount > 0 )
    {
      while( startCount < maxAdaptationLevel )
				{
          estimator_.estimateAndMark(tolerance_);
          writeData( eocDataOutput, timeProvider, eocDataOutput.willWrite( timeProvider ) );

          adaptManager.adapt();
					
					initializeStep(timeProvider );

					if( verbose )
						std::cout << "start: " << startCount << " grid size: " << grid_.size(0)<<std::endl;
          ++startCount;
    
            
        }
      }
   timeProvider.provideTimeStepEstimate(1e-4);

			if ( fixedTimeStep_ > 1e-20 )
			  timeProvider.init( fixedTimeStep_ );
		  else
			  timeProvider.init();

     

//   timeProvider.provideTimeStepEstimate(1e-4);                                         
    std::cout<<"deltaT "<<timeProvider.deltaT()<<" estimate "<<dgOperator_.timeStepEstimate()<<std::endl;
 
     writeEnergy( timeProvider ,std::cout);    
 		  writeData( eocDataOutput, timeProvider, eocDataOutput.willWrite( timeProvider ) );
  
    if(calcresidual_)
    {
      std::cout<<"Residual\n";
      Uold->clear();
      dgOperator_(U,*Uold);
      U.assign(*Uold);
      writeData(eocDataOutput ,timeProvider , eocDataOutput.willWrite( timeProvider ));
    }
    else
    for( ; timeProvider.time() < endTime && timeProvider.timeStep() < maximalTimeSteps ;  )   
		{ 
     timeProvider.provideTimeStepEstimate(maxTimeStep);                                         
			const double tnow  = timeProvider.time();
			const double ldt   = timeProvider.deltaT();
      int newton_iterations;
      int ils_iterations; 
      counter  = timeProvider.timeStep();
					
			Dune::FemTimer::start(timeStepTimer_);
			
      // grid adaptation (including marking of elements)
			if( (adaptCount > 0) && (counter % adaptCount) == 0 )
					estimateMarkAdapt( adaptManager );
      if(Uold!=nullptr)
      {
        Uold->clear();
        dgOperator_(U,*Uold);
      }
			step( timeProvider, newton_iterations, ils_iterations, 
            max_newton_iterations, max_ils_iterations);
      // Check that no NAN have been generated
			if (! U.dofsValid()) 
			{
  			std::cout << "Loop(" << loop_ << "): Invalid DOFs" << std::endl;
	   		eocDataOutput.write(timeProvider);
        energyfile.close();
        abort();
			}


      if(false)
        {
        Uold->clear();
        dgOperator_(U,*Uold);
       //Dune::Fem::DGL2ProjectionImpl::project(problem().fixedTimeFunction( timeProvider.time()), *Uold);
 
        if(false)
            {
          U.assign(*Uold);
          timeStepError = stepError(U,*Uold);
              timeStepError/=ldt;
          Uold->assign(U);
        }
      }
 
      double timeStepEstimate=dgOperator_.timeStepEstimate();	
  //    double diffTimeStep=dgOperator_.maxDiffusionTimeStep();
    //  double advTimeStep=dgOperator_.maxAdvectionTimeStep();
     
      if( (printCount > 0) && (counter % printCount == 0))
        {
          if( grid_.comm().rank() == 0 )
          {
    
            std::cout <<"step: " << counter << "  time = " << tnow << ", dt = " << ldt<<" ,timeStepEstimate " <<timeStepEstimate<<std::endl;
          if(false)
            {
              if(Uold!=nullptr)
              std::cout<< " ,Error between timesteps="<< timeStepError;
              std::cout<<std::endl;
            }   
          }
         
        writeEnergy( timeProvider , energyfile);

 
        }

    writeData( eocDataOutput , timeProvider , eocDataOutput.willWrite( timeProvider ) );
			writeCheckPoint( timeProvider , adaptManager );
    //statistics
    mindt = (ldt<mindt) ? ldt : mindt;
    maxdt = (ldt>maxdt) ? ldt : maxdt;
    averagedt += ldt;
    total_newton_iterations+=newton_iterations;
    total_ils_iterations+=ils_iterations;
  
      if(eocDataOutput.willWrite(timeProvider ) )    
        writeEnergy(timeProvider , energyfile);
  
 
    // next time step is prescribed by fixedTimeStep
    // it fixedTimeStep is not 0
    if ( fixedTimeStep_ > 1e-20 )
      timeProvider.next( fixedTimeStep_ );
    else
      timeProvider.next(); 
  
    }		/*end of timeloop*/
 
    // write last time step  
		writeData( eocDataOutput, timeProvider, true );

    energyfile.close();
	  averagedt /= double(timeProvider.timeStep());
	
		finalizeStep( timeProvider );                                  
    
		++loop_;
	
		// prepare the fixed time step for the next eoc loop
		fixedTimeStep_ /= fixedTimeStepEocLoopFactor_; 
	}
	
  inline double error ( TimeProviderType& timeProvider , DiscreteFunctionType& u )
	{
		typedef typename DiscreteFunctionType :: RangeType RangeType;
    Fem::L2Norm< GridPartType > l2norm(gridPart_);
    double error = l2norm.distance(problem().fixedTimeFunction(timeProvider.time()),u);
    return error;
	}
  
  //compute Error between old and NewTimeStep
  inline double stepError(DiscreteFunctionType uOld, DiscreteFunctionType& uNew)
	{
		typedef typename DiscreteFunctionType :: RangeType RangeType;
		Fem::L2Norm<GridPartType> l2norm(gridPart_);
	  return l2norm.distance(uOld, uNew);
	}

	virtual void finalizeStep(TimeProviderType& timeProvider)
	{ 
		DiscreteFunctionType& u = solution();
		bool doFemEoc = problem().calculateEOC( timeProvider, u, eocId_ );

		// ... and print the statistics out to a file
		if( doFemEoc )
			Fem::FemEoc::setErrors(eocId_, error(timeProvider, u));
	
		delete odeSolver_;
		odeSolver_ = 0;
		delete adaptationHandler_;
		adaptationHandler_ = 0;
	}

	//! restore all data from check point (overload to do something)
	virtual void restoreFromCheckPoint(TimeProviderType& timeProvider) {} 

	//! write a check point (overload to do something)
	virtual void writeCheckPoint(TimeProviderType& timeProvider,
															 AdaptationManagerType& am ) const {}

	virtual DiscreteFunctionType& solution () { return solution_; }



	virtual OdeSolverType* createOdeSolver(TimeProviderType&timeProvider) 
	{
		typedef PhaseFieldOdeSolver< DiscreteOperatorType > OdeSolverImpl;
		return new OdeSolverImpl(timeProvider, dgOperator_ );
	}

	virtual void finalize ( const int eocloop ) 
	{
		DiscreteFunctionType& U = solution(); 
  	DiscreteThetaType* theta1 = theta(); 
    DiscreteScalarType* totalenergy= energy(); 
 		DiscreteFunctionType* addVars = additionalVariables();
    // DiscreteSigmaType* sig= sigma();
	
    if( eocLoopData_ == 0 ) 
			{
			
        eocDataTup_ = IOTupleType( &U,addVars,totalenergy,theta1 ); 
        //eocDataTup_ = IOTupleType( &U,sig ); 
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

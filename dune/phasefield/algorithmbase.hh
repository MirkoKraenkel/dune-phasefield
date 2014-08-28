#ifndef ALGORITHMBASE_HH 
#define ALGORITHMBASE_HH 
// system includes
#include <sstream>
#include <dune/common/nullptr.hh>

// dune-fem includes
#include <dune/fem/misc/bartonnackmaninterface.hh>
#include <dune/fem/io/file/datawriter.hh>
#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/base/base.hh>
#include <dune/fem/operator/common/automaticdifferenceoperator.hh>
//Note Problen should be independent of Operator/Scheme
//#include<dune/phasefield/assembledtraits.hh>


// include std libs
#include <iostream>
#include <string>
// fem includes
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/componentl2norm.hh>

#include <dune/fem/operator/projection/l2projection.hh>
#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>

//adaption
#include <dune/fem-dg/operator/adaptation/estimatorbase.hh>
#include <dune/fem/space/discontinuousgalerkin/localrestrictprolong.hh>
#include <dune/fem-dg/operator/adaptation/adaptation.hh>
//post processing




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
          class Impl>
class PhasefieldAlgorithmBase
: public Fem::BartonNackmanInterface< PhasefieldAlgorithmBase< GridImp , AlgorithmTraits, Impl>, Impl >
{

public:
	//type of Grid
	typedef GridImp GridType;
  typedef Fem::BartonNackmanInterface< PhasefieldAlgorithmBase< GridImp , AlgorithmTraits, Impl>, Impl > BaseType;

  enum{ dimDomain= GridType::dimensionworld };

  //traits class gathers types depending on the problem and the operator which is specified there 
	typedef AlgorithmTraits Traits;
	//problem depending types
	
  typedef typename Traits::ProblemGeneratorType ProblemGeneratorType;
	typedef typename Traits::GridPartType GridPartType;
  //for interpolation of initial Data 
	typedef typename Traits::LagrangeGridPartType LagrangeGridPartType;
	
  //Problem Dependent Type
	typedef typename Traits :: InitialDataType             InitialDataType;
	typedef typename Traits :: ModelType                   ModelType;
	typedef typename Traits :: FluxType                    FluxType;
  
  //discrete spaces
  typedef typename Traits::DiscreteSpaceType       DiscreteSpaceType;
  typedef typename Traits::DiscreteScalarSpaceType DiscreteScalarSpaceType;
  
  //discrete functions
	typedef typename Traits::DiscreteFunctionType DiscreteFunctionType;
	typedef typename Traits::DiscreteScalarType   DiscreteScalarType;
 
  
  
  typedef typename DiscreteSpaceType::FunctionSpaceType FunctionSpaceType;

  //Operator
  typedef typename Traits :: DiscreteOperatorType DiscreteOperatorType;
 
  //Adaptation
  typedef typename Traits :: RestrictionProlongationType RestrictionProlongationType;
 	// type of adaptation manager 
	typedef Dune::Fem::AdaptationManager< GridType, RestrictionProlongationType > AdaptationManagerType;
  typedef AdaptationHandler< GridType,typename DiscreteSpaceType::FunctionSpaceType >  AdaptationHandlerType;
  
  typedef typename Dune::Fem::LagrangeDiscreteFunctionSpace< FunctionSpaceType, LagrangeGridPartType, 2> InterpolationSpaceType;
  typedef typename Dune::Fem::AdaptiveDiscreteFunction< InterpolationSpaceType > InterpolationFunctionType;


	////type of IOTuple 
  typedef typename Traits::IOTupleType IOTupleType;
  // type of data 
  // writer 
  typedef Dune::Fem::DataOutput< GridType,IOTupleType >    DataWriterType;
  typedef Dune::Fem::CheckPointer< GridType >   CheckPointerType;

	// type of ime provider organizing time for time loops 
	typedef Dune::Fem::GridTimeProvider< GridType >                 TimeProviderType;

	//MemberVaribles
protected:
	GridType&               grid_;
	GridPartType            gridPart_;
	DiscreteSpaceType       space_;
	DiscreteScalarSpaceType energySpace_;
	Dune::Fem::IOInterface* eocLoopData_;
	IOTupleType             eocDataTup_; 
	unsigned int            timeStepTimer_; 
	unsigned int            loop_ ; 
	double                  fixedTimeStep_;        
	double                  fixedTimeStepEocLoopFactor_;       
  std::string             energyFilename_;
  DiscreteFunctionType    solution_;
	DiscreteFunctionType    oldsolution_; 
  DiscreteScalarType*     energy_;
  const InitialDataType*  problem_;
  ModelType*              model_;
  AdaptationHandlerType*  adaptationHandler_;
  Timer                   overallTimer_;
  const unsigned int      eocId_;
  double tolerance_;
  bool interpolateInitialData_;
  bool computeResidual_;
  double timeStepTolerance_;
  ////////////////////////////////////////
public:
	//Constructor
	PhasefieldAlgorithmBase(GridType& grid):
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
		oldsolution_( "oldsolution", space() ),
	  energy_( Fem :: Parameter :: getValue< bool >("phasefield.energy", false) ? new DiscreteScalarType("energy",energyspace()) : nullptr ),
		problem_( ProblemGeneratorType::problem() ),
    model_( new ModelType( problem() ) ),
    adaptationHandler_( 0 ),
    overallTimer_(),
    eocId_( Fem::FemEoc::addEntry(std::string("L2error")) ),
    tolerance_(Fem::Parameter :: getValue< double >("phasefield.adaptTol", 100)),
    interpolateInitialData_( Fem :: Parameter :: getValue< bool >("phasefield.interpolinitial" , false ) ),
    computeResidual_( Fem :: Parameter :: getValue< bool >("phasefield.calcresidual" , false ) ),
    timeStepTolerance_( Fem :: Parameter :: getValue< double >( "phasefield.timesteptolerance",-1. ) )
    {
    }

  //! destructor 
  /*virtual*/ ~PhasefieldAlgorithmBase()
  {
		delete energy_;
		energy_=0;
		delete problem_ ;
		problem_ = 0;
		delete model_;
    model_=0;
    delete adaptationHandler_ ;
		adaptationHandler_ = 0;
  }

	//some acces methods
	//space
	DiscreteSpaceType& space () { return space_; }	

	DiscreteScalarSpaceType& energyspace () { return energySpace_; }
	
  size_t gridSize() const
	{
		size_t grSize=grid_.size(0);
		return grid_.comm().sum(grSize);
	}

  DiscreteFunctionType& oldsolution () {return oldsolution_;}
  DiscreteScalarType* energy ()        { return energy_;}

 
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

	void initializeStep ( TimeProviderType& timeprovider )
	{
		DiscreteFunctionType& U = solution();
		DiscreteFunctionType& Uold = oldsolution();
    //Create OdeSolver if necessary
    
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
       Dune::Fem::DGL2ProjectionImpl::project(problem().fixedTimeFunction(timeprovider.time()), U);
      }

     Uold.assign(U); 
  }

  void step ( TimeProviderType& timeProvider,
						  int& newton_iterations,
						  int& ils_iterations,
						  int& max_newton_iterations,
						  int& max_ils_iterations)
	{
 	  BaseType::asImp().step( timeProvider,
                            newton_iterations,
                            ils_iterations,
                            max_newton_iterations,
                            max_ils_iterations);
  }
	
	//! estimate and mark solution 
  /*virtual*/ void estimateMarkAdapt( AdaptationManagerType& am ){};

  template<class Stream>
  void writeEnergy( TimeProviderType& timeProvider,Stream& str)
  { 
    BaseType::asImp().writeEnergy( timeProvider, str);
  }
  
  //! write data, if pointer to additionalVariables is true, they are calculated first 
  virtual void writeData( DataWriterType& eocDataOutput,
									TimeProviderType& timeProvider,
                  const bool reallyWrite )
	{
    if( reallyWrite )
		{
	   //setup additionals variables if needed			 
    }
	  // write the data 
		eocDataOutput.write( timeProvider );
	}
	
  IOTupleType getDataTuple()
  {
   return BaseType::asImp().getDataTuple();
	}
	
  void initializeSolver( TimeProviderType& timeProvider)
  {
    BaseType::asImp().initializeSolver( timeProvider);
  }
  /*virtual*/ void operator()(int loopNumber,double& averagedt, double& mindt, double& maxdt,
                          size_t& counter, int& total_newton_iterations, int& total_ils_iterations,
                          int& max_newton_iterations, int& max_ils_iterations)
	{	

    double timeStepError=std::numeric_limits<double>::max();
    //some setup stuff
		const bool verbose = Dune::Fem::Parameter :: verbose ();
  	int printCount = Dune::Fem::Parameter::getValue<int>("phasefield.printCount", -1);
    printCount+=loopNumber*printCount; 
		int adaptCount = 0;
		int maxAdaptationLevel = 0;
		int startLevel = 0 ;
        
     Dune::Fem::AdaptationMethod< GridType > am( grid_ );
		
     if( am.adaptive() )
			{
    		adaptCount = Dune::Fem::Parameter::getValue<int>("fem.adaptation.adaptcount");
	    	maxAdaptationLevel = Dune::Fem::Parameter::getValue<int>("fem.adaptation.finestLevel");
 		    startLevel = Dune::Fem::Parameter::getValue< int >("phasefield.startLevel",0);
        maxAdaptationLevel-=startLevel;
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

    CheckPointerType checkPointer( grid_ , timeProvider );


    DiscreteFunctionType& U = solution();
    DiscreteFunctionType& Uold = oldsolution(); 

    Dune::Fem::persistenceManager << solution(); 
    RestrictionProlongationType rp( U );
    
    rp.setFatherChildWeight( Dune::DGFGridInfo<GridType> :: refineWeight() );
 
    // create adaptation manager 
    AdaptationManagerType adaptManager(grid_,rp);
		

    // restoreData if checkpointing is enabled (default is disabled)
    IOTupleType dataTuple=getDataTuple();//( &solution() ,nullptr,  this->energy());
	
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
    	initializeStep( timeProvider );
	 		initializeSolver( timeProvider);
      writeData( eocDataOutput, timeProvider, eocDataOutput.willWrite( timeProvider ) );
 
  		// adapt the grid to the initial data
	  	int startCount = 0;
		  
      if( adaptCount > 0 )
      {
        while( startCount < maxAdaptationLevel )
				{
         // estimator_.estimateAndMark(tolerance_);
          //writeData( eocDataOutput, timeProvider, eocDataOutput.willWrite( timeProvider ) );
          adaptManager.adapt();
					if( verbose )
						std::cout << "start: " << startCount << " grid size: " << grid_.size(0)<<std::endl;
          ++startCount;
        }
      }
      timeProvider.provideTimeStepEstimate(maxTimeStep);
    	// start first time step with prescribed fixed time step 
	    // if it is not 0 otherwise use the internal estimate
			if ( fixedTimeStep_ > 1e-20 )
			  timeProvider.init( fixedTimeStep_ );
		  else
			  timeProvider.init();
 		  writeData( eocDataOutput, timeProvider, eocDataOutput.willWrite( timeProvider ) );
    }
    else
    {
      checkPointer.restoreData( grid_, "checkpoint" ); 
    }
    
    if(computeResidual_)
    {
      computeResidual( U,Uold,timeProvider,eocDataOutput);      
    }
    else
    for( ; timeProvider.time() < endTime && timeProvider.timeStep() < maximalTimeSteps ;  )   
		{ 
			const double tnow  = timeProvider.time();
			const double ldt   = timeProvider.deltaT();
      int newton_iterations;
      int ils_iterations; 
      counter  = timeProvider.timeStep();
					
			Dune::FemTimer::start(timeStepTimer_);
			
      // grid adaptation (including marking of elements)
			if( (adaptCount > 0) && (counter % adaptCount) == 0 )
					estimateMarkAdapt( adaptManager );
		
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

     double timeStepEstimate=0;//dgOperator_.timeStepEstimate();	
     
     if( (printCount > 0) && (counter % printCount == 0))
			{
       // if( grid_.comm().rank() == 0 )
        {
          std::cout <<"step: " << counter << "  time = " << tnow << ", dt = " << ldt<<" ,timeStepEstimate " <<timeStepEstimate;
          timeStepError=stepError(U,Uold);	
          timeStepError/=ldt;
          std::cout<< " ,Error between timesteps="<< timeStepError;
          std::cout<<std::endl;
             
          std::cout<<" linearIterations: "<<ils_iterations<<" newtonIterations: "<<newton_iterations<<std::endl;
        }
         
        writeEnergy( timeProvider , energyfile);
     }
    writeData( eocDataOutput , timeProvider , eocDataOutput.willWrite( timeProvider ) );
    checkPointer.write(timeProvider);
    //statistics
    mindt = (ldt<mindt) ? ldt : mindt;
    maxdt = (ldt>maxdt) ? ldt : maxdt;
    averagedt += ldt;
    total_newton_iterations+=newton_iterations;
    total_ils_iterations+=ils_iterations;
  
    // next time step is prescribed by fixedTimeStep
    // it fixedTimeStep is not 0
    if ( fixedTimeStep_ > 1e-20 )
      timeProvider.next( fixedTimeStep_ );
    else
      timeProvider.next(); 
    Uold.assign(U); 
     
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
  
  void computeResidual ( DiscreteFunctionType& U, 
                         DiscreteFunctionType& Uold,
                         TimeProviderType& timeProvider,
                         DataWriterType& eocDataOutput)
  {
  }
  inline double error ( TimeProviderType& timeProvider , DiscreteFunctionType& u )
	{
		typedef typename DiscreteFunctionType :: RangeType RangeType;
    Fem::L2Norm< GridPartType > l2norm(gridPart_);
    double error = l2norm.distance(problem().fixedTimeFunction(timeProvider.time()),u);
    return error;
	}
  
  inline double densityError ( TimeProviderType& timeProvider , DiscreteFunctionType& u )
	{
		typedef typename DiscreteFunctionType :: RangeType RangeType;
    std::vector<unsigned int> comp{0.};
    Fem::ComponentL2Norm< GridPartType > l2norm(gridPart_, comp );
    double error = l2norm.distance(problem().fixedTimeFunction(timeProvider.time()),u);
    return error;
	}
  inline double phasefieldError ( TimeProviderType& timeProvider , DiscreteFunctionType& u )
	{
		typedef typename DiscreteFunctionType :: RangeType RangeType;
    std::vector<unsigned int> comp{ dimDomain +1 };
    Fem::ComponentL2Norm< GridPartType > l2norm(gridPart_, comp );
    double error = l2norm.distance(problem().fixedTimeFunction(timeProvider.time()),u);
    return error;
	}
  inline double velocityError ( TimeProviderType& timeProvider , DiscreteFunctionType& u )
	{
		typedef typename DiscreteFunctionType :: RangeType RangeType;
    std::vector<unsigned int> comp;
    for (int i = 0 ; i < dimDomain ; ++i)
      comp.push_back(1+i);
    Fem::ComponentL2Norm< GridPartType > l2norm(gridPart_, comp );
    double error = l2norm.distance(problem().fixedTimeFunction(timeProvider.time()),u);
    return error;
	}
  
  //compute Error between old and NewTimeStep
  inline double stepError(DiscreteFunctionType uOld, DiscreteFunctionType& uNew)
	{
		typedef typename DiscreteFunctionType :: RangeType RangeType;
//		Fem::L2Norm<GridPartType> l2norm(gridPart_);
    std::vector<unsigned int> comp{ dimDomain +2 , dimDomain+3 };
    Fem::ComponentL2Norm< GridPartType > l2norm(gridPart_, comp ,POLORDER);
    
    return l2norm.distance(uOld, uNew);
	}

	/*virtual*/ void finalizeStep(TimeProviderType& timeProvider)
	{ 
		DiscreteFunctionType& u = solution();
		bool doFemEoc = problem().calculateEOC( timeProvider, u, eocId_ );

		// ... and print the statistics out to a file
		if( doFemEoc )
			Fem::FemEoc::setErrors(eocId_, phasefieldError(timeProvider, u));
	
		delete adaptationHandler_;
		adaptationHandler_ = 0;
	}

	//! restore all data from check point (overload to do something)
	/*virtual*/ void restoreFromCheckPoint(TimeProviderType& timeProvider) {} 

	//! write a check point (overload to do something)
	/*virtual*/ void writeCheckPoint(TimeProviderType& timeProvider,
															 AdaptationManagerType& am ) const {}

	/*virtual*/ DiscreteFunctionType& solution () { return solution_; }


	/*virtual*/ void finalize ( const int eocloop ) 
	{
		DiscreteFunctionType& U = solution(); 
	
    if( eocLoopData_ == 0 ) 
			{
       eocDataTup_ = getDataTuple(); 
       eocLoopData_ = new DataWriterType( grid_, eocDataTup_ );
      }

			{
				problem().finalizeSimulation( U,  eocloop );
			}
		eocLoopData_->writeData( eocloop );
	}
};



}//end namespace Dune








#endif

#ifndef FEMHOWTO_BASEEVOL_HH
#define FEMHOWTO_BASEEVOL_HH

// system includes
#include <sstream>

// dune-fem includes
#include <dune/fem/io/file/datawriter.hh>
#include <dune/fem/space/common/adaptmanager.hh>

// local includes
#include "base.hh"
#include "../operator/cons2prim.hh"
#include "../operator/energyconverter.hh"
/////////////////////////////////////////////////////////////////////////////
//
//  EOC output parameter class 
//
/////////////////////////////////////////////////////////////////////////////

struct EocDataOutputParameters :   /*@LST1S@*/
       public Dune::LocalParameter<Dune::DataWriterParameters,EocDataOutputParameters> 
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
};                              /*@LST1E@*/

/////////////////////////////////////////////////////////////////////////////
//
//  Base class for evolutionary problems 
//
/////////////////////////////////////////////////////////////////////////////
template <class TraitsImp> 
class AlgorithmBase
{
  typedef TraitsImp Traits ;
public:
  // type of Grid
  typedef typename Traits :: GridType                  GridType;

  // type of the problem
  typedef typename Traits :: InitialDataType           ProblemType;

  // type of the model 
  typedef typename Traits :: ModelType                 ModelType;

  // choose a suitable GridView
  typedef typename Traits :: GridPartType              GridPartType;

  // the space type
  typedef typename Traits :: DiscreteSpaceType         DiscreteSpaceType;

  // the discrete function type
  typedef typename Traits :: DiscreteFunctionType      DiscreteFunctionType;
 // the discrete function type
  typedef typename Traits :: DiscreteGradientType      DiscreteGradientType;
//     // the space type
// //   typedef typename Traits :: ScalarDiscreteSpaceType         ScalarDiscreteSpaceType;

//   // the discrete function type
//   typedef typename Traits :: ScalarDFType      ScalarDFType;



  // the indicator function type (for limiting only)
  typedef typename Traits :: IndicatorType             IndicatorType;

  // the type of domain and range
  typedef typename DiscreteFunctionType::DomainType    DomainType;
  typedef typename DiscreteFunctionType::RangeType     RangeType;

  // the restriction and prolongation projection 
  typedef typename Traits :: RestrictionProlongationType RestrictionProlongationType;

  // type of adaptation manager 
  typedef Dune::AdaptationManager< GridType, RestrictionProlongationType > AdaptationManagerType;

  // type of IOTuple 
  typedef Dune::tuple< DiscreteFunctionType*,  DiscreteFunctionType*,IndicatorType* > IOTupleType; 
  // type of data writer 
  typedef Dune::DataWriter< GridType, IOTupleType >    DataWriterType;

  // type of time provider organizing time for time loops 
  typedef Dune::GridTimeProvider< GridType >                 TimeProviderType;

  //! constructor 
  AlgorithmBase(GridType& grid) 
   : grid_( grid ),
     gridPart_( grid_ ),
     space_( gridPart_ ),
     eocLoopData_( 0 ),
     eocDataTup_(),
     // Initialize Timer for CPU time measurements
     timeStepTimer_( Dune::FemTimer::addTo("max time/timestep") ),
     loop_( 0 ),
     fixedTimeStep_( Dune::Parameter::getValue<double>("fixedTimeStep",0) ),
     fixedTimeStepEocLoopFactor_( Dune::Parameter::getValue<double>("fixedTimeStepEocLoopFactor",1.) )
  {
  }

  //! also have virtual destructor 
  virtual ~AlgorithmBase() 
  {
    delete eocLoopData_; 
    eocLoopData_ = 0 ;
  } 

  //! return reference to space 
  DiscreteSpaceType& space()
  {
    return space_;
  }
//   ScalarDiscreteSpaceType& energyspace()
//   {
//     return energyspace_;
//   }

  // return size of grid 
  virtual size_t gridSize() const 
  { 
    size_t grSize = grid_.size( 0 );
    return grid_.comm().sum( grSize );
  }

  virtual std::string dataPrefix () const = 0;

  virtual const ProblemType& problem() const = 0;
  virtual const ModelType& model() const = 0;

  // return reference to solution 
  virtual DiscreteFunctionType& solution () = 0;


 //  virtual DiscreteGradientType& theta () = 0;
  
  
//   virtual ScalarDFType& energy () = 0;
 
  // return reference to additional variables 
  virtual DiscreteFunctionType* additionalVariables () = 0;
  // return reference to additional variables 
  virtual IndicatorType* indicator()  { return 0; }

  //! initialize method for time loop, i.e. L2-project initial data 
  virtual void initializeStep(TimeProviderType& tp) = 0;

  //! solve one time step 
  virtual void step(TimeProviderType& tp,
                    int& newton_iterations,
                    int& ils_iterations,
                    int& max_newton_iterations,
                    int& max_ils_iterations) = 0;

  //! call limiter if necessary 
  virtual void limitSolution () {} 

  //! finalize problem, i.e. calculate EOC ...
  virtual void finalizeStep(TimeProviderType& tp) = 0;

  //! restore all data from check point (overload to do something)
  virtual void restoreFromCheckPoint(TimeProviderType& tp) {} 

  //! write a check point (overload to do something)
  virtual void writeCheckPoint(TimeProviderType& tp,
                               AdaptationManagerType& am ) const {}

  //! estimate and mark cell for refinement/coarsening
  virtual void estimateMarkAdapt( AdaptationManagerType& am ) = 0;

  //! write data, if pointer to additionalVariables is true, they are calculated first 
  void writeData( DataWriterType& eocDataOutput, 
                  TimeProviderType& tp, 
                  const bool reallyWrite )
  {

    if( reallyWrite ) 
    {
      DiscreteFunctionType* additionalVariables = this->additionalVariables();
      // calculate DG-projection of additional variables
      if ( additionalVariables )
      {
	// calculate additional variables from the current num. solution
        setupAdditionalVariables( solution(), model(), *additionalVariables );
      }
    
      // energyconverter(solution(),theta(),model(),energy());  
    }

    // write the data 
    eocDataOutput.write( tp );    
  }


  //! default time loop implementation, overload for changes 
  virtual void operator()(double& averagedt, double& mindt, double& maxdt, size_t& counter,
                          int& total_newton_iterations, int& total_ils_iterations,
                          int& max_newton_iterations, int& max_ils_iterations)
  {
    const bool verbose = Dune::Parameter :: verbose ();
    int printCount = Dune::Parameter::getValue<int>("femhowto.printCount", -1);

    // if adaptCount is 0 then no dynamics grid adaptation
    int adaptCount = 0;
    int maxAdaptationLevel = 0;
    Dune::AdaptationMethod< GridType > am( grid_ );
    if( am.adaptive() )
      {
	adaptCount = Dune::Parameter::getValue<int>("fem.adaptation.adaptcount");
	maxAdaptationLevel = Dune::Parameter::getValue<int>("fem.adaptation.finestLevel");
      }

    
    double maxTimeStep =
      Dune::Parameter::getValue("femhowto.maxTimeStep", std::numeric_limits<double>::max());
    const double startTime = Dune::Parameter::getValue<double>("femhowto.startTime", 0.0);
    const double endTime   = Dune::Parameter::getValue<double>("femhowto.endTime");

    const int maximalTimeSteps =
      Dune::Parameter::getValue("femhowto.maximaltimesteps", std::numeric_limits<int>::max());

    // for statistics
    maxdt     = 0.;
    mindt     = std::numeric_limits<double>::max();
    averagedt = 0.;
    total_newton_iterations = 0;
    total_ils_iterations = 0;
    max_newton_iterations = 0;
    max_ils_iterations = 0;

    // Initialize TimeProvider
    TimeProviderType tp(startTime, grid_ );  /*@\label{fv:timeprovider}@*/

    DiscreteFunctionType& solution = this->solution(); 
  //   DiscreteGradientType& theta1 = this->theta(); 
//     ScalarDFType& en=this->energy();
    RestrictionProlongationType rp( solution );

    // set refine weight 
    rp.setFatherChildWeight( Dune::DGFGridInfo<GridType> :: refineWeight() );

    // create adaptation manager 
    AdaptationManagerType adaptManager( grid_ , rp );

    // restoreData if checkpointing is enabled (default is disabled)
    restoreFromCheckPoint( tp );

    // tuple with additionalVariables 
    IOTupleType dataTuple( &solution,   this->additionalVariables(), indicator() );

    // type of the data writer
    DataWriterType eocDataOutput( grid_, dataTuple, tp, EocDataOutputParameters( loop_ , dataPrefix() ) );

    // set initial data (and create ode solver)
    initializeStep( tp ); /*@LST0E@*/ /*@\label{fv:stepperInitialize}@*/

    // start first time step with prescribed fixed time step 
    // if it is not 0 otherwise use the internal estimate
    tp.provideTimeStepEstimate(maxTimeStep);
    if ( fixedTimeStep_ > 1e-20 )
      tp.init( fixedTimeStep_ );
    else
      tp.init();

    // adapt tshe grid to the initial data
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

    // write data 
     writeData( eocDataOutput, tp, eocDataOutput.willWrite( tp ) );

    //**********************************************                    /*@LST0S@*/
    //* Time Loop                                  *
    //**********************************************
    for( ; tp.time() < endTime; )    /*@\label{fv:timeLoop}@*/
    {
      tp.provideTimeStepEstimate(maxTimeStep);                                          /*@LST0E@*/
      const double tnow  = tp.time();
      const double ldt   = tp.deltaT();
      int newton_iterations;
      int ils_iterations;
      counter  = tp.timeStep();

      //************************************************        /*@LST0S@*/
      //* Compute an ODE timestep                      *
      //************************************************
      Dune::FemTimer::start(timeStepTimer_);

      // grid adaptation (including marking of elements)
      if( (adaptCount > 0) && (counter % adaptCount) == 0 )
        estimateMarkAdapt( adaptManager );

      // maximality of max_newton_iterations, max_ils_iterations
      // is dealt within this step(...) function
      step( tp, newton_iterations, ils_iterations,
            max_newton_iterations, max_ils_iterations ); /*@\label{fv:step}@*/

      Dune::FemTimer::stop(timeStepTimer_,Dune::FemTimer::max);              /*@LST0E@*/
      
  
      
      // Check that no NAN have been generated
      if (! solution.dofsValid()) {
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
      //theta1.print(cout);
      // write data 
      writeData( eocDataOutput, tp, eocDataOutput.willWrite( tp ) );

      // possibly write check point (default is disabled)
      writeCheckPoint( tp, adaptManager );

      // some statistics
      mindt = (ldt<mindt) ? ldt : mindt;
      maxdt = (ldt>maxdt) ? ldt : maxdt;
      averagedt += ldt;
      total_newton_iterations += newton_iterations;
      total_ils_iterations += ils_iterations;

      // next time step is prescribed by fixedTimeStep
      // it fixedTimeStep is not 0
      if ( fixedTimeStep_ > 1e-20 )
        tp.next( fixedTimeStep_ );
      else
        tp.next();

      // for debugging only 
      if( tp.timeStep() >= maximalTimeSteps ) break ;

      // possibly write a grid solution
      problem().postProcessTimeStep( grid_, solution, tp.time() );

    } /****** END of time loop *****/ /*@LST0S@*/

    // write last time step  
    writeData( eocDataOutput, tp, true );

    averagedt /= double(tp.timeStep());
    if(verbose)
    {
      std::cout << "Minimum dt: " << mindt
         << "\nMaximum dt: " << maxdt
         << "\nAverage dt: " << averagedt << std::endl;
    }

    // finalize eoc step 
    finalizeStep( tp );                                   /*@\label{fv:finalize}@*/
    
    // increase loop counter
    ++loop_;

    // prepare the fixed time step for the next eoc loop
    fixedTimeStep_ /= fixedTimeStepEocLoopFactor_; 
  }


  //! finalize problem, i.e. calculated EOC ...
  virtual void finalize( const int eocloop ) 
  {
    DiscreteFunctionType& U = solution(); 
 //    DiscreteGradientType& theta1 = theta(); 
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

protected:      
  GridType&            grid_;
  GridPartType         gridPart_;
  DiscreteSpaceType    space_;
  // ScalarDiscreteSpaceType    energyspace_;
  Dune::IOInterface*   eocLoopData_;
  IOTupleType          eocDataTup_; 
  unsigned int timeStepTimer_; 
  unsigned int loop_ ; 

  // use fixed time step if fixedTimeStep>0
  double fixedTimeStep_;        
  double fixedTimeStepEocLoopFactor_;        
};

#endif

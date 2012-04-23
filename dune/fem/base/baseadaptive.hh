#ifndef FEMHOWTO_BASEADAPTIVE_HH
#define FEMHOWTO_BASEADAPTIVE_HH

#include <dune/fem/space/common/adaptmanager.hh>

//using namespace Dune;

#include "base.hh"

// return true if adaptation is used 
template <class HGridType>
bool checkAdaptive(HGridType& grid)
{
  Dune::AdaptationMethod<HGridType> am( grid );
  return am.adaptive();
}


template <class Algorithm>
static void computeAdaptive(Algorithm & algorithm) 
{
  typedef typename Algorithm::DiscreteSpaceType         DiscreteSpaceType;
  typedef typename Algorithm::DiscreteFunctionType      DiscreteFunctionType;
  typedef typename Algorithm::GridType                  GridType;
  typedef typename Algorithm::GridPartType              GridPartType;
  typedef typename Algorithm::IOTupleType               IOTupleType;
  DiscreteSpaceType& space = algorithm.space();
  GridPartType& gridPart = space.gridPart();
  GridType& grid = gridPart.grid();

  // if adaptation is disabled fallback to normal compute 
  if ( ! checkAdaptive( grid )) 
  {
    std::cout << "Adaptivity disabled, fallback to EOC calculation! " << std::endl;
    compute(algorithm);
    return ;
  }

  // solution function
  DiscreteFunctionType u( "solution", space );
#ifdef WRITE_PRIM_VARIABLES
  DiscreteFunctionType additionalVariables( "additional", space );
#endif

  // get tolerance for adaptive problem 
  double tolerance = 1.0;  // default value 
  tolerance = Dune::Parameter::getValue("femhowto.tolerance", tolerance);
  // get some parameters
  const int eocSteps   = Dune::Parameter::getValue<int>("femhowto.eocSteps", 1);

  // Initialize the DataOutput that writes the solution on the harddisk in a
  // format readable by e.g. Paraview in each adaption step
#ifdef WRITE_PRIM_VARIABLES
  // tuple carrying PDE unknows and additional variables
  IOTupleType dataTup ( &u, &additionalVariables );
#else
  IOTupleType dataTup ( &u );
#endif
  typedef Dune::DataOutput<GridType, IOTupleType> DataOutputType;
  DataOutputType dataOutput( grid, dataTup );

  const unsigned int adaptStepTimer = Dune::FemTimer::addTo("totaltime (time of first step, max time of steps)", 2);
  double avgTimeStep = 0;
  double minTimeStep = 0;
  double maxTimeStep = 0;
  size_t time_step_counter = 0;
  int total_newton_iterations = 0;
  int total_ils_iterations = 0;
  int max_newton_iterations = 0;
  int max_ils_iterations = 0;

  // do one time step for initial adaptation
  Dune::FemTimer::start(adaptStepTimer);
  Dune::FemTimer::start(adaptStepTimer, 1);
#ifdef WRITE_PRIM_VARIABLES
  algorithm( u, additionalVariables, avgTimeStep, minTimeStep, maxTimeStep,
             time_step_counter, total_newton_iterations, total_ils_iterations,
             max_newton_iterations, max_ils_iterations );
#else
  algorithm( u, u, avgTimeStep, minTimeStep, maxTimeStep,
             time_step_counter, total_newton_iterations, total_ils_iterations,
             max_newton_iterations, max_ils_iterations );
#endif
  Dune::FemTimer::stop(adaptStepTimer, 1);
  Dune::FemTimer::stop(adaptStepTimer, Dune::FemTimer::sum);

#ifdef WRITE_PRIM_VARIABLES
  // finalize it with additional variables, which are of more interest anyway
  algorithm.finalize( additionalVariables, 1 );
#else
  algorithm.finalize( u, 1 );
#endif
  dataOutput.writeData(0);

  typedef typename Algorithm :: RestrictionProlongationType RestrictionProlongationType;
  typedef typename Algorithm :: AdaptationManagerType AdaptationManagerType;

  // create Restriction and Prolongation Operator 
  RestrictionProlongationType rpSolution( u );
  // create Adaptation Manager 
  AdaptationManagerType adaptManager( grid, rpSolution );

  int counter=1;
  while ( algorithm.estimateAndMark(u, tolerance) )
  {
    // adapt grid and restrict and prolong solution (see AdaptionManager)
    adaptManager.adapt();

    // re-calculate solution 
    Dune::FemTimer::start(adaptStepTimer);
    Dune::FemTimer::start(adaptStepTimer, 2);
#ifdef WRITE_PRIM_VARIABLES
    algorithm( u, additionalVariables, avgTimeStep, minTimeStep, maxTimeStep,
               time_step_counter, total_newton_iterations, total_ils_iterations,
               max_newton_iterations, max_ils_iterations );
#else
    algorithm( u, u, avgTimeStep, minTimeStep, maxTimeStep,
               time_step_counter, total_newton_iterations, total_ils_iterations,
               max_newton_iterations, max_ils_iterations );
#endif
    Dune::FemTimer::stop(adaptStepTimer, 2, Dune::FemTimer::max);
    Dune::FemTimer::stop(adaptStepTimer, Dune::FemTimer::sum);

#ifdef WRITE_PRIM_VARIABLES
    // finalize it with additional variables, which are of more interest anyway
    algorithm.finalize( additionalVariables, 1 );
#else
    algorithm.finalize( u, 1 );
#endif
    // write data file for time and timestep (internal step number is increased)
    dataOutput.writeData(counter);
    ++counter;
  }
  Dune::FemTimer::printFile("./timer.out");
  Dune::FemTimer::removeAll();
}
#endif

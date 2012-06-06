#ifndef FEMHOWTO_BASE_HH
#define FEMHOWTO_BASE_HH

#include <sstream>

#include <dune/common/misc.hh>
#include <dune/common/fvector.hh>
#include <dune/common/version.hh>

#if ! DUNE_VERSION_NEWER_REV(DUNE_GRID,2,1,0)
#include <dune/grid/io/file/dgfparser/dgfgridtype.hh>
#endif

#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/misc/utility.hh>
#include <dune/fem/misc/l2error.hh>
#include <dune/fem/misc/femeoc.hh>
#include <dune/fem/misc/femtimer.hh>
#include <dune/fem/misc/gridwidth.hh>
#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/io/file/datawriter.hh>


//! get memory in MB 
inline double getMemoryUsage()
{
  struct rusage info;
  getrusage( RUSAGE_SELF, &info );
  return (info.ru_maxrss / 1024.0);
}

template< class HGridType >
Dune::GridPtr< HGridType > initialize( const std::string& problemDescription )
{ 
  // ----- read in runtime parameters ------
  const std::string filekey = Dune::IOInterface::defaultGridKey( HGridType::dimension );
  const std::string filename = Dune::Parameter::getValue< std::string >( filekey ); /*@\label{base:param0}@*/

  // initialize grid
  Dune::GridPtr< HGridType > gridptr(filename);
  Dune::Parameter::appendDGF( filename );

  // load balance grid in case of parallel runs 
  gridptr->loadBalance();

  // output of error and eoc information
  std::string eocOutPath = Dune::Parameter::getValue<std::string>("femhowto.eocOutputPath",  /*@\label{fv:param1}@*/
                                                            std::string("./"));

  Dune::FemEoc::initialize(eocOutPath, "eoc", problemDescription); /*@\label{fv:femeocInit}@*/

  // and refine the grid until the startLevel is reached
  const int startLevel = Dune::Parameter::getValue<int>("femhowto.startLevel", 0);
  for(int level=0; level < startLevel ; ++level)
    Dune::GlobalRefine::apply(*gridptr, 1 ); /*@\label{fv:globalRefine1}@*/
  return gridptr;
} /*@LST0E@*/


/////////////////////////////////////////////////////////////////////////////
//
//  compute algorithm 
//
/////////////////////////////////////////////////////////////////////////////
template <class Algorithm>
void compute(Algorithm& algorithm)
{
  typedef typename Algorithm::DiscreteFunctionType DiscreteFunctionType;
  typename Algorithm::DiscreteSpaceType& space = algorithm.space();
  typename Algorithm::GridPartType& gridPart = space.gridPart();
  typedef typename Algorithm::GridPartType::GridType GridType;
  GridType& grid = gridPart.grid();

  // get some parameters
  const int eocSteps   = Dune::Parameter::getValue<int>("femhowto.eocSteps", 1);

  typename Algorithm::IOTupleType dataTup ( &algorithm.solution() );
  typedef Dune::DataOutput<GridType, typename Algorithm::IOTupleType> DataOutputType;
  DataOutputType dataOutput( grid, dataTup );

  const unsigned int femTimerId = Dune::FemTimer::addTo("timestep");
  for(int eocloop=0; eocloop < eocSteps; ++eocloop)
  {
    Dune::FemTimer :: start(femTimerId);

    // do one step
    double avgTimeStep = 0;
    double minTimeStep = 0;
    double maxTimeStep = 0;
    size_t counter = 0;
    int total_newton_iterations = 0;
    int total_ils_iterations = 0;
    int max_newton_iterations = 0;
    int max_ils_iterations = 0;
#if MYALGO
		algorithm();

#else
    algorithm( avgTimeStep, minTimeStep, maxTimeStep,
               counter, total_newton_iterations, total_ils_iterations,
               max_newton_iterations, max_ils_iterations );
#endif
    double runTime = Dune::FemTimer::stop(femTimerId);

    Dune::FemTimer::printFile("./timer.out");
    Dune::FemTimer::reset(femTimerId);

    // calculate grid width
    const double h = Dune::GridWidth::calcGridWidth(gridPart);

    algorithm.finalize( eocloop );

    if( Dune::Parameter :: verbose() )
      Dune::FemEoc::write(h,grid.size(0), runTime, counter,avgTimeStep,minTimeStep,
                    maxTimeStep, total_newton_iterations, total_ils_iterations, 
                    max_newton_iterations, max_ils_iterations, std::cout);
    else
      Dune::FemEoc::write(h,grid.size(0),runTime,counter,avgTimeStep,minTimeStep,
                    maxTimeStep,total_newton_iterations,total_ils_iterations,
                    max_newton_iterations, max_ils_iterations);

    dataOutput.writeData(eocloop);

    // Refine the grid for the next EOC Step. If the scheme uses adaptation,
    // the refinement level needs to be set in the algorithms' initialize method.
    if(eocloop < eocSteps-1)
    {
      Dune::GlobalRefine::apply(grid,Dune::DGFGridInfo<GridType>::refineStepsForHalf());
      grid.loadBalance();
    }
  } /***** END of EOC Loop *****/
  Dune::FemTimer::removeAll();
}

#endif

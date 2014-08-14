#ifndef FEMHOWTO_BASE_HH
#define FEMHOWTO_BASE_HH

#include <sstream>

#include <dune/common/fvector.hh>
#include <dune/common/version.hh>

#if ! DUNE_VERSION_NEWER_REV(DUNE_GRID,2,1,0)
#include <dune/grid/io/file/dgfparser/dgfgridtype.hh>
#endif

#include <dune/fem/misc/mpimanager.hh>
//#include <dune/fem/misc/utility.hh>
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
  bool restart = Dune::Fem::Parameter::getValue<bool>("phasefield.restart",false);

  Dune::GridPtr< HGridType > gridptr;
  
  if(!restart)
    {
      // ----- read in runtime parameters ------
      const std::string filekey = Dune::Fem::IOInterface::defaultGridKey( HGridType::dimension );
      const std::string filename = Dune::Fem::Parameter::getValue< std::string >( filekey ); /*@\label{base:param0}@*/

      // initialize grid
      gridptr = Dune::GridPtr< HGridType >(filename);
      Dune::Fem::Parameter::appendDGF( filename );

      // load balance grid in case of parallel runs 
      gridptr->loadBalance();

      // output of error and eoc information
      std::string eocOutPath = Dune::Fem::Parameter::getValue<std::string>("fem.prefix", std::string("."));

      Dune::Fem::FemEoc::initialize(eocOutPath, "eoc", problemDescription); 

      // and refine the grid until the startLevel is reached
      const int startLevel = Dune::Fem::Parameter::getValue<int>("phasefield.startLevel", 0);
      for(int level=0; level < startLevel ; ++level)
        Dune::Fem::GlobalRefine::apply(*gridptr, 1 ); 
  
      return gridptr;
    }
  else
    {
      std::cout<<"restore the grid\n";
      gridptr =Dune::Fem::CheckPointer< HGridType > :: restoreGrid( "checkpoint" );
      return gridptr;
    }
} 


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
  const int eocSteps   = Dune::Fem::Parameter::getValue<int>("phasefield.eocSteps", 1);
  typename Algorithm::IOTupleType dataTup =algorithm.getDataTuple();
  typedef Dune::Fem::DataOutput<GridType, typename Algorithm::IOTupleType> DataOutputType;
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
		
    algorithm(eocloop , 
              avgTimeStep, 
              minTimeStep, 
              maxTimeStep,
              counter, 
              total_newton_iterations, 
              total_ils_iterations,
              max_newton_iterations, 
              max_ils_iterations );


    double runTime = Dune::FemTimer::stop(femTimerId);

    Dune::FemTimer::printFile("./timer.out");
    Dune::FemTimer::reset(femTimerId);

    // calculate grid width
    const double h = Dune::Fem::GridWidth::calcGridWidth(algorithm.space().gridPart());

    algorithm.finalize( eocloop );
    if( Dune::Fem::Parameter :: verbose() )
      Dune::Fem::FemEoc::write(h,grid.size(0), runTime, counter,avgTimeStep,minTimeStep,
                    maxTimeStep, total_newton_iterations, total_ils_iterations, 
                    max_newton_iterations, max_ils_iterations, std::cout);
    else
      Dune::Fem::FemEoc::write(h,grid.size(0),runTime,counter,avgTimeStep,minTimeStep,
                    maxTimeStep,total_newton_iterations,total_ils_iterations,
                    max_newton_iterations, max_ils_iterations);
   
    DataOutputType dataOutput( grid, dataTup );

 //       dataOutput.writeData(eocloop);
    // Refine the grid for the next EOC Step. If the scheme uses adaptation,
    // the refinement level needs to be set in the algorithms' initialize method.
    if(eocloop < eocSteps-1)
    {
      Dune::Fem::GlobalRefine::apply(grid,Dune::DGFGridInfo<GridType>::refineStepsForHalf());
    }


    } /***** END of EOC Loop *****/
  Dune::FemTimer::removeAll();
}

#endif

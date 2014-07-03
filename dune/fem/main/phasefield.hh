
// in dbug mode also enable FieldVector checking and dune devel mode
#ifndef DNDEBUG 
#define DUNE_ISTL_WITH_CHECKING
#define DUNE_DEVEL_MODE
#endif


#include <config.h>

#include <dune/common/version.hh>
#include <dune/fem/misc/threads/threadmanager.hh>

#if defined USE_SMP_PARALLEL 
#if HAVE_DUNE_FEM_DG 
#define NSMOD_USE_SMP_PARALLEL
#endif
#endif

<<<<<<< Updated upstream
<<<<<<< Updated upstream
#if PARAGRID
#include <dune/grid/parallelgrid.hh>
#include <dune/grid/parallelgrid/dgfparser.hh>
=======
#if HAVE_DUNE_METAGRID
#include <dune/grid/parallelgrid.hh>
>>>>>>> Stashed changes
=======
#if HAVE_DUNE_METAGRID
#include <dune/grid/parallelgrid.hh>
>>>>>>> Stashed changes
#endif

#include <dune/fem/base/base.hh>

// problem dependent
#include <problemcreator.hh>
#if MIXED
#include <dune/phasefield/assembledalgorithm.hh>
#else
#include <dune/phasefield/phasefieldalgorithm.hh> 
#endif
#include <dune/fem/space/common/allgeomtypes.hh>


namespace simulation{

  void simulate()
  {
    Fem::FemEoc::clear();

#if PARAGRID 
   typedef Dune::GridSelector :: GridType HostGridType;
   typedef ProblemGenerator< HostGridType > ProblemGeneratorType;

    // use problem specific initialize method since some problems do different things
    // there, e.g. poisson 
		const std::string Flux="Phasefieldflux";
    Dune::GridPtr<HostGridType> gridptr = ProblemGeneratorType :: initializeGrid( Flux );

    // get grid reference 
    typedef ParallelGrid< HostGridType > GridType;
    GridType grid( *gridptr );
    grid.loadBalance();
#else
    typedef Dune::GridSelector :: GridType GridType;
    typedef ProblemGenerator< GridType > ProblemGeneratorType;

    // use problem specific initialize method since some problems do different things
    // there, e.g. poisson 
		const std::string Flux="Phasefieldflux";
    Dune::GridPtr<GridType> gridptr = ProblemGeneratorType :: initializeGrid( Flux );

    // get grid reference 
    GridType& grid = *gridptr;

#endif 
    typedef AlgorithmTraits<GridType,ProblemGeneratorType,POLORDER> AlgoTraits;
		PhasefieldAlgorithm<GridType,AlgoTraits,POLORDER> stepper(grid);
   
    //defined in base.hh
    compute( stepper );
    //  compute( grid );
  } 

} // end namespace LOOPSPACE

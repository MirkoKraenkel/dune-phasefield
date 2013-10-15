
// in dbug mode also enable FieldVector checking and dune devel mode
#ifndef DNDEBUG 
#define DUNE_ISTL_WITH_CHECKING
#define DUNE_DEVEL_MODE
#endif

// -1 means higher order FV 
#if POLORDER == -1 
#define HIGHER_ORDER_FV 
#undef POLORDER 
#define POLORDER 0
#endif

#include <config.h>

#include <dune/common/version.hh>
#include <dune/fem/misc/threads/threadmanager.hh>

#if defined USE_SMP_PARALLEL 
#if HAVE_DUNE_FEM_DG 
#define NSMOD_USE_SMP_PARALLEL
#endif
#endif

#warning "AC_HH"

#include <dune/fem/base/base.hh>

// problem dependent
#include <acproblemcreator.hh>
#include <dune/phasefield/allencahnalgorithm.hh> 
#include <dune/fem/space/common/allgeomtypes.hh>

#include <dune/grid/io/visual/grapedatadisplay.hh>

namespace simulation{

  void simulate()
  {
    Fem::FemEoc::clear();

    typedef Dune::GridSelector :: GridType GridType;
    typedef ProblemGenerator< GridType > ProblemGeneratorType;


    // use problem specific initialize method since some problems do different things
    // there, e.g. poisson 
		const std::string Flux="Phasefieldflux";
    Dune::GridPtr<GridType> gridptr = ProblemGeneratorType :: initializeGrid( Flux );

    // get grid reference 
    GridType & grid = *gridptr;

		typedef AlgorithmTraits<GridType,ProblemGeneratorType,POLORDER> AlgoTraits;
		AllenCahnAlgorithm<GridType,AlgoTraits,POLORDER> stepper(grid);
   
    //defined in base.hh
    compute( stepper );
  } 

} // end namespace LOOPSPACE

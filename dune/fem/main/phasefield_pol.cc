#if defined GRIDDIM
#ifndef CODEDIM
#define CODEDIM GRIDDIM
#endif
#endif

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
#include <dune/fem/misc/threadmanager.hh>

#if defined USE_SMP_PARALLEL 
#if HAVE_DUNE_FEM_DG 
#define NSMOD_USE_SMP_PARALLEL
#endif
#endif

//#include "codegen.hh"

#include <dune/fem/base/base.hh>

// problem dependent
#include <problemcreator.hh>

#include "stepper.hh"

#include <dune/fem/space/common/allgeomtypes.hh>

#include <dune/grid/io/visual/grapedatadisplay.hh>

namespace simulation{

  void simulate()
  {
    FemEoc::clear();

    typedef Dune::GridSelector :: GridType GridType;
    typedef ProblemGenerator< GridType > ProblemTraits;

    // ProblemType is a Dune::Function that evaluates to $u_0$ and also has a
    // method that gives you the exact solution.
    //typedef NSProblemType< GridType > ProblemType;
    //ProblemType problem;

    // Note to me: problem description is for FemEOC
    const std::string advFlux  = ProblemTraits :: advectionFluxName();
    const std::string diffFlux = ProblemTraits :: diffusionFluxName();

    // use problem specific initialize method since some problems do different things
    // there, e.g. poisson 
    Dune::GridPtr<GridType> gridptr = ProblemTraits :: initializeGrid( advFlux + diffFlux );

    // get grid reference 
    GridType & grid = *gridptr;
		
    Stepper<GridType, 
            ProblemTraits, 
            POLORDER> stepper( grid );


    compute( stepper );
  } 

} // end namespace LOOPSPACE


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

#if PARAGRID
#include <dune/grid/parallelgrid.hh>
#include <dune/grid/parallelgrid/dgfparser.hh>
#endif

#include <dune/fem/base/base.hh>

// problem dependent
#if NSK
#include <src/problems/korteweg/nskproblemcreator.hh>
#include <src/problems/kortewegmixed/mixednskproblemcreator.hh>
#else
#if MIXED
#if SPLIT
#include <src/problems/mixedscheme/splitproblemcreator.hh>
#include <dune/phasefield/splitalgoderived.hh>
#else
#include <src/problems/mixedscheme/mixedproblemcreator.hh>
#include <dune/phasefield/assembledalgoderived.hh>
#endif
#else
#include <src/problems/passscheme/problemcreator.hh>
#include <dune/phasefield/passalgoderived.hh> 
#endif
#endif
#include <dune/fem/space/common/allgeomtypes.hh>


namespace simulation{

template< bool mixed,int Polorder,class GridImp>
struct SchemeTraits
{
};
#if MIXED
template< int Polorder,class GridImp>
struct SchemeTraits< true, Polorder, GridImp >
{
  typedef GridImp GridType;
  typedef MixedProblemGenerator<GridType> ProblemGeneratorType;
  typedef MixedAlgorithmTraits< GridType , ProblemGeneratorType ,Polorder> AlgorithmTraitsType;
  typedef AssembledAlgorithm< GridType,AlgorithmTraitsType> AlgorithmType;
};
#else
template< int Polorder,class GridImp>
struct SchemeTraits< false ,Polorder, GridImp >
{
  typedef GridImp GridType;
  typedef ProblemGenerator<GridType> ProblemGeneratorType;
  typedef AlgorithmTraits< GridType , ProblemGeneratorType , Polorder > AlgorithmTraitsType;
  typedef PassAlgorithm< GridType , AlgorithmTraitsType > AlgorithmType;
};

#endif

  void simulate()
  {
    Fem::FemEoc::clear();

#if MIXED
    const bool mixed = true;
#else
    const bool mixed = false;
#endif
#if PARAGRID 
   typedef Dune::GridSelector :: GridType HostGridType;
   typedef MixedProblemGenerator< HostGridType > ProblemGeneratorType;
      // use problem specific initialize method since some problems do different things
    // there, e.g. poisson 
		const std::string Flux="Phasefieldflux";
    Dune::GridPtr<HostGridType> gridptr = ProblemGeneratorType :: initializeGrid( Flux );
    
    // get grid reference 
    typedef ParallelGrid< HostGridType > GridType;
    typedef MixedAlgorithmTraits< GridType , ProblemGeneratorType ,POLORDER> AlgorithmTraitsType;
    typedef AssembledAlgorithm< GridType,AlgorithmTraitsType> AlgorithmType;

    GridType grid( *gridptr );
    grid.loadBalance();
#else
    typedef Dune::GridSelector :: GridType GridType;
    
    typedef SchemeTraits< mixed , POLORDER, GridType> SchemeTraitsType; 
  
    typedef typename SchemeTraitsType::ProblemGeneratorType ProblemGeneratorType;
//    typedef typename SchemeTraitsType::AlgorithmTraitsType AlgoTraits;
    typedef typename SchemeTraitsType::AlgorithmType AlgorithmType;

    // typedef MixedProblemGenerator< GridType > ProblemGeneratorType;

    // use problem specific initialize method since some problems do different things
    // there, e.g. poisson 
		const std::string Flux="Phasefieldflux";
    Dune::GridPtr<GridType> gridptr = ProblemGeneratorType :: initializeGrid( Flux );

    // get grid reference 
    GridType& grid = *gridptr;
#endif 
    AlgorithmType algorithm( grid );  
#if 0
typedef AlgorithmTraits<GridType,ProblemGeneratorType,POLORDER> AlgoTraits;
#if MIXED
	  AssembledAlgorithm< GridType,AlgoTraits >  stepper( grid );
#else
    PassAlgorithm< GridType, AlgoTraits > stepper( grid );  
#endif
#endif
    //defined in base.hh
    compute( algorithm );
  } 

} // end namespace LOOPSPACE

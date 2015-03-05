#ifndef PHASEFIELD_PROBLEMTYPE_HH
#define PHASEFIELD_PROBLEMTYPE_HH
#include <dune/common/version.hh>
#include <dune/fem/io/parameter.hh>

#if DUNE_VERSION_NEWER_REV(DUNE_GRID,2,1,0)
// DGF gridtype 
#else
#include <dune/grid/io/file/dgfparser/dgfgridtype.hh>
#endif

template< int dimension, bool mixed>
struct RangeTypeProvider;

template< int dimension>
struct RangeTypeProvider<dimension,true>
{
#if LAMBDASCHEME 
  enum{ rangeDim=3*dimension+4 };
#else
enum{ rangeDim=2*dimension+4 };
#endif

#if 0
#warning "FEMSCHEME"
  enum{ rangeDim=dimension+4 };
#endif
};


template< int dimension>
struct RangeTypeProvider<dimension,false>
{
  enum{ rangeDim=dimension+2};
};






///////////////////////////////////////
// AVAILABLE PROBLEMS
///////////////////////////////////////
#if PROBLEM==1
#include "travellingwave.hh"
#ifdef MIXED
typedef TravelProblem< GridSelector :: GridType, 
 RangeTypeProvider<GridSelector::GridType::dimensionworld,true> 
 >  PhaseProblemType;
#else
typedef TravelProblem< GridSelector :: GridType, 
 RangeTypeProvider<GridSelector::GridType::dimensionworld,false> 
 >  PhaseProblemType;
#endif
#elif PROBLEM==2
#if MIXED
#include "../mixedscheme/bubbleensemble2.hh"
typedef BubbleEnsemble< GridSelector :: GridType,
        RangeTypeProvider< GridSelector::GridType::dimensionworld,true>
        >PhaseProblemType;
#else
#include "../mixedscheme/bubbleensemble2.hh"
typedef BubbleEnsemble< GridSelector :: GridType,
        RangeTypeProvider< GridSelector::GridType::dimensionworld,false>
        >PhaseProblemType;
#endif
#elif PROBLEM==3
#include "../mixedscheme/mixproblem2d.hh"
#if MIXED
typedef MixProblem< GridSelector :: GridType,
        RangeTypeProvider< GridSelector::GridType::dimensionworld,true>
        >PhaseProblemType;
#else
typedef MixProblem< GridSelector :: GridType,
        RangeTypeProvider< GridSelector::GridType::dimensionworld,false>
        >PhaseProblemType;
#endif
#elif PROBLEM==4
#include "../mixedscheme/tanhproblem.hh"
#if MIXED
typedef TanhProblem< GridSelector :: GridType,
        RangeTypeProvider< GridSelector::GridType::dimensionworld,true>
        >PhaseProblemType;
#else
typedef TanhProblem< GridSelector :: GridType,
        RangeTypeProvider< GridSelector::GridType::dimensionworld,false>
        >PhaseProblemType;
#endif
#elif PROBLEM==5
#include "../mixedscheme/datareadproblem.hh"
#if MIXED
typedef DataReadProblem< GridSelector :: GridType,
        RangeTypeProvider< GridSelector::GridType::dimensionworld,true>
        >PhaseProblemType;
#else
typedef DatareadProblem< GridSelector :: GridType,
        RangeTypeProvider< GridSelector::GridType::dimensionworld,false>
        >PhaseProblemType;
#endif
#elif PROBLEM==6
#include "../mixedscheme/sourceproblem.hh"
#if MIXED
typedef SourceProblem< GridSelector :: GridType,
        RangeTypeProvider< GridSelector::GridType::dimensionworld,true>
        >PhaseProblemType;
#else
typedef SourceProblem< GridSelector :: GridType,
        RangeTypeProvider< GridSelector::GridType::dimensionworld,false>
        >PhaseProblemType;
#endif



#else      
#error "No valid problem number specified"
#endif




#define PROBLEM_HAS_SOLUTION



#endif

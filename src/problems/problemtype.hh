#ifndef PHASEFIELD_PROBLEMTYPE_HH
#define PHASEFIELD_PROBLEMTYPE_HH
#include <dune/common/version.hh>
#include <dune/fem/io/parameter.hh>


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
#if GRIDDIM==2
#include "../mixedscheme/bubbleensemble2.hh"
#else
#include "../mixedscheme/bubbleensemble1d.hh"
#endif
#if MIXED
typedef BubbleEnsemble< GridSelector :: GridType,
        RangeTypeProvider< GridSelector::GridType::dimensionworld,true>
        >PhaseProblemType;
#else
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
#if SPLIT
typedef NvStTanhProblem< GridSelector :: GridType,
        RangeTypeProvider< GridSelector::GridType::dimensionworld,true>
        > NvStPhaseProblemType;
typedef AcTanhProblem< GridSelector :: GridType,
        RangeTypeProvider< GridSelector::GridType::dimensionworld,true>
        > AcPhaseProblemType;
#endif
#else
typedef TanhProblem< GridSelector :: GridType,
        RangeTypeProvider< GridSelector::GridType::dimensionworld,false>
        >PhaseProblemType;
#endif
#elif PROBLEM==5
#include "../mixedscheme/sharpproblem.hh"
#if MIXED
typedef SharpProblem< GridSelector :: GridType,
        RangeTypeProvider< GridSelector::GridType::dimensionworld,true>
        >PhaseProblemType;
#else
typedef SharpProblem< GridSelector :: GridType,
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

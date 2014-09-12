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
struct NSKRangeTypeProvider;

template< int dimension>
struct NSKRangeTypeProvider<dimension,true>
{
  enum{ rangeDim=2*dimension+2 };
};


template< int dimension>
struct NSKRangeTypeProvider<dimension,false>
{
  enum{ rangeDim=dimension+1};
};







///////////////////////////////////////
// AVAILABLE PROBLEMS
///////////////////////////////////////
#if PROBLEM==1
#include "tanhproblem.hh"
#if MIXED
typedef TanhBubbleProblem< GridSelector:: GridType,
        NSKRangeTypeProvider< GridSelector::GridType::dimensionworld , true >
        >PhaseProblemType;
#else
typedef TanhBubbleProblem< GridSelector:: GridType,
        NSKRangeTypeProvider< GridSelector::GridType::dimensionworld , false >
        >  PhaseProblemType;
#endif
#else      
#error "No valid problem number specified"
#endif




#define PROBLEM_HAS_SOLUTION



#endif

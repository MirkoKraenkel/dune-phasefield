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
#if DGSCHEME
#warning "DGSCHEME"
  enum{ rangeDim=2*dimension+4 };
#elif FEMSCHEME
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
#include "tanhproblembalanced.hh"
typedef TanhProblem< GridSelector:: GridType>  PhaseProblemType;
#elif PROBLEM==2
#include "perfectgas_problem.hh"
typedef PhaseProblem< GridSelector :: GridType >  PhaseProblemType;
#elif PROBLEM==3
#include "constantproblem.hh"
typedef ConstantProblem< GridSelector:: GridType> PhaseProblemType;
#elif PROBLEM==4
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
#elif PROBLEM==5
#include "bubbleproblem.hh"
typedef BubbleProblem< GridSelector:: GridType>  PhaseProblemType;
#elif PROBLEM==7
#include "../mixedscheme/heatproblem.hh"
typedef HeatProblem< GridSelector :: GridType,
        RangeTypeProvider< GridSelector::GridType::dimensionworld,true>
        >PhaseProblemType;
#else
#error "No valid problem number specified"
#endif




#define PROBLEM_HAS_SOLUTION



#endif

#ifndef PHASEFIELD_PROBLEMTYPE_HH
#define PHASEFIELD_PROBLEMTYPE_HH
#include <dune/common/version.hh>
#include <dune/fem/io/parameter.hh>

#if DUNE_VERSION_NEWER_REV(DUNE_GRID,2,1,0)
// DGF gridtype 
#else
#include <dune/grid/io/file/dgfparser/dgfgridtype.hh>
#endif



///////////////////////////////////////
// AVAILABLE PROBLEMS
///////////////////////////////////////
#if PROBLEM==1
#include "tanhproblem.hh"
typedef TanhProblem< GridSelector:: GridType>  PhaseProblemType;
 #else      
#error "No valid problem number specified"
#endif




#define PROBLEM_HAS_SOLUTION



#endif

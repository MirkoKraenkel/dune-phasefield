#ifndef FEMHOWTO_NSEQ_PROBLEMTYPE_HH
#define FEMHOWTO_NSEQ_PROBLEMTYPE_HH

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
#include "problem1d.hh"

typedef PhaseWaves< GridSelector :: GridType >  PhaseProblemType;

#define PROBLEM_HAS_SOLUTION



#endif

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

#if PROBLEM==1
   // a nonstationary exact solution 
   // for Navier-Stokes equations
#warning "Phaseproblen choosen":
#if GRIDDIM==2
#include "problem2d.hh"
#else
#include "problem1d.hh"
#endif
   typedef PhaseWaves< GridSelector :: GridType >  PhaseProblemType;

#define PROBLEM_HAS_SOLUTION

#else
   #error "No valid problem number specifiedxx"
#endif

#endif

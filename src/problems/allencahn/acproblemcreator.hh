#ifndef ACRPOBLEMCREATOR_HHC
#define ACRPOBLEMCREATOR_HH
#include <config.h>
#warning "ACCREATING"
// system includes
#include <string>

// dune includes
#include <dune/common/version.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>

#if DUNE_VERSION_NEWER_REV(DUNE_GRID,2,1,0)
#else
#include <dune/grid/io/file/dgfparser/dgfgridtype.hh>
#endif

// local includes
#include <dune/fem/fluxes/eulerfluxes.hh>
#include <dune/fem-dg/operator/fluxes/diffusionflux.hh>

#include <dune/fem/operator/discretemodelcommon.hh>

#include "problemtype.hh"

#include <dune/phasefield/modelling/allencahnmodel.hh>


template< class GridType > 
struct ProblemGenerator 
{
  typedef PhaseProblemType ProblemType;

  template< class GridPart >
  struct Traits
  {
    typedef ProblemType InitialDataType;

    typedef AllenCahnModel<GridPart,InitialDataType> ModelType;  
    // choice of diffusion flux (see diffusionflux.hh for methods)
    static const Dune :: DGDiffusionFluxIdentifier PrimalDiffusionFluxId 
           = Dune :: method_general ;


    typedef LLFFlux<ModelType> FluxType;
 };

  static inline std::string advectionFluxName()
  {
		return "Non";
 }

  static inline std::string diffusionFluxName()
  {
    return "LDG";
  }

  static inline Dune::GridPtr<GridType>               
  initializeGrid( const std::string description )
  {
    // use default implementation 
    return initialize< GridType > ( description );
  }

  static ProblemType* problem()
  {
    return new ProblemType ();
  }

};
//#include<dune/phasefield/phasefieldalgorithm_eq.hh>


#endif // PHASE_PROBLEMCREATOR_HH

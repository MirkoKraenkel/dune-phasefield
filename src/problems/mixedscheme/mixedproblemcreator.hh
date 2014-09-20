#ifndef PHASE_MIXERDPOBLEMCREATOR_HH
#define PHASE_MIXERDPOBLEMCREATOR_HH
#include <config.h>
#warning "MixedProblemCreator"
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



#include "../passscheme/problemtype.hh"

//#include <dune/fem/operator/assembled/heatmodel.hh>
#include <dune/fem/operator/assembled/models/phasefieldmodel_sources.hh>
#if IMPLICITTAU
#include <dune/fem/operator/assembled/fluxes/fluximpl.hh>
#else
#include <dune/fem/operator/assembled/fluxes/flux.hh>
#endif
template< class GridType > 
struct MixedProblemGenerator 
{
  typedef PhaseProblemType ProblemType;
  
  enum{ dimRange=ProblemType::dimRange};
  template< class GridPart >
  struct Traits
  {
    typedef PhaseProblemType  InitialDataType;

    typedef PhasefieldModel< typename GridPart::GridType, InitialDataType > ModelType;
    // choice of diffusion flux (see diffusionflux.hh for methods)

    // ******************************** NUMERICAL FLUX *****************************
	  typedef MixedFlux<ModelType > FluxType;
    // ****************************** END NUMERICAL FLUX ***************************
  };

  static inline std::string advectionFluxName()
  {
    return "mixedFlux\n";
  }

  static inline std::string diffusionFluxName()
  {
    return "IP";
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


#endif // PHASE_MIXERDPOBLEMCREATOR_HH


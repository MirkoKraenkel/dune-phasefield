#ifndef PHASE_MIXERDNSKPOBLEMCREATOR_HH
#define PHASE_MIXERDMSKPOBLEMCREATOR_HH
#include <config.h>
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



#include "../korteweg/problemtype.hh"

#include <dune/fem/operator/assembled/models/nskmodel.hh>
#if IMPLICITTAU
#include <dune/fem/operator/assembled/fluxes/nskfluximpl.hh>
#else
#include <dune/fem/operator/assembled/fluxes/nskflux.hh>
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

    typedef NSKModel< typename GridPart::GridType, InitialDataType > ModelType;

    // ******************************** NUMERICAL FLUX *****************************
	  typedef NSKMixedFlux<ModelType > FluxType;
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


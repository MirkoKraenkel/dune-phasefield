#ifndef PHASE_MIXERDPOBLEMCREATOR_HH
#define PHASE_MIXERDPOBLEMCREATOR_HH
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



#include "../problemtype.hh"

#include <dune/fem/operator/assembled/models/phasefieldmodel.hh>
#include <dune/fem/operator/assembled/splitop/fluxes/splitflux.hh>

template< class GridType > 
struct MixedProblemGenerator 
{
  typedef PhaseProblemType ProblemType;
  typedef NvStPhaseProblemType NvStProblemType;
  typedef AcPhaseProblemType AcProblemType;
  enum{ dimRange=ProblemType::dimRange};
  template< class GridPart >
  struct Traits
  {
    typedef PhaseProblemType InitialDataType;
    typedef NvStPhaseProblemType  NvStInitialDataType;
    typedef AcPhaseProblemType AcInitialDataType;
    typedef PhasefieldModel< typename GridPart::GridType,InitialDataType > ModelType;

    // ******************************** NUMERICAL FLUX *****************************
	  typedef NvStMixedFlux<ModelType > NvStFluxType;
    typedef AcMixedFlux<ModelType> AcFluxType;
    // ****************************** END NUMERICAL FLUX ***************************
  }
  
  ;

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
    return new ProblemType();
  }


  static NvStProblemType* nvstproblem()
  {
    return new NvStProblemType ();
  }

  static AcProblemType* acproblem()
  {
    return new AcProblemType ();
  }

};


#endif // PHASE_MIXERDPOBLEMCREATOR_HH


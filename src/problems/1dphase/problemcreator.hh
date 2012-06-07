#ifndef PHASE_RPOBLEMCREATOR_HH
#define PHASE_RPOBLEMCREATOR_HH
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

// local includes
#include <dune/fem/fluxes/eulerfluxes.hh>
#include <dune/fem-dg/operator/fluxes/diffusionflux.hh>
#if NONCON
#include <dune/fem/operator/discretemodelcommon.hh>
#else
#include <dune/fem/operator/projectiondiscretemodelcommon.hh>
#endif

#include "problemtype.hh"
#if WELLBALANCED
#include "wellbalancedmodel.hh"
#else
#include "model.hh"
#endif



template< class GridType > 
struct ProblemGenerator 
{
  typedef PhaseProblemType ProblemType;

  template< class GridPart >
  struct Traits
  {
    typedef ProblemType InitialDataType;
    typedef Dune::PhaseModel< GridPart, InitialDataType > ModelType;
    // choice of diffusion flux (see diffusionflux.hh for methods)
    static const Dune :: DGDiffusionFluxIdentifier PrimalDiffusionFluxId 
           = Dune :: method_general ;

// ******************************** NUMERICAL FLUX *****************************
#if WELLBALANCED
#warning "FLUX: WB"
		typedef WBFlux< ModelType > FluxType;
#else
#if (FLUX==1)
#warning "FLUX: LLF"
    typedef LLFFlux< ModelType > FluxType;
#elif (FLUX==2)

#else
#error "Set the flag FLUX! See Makefile.am for details!"
#endif
#endif    
// ****************************** END NUMERICAL FLUX ***************************
  };

  static inline std::string advectionFluxName()
  {
#if (FLUX==1)
    return "LLF";
#elif (FLUX==2)
		return "WB";
#endif
 }

  static inline std::string diffusionFluxName()
  {
#ifdef EULER
    return "";
#elif (defined PRIMALDG)
    return Dune::Parameter::getValue< std::string >("dgdiffusionflux.method");
#else
    return "LDG";
#endif
  }

  static inline Dune::GridPtr<GridType>               
  initializeGrid( const std::string description )
  {
    // use default implementation 
    return initialize< GridType > ( description );
  }

  static ProblemType* problem()
  {
    // choice of explicit or implicit ode solver
    return new ProblemType ();
  }

};
#if MYALGO
#include<dune/phasefield/phasefieldalgorithm.hh>
#else
#include <dune/fem/main/stepper.hh>
#include <dune/fem/main/steppertraits.hh>
#endif


#endif // FEMHOWTO_NSEQ_RPOBLEMCREATOR_HH

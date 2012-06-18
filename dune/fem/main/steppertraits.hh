#ifndef DUNE_STEPPERTRAITS_HH
#define DUNE_STEPPERTRAITS_HH

#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/space/dgspace/localrestrictprolong.hh>
#include <dune/fem/solver/odesolver.hh>

template <class GridImp,
          class ProblemTraits, 
          int polOrd>             
struct StepperTraits 
{
  enum { polynomialOrder = polOrd };

  // type of Grid
  typedef GridImp                                   GridType;

  // Choose a suitable GridView
  typedef Dune :: DGAdaptiveLeafGridPart< GridType >       GridPartType;

  // problem dependent types 
  typedef typename ProblemTraits :: template Traits< GridPartType > :: InitialDataType  InitialDataType;
  typedef typename ProblemTraits :: template Traits< GridPartType > :: ModelType        ModelType;
  typedef typename ProblemTraits :: template Traits< GridPartType > :: FluxType         FluxType;
  static const Dune::DGDiffusionFluxIdentifier DiffusionFluxId 
    = ProblemTraits :: template Traits< GridPartType > ::PrimalDiffusionFluxId ;

  // The DG Operator (using 2 Passes)
  // The first operator is sum of the other two
  // The other two are needed for semi-implicit time discretization

#warning "No limiter is applied to the numerical solution !!"
  typedef Dune :: DGAdvectionDiffusionOperator< ModelType, FluxType,
						DiffusionFluxId,  polynomialOrder >            DgType; /*@LST1E@*/
  typedef Dune :: DGAdvectionDiffusionOperator< ModelType, FluxType,
						DiffusionFluxId, polynomialOrder >      DgAdvectionType; /*@LST1E@*/
  typedef Dune :: DGDiffusionOperator< ModelType, FluxType,
                               DiffusionFluxId, polynomialOrder >      DgDiffusionType; /*@LST1E@*/

  // The discrete function for the unknown solution is defined in the DgOperator
  typedef typename DgType :: DestinationType                         DiscreteFunctionType;
  typedef typename DgType :: Destination1Type                        DiscreteGradientType;
  typedef typename DgType :: ScalarDFType                            ScalarDFType;

  // The indicator function in case of limiting 
  typedef typename DgAdvectionType :: IndicatorType                  IndicatorType;

  // ... as well as the Space type
  typedef typename DgType :: SpaceType                               DiscreteSpaceType;
  typedef typename DgType :: ScalarDiscreteFunctionSpaceType         ScalarDiscreteSpaceType;
  
  // The ODE Solvers                                                         /*@LST1S@*/
  typedef DuneODE :: OdeSolverInterface< DiscreteFunctionType > OdeSolverType ;

  // type of restriction/prolongation projection for adaptive simulations 
  typedef Dune :: RestrictProlongDefault< DiscreteFunctionType > RestrictionProlongationType;
};

#endif

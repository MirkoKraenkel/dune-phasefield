#ifndef ALGORITHM_TRAITS_HH
#define ALGORITHM_TRAITS_HH

#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/space/dgspace/localrestrictprolong.hh>
#include <dune/fem/solver/odesolver.hh>



template <class GridImp,
          class ProblemGeneratorImp, 
          int polOrd>             
struct AlgorithmTraitsTraits 
{
  enum { polynomialOrder = polOrd };

  // type of Grid
  typedef GridImp                                  GridType;
	typedef ProblemGeneratorImp                      ProblemGeneratorType;

  // Choose a suitable GridView
  typedef Dune :: DGAdaptiveLeafGridPart< GridType >       GridPartType;

  // problem dependent types 
  typedef typename ProblemGeneratorType :: template Traits< GridPartType > :: InitialDataType  InitialDataType;
  typedef typename ProblemGeneratorType :: template Traits< GridPartType > :: ModelType        ModelType;
  typedef typename ProblemGeneratorType :: template Traits< GridPartType > :: FluxType         FluxType;

	typedef Dune :: WellBalancedPhasefieldOperatorr< ModelType, FluxType,DiffusionFluxId,  polynomialOrder >         DiscreteOpType;
	
  typedef typename DiscreteOpType :: DestinationType                         DiscreteFunctionType;
  typedef typename DiscreteOpType :: Destination1Type                        DiscreteGradientType;
	typedef typename DiscreteOpType :: Destination2Type                        DiscreteThetaType;
	typedef typename DiscreteOpType :: ScalarDFType                            ScalarDFType;
	// ... as well as the Space type
  typedef typename DgType :: SpaceType                               DiscreteSpaceType;
  typedef typename DgType :: Space2Type                              DiscreteThetaSpaceType;

  typedef typename DgType :: ScalarDiscreteFunctionSpaceType         ScalarDiscreteSpaceType;
	

	// The ODE Solvers                                                         /*@LST1S@*/
  typedef DuneODE :: OdeSolverInterface< DiscreteFunctionType > OdeSolverType ;
  // type of restriction/prolongation projection for adaptive simulations 
  typedef Dune :: RestrictProlongDefault< DiscreteFunctionType > RestrictionProlongationType;
};




#endif

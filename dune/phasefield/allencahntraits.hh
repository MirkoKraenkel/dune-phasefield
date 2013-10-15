#ifndef ALGORITHM_TRAITS_HH
#define ALGORITHM_TRAITS_HH

#include <dune/fem/gridpart/adaptiveleafgridpart.hh>

#include <dune/fem/space/discontinuousgalerkin/localrestrictprolong.hh>

#include <dune/fem/solver/odesolver.hh>

#include <dune/fem/operator/allencahnoperator.hh>
template <class GridImp,
          class ProblemGeneratorImp, int polOrd>             
struct AlgorithmTraits 
{
  enum { polynomialOrder = polOrd };

  // type of Grid
  typedef GridImp                                  GridType;
	typedef ProblemGeneratorImp                      ProblemGeneratorType;

  // Choose a suitable GridView
  typedef Dune :: Fem::DGAdaptiveLeafGridPart< GridType >       GridPartType;
  typedef Dune :: Fem::AdaptiveLeafGridPart< GridType >         LagrangeGridPartType;


  // problem dependent types 
  typedef typename ProblemGeneratorType :: template Traits< GridPartType > :: InitialDataType  InitialDataType;
  typedef typename ProblemGeneratorType :: template Traits< GridPartType > :: ModelType        ModelType;
  typedef typename ProblemGeneratorType :: template Traits< GridPartType > :: FluxType         FluxType;
 
  
  typedef Dune::Fem::GridFunctionAdapter<Velocity,GridPartType>     VelocityType;	
	//typedef Dune :: WellBalancedPhasefieldOperator< ModelType, FluxType,DiffusionFluxId,  polynomialOrder >    DiscreteOperatorType;
	typedef Dune :: DGAllenCahnOperator< ModelType, FluxType, VelocityType,polynomialOrder >    DiscreteOperatorType;
	
	// Types of all Discretefuntions used in the simulation: Destinations of the passes, Scalarspace for energy calculation
  typedef typename DiscreteOperatorType :: DestinationType                         DiscreteFunctionType;
	typedef typename DiscreteOperatorType :: Destination1Type                        DiscreteSigmaType;
//	typedef typename DiscreteOperatorType :: DiscreteScalarType                      DiscreteScalarType;

	// ... as well as the Space type
	typedef typename  DiscreteOperatorType :: SpaceType                              DiscreteSpaceType;
  typedef typename  DiscreteOperatorType :: Space1Type                             SigmaDiscreteSpaceType;   
 // typedef typename  DiscreteOperatorType :: ScalarDiscreteFunctionSpaceType        ScalarDiscreteSpaceType;
	

	// The ODE Solvers                                                         /*@LST1S@*/
  typedef DuneODE :: OdeSolverInterface< DiscreteFunctionType > OdeSolverType ;
  // type of restriction/prolongation projection for adaptive simulations 
  typedef Dune :: Fem::RestrictProlongDefault< DiscreteFunctionType > RestrictionProlongationType;
};




#endif

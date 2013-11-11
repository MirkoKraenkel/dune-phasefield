#ifndef ALGORITHM_TRAITS_HH
#define ALGORITHM_TRAITS_HH


#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/space/discontinuousgalerkin.hh>
#include <dune/fem/space/discontinuousgalerkin/localrestrictprolong.hh>
#include <dune/fem/function/adaptivefunction.hh>

#include <dune/fem/operator/assembled/mixedoperator.hh>
template <class GridImp,
          class ProblemGeneratorImp, int polOrd>             
struct AlgorithmTraits 
{
  enum { polynomialOrder = polOrd };

 

  // type of Grid
  typedef GridImp                                  GridType;
	typedef typename GridType::ctype                 ctype;
  
  
  
  typedef ProblemGeneratorImp                      ProblemGeneratorType;

 // Choose a suitable GridView
  typedef Dune :: Fem::DGAdaptiveLeafGridPart< GridType >       GridPartType;
  typedef Dune :: Fem::AdaptiveLeafGridPart< GridType >         LagrangeGridPartType;
  
  enum{ dimDomain = GridType::dimensionworld };
  enum{ dimRange = 2*dimDomain+4 };

 
  // problem dependent types 
  typedef typename ProblemGeneratorType :: template Traits< GridPartType > :: InitialDataType  InitialDataType;
  typedef typename ProblemGeneratorType :: template Traits< GridPartType > :: ModelType        ModelType;
  typedef typename ProblemGeneratorType :: template Traits< GridPartType > :: FluxType         FluxType;
	
  // FunctionSpaces
  typedef typename Dune::Fem::FunctionSpace<ctype, double, dimDomain,dimRange> FunctionSpaceType;
  typedef typename Dune::Fem::DiscontinuousGalerkinSpace<FunctionSpaceType,GridPartType,polOrd,Dune::Fem::CachingStorage> DiscreteSpaceType;
  typedef typename Dune::Fem::FunctionSpace<ctype, double, dimDomain,1>  EnergySpaceType;
  typedef typename Dune::Fem::DiscontinuousGalerkinSpace<EnergySpaceType,GridPartType,polOrd,Dune::Fem::CachingStorage> DiscreteEnergySpaceType;


  // DiscreteFunctions
  typedef typename Dune::Fem::AdaptiveDiscreteFunction<DiscreteSpaceType> DiscreteFunctionType;


  typedef DGPhasefieldOperator<DiscreteFunctionType,ModelType,FluxType> DiscreteOperatorType;


	

  // type of restriction/prolongation projection for adaptive simulations 
  typedef Dune :: Fem::RestrictProlongDefault< DiscreteFunctionType > RestrictionProlongationType;
};




#endif

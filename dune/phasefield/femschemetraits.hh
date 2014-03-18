#ifndef ALGORITHM_TRAITS_HH
#define ALGORITHM_TRAITS_HH


#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/space/discontinuousgalerkin.hh>
#include <dune/fem/space/discontinuousgalerkin/localrestrictprolong.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/operator/linear/spoperator.hh>


//include solvers
#include <dune/fem/solver/oemsolver.hh>
#include <dune/fem/util/oemwrapper.hh>

#if HAVE_DUNE_ISTL
#include <dune/fem/function/blockvectorfunction.hh>
#include <dune/fem/operator/linear/istloperator.hh>
#include <dune/fem/solver/istlsolver.hh>
#endif

#if HAVE_UMFPACK
#include <dune/fem/solver/umfpacksolver.hh>
#endif

#if HAVE_PETSC
#include <dune/fem/function/petscdiscretefunction/petscdiscretefunction.hh>
#include <dune/fem/operator/linear/petscoperator.hh>
#include <dune/fem/solver/petscsolver.hh>
#endif




#include <dune/fem/operator/assembled/femoperator.hh>



template <class GridImp,
          class ProblemGeneratorImp,int polOrd>             
struct AlgorithmTraits 
{
  
  // type of Grid
  typedef GridImp                                             GridType;
	typedef typename GridType::ctype                            ctype;
  
  //Type for choosing the Problem 
  typedef ProblemGeneratorImp                                  ProblemGeneratorType;

 // Choose a suitable GridView
  typedef Dune :: Fem::AdaptiveLeafGridPart< GridType >        GridPartType;

  typedef GridPartType                                         LagrangeGridPartType;

  enum { polynomialOrder = polOrd };

  enum { dimDomain = GridType::dimensionworld };
  
  //( rho , v_1...v_n , phi , mu , tau ) 
  enum { dimRange = dimDomain+4 };


    // problem dependent types 
  typedef typename ProblemGeneratorType :: template Traits< GridPartType > :: InitialDataType  InitialDataType;
  typedef typename ProblemGeneratorType :: template Traits< GridPartType > :: ModelType        ModelType;
  typedef typename ProblemGeneratorType :: template Traits< GridPartType > :: FluxType         FluxType;
	
  // FunctionSpaces
  typedef typename Dune::Fem::FunctionSpace<ctype, double, dimDomain,dimRange> FunctionSpaceType;
  typedef typename Dune::Fem::LagrangeDiscreteFunctionSpace<FunctionSpaceType,GridPartType,polOrd,Dune::Fem::CachingStorage> DiscreteSpaceType;

  typedef typename Dune::Fem::FunctionSpace<ctype, double, dimDomain,1>  EnergySpaceType;
  typedef typename Dune::Fem::LagrangeDiscreteFunctionSpace<EnergySpaceType  ,GridPartType,polOrd,Dune::Fem::CachingStorage> DiscreteEnergySpaceType;


  // DiscreteFunctions
  typedef typename Dune::Fem::AdaptiveDiscreteFunction<DiscreteSpaceType>       DiscreteFunctionType;
  
  typedef typename Dune::Fem::AdaptiveDiscreteFunction<DiscreteEnergySpaceType> DiscreteScalarType;

  
  typedef OEMWrapper<DiscreteFunctionType> JacobianOperatorType;
//  typedef Dune::Fem::SparseRowLinearOperator< DiscreteFunctionType, DiscreteFunctionType > JacobianOperatorType; 

//  typedef FDJacobianDGPhasefieldOperator<DiscreteFunctionType,ModelType,FluxType,JacobianType> DiscreteOperatorType;
  typedef FDJacobianFemPhasefieldOperator< DiscreteFunctionType, ModelType, JacobianOperatorType > DiscreteOperatorType;
//  typedef LocalFDOperator<DiscreteFunctionType,ModelType,FluxType,JacobianOperatorType>   DiscreteOperatorType;


//  typedef Dune::Fem::UMFPACKOp<DiscreteFunctionType,JacobianOperatorType> LinearSolverType;
  typedef Dune::Fem::OEMGMRESOp< DiscreteFunctionType, JacobianOperatorType >  LinearSolverType;

  // type of restriction/prolongation projection for adaptive simulations 
  typedef Dune :: Fem::RestrictProlongDefault< DiscreteFunctionType > RestrictionProlongationType;
};




#endif

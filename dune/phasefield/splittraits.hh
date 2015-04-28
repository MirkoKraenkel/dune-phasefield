#ifndef SPLITSCHEME_TRAITS_HH
#define SPLITSCHEME_TRAITS_HH


#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/space/discontinuousgalerkin.hh>
#include <dune/fem/space/discontinuousgalerkin/localrestrictprolong.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/operator/linear/spoperator.hh>


//include solvers
#include <dune/fem/solver/oemsolver.hh>
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


#include <dune/fem/operator/assembled/splitop/navierstokesintegrator.hh>
#include <dune/fem/operator/assembled/splitop/allencahnintegrator.hh>
#include <dune/fem/operator/assembled/mixedoperator.hh>

#if FD
#include <dune/fem/operator/assembled/localfdoperator.hh>
#elif COUPLING
#include <dune/fem/operator/assembled/splitop/fluxes/jacobiansplitflux.hh>
#include <dune/fem/operator/assembled/splitop/allencahntensor.hh>
#include <dune/fem/operator/assembled/splitop/navierstokestensor.hh>
#include <dune/fem/operator/assembled/splitop/matrixoperator.hh>
#endif

//adaptation
//#include <dune/fem/adaptation/jumpestimator.hh>
//#include <dune/fem/adaptation/splitestimator.hh>
#include <dune/fem/adaptation/mixedestimator.hh>

template <class GridImp,
          class ProblemGeneratorImp,int polOrd>             
struct MixedAlgorithmTraits 
{
  enum { polynomialOrder = polOrd };

  // type of Grid
  typedef GridImp                                  GridType;
	typedef typename GridType::ctype                 ctype;
  
  //Type for choosing the Problem 
  typedef ProblemGeneratorImp                      ProblemGeneratorType;

 // Choose a suitable GridView
  typedef Dune :: Fem::DGAdaptiveLeafGridPart< GridType, All_Partition >       GridPartType;
  typedef Dune :: Fem::AdaptiveLeafGridPart< GridType >         LagrangeGridPartType;
  
  enum{ dimDomain = GridType::dimensionworld };
  
  //(rho,v_1...v_n,phi,mu,tau,sigma_1...sigma_n)
  enum{ dimRange=ProblemGeneratorType::dimRange/2};
  
  // problem dependent types 
  typedef typename ProblemGeneratorType :: template Traits< GridPartType > :: InitialDataType  InitialDataType;
  typedef typename ProblemGeneratorType :: template Traits< GridPartType > :: NvStInitialDataType NvStInitialDataType;
  typedef typename ProblemGeneratorType :: template Traits< GridPartType > :: AcInitialDataType AcInitialDataType;
  typedef typename ProblemGeneratorType :: template Traits< GridPartType > :: ModelType        ModelType;
  typedef typename ProblemGeneratorType :: template Traits< GridPartType > :: AcFluxType       AcFluxType;
	typedef typename ProblemGeneratorType :: template Traits< GridPartType > :: NvStFluxType     NvStFluxType;
	
  // FunctionSpaces
  typedef typename Dune::Fem::FunctionSpace<ctype, double, dimDomain,dimRange> FunctionSpaceType;
  typedef typename Dune::Fem::DiscontinuousGalerkinSpace<FunctionSpaceType,GridPartType,polOrd,Dune::Fem::CachingStorage> DiscreteSpaceType;

  typedef typename Dune::Fem::FunctionSpace<ctype, double, dimDomain,1>  ScalarSpaceType;
  typedef typename Dune::Fem::DiscontinuousGalerkinSpace<ScalarSpaceType  ,GridPartType,polOrd,Dune::Fem::CachingStorage> DiscreteScalarSpaceType;


  // DiscreteFunctions
#if WANT_ISTL
  typedef typename Dune::Fem::ISTLBlockVectorDiscreteFunction<DiscreteSpaceType> DiscreteFunctionType;
  typedef typename Dune::Fem::ISTLBlockVectorDiscreteFunction<DiscreteScalarSpaceType> DiscreteScalarType;
  typedef Dune::Fem::ISTLLinearOperator< DiscreteFunctionType, DiscreteFunctionType > JacobianOperatorType;
#if FD 
  typedef PhasefieldAllenCahnIntegrator< DiscreteFunctionType,DiscreteFunctionType,ModelType,AcFluxType> ACIntegratorType;
  typedef DGOperator< DiscreteFunctionType, ACIntegratorType> ACOperatorType;
  typedef PhasefieldNavierStokesIntegrator< DiscreteFunctionType,DiscreteFunctionType,ModelType,NvStFluxType> NvStIntegratorType;
  typedef DGOperator< DiscreteFunctionType, NvStIntegratorType> NvStOperatorType;
  typedef LocalFDOperator< ACOperatorType, JacobianOperatorType> DiscreteAllenCahnOperatorType;
  typedef LocalFDOperator< NvStOperatorType, JacobianOperatorType> DiscreteNavierStokesOperatorType;
#elif COUPLING
  typedef AllenCahnJacobianFlux<ModelType> AcJacFluxType; 
  typedef PhasefieldAllenCahnTensor< DiscreteFunctionType,DiscreteFunctionType,ModelType,AcFluxType,AcJacFluxType> ACIntegratorType;
  typedef DGOperator< DiscreteFunctionType, ACIntegratorType> ACOperatorType;
  typedef MatrixOperator< ACOperatorType, ACIntegratorType, JacobianOperatorType> DiscreteAllenCahnOperatorType;
  
  typedef NavierStokesJacobianFlux<ModelType> NvStkJacFluxType;
  typedef PhasefieldNavierStokesTensor< DiscreteFunctionType,DiscreteFunctionType,ModelType,NvStFluxType,NvStkJacFluxType>   NvStIntegratorType;
  typedef DGOperator< DiscreteFunctionType, NvStIntegratorType> NvStOperatorType;
  typedef MatrixOperator< NvStOperatorType, NvStIntegratorType, JacobianOperatorType> DiscreteNavierStokesOperatorType; 
#endif
#endif

#if BICG
  typedef typename Dune::Fem::ISTLBICGSTABOp< DiscreteFunctionType, JacobianOperatorType > LinearSolverType; 
#else
  typedef typename Dune::Fem::ISTLGMResOp< DiscreteFunctionType , JacobianOperatorType > LinearSolverType;
#endif


  


  //typedef JumpEstimator<DiscreteFunctionType,ModelType> EstimatorType;
	typedef MixedEstimator<DiscreteFunctionType,ModelType> EstimatorType;
  typedef Dune::Fem::LocalFunctionAdapter<EstimatorType> EstimatorDataType;
  //Pointers for (rho,v,phi,mu,tau,sigma) and (pressure,totalenergy)
  typedef Dune::tuple< DiscreteFunctionType*,DiscreteFunctionType*,EstimatorDataType*,DiscreteScalarType*, DiscreteScalarType*> IOTupleType; 

  // type of restriction/prolongation projection for adaptive simulations 
  typedef Dune :: Fem::RestrictProlongDefault< DiscreteFunctionType > RestrictionProlongationType;
};




#endif

#ifndef FEMSCHEME_TRAITS_HH
#define FEMSCHEME_TRAITS_HH


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

#if NSK
#include <dune/fem/operator/assembled/nskintegrator.hh>
#else
#include <dune/fem/operator/assembled/integrator.hh>
#endif

#include <dune/fem/operator/assembled/mixedoperator.hh>

#if MATRIXFREE
#include <dune/fem/util/oemwrapper.hh>
#elif FD 
#include <dune/fem/operator/assembled/localfdoperator.hh>
#elif COUPLING 
#if NSK
#include <dune/fem/operator/assembled/fluxes/nskjacobianflux.hh>
#include <dune/fem/operator/assembled/nsktensor.hh>
#else
#include <dune/fem/operator/assembled/fluxes/jacobianflux.hh>
#include <dune/fem/operator/assembled/phasefieldtensor.hh>
#endif

#include <dune/fem/operator/assembled/matrixoperator.hh>
#endif

//adaptation

#include <dune/fem/adaptation/rhosigmaestimator.hh>


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
  enum{ dimRange=ProblemGeneratorType::dimRange};
  // problem dependent types 
  typedef typename ProblemGeneratorType :: template Traits< GridPartType > :: InitialDataType  InitialDataType;
  typedef typename ProblemGeneratorType :: template Traits< GridPartType > :: ModelType        ModelType;
  typedef typename ProblemGeneratorType :: template Traits< GridPartType > :: FluxType         FluxType;
	
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
#if NSK
  using PhasefieldIntegratorType=NSKIntegrator< DiscreteFunctionType, ModelType, FluxType >;
#else
  using PhasefieldIntegratorType=PhasefieldMixedIntegrator< DiscreteFunctionType, ModelType, FluxType >;
#endif 
  using PhasefieldOperatorType=DGOperator< DiscreteFunctionType ,PhasefieldIntegratorType >;
  using DiscreteOperatorType= LocalFDOperator<PhasefieldOperatorType , JacobianOperatorType>;
#else
#if NSK
  using JacFluxType=NskJacobianFlux<ModelType>;
  using PhasefieldIntegratorType=NskMixedTensor<DiscreteFunctionType,ModelType,FluxType,JacFluxType>;
#else
  using JacFluxType=JacobianFlux<ModelType>;
  using PhasefieldIntegratorType=PhasefieldMixedTensor<DiscreteFunctionType,ModelType,FluxType,JacFluxType>;
#endif
  using PhasefieldOperatorType=DGOperator< DiscreteFunctionType,PhasefieldIntegratorType> ;
  using DiscreteOperatorType=MatrixOperator<PhasefieldOperatorType,PhasefieldIntegratorType,JacobianOperatorType>;
#endif


#if BICG
  typedef typename Dune::Fem::ISTLBICGSTABOp< DiscreteFunctionType, JacobianOperatorType > LinearSolverType; 
#else
  typedef typename Dune::Fem::ISTLGMResOp< DiscreteFunctionType , JacobianOperatorType > LinearSolverType;
#endif

#else  
  typedef typename Dune::Fem::AdaptiveDiscreteFunction<DiscreteSpaceType> DiscreteFunctionType;
  typedef typename Dune::Fem::AdaptiveDiscreteFunction<DiscreteScalarSpaceType> DiscreteScalarType;

  
#if MATRIXFREE 
  typedef OEMWrapper<DiscreteFunctionType> JacobianOperatorType;
  typedef FDJacobianDGPhasefieldOperator<DiscreteFunctionType,ModelType,FluxType,JacobianOperatorType> DiscreteOperatorType;
  typedef Dune::Fem::OEMGMRESOp<DiscreteFunctionType,JacobianOperatorType> LinearSolverType;
#else 
  typedef Dune::Fem::SparseRowLinearOperator< DiscreteFunctionType, DiscreteFunctionType> JacobianOperatorType; 
  #if NSK
typedef LocalFDOperator<DiscreteFunctionType,ModelType,FluxType,JacobianOperatorType>  DiscreteOperatorType;
  #else
typedef PhasefieldJacobianOperator<DiscreteFunctionType,ModelType,FluxType,JacobianOperatorType>  DiscreteOperatorType;
#endif
//typedef Dune::Fem::UMFPACKOp<DiscreteFunctionType,JacobianOperatorType> LinearSolverType;
  typedef Dune::Fem::OEMGMRESOp<DiscreteFunctionType,JacobianOperatorType> LinearSolverType;
#endif
#endif

  //typedef JumpEstimator<DiscreteFunctionType,ModelType> EstimatorType;
	typedef MixedEstimator<DiscreteFunctionType,ModelType> EstimatorType;
  typedef Dune::Fem::LocalFunctionAdapter<EstimatorType> EstimatorDataType;
  //Pointers for (rho,v,phi,mu,tau,sigma) and (pressure,totalenergy)
  //typedef Dune::tuple< DiscreteFunctionType*,EstimatorDataType*,DiscreteScalarType*, DiscreteScalarType*> IOTupleType; 
  typedef Dune::tuple< DiscreteFunctionType*,DiscreteFunctionType*,DiscreteScalarType*, DiscreteScalarType*> IOTupleType; 


  // type of restriction/prolongation projection for adaptive simulations 
  typedef Dune :: Fem::RestrictProlongDefault< DiscreteFunctionType > RestrictionProlongationType;
};




#endif

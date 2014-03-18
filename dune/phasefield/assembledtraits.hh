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





//#include <dune/fem/operator/assembled/mixedoperator.hh>
#if MATRIXFREE
#include <dune/fem/operator/assembled/newheatoperator.hh>
#else
#if VISIT
#include <dune/fem/operator/assembled/mixedopjacobian.hh>
#else
#include <dune/fem/operator/assembled/localfdoperator.hh>
#endif
#endif


template <class GridImp,
          class ProblemGeneratorImp,int polOrd>             
struct AlgorithmTraits 
{
  enum { polynomialOrder = polOrd };

  // type of Grid
  typedef GridImp                                  GridType;
	typedef typename GridType::ctype                 ctype;
  
  //Type for choosing the Problem 
  typedef ProblemGeneratorImp                      ProblemGeneratorType;

 // Choose a suitable GridView
  typedef Dune :: Fem::DGAdaptiveLeafGridPart< GridType >       GridPartType;
  typedef Dune :: Fem::AdaptiveLeafGridPart< GridType >         LagrangeGridPartType;
  
  enum{ dimDomain = GridType::dimensionworld };
  
  //(rho,v_1...v_n,phi,mu,tau,sigma_1...sigma_n)
#if RHOMODEL
  enum{ dimRange = 2*dimDomain+5 };
#else
  enum{ dimRange = 2*dimDomain+4 };
#endif
  // problem dependent types 
  typedef typename ProblemGeneratorType :: template Traits< GridPartType > :: InitialDataType  InitialDataType;
  typedef typename ProblemGeneratorType :: template Traits< GridPartType > :: ModelType        ModelType;
  typedef typename ProblemGeneratorType :: template Traits< GridPartType > :: FluxType         FluxType;
	
  // FunctionSpaces
  typedef typename Dune::Fem::FunctionSpace<ctype, double, dimDomain,dimRange> FunctionSpaceType;
  typedef typename Dune::Fem::DiscontinuousGalerkinSpace<FunctionSpaceType,GridPartType,polOrd,Dune::Fem::CachingStorage> DiscreteSpaceType;

  typedef typename Dune::Fem::FunctionSpace<ctype, double, dimDomain,1>  EnergySpaceType;
  typedef typename Dune::Fem::DiscontinuousGalerkinSpace<EnergySpaceType  ,GridPartType,polOrd,Dune::Fem::CachingStorage> DiscreteEnergySpaceType;


  // DiscreteFunctions
#if WANT_ISTL
  typedef typename Dune::Fem::ISTLBlockVectorDiscreteFunction<DiscreteSpaceType> DiscreteFunctionType;
  typedef typename Dune::Fem::ISTLBlockVectorDiscreteFunction<DiscreteEnergySpaceType> DiscreteScalarType;
  typedef Dune::Fem::ISTLLinearOperator< DiscreteFunctionType, DiscreteFunctionType > JacobianOperatorType;
#if VISIT
  typedef PhasefieldJacobianOperator<DiscreteFunctionType,ModelType,FluxType,JacobianOperatorType>  DiscreteOperatorType;
#else
  typedef LocalFDOperator<DiscreteFunctionType,ModelType,FluxType,JacobianOperatorType>  DiscreteOperatorType;
#endif
#if BICG
  typedef typename Dune::Fem::ISTLBICGSTABOp< DiscreteFunctionType, JacobianOperatorType > LinearSolverType; 
#else
 typedef typename Dune::Fem::ISTLGMResOp< DiscreteFunctionType , JacobianOperatorType > LinearSolverType;
#endif
#else  
  typedef typename Dune::Fem::AdaptiveDiscreteFunction<DiscreteSpaceType> DiscreteFunctionType;
  typedef typename Dune::Fem::AdaptiveDiscreteFunction<DiscreteEnergySpaceType> DiscreteScalarType;

  

#if MATRIXFREE 
  typedef OEMWrapper<DiscreteFunctionType> JacobianOperatorType;
  typedef FDJacobianDGPhasefieldOperator<DiscreteFunctionType,ModelType,FluxType,JacobianOperatorType> DiscreteOperatorType;
  typedef Dune::Fem::OEMGMRESOp<DiscreteFunctionType,JacobianOperatorType> LinearSolverType;
#else 
  typedef Dune::Fem::SparseRowLinearOperator< DiscreteFunctionType, DiscreteFunctionType> JacobianOperatorType; 
  typedef LocalFDOperator<DiscreteFunctionType,ModelType,FluxType,JacobianOperatorType>  DiscreteOperatorType;
  typedef Dune::Fem::UMFPACKOp<DiscreteFunctionType,JacobianOperatorType> LinearSolverType;
//  typedef Dune::Fem::OEMGMRESOp<DiscreteFunctionType,JacobianOperatorType> LinearSolverType;

#endif
#endif
  // type of restriction/prolongation projection for adaptive simulations 
  typedef Dune :: Fem::RestrictProlongDefault< DiscreteFunctionType > RestrictionProlongationType;
};




#endif

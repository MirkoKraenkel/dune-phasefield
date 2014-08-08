#ifndef ASSEMBLEDALGORITHM_HH
#define ASSEMBLEDALGORITHM_HH
#include "algorithmbase.hh"
//solver
#include <dune/fem/solver/cginverseoperator.hh>
#include <dune/fem/solver/oemsolver.hh>
#include <dune/fem/solver/newtoninverseoperator.hh>
//adaptation
#include <dune/fem/adaptation/mixedestimator.hh>
//post processing
#include <dune/fem/operator/assembled/energyconverter.hh>
#include <dune/phasefield/util/cons2prim.hh>
#include <dune/fem/operator/assembled/boundary.hh>

namespace Dune{


template< class GridImp,
          class AlgorithmTraits>
class AssembledAlgorithm: public PhasefieldAlgorithmBase<GridImp,AlgorithmTraits>
{
 public:
 typedef PhasefieldAlgorithmBase<GridImp,AlgorithmTraits> BaseType;
  typedef typename BaseType::Traits Traits;
  typedef typename BaseType::GridType GridType;
  typedef typename BaseType::DiscreteFunctionType DiscreteFunctionType;
  typedef typename BaseType::DiscreteScalarType DiscreteScalarType;
  typedef typename BaseType::ModelType ModelType;
  typedef typename BaseType::TimeProviderType TimeProviderType; 
  typedef typename BaseType::AdaptationManagerType AdaptationManagerType;
  //Operator
  typedef typename Traits :: DiscreteOperatorType DiscreteOperatorType;
  typedef typename Traits :: JacobianOperatorType JacobianOperatorType;
 	typedef PhasefieldBoundaryCorrection<DiscreteFunctionType, ModelType> BoundaryCorrectionType;
  //Solver
  typedef typename Traits :: LinearSolverType LinearSolverType;
  typedef typename Dune::Fem::NewtonInverseOperator< JacobianOperatorType, LinearSolverType > NewtonSolverType; 
  //adaptation
  typedef MixedEstimator<DiscreteFunctionType,ModelType> EstimatorType;
	typedef Dune::Fem::LocalFunctionAdapter<EstimatorType> EstimatorDataType;

  using BaseType::grid_;
  using BaseType::gridPart_;
  using BaseType::space_;
  using BaseType::solution_;
  using BaseType::model_;
  using BaseType::overallTimer_;
  using BaseType::tolerance_;
  private:
    DiscreteOperatorType    dgOperator_;
    BoundaryCorrectionType  boundaryCorrection_;
    DiscreteFunctionType    start_;
    EstimatorType           estimator_;
    EstimatorDataType       estimatorData_;

  public:
  //Constructor
  AssembledAlgorithm( GridType& gridImp):
      BaseType( gridImp),
      dgOperator_( *model_ , space()),
      boundaryCorrection_( *model_ , space()),
      start_( "start", space() ),
      estimator_( solution_, grid_, model()),
      estimatorData_( "estimator", estimator_, gridPart_, space_.order() )
      {
        start_.clear();
      }

    void step ( TimeProviderType& timeProvider,
		  				  int& newton_iterations,
			  			  int& ils_iterations,
				  		  int& max_newton_iterations,
					  	  int& max_ils_iterations)
	  {
      const double time=timeProvider.time();
      const double deltaT=timeProvider.deltaT();
    
      DiscreteFunctionType& U=solution();

      dgOperator_.setPreviousTimeStep(U);
      dgOperator_.setTime(time);
      dgOperator_.setDeltaT(deltaT);
     
      NewtonSolverType invOp(dgOperator_); 
      start_.clear();
   
      invOp(start_,U);
      //boundaryCorrection_(U);
      //reset overall timer
      overallTimer_.reset(); 
     
		  newton_iterations     = invOp.iterations();
      ils_iterations        = invOp.linearIterations();;
      max_newton_iterations = std::max(max_newton_iterations,newton_iterations);
      max_ils_iterations    = std::max(max_ils_iterations,ils_iterations);
	}
	  
  virtual void estimateMarkAdapt( AdaptationManagerType& am)
  {
    estimator_.estimateAndMark(tolerance_);
    am.adapt();
  }
  using BaseType::space;
  using BaseType::energyspace;
  using BaseType::gridSize;
  using BaseType::oldsolution;
  using BaseType::solution;
  using BaseType::energy;
  using BaseType::dataPrefix;
  using BaseType::problem;
  using BaseType::model;

  template<class Stream>
  void writeEnergy( TimeProviderType& timeProvider,
                    Stream& str)
  {
    DiscreteScalarType* totalenergy = energy();
    { 
      double kineticEnergy;
      double chemicalEnergy; 
      double energyIntegral =energyconverter(solution(),model(),*totalenergy,kineticEnergy,chemicalEnergy);
      str<<std::setprecision( 20 )<<timeProvider.time()<<"\t"<<energyIntegral<<"\t"<<chemicalEnergy<<"\t"<<kineticEnergy<<"\n";
    } 
  }
 




};
}//end namespace Dune



#endif

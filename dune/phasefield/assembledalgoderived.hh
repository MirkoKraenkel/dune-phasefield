#ifndef ASSEMBLEDALGORITHM_HH
#define ASSEMBLEDALGORITHM_HH
#include "assembledtraits.hh"
#include "algorithmbase.hh"
//solver
#include <dune/fem/solver/cginverseoperator.hh>
#include <dune/fem/solver/oemsolver.hh>
#include <dune/fem/solver/newtoninverseoperator.hh>
//post processing
#include <dune/fem/operator/assembled/energyconverter.hh>

#include <dune/phasefield/util/cons2prim.hh>
#include <dune/fem/operator/assembled/boundary.hh>

namespace Dune{


template< class GridImp,
          class AlgorithmTraits>
class AssembledAlgorithm: public PhasefieldAlgorithmBase< GridImp,AlgorithmTraits, AssembledAlgorithm< GridImp,AlgorithmTraits > > 
{
 public:
 typedef AssembledAlgorithm< GridImp, AlgorithmTraits> ThisType;
 typedef PhasefieldAlgorithmBase< GridImp , AlgorithmTraits , ThisType > BaseType;

  typedef typename BaseType::Traits Traits;
  typedef typename BaseType::GridType GridType;
  typedef typename BaseType::DiscreteFunctionType DiscreteFunctionType;
  typedef typename BaseType::DiscreteScalarType DiscreteScalarType;
  typedef typename BaseType::ModelType ModelType;
  typedef typename BaseType::TimeProviderType TimeProviderType; 
  typedef typename BaseType::AdaptationManagerType AdaptationManagerType;
  //Operator
  typedef typename Traits :: PhasefieldIntegratorType IntegratorType;
  typedef typename Traits :: PhasefieldOperatorType PhasefieldOperatorType;
  typedef typename Traits :: DiscreteOperatorType DiscreteOperatorType;
  typedef typename Traits :: JacobianOperatorType JacobianOperatorType;
 	typedef PhasefieldBoundaryCorrection<DiscreteFunctionType, ModelType> BoundaryCorrectionType;
  //Solver
  typedef typename Traits :: LinearSolverType LinearSolverType;
  typedef typename Dune::Fem::NewtonInverseOperator< JacobianOperatorType, LinearSolverType > NewtonSolverType; 
  //adaptation
  typedef typename Traits::EstimatorType EstimatorType;
  typedef typename Traits::EstimatorDataType EstimatorDataType;
  typedef typename Traits::IOTupleType IOTupleType;
  typedef typename BaseType::DataWriterType DataWriterType;

  using BaseType::grid_;
  using BaseType::gridPart_;
  using BaseType::space_;
  using BaseType::fixedTimeStep_;
  using BaseType::solution_;
  using BaseType::oldsolution_;
  using BaseType::energy_;
  using BaseType::model_;
  using BaseType::overallTimer_;
  using BaseType::tolerance_;
  private:
  IntegratorType          integrator_;
  PhasefieldOperatorType  pfOperator_;
  DiscreteOperatorType    dgOperator_;
  BoundaryCorrectionType  boundaryCorrection_;
  DiscreteFunctionType    start_;
  DiscreteScalarType      pressure_;
  EstimatorType           estimator_;
  EstimatorDataType       estimatorData_;
  double                  surfaceEnergy_;
  int                     maxNewtonIter_;
  int                     minNewtonIter_;
  int                     maxLinearIter_;
  double                  timestepfactor_;
  double                  maxTimeStep_;

  public:
  //Constructor
  AssembledAlgorithm( GridType& gridImp):
    BaseType( gridImp ),
    integrator_( *model_ ,space()),
    pfOperator_( integrator_, space()),
    dgOperator_( pfOperator_, space()),
    boundaryCorrection_( *model_ , space()),
    start_( "start", space() ),
    pressure_( "pressure", energyspace()),
    estimator_( solution_, grid_, model()),
    estimatorData_( "estimator", estimator_, gridPart_, space_.order() ),
    surfaceEnergy_(0.),
    maxNewtonIter_( Dune::Fem::Parameter::getValue<int>("phasefield.maxNewtonIterations") ),
    minNewtonIter_( Dune::Fem::Parameter::getValue<int>("phasefield.minNewtonIterations") ),
    maxLinearIter_( Dune::Fem::Parameter::getValue<int>("phasefield.maxLinearIterations") ),
    timestepfactor_( Dune::Fem::Parameter::getValue<double>("phasefield.timestepfactor")),
    maxTimeStep_( Dune::Fem::Parameter::getValue<double>("phasefield.maxTimeStep"))
    {
      start_.clear();
    }

//  IOTupleType getDataTuple (){ return IOTupleType( &solution(), &estimatorData_,&pressure_, this->energy());}
    IOTupleType getDataTuple (){ return IOTupleType( &solution(), &oldsolution(),&pressure_, this->energy());}


    DiscreteScalarType&  getpressure(){ return pressure_;}

    void initializeSolver ( typename BaseType::TimeProviderType& timeProvider){}

    template<class stream>
    void step ( TimeProviderType& timeProvider,
                int& newton_iterations,
                int& ils_iterations,
                int& max_newton_iterations,
                int& max_ils_iterations,
                stream& str)
	  {
      const double time=timeProvider.time();
      const double deltaT=timeProvider.deltaT();
      DiscreteFunctionType& U=solution();
      DiscreteFunctionType& Uold=oldsolution();
      Uold.assign(U);
      dgOperator_.setPreviousTimeStep(U);
      dgOperator_.setTime(time);
      dgOperator_.setDeltaT(deltaT);
      NewtonSolverType invOp(dgOperator_); 
      start_.clear();

      invOp(start_,U);
      
      //reset overall timer
      overallTimer_.reset();

      newton_iterations     = invOp.iterations();
      ils_iterations        = invOp.linearIterations();;
      max_newton_iterations = std::max(max_newton_iterations,newton_iterations);
      max_ils_iterations    = std::max(max_ils_iterations,ils_iterations);
          
      if( !invOp.converged() )
      {
        std::cout<<"no convergence!"<<"\n";
        std::cout<<"NewtonIterations: "<<newton_iterations<<" ( "<<ils_iterations<<" )  with deltaT="<<deltaT;
        str<<" aborted at T = "<<time<<"\n";
        str.close();
        abort();
        //timeProvider.invalidateTimeStep();
        //reduceTimeStep( timeProvider, deltaT);
        //U.assign(Uold);
      }


      if (fixedTimeStep_>1e-20)
      {
      }
      else
      {
        if( newton_iterations > maxNewtonIter_ || ils_iterations > maxLinearIter_ )
        {
          timestepfactor_*=0.5;
          std::cout<<"ReducedTimeStep= "<<timestepfactor_<<std::endl;
        }
        else if( newton_iterations < minNewtonIter_ && deltaT < maxTimeStep_ )
        {
          if(timestepfactor_*1.5<1)
            timestepfactor_*=1.5;
          std::cout<<"IncreaseTimeStep= "<<timestepfactor_<<std::endl;
        }
        double timeesti=timeStepEstimate();
        timeProvider.provideTimeStepEstimate( timeesti );
      }
    }
    
    void estimateMarkAdapt( AdaptationManagerType& am)
    {
      estimate();
      adapt(am);
    }
  
    void estimate()
    {
      estimator_.estimateAndMark(tolerance_);
    }

    void adapt( AdaptationManagerType& am)
    {
      am.adapt();
    }

    double timeStepEstimate ()
    {
      return std::min( maxTimeStep_, timestepfactor_*dgOperator_.timeStepEstimate());
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
    void writeEnergy ( TimeProviderType& timeProvider,
                       Stream& str,
                       int iter,
                       double dt)
    {

      DiscreteScalarType* totalenergy = energy();
      DiscreteScalarType&  pressure=getpressure();

      double kineticEnergy;
      double chemicalEnergy; 
      double surfaceEnergy;
      double energyIntegral =energyconverter(solution(),model(),*totalenergy,pressure,kineticEnergy,chemicalEnergy,surfaceEnergy);
      str<<std::setprecision( 20 )<<timeProvider.time()<<"\t"
        <<energyIntegral<<"\t"
        <<chemicalEnergy<<"\t"
        <<kineticEnergy<<"\t"
        <<surfaceEnergy<<"\t"
        <<iter<<"\t"
        <<dt<<"\n";
    }
 
    virtual void writeData( DataWriterType& eocDataOutput,
                            TimeProviderType& timeProvider,
                            const bool reallyWrite )
	  {
      DiscreteScalarType* totalenergy = energy();
      DiscreteScalarType&  pressure=getpressure();
      if( reallyWrite )
		  {

	      double kineticEnergy;
        double chemicalEnergy;
        double surfaceEnergy;
        energyconverter(solution(),model(),*totalenergy,pressure,kineticEnergy,chemicalEnergy,surfaceEnergy);
      }

      // write the data
      eocDataOutput.write( timeProvider );
	  }

};
}//end namespace Dune



#endif

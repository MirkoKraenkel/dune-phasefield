#ifndef SPLITALGORITHM_HH
#define SPLITALGORITHM_HH
#include "splittraits.hh"
#include "algorithmbase.hh"
//solver
#include <dune/fem/solver/cginverseoperator.hh>
#include <dune/fem/solver/oemsolver.hh>
#include <dune/fem/solver/newtoninverseoperator.hh>
//post processing
#include <dune/fem/operator/assembled/splitop/splitenergyconverter.hh>
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
  typedef typename BaseType::GridPartType GridPartType;
  //Operator
  typedef typename Traits :: ACIntegratorType ACIntegratorType;
  typedef typename Traits :: NvStIntegratorType NvStIntegratorType;
  typedef typename Traits :: ACOperatorType ACOperatorType;
  typedef typename Traits :: NvStOperatorType NvStOperatorType;
  typedef typename Traits :: DiscreteNavierStokesOperatorType DiscreteNvStkOperatorType;
  typedef typename Traits :: DiscreteAllenCahnOperatorType DiscreteACOperatorType;
  
  //Problem dependent type
  typedef typename Traits :: AcInitialDataType AcInitialDataType;
  typedef typename Traits :: NvStInitialDataType NvStInitialDataType;
  
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
    NvStIntegratorType nvstintegrator_;
    NvStOperatorType nvstOp_;
    DiscreteNvStkOperatorType dgNavStkOperator_;
    ACIntegratorType acintegrator_;
    ACOperatorType acOp_; 
    DiscreteACOperatorType  dgACOperator_; 
    BoundaryCorrectionType  boundaryCorrection_;
    DiscreteFunctionType    start_;
    DiscreteFunctionType    acsolution_;
    DiscreteFunctionType    acoldsolution_;
    const NvStInitialDataType*         nvstproblem_;
    const AcInitialDataType*           acproblem_;
    DiscreteScalarType      pressure_;
    EstimatorType           estimator_;
    EstimatorDataType       estimatorData_;
    double                  surfaceEnergy_;
    bool                    stepconverged_;
    int                     maxNewtonIter_;
    double                  timestepfactor_;
    double                  itertol_;
  public:
  //Constructor
  AssembledAlgorithm( GridType& gridImp):
      BaseType( gridImp),
      nvstintegrator_( *model_, space()),
      nvstOp_( nvstintegrator_,space()),
      dgNavStkOperator_( nvstOp_, space()),
      acintegrator_( *model_,space()),
      acOp_( acintegrator_,space()),
      dgACOperator_( acOp_, space()),
      boundaryCorrection_( *model_ , space()),
      start_( "start", space() ),
      acsolution_( "acsolution", space()),
      acoldsolution_( "acoldsolution", space()),
      nvstproblem_( Traits::ProblemGeneratorType::nvstproblem()),
      acproblem_( Traits::ProblemGeneratorType::acproblem()),
      pressure_( "pressure", energyspace()),
      estimator_( solution_, grid_, model()),
      estimatorData_( "estimator", estimator_, gridPart_, space_.order() ),
      surfaceEnergy_(0.),
      stepconverged_(false),
      maxNewtonIter_( Dune::Fem::Parameter::getValue<int>("phasefield.maxNewtonIterations") ),
      timestepfactor_( Dune::Fem::Parameter::getValue<double>("phasefield.timestepfactor") ),
      itertol_( Dune::Fem::Parameter::getValue<double>("phasefield.iterationtolerance"))
      {
        start_.clear();
      }

    IOTupleType getDataTuple (){ return IOTupleType( &solution(), &acsolution(),&estimatorData_,&pressure_, this->energy());}

    DiscreteScalarType&  getpressure(){ return pressure_;}

    DiscreteFunctionType& acsolution(){ return acsolution_;}
    DiscreteFunctionType& acoldsolution(){ return acoldsolution_;}
    void initializeSolver ( typename BaseType::TimeProviderType& timeProvider){}

    void reduceTimeStep( TimeProviderType& timeProvider, double ldt)
      {
        double factor=0.5;
        double newtime=std::min(factor*ldt,timeStepEstimate());
        timeProvider.provideTimeStepEstimate( newtime);
        //timeProvider.invalidateTimeStep();
      }
    
    void initializeStep( TimeProviderType& timeProvider)
    {

      DiscreteFunctionType& U=solution();
      DiscreteFunctionType& Uold=oldsolution();
      DiscreteFunctionType& acU=acsolution();
      DiscreteFunctionType& acUold=acoldsolution();
      Dune::Fem::DGL2ProjectionImpl::project( *nvstproblem_, U );
      Dune::Fem::DGL2ProjectionImpl::project( *acproblem_,acU );
      Uold.assign(U);
      acUold.assign(acU);
    }

    void step ( TimeProviderType& timeProvider,
		  				  int& newton_iterations,
			  			  int& ils_iterations,
				  		  int& max_newton_iterations,
					  	  int& max_ils_iterations)
	  {
      const double time=timeProvider.time();
      const double deltaT=timeProvider.deltaT();
      Fem::L2Norm<GridPartType> l2norm(gridPart_);      
      
      DiscreteFunctionType& U=solution();
      DiscreteFunctionType& UAc=acsolution();
      DiscreteFunctionType& Uold=oldsolution();
      DiscreteFunctionType& UAcOld=acoldsolution();
      dgNavStkOperator_.setPreviousTimeStep(U);
      dgNavStkOperator_.setTime(time);
      dgNavStkOperator_.setDeltaT(deltaT);
      
      dgACOperator_.setPreviousTimeStep(UAc);
      dgACOperator_.setTime(time);
      dgACOperator_.setDeltaT(deltaT);
      
      NewtonSolverType invOp(dgNavStkOperator_); 
      NewtonSolverType invOp2( dgACOperator_);  
      double errorNvSt,errorAc;    
      int counter=0; 
      do 
        {
          dgNavStkOperator_.integrator().setAddVariables(UAc);
      
          start_.clear();
          invOp(start_,U);
          errorNvSt=l2norm.distance(U, Uold);
          dgACOperator_.integrator().setAddVariables(U);
          start_.clear();
          invOp2(start_,UAc);
          errorAc=l2norm.distance(UAc,UAcOld);
          Uold.assign(U);
          UAcOld.assign(UAc);
          ++counter;
      }
      while( errorAc>itertol_ || errorNvSt> itertol_);
      //std::cout<<"counter="<<counter<<"\n";
      //reset overall timer
      overallTimer_.reset();

      newton_iterations     = counter;//invOp.iterations();
      ils_iterations        = invOp.linearIterations();;
      max_newton_iterations = std::max(max_newton_iterations,newton_iterations);
      max_ils_iterations    = std::max(max_ils_iterations,ils_iterations);
          
      if( !invOp.converged() )
        {
           std::cout<<"no convergence!"<<"\n";
           std::cout<<"NewtonIterations: "<<newton_iterations<<" ( "<<ils_iterations<<" )  with deltaT="<<deltaT;
           abort(); 
        }
     if (fixedTimeStep_>1e-20)
      {
     
      }
     else if( newton_iterations >  3 )
        {
          std::cout<<"NewtonIterations: "<<newton_iterations<<" ( "<<ils_iterations<<" )  with deltaT="<<deltaT;
          reduceTimeStep( timeProvider , deltaT );
         }
     else if( newton_iterations < 2 && timeProvider.deltaT() < timeStepEstimate())
        {
          timeProvider.provideTimeStepEstimate( 1.2*deltaT);
        }
      else
        {
          timeProvider.provideTimeStepEstimate( timeStepEstimate());
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

  double timeStepEstimate()
  {
    return timestepfactor_*std::min( dgNavStkOperator_.timeStepEstimate(),dgACOperator_.timeStepEstimate());
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
                    Stream& str,
                    int iter,
                    double dt)
  {
    DiscreteScalarType* totalenergy = energy();
    DiscreteScalarType&  pressure=getpressure();
    { 
      double kineticEnergy;
      double chemicalEnergy; 
      double surfaceEnergy;
      double energyIntegral =splitenergyconverter(solution(),acsolution(),model(),*totalenergy,pressure,kineticEnergy,chemicalEnergy,surfaceEnergy);
      str<<std::setprecision( 20 )<<timeProvider.time()<<"\t"
        <<energyIntegral<<"\t"
        <<chemicalEnergy<<"\t"
        <<kineticEnergy<<"\t"
        <<surfaceEnergy<<"\t"
        <<iter<<"\t"
        <<dt<<"\n";
    } 
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
      double energyIntegral = splitenergyconverter(solution(),acsolution(),model(),*totalenergy,pressure,kineticEnergy,chemicalEnergy,surfaceEnergy);
      std::cout<<"Energy="<<energyIntegral<<"\n";   
   }
	 // write the data
   eocDataOutput.write( timeProvider );
	}
};
}//end namespace Dune



#endif

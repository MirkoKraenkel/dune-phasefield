#ifndef DUNE_PHASEFIELD_ODESOLVER_HH
#define DUNE_PHASEFIELD_ODESOLVER_HH

#include <limits>
#include <dune/fem/solver/odesolver.hh>

namespace Dune {


struct PhaseFieldOdeSolverParameters : public DuneODE :: ODEParameters 
{
	
  PhaseFieldOdeSolverParameters* clone() const 
  {
    return new PhaseFieldOdeSolverParameters( *this );
  }
};

template <class Operator> 
class PhaseFieldOdeSolver : 
  public DuneODE :: OdeSolverInterface< typename Operator :: DestinationType >
{
  typedef Operator           OperatorType;

public:
  typedef typename OperatorType :: DestinationType DestinationType ;
  typedef DestinationType  DiscreteFunctionType;

  // The ODE Solvers                                                  
  typedef DuneODE :: OdeSolverInterface< DestinationType >           OdeSolverInterfaceType;
  typedef DuneODE :: ImplicitOdeSolver< DiscreteFunctionType >       ImplicitOdeSolverType; 
  typedef DuneODE :: ExplicitOdeSolver< DiscreteFunctionType >       ExplicitOdeSolverType; 

  typedef typename OdeSolverInterfaceType :: MonitorType MonitorType;

protected:
  OperatorType&          operator_;
  Fem::TimeProviderBase&      timeProvider_; 
  const PhaseFieldOdeSolverParameters* param_;
	
	OdeSolverInterfaceType* odeSolver_;
	
  const int verbose_ ;
  const int rkSteps_ ; 
  const int odeSolverType_;
  int imexCounter_ , exCounter_;
  int minIterationSteps_, maxIterationSteps_ ;
  bool imex_ ;
public:
  PhaseFieldOdeSolver( Fem::TimeProviderBase& tp,
											 OperatorType& op, 
											 const PhaseFieldOdeSolverParameters& parameter= PhaseFieldOdeSolverParameters() )
   : operator_( op ), 
     timeProvider_( tp ),
     param_( parameter.clone() ),
     odeSolver_( 0 ),
     verbose_( param_->verbose() ),
     rkSteps_( obtainRungeKuttaSteps() ),
     odeSolverType_( obtainOdeSolverType() ),
     imexCounter_( 0 ), exCounter_ ( 0 ),
     minIterationSteps_( std::numeric_limits< int > :: max() ),
     maxIterationSteps_( 0 ),
     imex_( odeSolverType_ > 1 )
  {
    // create implicit or explicit ode solver
    if( odeSolverType_ == 0 )
    {
      odeSolver_ = new ExplicitOdeSolverType( operator_, tp, rkSteps_);
    }
    else if (odeSolverType_ == 1)
    {
      odeSolver_ = new ImplicitOdeSolverType( operator_, tp, rkSteps_);
    }
		else 
    {
      DUNE_THROW(NotImplemented,"Wrong ODE solver selected");
    }
  }

  //! destructor 
  ~PhaseFieldOdeSolver() 
  {
    delete param_;           param_ = 0;
    delete odeSolver_;       odeSolver_ = 0;
  
  }

  //! initialize method 
  void initialize( const DestinationType& U )
  {
    assert( odeSolver_ );
    odeSolver_->initialize( U );

  }

  //! solver the ODE 
  void solve( DestinationType& U , 
              MonitorType& monitor ) 
  {
    // take CPU time of solution process 
    Timer timer ;
		
    operator_.switchupwind();
  
    // reset compute time counter 
    resetComputeTime();

    double maxAdvStep  = 0;
    double maxDiffStep = 0;

   
    {
      assert( odeSolver_ );
      odeSolver_->solve( U, monitor );

     
      const int iterationSteps = monitor.newtonIterations_ * monitor.linearSolverIterations_ ;
      minIterationSteps_ = std::min( minIterationSteps_, iterationSteps );
      maxIterationSteps_ = std::max( maxIterationSteps_, iterationSteps );
    }

   
		
    // store needed time 
    monitor.odeSolveTime_     = timer.elapsed();
    monitor.operatorTime_     = operatorTime();
    monitor.numberOfElements_ = numberOfElements();
  }

  //! return CPU time needed for the operator evaluation
  double operatorTime() const 
  {
   return operator_.computeTime();
  }

  //! return number of elements meat during operator evaluation 
  size_t numberOfElements() const  
  {
		return operator_.numberOfElements();
  }

  void description(std::ostream&) const {}

  // gather information from the space operator, the time integratior
  // and the problem to output before each table in tex file
  std::string description() const
  {
    std::string latexInfo;

    if ((odeSolverType_==0) || (odeSolverType_==1))
      latexInfo = operator_.description();

    std::stringstream odeInfo; 
    odeSolver_->description( odeInfo );

    latexInfo += odeInfo.str() + "\n";
    return latexInfo;
  }

protected:
  int obtainOdeSolverType () const 
  {
    // we need this choice of explicit or implicit ode solver
    // defined here, so that it can be used later in two different
    // methods
    static const std::string odeSolver[]  = { "EX", "IM" };
    std::string key( "fem.ode.odesolver" );
    if( Fem::Parameter :: exists( key ) )
      return Fem::Parameter::getEnum( "fem.ode.odesolver", odeSolver, 0 );
    else 
    {
      std::cerr << "WARNING: deprecated key, use `fem.ode.odesolver' instread!" << std::endl;
      return Fem::Parameter::getEnum( "femhowto.odesolver", odeSolver, 0 );
    }
  }

  int obtainRungeKuttaSteps() const
  {
    std::string key("fem.ode.order");
    if ( Fem::Parameter :: exists( key ) )
      return Fem::Parameter::getValue< int > ( key );
    else
      return operator_.space().order() + 1;
  }

  void resetComputeTime() const 
  {
    // this will reset the internal time counters 
    operator_.computeTime() ;
	}
}; // end PhaseFieldOdeSolver 

} // end namespace Dune 
#endif

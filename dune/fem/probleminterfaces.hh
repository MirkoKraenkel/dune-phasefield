#ifndef DUNE_PROBLEMINTERFACE_HH
#define DUNE_PROBLEMINTERFACE_HH

#include <dune/common/version.hh>
#include <dune/fem/function/common/timedependentfunction.hh>
#include <dune/fem/misc/gridsolution.hh>

namespace Dune {

/**
 * @brief describes the interface for
 * initial and exact solution of the advection-diffusion model
 */
template< class FunctionSpaceImp, bool constantVelocity >
class EvolutionProblemInterface
{
  typedef EvolutionProblemInterface< FunctionSpaceImp,
                                     constantVelocity >              ThisType;

public:
  typedef FunctionSpaceImp                                           FunctionSpaceType;

  enum { ConstantVelocity = constantVelocity };
  enum { dimDomain = FunctionSpaceType :: dimDomain };
  enum { dimRange  = FunctionSpaceType :: dimRange  };

  typedef typename FunctionSpaceType :: DomainType                   DomainType;
  typedef typename FunctionSpaceType :: RangeType                    RangeType;
  typedef typename FunctionSpaceType :: DomainFieldType              DomainFieldType;
  typedef typename FunctionSpaceType :: RangeFieldType               RangeFieldType;
  typedef typename FunctionSpaceType :: JacobianRangeType            JacobianRangeType;

  typedef Fem :: Parameter ParameterType ;

  /**
   * @brief define problem parameters
   */
protected:
  EvolutionProblemInterface()
  : writeGridSolution_( ParameterType::getValue<bool>("gridsol.writesolution", false) ),
    saveStep_( ParameterType :: getValue< double >("gridsol.firstwrite") ),
    saveInterval_( ParameterType :: getValue< double >("gridsol.savestep") ),
    writeCounter_( 0 )
  {
  }

public:
  typedef Fem :: TimeDependentFunction< ThisType > TimeDependentFunctionType;

  //! turn timedependent function into function by fixing time 
  TimeDependentFunctionType fixedTimeFunction( const double time ) const 
  {
    return TimeDependentFunctionType( *this, time );
  }

  // return prefix for data loops 
  virtual std::string dataPrefix() const 
  {
    return std::string( "loop" );
  }

  //! destructor
  virtual ~EvolutionProblemInterface() {}

  //! return true is source term is available 
  virtual inline bool hasStiffSource() const { return true ; }
  virtual inline bool hasNonStiffSource() const { return false; }

  //! stiff source term
  virtual inline double stiffSource( const DomainType& arg, 
                              const double time, 
                              const RangeType& u, 
                              RangeType& res) const 
  {
    abort();
    res = 0;
    return 0.0;
  }

  //! non stiff source term
  virtual inline double nonStiffSource( const DomainType& arg, 
                                 const double time, 
                                 const RangeType& u, 
                                 RangeType& res) const 
  {
    abort();
    res = 0;
    return 0.0;
  }

  /** \brief decide if the refinement/coarsening should be allowed in certain regions of the
   *    computational domain
   *  \param[in] x Global coordinates
   *  \return whether or not the refinement is allowed
   */
  bool allowsRefinement( const DomainType& x ) const
  {
    // by default refinement is allowed in the whole computational domain
    return true;
  }

  //! use both indicator for a grid adaptation
  inline bool twoIndicators() const
  {
    return false;
  }

  /** \brief return a first indicator for an adaptation of the grid
   *  \param[in] u Value of the numerical solution in a point
   *
   *  \note NaN does no adaptation
   *
   *  \return scalar indicator of importance for a grid adaptation
   *    like pot. temperature, density etc
   */
  inline double indicator1( const DomainType& x, const RangeType& u ) const
  {
    return std::numeric_limits<float>::quiet_NaN();
  }

  /** \brief return a second indicator for an adaptation of the grid
   *  \param[in] u Value of the numerical solution in a point
   *
   *  \note NaN does no adaptation
   *
   *  \return scalar indicator of importance for a grid adaptation
   *    like pot. temperature, density etc
   */
  inline double indicator2( const DomainType& x, const RangeType& u ) const
  {
    return std::numeric_limits<float>::quiet_NaN();
  }


  //! return diffusion coefficient (default returns epsilon)
  virtual inline double diffusion ( const RangeType& u, const JacobianRangeType& gradU ) const 
  {
    return epsilon();
  }

  /** \brief start time of problem */
  virtual double startTime() const { return 0.0; }

  /** \brief diffusion coefficient of problem */
  virtual double epsilon() const { return 0.0; }

  /**
   * @brief getter for the velocity
   */
  virtual void velocity(const DomainType& x, DomainType& v) const {}

  /**
   * @brief old version of the exact solution
   *
   * old version of evaluate(const DomainType& arg, double t, RangeType& res),
   * which is still needed by the DataWriter
   */
  virtual inline void evaluate(const double t,
                               const DomainType& arg, RangeType& res) const 
  {
    evaluate(arg, t, res);
  }

  /**
   * @brief evaluate exact solution, to be implemented in derived classes
   */
  virtual void evaluate(const DomainType& arg,
                        const double t, RangeType& res) const = 0 ;

  /**
   * @brief latex output for EocOutput, default is empty
   */
  virtual std::string description() const
  {
    return std::string("");
  }

  /*  \brief calculate EOC between exact and numerical solution 
   *  
   *  \param[in] tp     time provider 
   *  \param[in] u      numerical solution 
   *  \param[in] eocId  for FemEoc
   * 
   *  \return true if no EOC was calculated 
   */
  template< class TimeProviderType, class DiscreteFunctionType >
  bool calculateEOC( TimeProviderType& tp, DiscreteFunctionType& u,
                     const unsigned int eocId ) const
  {
    // return true means that EOC is calculated in Stepper
    return true ;
  }

  /* \brief 
   */
  template< class GridType, class DiscreteFunctionType >
  inline void postProcessTimeStep( const GridType& grid, 
                                   const DiscreteFunctionType& solution, 
                                   const double time ) const 
  {
    if( writeGridSolution_ && time > saveStep_ )
    {
      typedef Dune :: Fem :: GridSolutionVector< GridType, DiscreteFunctionType > ExSolGrid;
      ExSolGrid :: writeDiscreteFunction( grid, solution, time, writeCounter_ );
      ++writeCounter_;
      saveStep_ += saveInterval_;
    }
  }


  /*  \brief finalize the simulation using the calculated numerical
   *  solution u for this problem
   *  
   *  \param[in] variablesToOutput Numerical solution in the suitably chosen variables
   *  \param[in] eocloop Specific EOC loop
   */
  template< class DiscreteFunctionType >
  void finalizeSimulation( DiscreteFunctionType& variablesToOutput,
                           const int eocloop) const
  {}


protected:
  const bool writeGridSolution_;
  mutable double saveStep_ ;
  const double saveInterval_ ;
  mutable int writeCounter_ ;
};



// ProblemInterface
//-----------------

template <class FunctionSpaceImp>
class ProblemInterface
{
public:
  typedef FunctionSpaceImp                                           FunctionSpaceType;
  typedef ProblemInterface< FunctionSpaceType >                      ThisType;

  enum { dimDomain = FunctionSpaceType :: dimDomain };

  typedef typename FunctionSpaceType :: DomainType                   DomainType;
  typedef typename FunctionSpaceType :: RangeType                    RangeType;
  typedef typename FunctionSpaceType :: JacobianRangeType            JacobianRangeType;
  typedef typename FunctionSpaceType :: DomainFieldType              DomainFieldType;
  typedef typename FunctionSpaceType :: RangeFieldType               RangeFieldType;

  typedef FieldMatrix< RangeFieldType, dimDomain, dimDomain >        DiffusionMatrixType;

public:

  //! destructor
  virtual ~ProblemInterface() {}

  //! evaluates the right hand side (Function like bahavior)
  inline void evaluate(const DomainType& x, RangeType& ret) const
  {
    u(x , ret);
  }

  //! the right hand side
  virtual void f(const DomainType& x, RangeType& ret) const = 0;

  //! the exact solution
  virtual void u(const DomainType& x, RangeType& ret) const = 0;

  //! the diffusion matrix
  virtual void K(const DomainType& x, DiffusionMatrixType& m) const = 0;

  //! returns true if diffusion coefficient is constant
  virtual bool constantK() const = 0;

  //! the Dirichlet boundary data function
  virtual void g(const DomainType& x, RangeType& ret) const
  {
    u( x, ret );
  }

  //! the Neumann boundary data function
  virtual void psi(const DomainType& x,
                   JacobianRangeType& dn) const
  {
    abort();
    // gradient( x, dn );
  }

  /**
   * @brief getter for the velocity
   */
  virtual void velocity(const DomainType& x, DomainType& v) const 
  {
    v = 0;  
  }

  //! the gradient of the exact solution
  //virtual void gradient(const DomainType& x,
  //                      JacobianRangeType& grad) const = 0;

  //! return whether boundary is Dirichlet (true) or Neumann (false)
  virtual bool dirichletBoundary(const int bndId, const DomainType& x) const
  {
    return true;
  }

  // return prefix for data loops 
  virtual std::string dataPrefix() const 
  {
    return std::string( "loop" );
  }

  /**
   * @brief latex output for EocOutput, default is empty
   */
  virtual std::string description() const
  {
    return std::string("");
  }

protected:
  ProblemInterface() : exactSolution_( *this ) {}

  //! the exact solution to the problem for EOC calculation
  class ExactSolution
  : public Fem:: Function< FunctionSpaceType, ExactSolution >
  {
  private:
    typedef Fem:: Function< FunctionSpaceType, ExactSolution >      BaseType;

    typedef ProblemInterface< FunctionSpaceType>   DataType;
  protected:
    FunctionSpaceType  functionSpace_;
    const DataType &data_;

  public:
    inline ExactSolution ( const ThisType& data )
    : BaseType( ),
      functionSpace_(),
      data_( data )
    {
    }

    inline void evaluate ( const DomainType &x, RangeType &ret ) const
    {
      data_.u( x, ret );
    }

    inline void jacobian ( const DomainType &x, JacobianRangeType &ret ) const
    {
      data_.gradient( x, ret );
    }

    inline void evaluate (const DomainType &x,
                          const double time, RangeType &phi ) const
    {
      evaluate( x, phi );
    }
  }; // end class ExactSolution

public:
  //! type of function converter for exact solution and gradient
  typedef ExactSolution ExactSolutionType;

protected:
  ExactSolutionType exactSolution_;

public:
  const ExactSolutionType& exactSolution() const { return exactSolution_; }
};

}
#endif  /*DUNE_PROBLEMINTERFACE_HH*/

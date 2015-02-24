// vim: set expandtab ts=2 sw=2 sts=2:
#ifndef DUNE_FEM_NEWTONINVERSEOPERATOR_HH
#define DUNE_FEM_NEWTONINVERSEOPERATOR_HH
#define TIMING 0
#include <cfloat>
#include <iostream>

#include <dune/fem/io/parameter.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/common/differentiableoperator.hh>
#include <dune/fem/misc/femtimer.hh>
namespace Dune
{

  namespace Fem
  {

    // NewtonParameter
    // ---------------

    struct NewtonParameter
#ifndef DOXYGEN
    : public LocalParameter< NewtonParameter, NewtonParameter >
#endif
    {
      NewtonParameter () {}

      virtual double toleranceParameter () const
      {
        return Parameter::getValue< double >( "fem.solver.newton.tolerance", 1e-6 );
      }

      virtual double linAbsTolParameter ( const double &tolerance )  const
      {
        return Parameter::getValue< double >( "fem.solver.newton.linabstol", tolerance / 8 );
      }

      virtual double linReductionParameter ( const double &tolerance ) const
      {
        return Parameter::getValue< double >( "fem.solver.newton.linreduction", tolerance / 8 );
      }

      virtual bool verbose () const
      {
        const bool v = Parameter::getValue< bool >( "fem.solver.verbose", false );
        return Parameter::getValue< bool >( "fem.solver.newton.verbose", v );
      }

      virtual bool linearSolverVerbose () const
      {
        const bool v = Parameter::getValue< bool >( "fem.solver.verbose", false );
        return Parameter::getValue< bool >( "fem.solver.newton.linear.verbose", v );
      }

      virtual bool matrixout() const
      {
        return Parameter::getValue< bool >( "fem.solver.matrixout", false );
      }
      virtual int maxIterationsParameter () const
      {
        return Parameter::getValue< int >( "fem.solver.newton.maxiterations", std::numeric_limits< int >::max() );
      }

      virtual int maxLinearIterationsParameter () const
      {
        return Parameter::getValue< int >( "fem.solver.newton.maxlineariterations", std::numeric_limits< int >::max() );
      }
    };

    class NewtonControlDefault
    {
      public:
       NewtonControlDefault()
       :assemblesteps_( Parameter::getValue< int >( "fem.solver.newton.assemblesteps",1)) 
       {}


       bool calcJacobian(int iterations, double  reduction) const
       {
         return (iterations%assemblesteps_==0);
       }

      private:
      int assemblesteps_;
    };
  
    // NewtonInverseOperator
    // ---------------------

    /** \class NewtonInverseOperator
     *  \brief inverse operator based on a newton scheme
     *
     *  \tparam  Op      operator to invert (must be a DifferentiableOperator)
     *  \tparam  LInvOp  linear inverse operator
     *
     *  \note Verbosity of the NewtonInverseOperator is controlled via the
     *        paramter <b>fem.solver.newton.verbose</b>; it defaults to
     *        <b>fem.solver.verbose</b>.
     */
    template< class JacobianOperator, class LInvOp, class Control=NewtonControlDefault >
    class NewtonInverseOperator
    : public Operator< typename JacobianOperator::RangeFunctionType, typename JacobianOperator::DomainFunctionType >
    {
      typedef NewtonInverseOperator< JacobianOperator, LInvOp, Control > ThisType;
      typedef Operator< typename JacobianOperator::RangeFunctionType, typename JacobianOperator::DomainFunctionType > BaseType;
      typedef Control ControlType;

    public:
      //! type of operator's Jacobian
      typedef JacobianOperator JacobianOperatorType;

      //! type of operator to invert
      typedef DifferentiableOperator< JacobianOperatorType > OperatorType;

      //! type of linear inverse operator
      typedef LInvOp LinearInverseOperatorType;

      typedef typename BaseType::DomainFunctionType DomainFunctionType;
      typedef typename BaseType::RangeFunctionType RangeFunctionType;

      typedef typename BaseType::DomainFieldType DomainFieldType;

      /** constructor
       *
       *  \param[in]  op       operator to invert
       *
       *  \note The tolerance is read from the paramter
       *        <b>fem.solver.newton.tolerance</b>
       */
      explicit NewtonInverseOperator ( const OperatorType &op,
                                       const NewtonParameter &parameter = NewtonParameter() )
      : op_( op ),
        tolerance_( parameter.toleranceParameter() ),
        linAbsTol_( parameter.linAbsTolParameter( tolerance_ ) ),
        linReduction_( parameter.linReductionParameter( tolerance_ ) ),
        verbose_( parameter.verbose() && MPIManager::rank () == 0 ),
        matrixout_( parameter.matrixout()),
        linVerbose_( parameter.linearSolverVerbose() ),
        maxIterations_( parameter.maxIterationsParameter() ),
        maxLinearIterations_( parameter.maxLinearIterationsParameter() )
      {}

      /** constructor
       *
       *  \param[in]  op       operator to invert
       *  \param[in]  epsilon  tolerance for norm of residual
       */
      NewtonInverseOperator ( const OperatorType &op, const DomainFieldType &epsilon,
                              const NewtonParameter &parameter = NewtonParameter() )
      : op_( op ),
        tolerance_( epsilon ),
        linAbsTol_( parameter.linAbsTolParameter( tolerance_ ) ),
        linReduction_( parameter.linReductionParameter( tolerance_ ) ),
        verbose_( parameter.verbose() ),
        matrixout_( parameter.matrixout()),
        linVerbose_( parameter.linearSolverVerbose() ),
        maxIterations_( parameter.maxIterationsParameter() ),
        maxLinearIterations_( parameter.maxLinearIterationsParameter() )
      {}

      virtual void operator() ( const DomainFunctionType &u, RangeFunctionType &w ) const;

      int iterations () const { return iterations_; }
      int linearIterations () const { return linearIterations_; }

      bool converged () const
      {
        // check for finite |residual| - this also works for -ffinite-math-only (gcc)
        const bool finite = (delta_ < std::numeric_limits< DomainFieldType >::max());
        return finite && (iterations_ < maxIterations_) && (linearIterations_ < maxLinearIterations_);
      }

    private:
      const OperatorType &op_;
      const double tolerance_, linAbsTol_, linReduction_;
      const bool verbose_;
      const bool linVerbose_;
      const bool matrixout_;
      const int maxIterations_;
      const int maxLinearIterations_;

      mutable DomainFieldType delta_;
      mutable int iterations_;
      mutable int linearIterations_;
      ControlType control_;
      
    };



    // Implementation of NewtonInverseOperator
    // ---------------------------------------

    template< class JacobianOperator, class LInvOp, class Control> 
    inline void NewtonInverseOperator< JacobianOperator, LInvOp , Control >
      ::operator() ( const DomainFunctionType &u, RangeFunctionType &w ) const
    {
      DomainFunctionType residual( u );
      RangeFunctionType dw( w );
      RangeFunctionType wbacktrack( w );
      JacobianOperatorType jOp( "jacobianOperator", dw.space(), u.space() );
      double oldDelta;
      // compute initial residual
      op_( w, residual );
      residual -= u;
      delta_ = std::sqrt( residual.scalarProductDofs( residual ) );
      oldDelta=delta_;
      for( iterations_ = 0, linearIterations_ = 0; converged() && (delta_ > tolerance_); ++iterations_ )
      {
        if( verbose_ )
          std::cerr << "Newton iteration " << iterations_ << ": |residual| = " << delta_ << std::endl;

        // evaluate operator's jacobian
       if( control_.calcJacobian(iterations_, delta_))
       {
#if TIMING
          Dune::FemTimer::start();
#endif
          op_.jacobian( w, jOp );
#if TIMING
          std::cout<<"Assemble ="<<Dune::FemTimer::stop()<<"\n";
#endif
       }
#if 1 
if( matrixout_ )
        { 
          jOp.matrix().print(std::cout);
        }
#endif
   // David: With this factor, the tolerance of CGInverseOp is the absolute
        //        rather than the relative error
        //        (see also dune-fem/dune/fem/solver/inverseoperators.hh)
        const int remLinearIts = maxLinearIterations_ - linearIterations_;
        const LinearInverseOperatorType jInv( jOp, linReduction_, linAbsTol_, remLinearIts, linVerbose_ );

        dw.clear();
#if TIMING
        Dune::FemTimer::start();
#endif
        jInv( residual, dw );
#if TIMING
        std::cout<<"Solve ="<<Dune::FemTimer::stop()<<"\n";
#endif
        linearIterations_ += jInv.iterations();
        w -= dw;
        double damping=0.5;
        op_( w, residual );
        residual -= u;
        delta_ = std::sqrt( residual.scalarProductDofs( residual ) );
#if   0
        if( iterations_ > 0)    
        {      
          while( delta_>oldDelta)
            {
              
              w.assign(wbacktrack);
              damping/=2.;
              w.axpy(-1*damping,dw); 
              op_( w, residual );
              residual -= u;
              delta_ = std::sqrt( residual.scalarProductDofs( residual ) );
              std::cout<<"deltaBacktrack="<<delta_<<"\n"; 
              abort();
        }
}

        oldDelta=delta_;
#endif
      }
      if( verbose_ )
        std::cerr << "Newton iteration " << iterations_ << ": |residual| = " << delta_ << std::endl;
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_NEWTONINVERSEOPERATOR_HH

// vim: set expandtab ts=2 sw=2 sts=2:
#ifndef DUNE_FEM_NEWTONINVERSEOPERATOR_HH
#define DUNE_FEM_NEWTONINVERSEOPERATOR_HH
#define TIMING 0
#include <cfloat>
#include <iostream>
#include <limits>
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
        return Parameter::getValue< double >( "fem.solver.newton.linabstol",-1 );
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

    class NewtonBackTrackDefault
    {
      public:
        NewtonBackTrackDefault():
        thetaMin_(0.1),
        thetaMax_(0.5)
        {}

      //
      void set(double a, double b,double c ){a_=a;b_=b;c_=c;}

      double btFunction(double theta)
        {
          return a_*theta*theta+b_*theta+c_;
        }
      
      double btFunctionCritPt()
        {
          return b_/(2*a_);
        }

      double newTheta()
        {
          double res;
          if (btFunction( thetaMin_ ) < btFunction(thetaMax_))
            res=thetaMin_;
          else
            res=thetaMax_;

          if( btFunction(btFunctionCritPt())<btFunction(res) && btFunctionCritPt() <thetaMax_ && btFunctionCritPt()>thetaMin_)
            return btFunctionCritPt();
          else
            return res;
        }
      
        
      private:
        double thetaMin_,thetaMax_;
        double a_,b_,c_;
    };




    class NewtonControlDefault
    {
    public:
      NewtonControlDefault(double newtonTol):
      newtonTol_(newtonTol),
      assemblesteps_( Parameter::getValue< int >( "fem.solver.newton.assemblesteps",1)),
      gamma_( Parameter::getValue<double>( "fem.solver.newton.gamma", 2)),
      fac_( Parameter::getValue<double>("fem.solver.newton.etafactor",1)),
      minEta_( Parameter::getValue<double>( "fem.solver.newton.minEta",newtonTol/8)),
      choice_( Parameter::getValue<int>("fem.solver.newton.choice",1)) 
      {}
       
      double etaStart() const
      {
        return Parameter::getValue< double >( "fem.solver.newton.linreduction",newtonTol_/8  );
      }

      // return forcing term for residual of inner solver)
      double eta( double deltaOld, double deltaNow, double reduction) const
      {
        if(choice_==1)
          {
            return eta1( deltaOld,deltaNow,reduction);
          }
        else if(choice_==2)
          {
            return eta2( deltaOld,deltaNow);
          }
        else
          {
            return minEta_;
          }
       } 
      bool calcJacobian(int iterations, double  reduction) const
      {
        return (iterations%assemblesteps_==0);
      }

     private:

      double eta1( double deltaOld,double deltaNow, double reduction) const
      {
        return std::max( std::abs(deltaNow-reduction)/deltaOld, minEta_);
      }

      double eta2( double deltaOld,double deltaNow) const
      {
        return std::max( fac_*std::pow(deltaNow/deltaOld,gamma_), minEta_);
      }


      double newtonTol_;
      int assemblesteps_;
      double gamma_;
      double fac_;
      double minEta_;
      int choice_;
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
        linVerbose_( parameter.linearSolverVerbose() ),
        matrixout_( parameter.matrixout()),
        maxIterations_( parameter.maxIterationsParameter() ),
        maxLinearIterations_( parameter.maxLinearIterationsParameter() ),
        control_( tolerance_ ),
        failed_( false )
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
        linVerbose_( parameter.linearSolverVerbose() ),
        matrixout_( parameter.matrixout()),
        maxIterations_( parameter.maxIterationsParameter() ),
        maxLinearIterations_( parameter.maxLinearIterationsParameter() )
        failed_( false ),
      {}

      virtual void operator() ( const DomainFunctionType &u, RangeFunctionType &w ) const;

      int iterations () const { return iterations_; }

      int linearIterations () const { return linearIterations_; }

      bool converged () const
      {
        // check for finite |residual| - this also works for -ffinite-math-only (gcc)
        const bool finite = (delta_ < std::numeric_limits< DomainFieldType >::max());
        return finite  && (iterations_ < maxIterations_) && (linearIterations_ < maxLinearIterations_) && !failed_;
      }

    private:
      const OperatorType &op_;
      mutable double tolerance_, linAbsTol_, linReduction_;
      const bool verbose_;
      const bool linVerbose_;
      const bool matrixout_;
      const int maxIterations_;
      const int maxLinearIterations_;
      mutable DomainFieldType delta_;
      mutable int iterations_;
      mutable int linearIterations_;
      ControlType control_;
      mutable boole failed_;
      mutable NewtonBackTrackDefault backtrack_;

    };



    // Implementation of NewtonInverseOperator
    // ---------------------------------------

    template< class JacobianOperator, class LInvOp, class Control> 
    inline void NewtonInverseOperator< JacobianOperator, LInvOp , Control >
      ::operator() ( const DomainFunctionType &u, RangeFunctionType &w ) const
    {
      DomainFunctionType residual( u );
      DomainFunctionType resOld(u);
      RangeFunctionType dw( w );
      RangeFunctionType dwOld( w );
      RangeFunctionType wbacktrack(w);
      JacobianOperatorType jOp( "jacobianOperator", dw.space(), u.space() );
      double oldDelta;
      double eta=control_.etaStart();

      // compute initial residual
      //F(x_0)      
      op_( w, residual );
      residual -= u;
      // remember residual
      resOld.assign(residual);
      // ||F(x_0)||
      delta_ = std::sqrt( residual.scalarProductDofs( residual ) );
      oldDelta=delta_;
      
      for( iterations_ = 0, linearIterations_ = 0; converged() && (delta_ > tolerance_); ++iterations_ )
      {
        if( verbose_ )
          std::cerr << "Newton iteration " << iterations_ << ": |residual| = " << delta_ << std::endl;

        // evaluate operator's jacobian
        if( control_.calcJacobian(iterations_, delta_))
        {
          op_.jacobian( w, jOp );
        }
        //for debugging 
        if( matrixout_ )
        { 
          jOp.matrix().print(std::cout);
        }
        
        
        const LinearInverseOperatorType jInv( jOp, eta ,linAbsTol_, maxLinearIterations_,linVerbose_ );
        dw.clear();

        //solve the system 
        //dw=-s_k
        jInv( residual, dw );

        //linred=||F(x_k)+F'(x_k)s_k||< eta_k
        double linred=jInv.achievedreduction();
  
        linearIterations_ += jInv.iterations();

        RangeFunctionType dFx(w);
        dFx.clear();
        //F'(x_k)s_k

        jOp(dw,dFx);
        dFx*=-1;
        double dFxdots=dFx.scalarProductDofs(residual);

        //x_k+s_k
        w -= dw;


        //F(x_k-sk))
        op_( w, residual );
        residual -= u;

        //||F(x_k+s_k))||
        delta_ = std::sqrt( residual.scalarProductDofs( residual ) );

        eta=std::min(control_.eta(oldDelta,delta_,linred),0.9);

#if  1
        int count=0;
        double damping=0.9;
        if( iterations_ > 0)
        {
          while( delta_>(1-(1-eta)*1e-4)*oldDelta )
            {
              if(count>10)
                {
                  std::cout<<"Backtracking Failed!\n";
                  failed_=true;
                }
              backtrack_.set(delta_-2*dFxdots ,2*dFxdots , oldDelta);
              
              w.assign(wbacktrack);
              damping=backtrack_.newTheta();

              w.axpy(-1*damping,dw);
              op_( w, residual );

              residual -= u;
              delta_ = std::sqrt( residual.scalarProductDofs( residual ) );
              eta=1-damping*(1-eta);
              count++;
        }

      }

#endif
      wbacktrack.assign(w);
      oldDelta=delta_;
     }
      if( verbose_ )
        std::cerr << "Newton iteration " << iterations_ << ": |residual| = " << delta_ << std::endl;
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_NEWTONINVERSEOPERATOR_HH

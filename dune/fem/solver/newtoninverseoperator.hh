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
        NewtonBackTrackDefault()
        {}
            
      double btFunction(double theta)
        {
          return a_*theta*theta+b_*theta+c_;
        }
      
      double btFunctionCritPt()
        {
          return -2*a_/b_;
        }

      double newTheta()
        {
          return std::min( btFunctionCritPt() , std::min( btFunction(thetaMin_) , btFunction(thetaMax_) ) );
        }
      
        
      private:
        double thetaMin_,thetaMax_;
        double a_,b_,c_,t_;
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
        control_( tolerance_ )
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
        
        const int remLinearIts = maxLinearIterations_ - linearIterations_;
        
        
        const LinearInverseOperatorType jInv( jOp, eta ,linAbsTol_, maxLinearIterations_,linVerbose_ );
        dw.clear();

        //solve the system 
        //dw=s_0 
        jInv( residual, dw );

        //||F(x_0)+F'(x_0)s_0||
        double linred=jInv.achievedreduction();
  
        DomainFunctionType residualtmp(residual);
        RangeFunctionType dFs(dw);
        jOp(dw,dFs);
        residualtmp-=dFs;
        double myRed=std::sqrt( residualtmp.scalarProductDofs(residualtmp ));
        
        linearIterations_ += jInv.iterations();

        //x_1
        w -= dw;

        //F(x_1)
        op_( w, residual );
        residual -= u;

        //||F(x_1)||
        delta_ = std::sqrt( residual.scalarProductDofs( residual ) );

        if( verbose_) 
          std::cout<<"delta= "<<delta_<<"\t oldDelta= "<<oldDelta<<"\t linred= "<<linred<<"\t eta= "<<eta<<std::endl;

        eta=std::min(control_.eta(oldDelta,delta_,linred),0.9);

 #if  0
        double damping=0.5;
        if( iterations_ > 1)    
        {      
          while( delta_>(1-(1-eta)*1e-4)*oldDelta )
            {
              
              w.assign(wbacktrack);
              damping/=2.;
              w.axpy(-1*damping,dw); 
              op_( w, residual );
              residual -= u;
              delta_ = std::sqrt( residual.scalarProductDofs( residual ) );
              std::cout<<"deltaBacktrack="<<delta_<<"\n"; 
              eta=1-damping*(1-eta); 
        }

      }

     #endif
      oldDelta=delta_;
     }
      if( verbose_ )
        std::cerr << "Newton iteration " << iterations_ << ": |residual| = " << delta_ << std::endl;
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_NEWTONINVERSEOPERATOR_HH

#ifndef DUNE_FEM_DG_WELLBALOPERATOR_HH
#define DUNE_FEM_DG_WELLBALOPERATOR_HH
#warning "WELLBALANCEDOP"
#include <string>

// dune-fem includes
#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/pass/dgdiscretemodel.hh>
#include <dune/fem/pass/selection.hh>
#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/space/dgspace.hh>
#include <dune/fem/operator/common/spaceoperatorif.hh>

// dune-fem-dg includes
#include <dune/fem-dg/operator/limiter/limitpass.hh>

// note to me: it doesn't make sense to include primaldiscretemodel.hh
//             but it was design in that way. to be removed later!!!!!!
// local includes
//#include <dune/fem/operator/primaldiscretemodel.hh>

#include "wellbalanceddiscretemodel.hh" 

#include <dune/fem-dg/operator/dg/operatorbase.hh>
#include <dune/fem-dg/pass/dgpass.hh>



//PassTraits are defined  in <dune-phasefield/dune/fem/opeartor/projdiscretemodelcommon.hh>


namespace Dune {  
	
  template< class Model, class NumFlux, 
						int polOrd, bool advection = true , bool diffusion = true >
  class  DGAdvectionDiffusionOperator : 
    public Fem::SpaceOperatorInterface 
    < typename PassTraits< Model, Model::Traits::dimRange, polOrd > :: DestinationType >
  {
    // Id's for the three Passes (including StartPass)
    enum PassIdType { u, gradPass, acPass,navstkPass };    
 
 public:
    enum { dimRange =  Model::dimRange };
    enum { dimDomain = Model::Traits::dimDomain };
    
    typedef NumFlux NumFluxType;

    typedef PassTraits< Model, Model::Traits::dimRange, polOrd >     PassTraitsType ;
    
    typedef typename PassTraitsType::ScalarDiscreteFunctionSpaceType ScalarDiscreteFunctionSpaceType;
    
		typedef typename PassTraitsType::DiscreteScalarType DiscreteScalarType;
		typedef DiscreteScalarType ScalarDFType;
    

 
    typedef typename PassTraitsType::IndicatorType                   IndicatorType;
    typedef typename IndicatorType::DiscreteFunctionSpaceType        IndicatorSpaceType;

    // Pass 3 Model (advection)
    typedef AdvectionDiffusionLDGModel< Model, NumFluxType, polOrd, u, acPass, gradPass, advection, diffusion >    DiscreteModel3Type;
   
		typedef typename DiscreteModel3Type :: DiffusionFluxType  DiffusionFluxType;
    
		//Pass  2 (allen-cahn and chemical potential)
    typedef ThetaModel< Model, polOrd, u, gradPass >  DiscreteModel2Type;
    
    // Pass 1 Model (gradient)
    typedef GradientModel<Model, DiffusionFluxType, polOrd, u >       DiscreteModel1Type;
    
                                                                       /*@LST0E@*/
    typedef typename DiscreteModel1Type :: Traits                      Traits1;
    typedef typename DiscreteModel2Type :: Traits                      Traits2;
    typedef typename DiscreteModel3Type :: Traits                      Traits3;
    
    typedef typename Model :: Traits :: GridType                       GridType;
    
    typedef typename Traits3 :: DomainType                             DomainType;
    typedef typename Traits3 :: DiscreteFunctionType                   DiscreteFunction3Type;
    
  //   typedef typename Traits3 :: DiscreteFunctionType                   DiscreteFunction3Type;
		//Space of Gradient
		typedef typename Traits1 :: DiscreteFunctionSpaceType              Space1Type;
		//Space of Theta
    typedef typename Traits2 :: DiscreteFunctionSpaceType              Space2Type;
    //Space of the solutionvector
		typedef typename Traits3 :: DiscreteFunctionSpaceType              Space3Type;
   
		//sigma
    typedef typename Traits1 :: DestinationType                        Destination1Type;
    //theta
		typedef typename Traits2 :: DestinationType                        Destination2Type;
    //solutionn
		typedef typename Traits3 :: DestinationType                        Destination3Type;
  
    typedef Destination3Type                                           DestinationType;
    typedef Space3Type                                                 SpaceType;
  
    typedef typename Traits1 :: GridPartType                           GridPartType;
    
    typedef Fem::StartPass< DiscreteFunction3Type, u >                      Pass0Type; /*@LST0S@*/
    

    typedef LocalCDGPass< DiscreteModel1Type, Pass0Type, gradPass >    Pass1Type; /*@\label{ad:typedefpass1}@*/
    typedef LocalCDGPass< DiscreteModel2Type, Pass1Type, acPass >      Pass2Type; /*@\label{ad:typedefpass1}@*/
    typedef LocalCDGPass< DiscreteModel3Type, Pass2Type, navstkPass >  Pass3Type; /*@\label{ad:typedefpass2}@*//*@LST0E@*/
    
  public:
		DGAdvectionDiffusionOperator(GridType& grid , const NumFluxType& numf ) :
      grid_( grid ),
      model_( numf.model() ),
      numflux_( numf ),
      gridPart_( grid_ ),
      space1_(gridPart_),
      space2_(gridPart_),
      space3_(gridPart_),
      diffFlux_( gridPart_, model_ ),
      discModel1_(model_, diffFlux_ ),
      discModel2_(model_),
      discModel3_(model_, numflux_, diffFlux_),
      pass0_ (),
      pass1_(discModel1_, pass0_, space1_),    /*@\label{ad:initialisepass1}@*/
      pass2_(discModel2_, pass1_, space2_),     /*@\label{ad:initialisepass2}@*/
      pass3_(discModel3_, pass2_, space3_)     /*@\label{ad:initialisepass2}@*/
    { }

    void setTime(const double time) {
	    pass3_.setTime( time );
    }

    double timeStepEstimate() const {
	    return pass3_.timeStepEstimate();
    }

    void operator()( const DestinationType& arg, DestinationType& dest ) const {
	    pass3_( arg, dest );
    }
    
    void theta( const DestinationType& arg, Destination2Type& dest ) const {
      pass2_(arg,dest);
    }
    void gradient( const DestinationType& arg, Destination1Type& dest ) const {
      pass1_(arg,dest);
    }
    
     void energy(const DestinationType &u,const Destination1Type& du,  ScalarDFType& en)const
     {
       abort();
       // du.print(cout); 
     }


    inline const SpaceType& space() const {
	    return space3_;
    } 

    inline void switchupwind() 
    { 
      diffFlux_.switchUpwind();
    }

    inline double maxAdvectionTimeStep() const 
    {
      return discModel3_.maxAdvectionTimeStep();
    } 
    inline double maxDiffusionTimeStep() const 
    {
      return discModel3_.maxDiffusionTimeStep();
    } 

    inline void limit(const DestinationType& arg,DestinationType& dest) const
    {}
    
    inline double computeTime() const 
    {
      return pass3_.computeTime();
    }

    inline size_t numberOfElements () const 
    {
      return pass3_.numberOfElements();
    }
    
    void printmyInfo(std::string filename) const {
	    std::ostringstream filestream;
            filestream << filename;
            std::ofstream ofs(filestream.str().c_str(), std::ios::app);
            ofs << "LDG Op., polynomial order: " << polOrd << "\\\\\n\n";
            ofs.close();
    }

    std::string description() const
    {
      std::stringstream stream;
      stream <<" {\\bf LDG Diff. Op.}, flux formulation, order: " << polOrd+1
             <<", $\\eta = ";
      diffFlux_.diffusionFluxPenalty( stream );
      stream <<"$, {\\bf Adv. Flux:} ";
      if (FLUX==1)
        stream <<"LLF";
      else if (FLUX==2)
        stream <<"HLL";
      stream <<",\\\\\n";
      return stream.str();
    }

  private:
    GridType&           grid_;
    const Model&        model_;
    const NumFluxType&  numflux_;
    GridPartType        gridPart_;
    Space1Type          space1_;
    Space2Type          space2_;
    Space3Type          space3_;

  protected:
    DiffusionFluxType   diffFlux_;
    
  private:
    DiscreteModel1Type  discModel1_;
    DiscreteModel2Type  discModel2_;
    DiscreteModel3Type  discModel3_;
   
    Pass0Type           pass0_;
    Pass1Type           pass1_;
    Pass2Type           pass2_;
    Pass3Type           pass3_;
  };


  // LDGAdvectionTraits
  //-------------------
#if 0
  template <class Mod, class NumFlux, 
            int pOrd,
            bool advection>
  struct LDGAdvectionTraits
  {
    enum PassIdType { u, cdgpass };    /*@\label{ad:passids}@*/

    typedef Mod  Model;
    enum { dimRange = Model::dimRange };
    typedef NumFlux   NumFluxType;
    enum { polOrd = pOrd };
    typedef AdvectionDiffusionDGPrimalModel
      // put a method_none here to avoid diffusion 
      < Model, NumFluxType, method_none, polOrd, u, advection, false> DiscreteModelType;
  };

  
  // DGAdvectionOperator
  //--------------------

  template< class Model, class NumFlux, 
            DGDiffusionFluxIdentifier diffFluxId, // dummy parameter 
            int polOrd >
  class DGAdvectionOperator : public
    DGAdvectionDiffusionOperatorBase< LDGAdvectionTraits<Model, NumFlux, polOrd, true> >
  {
    typedef LDGAdvectionTraits<Model, NumFlux, polOrd, true> Traits;
    typedef DGAdvectionDiffusionOperatorBase< Traits >  BaseType;
    typedef typename BaseType :: GridType  GridType;
    typedef typename BaseType :: NumFluxType  NumFluxType;
  public:
    DGAdvectionOperator( GridType& grid , const NumFluxType& numf )
      : BaseType( grid, numf )
    {}

    std::string description() const
    {
      std::stringstream stream;
      stream <<"{\\bf Adv. Op.}, flux formulation, order: " << polOrd+1
             <<", {\\bf Adv. Flux:} ";
      stream <<",\\\\\n";
      return stream.str();
    }
  };
#endif

  // DGDiffusionOperator
  //--------------------

  template< class Model, class NumFlux, 
            DGDiffusionFluxIdentifier diffFluxId, // dummy parameter 
            int polOrd >
  class DGDiffusionOperator : public
    DGAdvectionDiffusionOperator< Model, NumFlux, diffFluxId, polOrd, false > 
  {
    typedef DGAdvectionDiffusionOperator< Model, NumFlux, diffFluxId, polOrd, false >  BaseType;
    typedef typename BaseType :: GridType  GridType;
    typedef typename BaseType :: NumFluxType  NumFluxType;

  public:
    DGDiffusionOperator( GridType& grid , const NumFluxType& numf )
      : BaseType( grid, numf )
    {}

    std::string description() const
    {
      std::stringstream stream;
      stream <<"{\\bf LDG Diff. Op.}, flux formulation, order: " << polOrd+1
             <<", $\\eta = ";
      diffFlux_.diffusionFluxPenalty( stream );
      stream <<"$, {\\bf Adv. Flux:} ";
      stream <<"None";
      stream <<",\\\\\n";
      return stream.str();
    }

  private:
    using BaseType::diffFlux_;
  };

#if 0
  // DGLimitedAdvectionDiffusionOperator
  //------------------------------------

  /** \class DGLimitedAdvectionDiffusionOperator
   *  \brief Dual operator for NS equtions with a limiting
   *         of the numerical solution
   *
   *  \tparam Model Analytical model
   *  \tparam NumFlux Numerical flux
   *  \tparam polOrd Polynomial degree
   *  \tparam advection Advection
   */
  template< class Model, class NumFlux, 
            DGDiffusionFluxIdentifier diffFluxId, // dummy parameter 
            int polOrd, bool advection = true >
  class DGLimitedAdvectionDiffusionOperator :
    public SpaceOperatorInterface 
      < typename PassTraits< Model, Model::Traits::dimRange, polOrd > :: DestinationType >
  {
    enum PassIdType { u, limitPassId, gradPassId, navstkPassId };    /*@\label{ad:passids}@*/

  public:
    enum { dimRange = Model::dimRange };
    enum { dimDomain = Model::Traits::dimDomain };

    typedef NumFlux NumFluxType;
    typedef PassTraits< Model, dimRange, polOrd > PassTraitsType;

    // Pass 2 Model (navstkPassId)
    typedef AdvectionDiffusionLDGModel
      < Model, NumFluxType, polOrd, limitPassId, gradPassId, advection, true > DiscreteModel3Type;

    // Pass 1 Model (gradPassId)
    typedef typename DiscreteModel3Type :: DiffusionFluxType  DiffusionFluxType;
    typedef GradientModel< Model, DiffusionFluxType, polOrd, limitPassId >
                                                                       DiscreteModel2Type;
    // The model of the limiter pass (limitPassId)
    typedef Fem :: StandardLimiterDiscreteModel< PassTraitsType, Model, u > LimiterDiscreteModelType;

    // Pass 0 Model (limitPassId)
    typedef LimiterDiscreteModelType                                   DiscreteModel1Type;


                                                                       /*@LST0E@*/
    typedef typename DiscreteModel1Type :: Traits                      Traits1;
    typedef typename DiscreteModel2Type :: Traits                      Traits2;
    typedef typename DiscreteModel3Type :: Traits                      Traits3;

    typedef typename Model :: Traits :: GridType                       GridType;

    typedef typename Traits3 :: DomainType                             DomainType;
    typedef typename Traits3 :: DiscreteFunctionType                   DiscreteFunction3Type;

    typedef typename Traits1 :: DiscreteFunctionSpaceType              Space1Type;
    typedef typename Traits2 :: DiscreteFunctionSpaceType              Space2Type;
    typedef typename Traits3 :: DiscreteFunctionSpaceType              Space3Type;
    typedef typename Traits2 :: DestinationType                        Destination2Type;
    typedef typename Traits3 :: DestinationType                        Destination3Type;
    typedef Destination3Type                                           DestinationType;
    typedef Space3Type                                                 SpaceType;

    typedef typename Traits2 :: GridPartType                           GridPartType;

    typedef StartPass< DiscreteFunction3Type, u >                      Pass0Type; /*@LST0S@*/
    typedef LimitDGPass< DiscreteModel1Type, Pass0Type, limitPassId >    Pass1Type; /*@\label{ad:typedefpass1}@*/
    typedef LocalCDGPass< DiscreteModel2Type, Pass1Type, gradPassId >    Pass2Type; /*@\label{ad:typedefpass1}@*/
    typedef LocalCDGPass< DiscreteModel3Type, Pass2Type, navstkPassId >  Pass3Type; /*@\label{ad:typedefpass2}@*//*@LST0E@*/

    typedef typename PassTraitsType::IndicatorType                     IndicatorType;
    typedef typename IndicatorType::DiscreteFunctionSpaceType          IndicatorSpaceType;

    template <class Limiter, int pOrd>
    struct LimiterCall
    {
      template <class ArgumentType, class DestinationType>
      static inline void limit(const Limiter& limiter,
                               ArgumentType* arg,
                               DestinationType& dest)
      {
        limiter.enableFirstCall();
        assert( arg );
        arg->assign(dest);
        limiter(*arg,dest);
        limiter.disableFirstCall();
      }
    };

    template <class Limiter>
    struct LimiterCall<Limiter,0>
    {
      template <class ArgumentType, class DestinationType>
      static inline void limit(const Limiter& limiter,
                               const ArgumentType* arg,
                               DestinationType& dest)
      {
      }
    }; 

  public:
    DGLimitedAdvectionDiffusionOperator( GridType& grid , const NumFluxType& numf ) 
      : grid_( grid )
      , model_( numf.model() )
      , numflux_( numf )
      , gridPart_( grid_ )
      , space1_(gridPart_)
      , space2_(gridPart_)
      , space3_(gridPart_)
      , uTmp_( (polOrd > 0) ? (new DestinationType("limitTmp", space3_)) : 0 )
      , fvSpc_( gridPart_ )
      , indicator_( "Indicator", fvSpc_ )
      , diffFlux_( gridPart_, model_ )
      , discModel1_( model_ )
      , discModel2_( model_, diffFlux_ )
      , discModel3_( model_, numflux_, diffFlux_ )
      , pass0_()
      , pass1_(discModel1_, pass0_, space1_)    /*@\label{ad:initialisepass1}@*/
      , pass2_(discModel2_, pass1_, space2_)    /*@\label{ad:initialisepass1}@*/
      , pass3_(discModel3_, pass2_, space3_)     /*@\label{ad:initialisepass2}@*/
    {
      discModel1_.setIndicator( &indicator_ );
    }

    ~DGLimitedAdvectionDiffusionOperator() { delete uTmp_; }

    void setTime(const double time) {
	    pass3_.setTime( time );
    }

    double timeStepEstimate() const {
	    return pass3_.timeStepEstimate();
    }

    void operator()( const DestinationType& arg, DestinationType& dest ) const {
	    pass3_( arg, dest );
      pass1_.enableFirstCall();
    }

    inline const SpaceType& space() const {
	    return space3_;
    } 

    inline void switchupwind() 
    { 
      diffFlux_.switchUpwind();
    }
    double limitTime() const
    {
      return pass1_.computeTime();
    }
    std::vector<double> limitSteps() const
    {
      return pass1_.computeTimeSteps();
    }
    const Pass1Type& limitPass() const
    {
      return pass1_;
    }

    inline void limit( const DestinationType& arg, DestinationType& dest) const
    {
      pass1_.enableFirstCall();
      LimiterCall< Pass1Type, polOrd >::limit( pass1_, uTmp_, dest );
    }
    
    void printmyInfo(std::string filename) const {
	    std::ostringstream filestream;
            filestream << filename;
            std::ofstream ofs(filestream.str().c_str(), std::ios::app);
            ofs << "Limited LDG Op., polynomial order: " << polOrd << "\\\\\n\n";
            ofs.close();
    }

    std::string description() const
    {
      std::stringstream stream;
      stream <<" {\\bf LDG Diff. Op.}, dual form, order: " << polOrd+1
             <<", penalty: ";
      diffFlux_.diffusionFluxPenalty( stream );
      stream <<", {\\bf Adv. Flux:} ";
      if (FLUX==1)
        stream <<"LLF";
      else if (FLUX==2)
        stream <<"HLL";
      stream <<",\\\\\n";
      return stream.str();
    }

  private:
    GridType&           grid_;
    const Model&        model_;
    const NumFluxType&  numflux_;
    GridPartType        gridPart_;
    Space1Type          space1_;
    Space2Type          space2_;
    Space3Type          space3_;
    mutable DestinationType* uTmp_;
    IndicatorSpaceType  fvSpc_;
    IndicatorType       indicator_;

  protected:
    DiffusionFluxType   diffFlux_;
    
  private:
    DiscreteModel1Type  discModel1_;
    DiscreteModel2Type  discModel2_;
    DiscreteModel3Type  discModel3_;
    Pass0Type           pass0_;
    Pass1Type           pass1_;
    Pass2Type           pass2_;
    Pass3Type           pass3_;
  };
#endif
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
#if 0
  template <class Mod, class NumFlux, 
            DGDiffusionFluxIdentifier diffFluxId,
            int pOrd, 
            bool advection, bool diffusion >
  struct AdaptationIndicatorTraits 
  {
    enum { u, cdgpass };
    
    typedef Mod  Model;
    enum { dimRange = Model::dimRange };
    typedef NumFlux NumFluxType;
    enum { polOrd = pOrd };
    typedef AdaptiveAdvectionDiffusionDGPrimalModel
      < Model, NumFluxType, diffFluxId, polOrd, u, advection, diffusion> DiscreteModelType;
  };

  // DGAdaptationIndicatorOperator 
  //------------------------------

  template< class Model, class NumFlux,
            DGDiffusionFluxIdentifier diffFluxId, int polOrd,
            bool advection, bool diffusion = false >
  struct DGAdaptationIndicatorOperator : public 
    DGAdvectionDiffusionOperatorBase< 
       AdaptationIndicatorTraits< Model, NumFlux, diffFluxId, polOrd, advection, diffusion > >
  {
    typedef AdaptationIndicatorTraits< Model, NumFlux, diffFluxId, polOrd, advection, diffusion > Traits ;
    typedef DGAdvectionDiffusionOperatorBase< Traits >  BaseType;
    typedef typename BaseType :: GridType  GridType;
    typedef typename BaseType :: NumFluxType  NumFluxType;

    DGAdaptationIndicatorOperator( GridType& grid , const NumFluxType& numf ) 
      : BaseType( grid, numf )
    {
      if (Fem::Parameter::verbose())
      {
        std::cerr <<"\nWARNING: The adaptation indicator based on Ohlberger's a-posteriori\n";
        std::cerr <<"         error estimator is not ment to be used with flux formulation.\n\n";
      }
    }

    std::string description() const
    {
      std::stringstream stream;
      discreteModel_.diffusionFlux().diffusionFluxName( stream );
      stream <<" {\\bf Adv. Op.} in primal formulation, order: " << polOrd+1
             <<", $\\eta = ";
      discreteModel_.diffusionFlux().diffusionFluxPenalty( stream );
      stream <<"$, $\\chi = ";
      discreteModel_.diffusionFlux().diffusionFluxLiftFactor( stream );
      stream << "$, {\\bf Adv. Flux:} " << numflux_.name() << ",\\\\" << std::endl;
      return stream.str();
    }

  protected:
    using BaseType::discreteModel_;
    using BaseType::numflux_;
  };
#endif
}
#endif

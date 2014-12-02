#ifndef DUNE_PHASEFIELD_FLUXPROJOPERATOR_HH
#define DUNE_PHASEFIELD_FLUXPROJOPERATOR_HH
#include <string>

// dune-fem includes
#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/pass/localdg/discretemodel.hh>
#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/space/discontinuousgalerkin.hh>
#include <dune/fem/operator/common/spaceoperatorif.hh>

// dune-fem-dg includes
#include <dune/fem-dg/operator/limiter/limitpass.hh>

// local includes
//#include <dune/fem-dg/operator/dg/primaldiscretemodel.hh>
//#include <dune/fem-dg/operator/dg/fluxdiscretemodel.hh>
#include "fluxprojdiscretemodel.hh" 
//#include <dune/fem-dg/operator/dg/operatorbase.hh>
#include <dune/fem-dg/pass/dgpass.hh>


namespace Dune {  

  // DGAdvectionDiffusionOperator
  //-----------------------------

  template< class Model, class NumFlux, 
            DGDiffusionFluxIdentifier diffFluxId,
            int polOrd, bool advection = true , bool diffusion = true >
  class DGAdvectionDiffusionOperator : 
  public Dune::Fem::SpaceOperatorInterface 
    < typename MyPassTraits< Model, Model::Traits::dimRange, polOrd > :: DestinationType >
  {
    // Id's for the three Passes (including StartPass)
    enum PassIdType { u, gradPass, projPass,advectPass };    /*@\label{ad:passids}@*/
    
  public:
    enum { dimRange =  Model::dimRange };
    enum { dimDomain = Model::Traits::dimDomain };
    
    typedef NumFlux NumFluxType;

    typedef MyPassTraits< Model, Model::Traits::dimRange, polOrd >     PassTraitsType ;
    
    typedef typename PassTraitsType::ScalarDiscreteFunctionSpaceType ScalarDiscreteFunctionSpaceType;
		typedef typename PassTraitsType::DiscreteScalarType DiscreteScalarType;   
		typedef typename PassTraitsType::ScalarDFType ScalarDFType;
    

 
    typedef typename PassTraitsType::IndicatorType                   IndicatorType;
    typedef typename IndicatorType::DiscreteFunctionSpaceType        IndicatorSpaceType;

    // Pass 3 Model (advection)
    typedef AdvectionDiffusionLDGModel< Model, NumFluxType, polOrd, u, projPass, gradPass, advection, diffusion >    DiscreteModel3Type;
    
    typedef ProjectionModel< Model, polOrd, u, gradPass >  DiscreteModel2Type;
    
    // Pass 1 Model (gradient)
    typedef typename DiscreteModel3Type :: DiffusionFluxType  DiffusionFluxType;
    typedef GradientModel<Model, DiffusionFluxType, polOrd, u >       DiscreteModel1Type;
    
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
   
    typedef typename Traits1 :: DestinationType                        Destination1Type;
    typedef typename Traits2 :: DestinationType                        Destination2Type;
    typedef typename Traits3 :: DestinationType                        Destination3Type;
  
    typedef Destination3Type                                           DestinationType;
    typedef Space3Type                                                 SpaceType;
  
    typedef typename Traits1 :: GridPartType                           GridPartType;
    
    typedef Dune::Fem::StartPass< DiscreteFunction3Type, u >                      Pass0Type;
    

    typedef LocalCDGPass< DiscreteModel1Type, Pass0Type, gradPass >    Pass1Type; 
    typedef LocalCDGPass< DiscreteModel2Type, Pass1Type, projPass >    Pass2Type;
    typedef LocalCDGPass< DiscreteModel3Type, Pass2Type, advectPass >  Pass3Type; 

  public:
    DGAdvectionDiffusionOperator( GridType& grid , const NumFluxType& numf ) :
      grid_( grid ),
      model_( numf.model() ),
      numflux_( numf ),
      gridPart_( grid_ ),
      space1_(gridPart_),
      space2_(gridPart_),
      space3_(gridPart_),
      diffFlux_( gridPart_, model_ ),
      problem1_(model_, diffFlux_ ),
      problem2_(model_),
      problem3_(model_, numflux_, diffFlux_),
      pass0_ (),
      pass1_(problem1_, pass0_, space1_),    /*@\label{ad:initialisepass1}@*/
      pass2_(problem2_, pass1_, space2_),     /*@\label{ad:initialisepass2}@*/
      pass3_(problem3_, pass2_, space3_)     /*@\label{ad:initialisepass2}@*/
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
    


    inline const SpaceType& space() const {
	    return space3_;
    } 

    inline void switchupwind() 
    { 
      diffFlux_.switchUpwind();
    }

    inline double maxAdvectionTimeStep() const 
    {
      return problem3_.maxAdvectionTimeStep();
    } 
    inline double maxDiffusionTimeStep() const 
    {
      return problem3_.maxDiffusionTimeStep();
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
    DiscreteModel1Type  problem1_;
    DiscreteModel2Type  problem2_;
    DiscreteModel3Type  problem3_;
   
    Pass0Type           pass0_;
    Pass1Type           pass1_;
    Pass2Type           pass2_;
    Pass3Type           pass3_;
  };


  // LDGAdvectionTraits
  //-------------------
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

}
#endif

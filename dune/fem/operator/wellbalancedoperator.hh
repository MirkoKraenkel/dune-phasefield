#ifndef DUNE_FEM_DG_WELLBALOPERATOR_HH
#define DUNE_FEM_DG_WELLBALOPERATOR_HH
#warning "WELLBALANCEDOP"
#include <string>

// dune-fem includes
#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/pass/localdg/discretemodel.hh>
#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/space/discontinuousgalerkin.hh>
#include <dune/fem/operator/common/spaceoperatorif.hh>



#include "wellbalanceddiscretemodel.hh" 

#include <dune/fem-dg/pass/dgpass.hh>



//PassTraits are defined  in <dune-phasefield/dune/fem/opeartor/projdiscretemodelcommon.hh>


namespace Dune {  
	
  template< class Model, class NumFlux,int polOrd, bool advection = true , bool diffusion = true >
  class  DGAdvectionDiffusionOperator : 
  public Fem::SpaceOperatorInterface< typename MyPassTraits< Model, Model::Traits::dimRange, polOrd > :: DestinationType >
  {
    // Id's for the three Passes (including StartPass)
    enum PassIdType { u, gradPass, acPass,navstkPass };    
 
   public:
    enum { dimRange =  Model::dimRange };
    enum { dimDomain = Model::Traits::dimDomain };
    
    typedef NumFlux NumFluxType;

    typedef MyPassTraits< Model, Model::Traits::dimRange, polOrd >     PassTraitsType ;
    
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
    
    typedef typename DiscreteModel1Type :: Traits                      Traits1;
    typedef typename DiscreteModel2Type :: Traits                      Traits2;
    typedef typename DiscreteModel3Type :: Traits                      Traits3;
    
    typedef typename Model :: Traits :: GridType                       GridType;
    
    typedef typename Traits3 :: DomainType                             DomainType;
    typedef typename Traits3 :: DiscreteFunctionType                   DiscreteFunction3Type;
    
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
    
    typedef Fem::StartPass< DiscreteFunction3Type, u >                 Pass0Type; 
    

    typedef LocalCDGPass< DiscreteModel1Type, Pass0Type, gradPass >    Pass1Type; 
    typedef LocalCDGPass< DiscreteModel2Type, Pass1Type, acPass >      Pass2Type; 
    typedef LocalCDGPass< DiscreteModel3Type, Pass2Type, navstkPass >  Pass3Type; 
    
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
      pass1_(discModel1_, pass0_, space1_),    
      pass2_(discModel2_, pass1_, space2_),   
      pass3_(discModel3_, pass2_, space3_)  
    { }

    void setTime(const double time) 
    {
	    pass3_.setTime( time );
    }

    double timeStepEstimate() const 
    {
      return pass3_.timeStepEstimate();
    }

    void operator()( const DestinationType& arg, DestinationType& dest ) const 
    {
	    pass3_( arg, dest );
    }
    
    void theta( const DestinationType& arg, Destination2Type& dest ) const 
    {
      pass2_(arg,dest);
    }
    
    void gradient( const DestinationType& arg, Destination1Type& dest ) const 
    {
      pass1_(arg,dest);
    }
    
    inline const SpaceType& space() const 
    {
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

#ifndef PROBLEM_HH
#define PROBLEM_HH

#include <dune/common/version.hh>

// dune-fem includes
#include <dune/fem/misc/linesegmentsampler.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/space/common/functionspace.hh>

// local includes
// #include "idealthermodynamics2interpol.hh"
#if 1
#include <dune/phasefield/modelling/combinedthermodynamics.hh>
#else
#include "thermoequal.hh"
#endif

#include <dune/fem/probleminterfaces.hh>




namespace Dune {

template <class GridType>
class PhaseProblem : public EvolutionProblemInterface<
                    Dune::Fem::FunctionSpace< double, double, GridType::dimension, GridType::dimension + 2 >,
                    true >
                  
{
 
  public:
  typedef Fem::FunctionSpace<typename GridType::ctype,
                        double, GridType::dimensionworld,
                        GridType::dimensionworld + 2 > FunctionSpaceType ;

  enum{ dimension = GridType::dimensionworld };
  enum { energyId = dimension + 1 };
  typedef typename FunctionSpaceType :: DomainFieldType   DomainFieldType;
  typedef typename FunctionSpaceType :: DomainType        DomainType;
  typedef typename FunctionSpaceType :: RangeFieldType    RangeFieldType;
  typedef typename FunctionSpaceType :: RangeType         RangeType;

  typedef CombinedThermodynamics<VapourDynamics,LiquidDynamics>  ThermodynamicsType;

  PhaseProblem() : 
		 myName_( "Problem1" ),
		 endTime_ ( Fem::Parameter::getValue<double>( "phasefield.endTime" )), 
		 mu_( Fem::Parameter :: getValue< double >( "mu" )),
		 delta_(Fem::Parameter::getValue<double>( "phasefield.delta" )),
		 c1_(Fem::Parameter :: getValue< double >( "c1")),
		 c2_(Fem::Parameter ::getValue< double> ("c2")),
		 smear_( Fem::Parameter::getValue<double> ("smear")),
		 phiscale_(Fem::Parameter::getValue<double> ("phiscale")),
		 gamma_(Fem::Parameter::getValue<double> ("gamma")),
     liqdyn_(),
     vapdyn_(),
     thermodyn_(vapdyn_,liqdyn_)
     {
      thermodyn_.init();
      std::cout<<"combinedProblem constructor\n";
     }


  // initialize A and B 
  double init(const bool returnA ) const ;

  // print info 
  void printInitInfo() const;

  // source implementations 
  inline bool hasStiffSource() { return true; }
  inline bool hasNonStiffSource() { return false; }
  inline double stiffSource( const double t, const DomainType& x,const RangeType& u,RangeType& res ) const;
  inline double nonStiffSource( const double t, const DomainType& x,const RangeType& u, RangeType& res ) const;

  // this is the initial data
  inline void evaluate( const DomainType& arg , RangeType& res ) const 
  {
    evaluate( 0., arg, res );
  }

  // evaluate function 
  inline void evaluate( const double t, const DomainType& x, RangeType& res ) const;

  // cloned method 
  inline void evaluate( const DomainType& x, const double t, RangeType& res ) const
  {
    evaluate( t, x, res );
  }



  template< class DiscreteFunctionType >
  void finalizeSimulation( DiscreteFunctionType& variablesToOutput,
                           const int eocloop) const
  {}


  const ThermodynamicsType& thermodynamics() const { return thermodyn_;}
  void printmyInfo( std::string filename ) const {}
  inline double endtime() const { return endTime_; }
  inline std::string myName() const { return myName_; }
  void paraview_conv2prim() const {}
  std::string description() const;
 
  inline double mu() const { return mu_; }
  inline double delta() const{return delta_;}
  inline double c1() const {return c1_;}
  inline double c2() const {return c2_;}
  inline double gamma() const {return gamma_;}
protected:
  const std::string myName_;
  const double endTime_;
  const double mu_;
  const double delta_;
  const double c1_,c2_;
  const double smear_;
  const double phiscale_;
  const double gamma_;
  const   LiquidDynamics  liqdyn_;
  const VapourDynamics vapdyn_;
  const ThermodynamicsType thermodyn_;
  
};


template <class GridType>
inline double PhaseProblem<GridType>
:: init(const bool returnA ) const 
{

  return 0;
}



template <class GridType>
inline void PhaseProblem<GridType>
:: printInitInfo() const
{}


template <class GridType>
inline double PhaseProblem<GridType>
:: stiffSource( const double t, const DomainType& x, const RangeType& u,RangeType& res ) const
{
  double p,Wphi;
  // thermodynamics().pressAndTempEnergyForm(u,p,Wphi);
  res[0]=0;
  res[1]=0;
  res[2]=Wphi;
  return 0.;
}


template <class GridType>
inline double PhaseProblem<GridType>
:: nonStiffSource( const double t, const DomainType& x, const RangeType& u,RangeType& res ) const
{
  abort();

  // time step restriction
  return 0.0;
}



template <class GridType>
inline void PhaseProblem<GridType>
:: evaluate( const double t, const DomainType& arg, RangeType& res ) const 
{
    
  double factor1,factor2,shift;
  factor2=c2_;
  factor2-=c1_;
  shift=factor2;
  shift*=0.5;
 
    { 
      double x=arg[0];
      factor1=smear_;
      
      double r=sqrt(x*x);
      
      res[0]=shift*(tanh((r-0.5)/factor1)+1);
      res[0]+=c1_;  
     
      res[1]=0;
  
      // res[2]=0.5;
      factor1=smear_;
      res[2]=phiscale_*0.5*tanh((r-0.5)/factor1)+0.5;
      res[2]*=res[0];
    
    }


}






template <class GridType>
inline std::string PhaseProblem<GridType>
:: description() const
{
  std::ostringstream stream;
  std::string returnString = stream.str();
  return returnString;
}

} // end namespace Dune
#endif

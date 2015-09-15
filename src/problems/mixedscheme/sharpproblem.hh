#ifndef HEATPROBLEM_HH
#define HEATPROBLEM_HH
#include <dune/common/version.hh>

// dune-fem includes
#include <dune/fem/misc/linesegmentsampler.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/space/common/functionspace.hh>

// local includes
#if SURFACE
#include <dune/phasefield/modelling/thermodsurface.hh>
#else
#include <dune/phasefield/modelling/thermodynamicsbalancedphases.hh>
#endif

#include <dune/fem-dg/models/defaultprobleminterfaces.hh>


namespace Dune {

template <class GridType, class RangeProvider>
class SharpProblem : public EvolutionProblemInterface<
                       Dune::Fem::FunctionSpace< double, double, GridType::dimension,RangeProvider::rangeDim >, true >
{
 
public:
  enum{ dimension = GridType::dimensionworld };
  enum{ dimDomain = dimension };
  enum{ phasefieldId = dimension + 1 };
  enum{ dimRange=RangeProvider::rangeDim};
  typedef Fem::FunctionSpace<typename GridType::ctype, double, GridType::dimensionworld,dimRange > FunctionSpaceType ;


  typedef typename FunctionSpaceType :: DomainFieldType   DomainFieldType;
  typedef typename FunctionSpaceType :: DomainType        DomainType;
  typedef typename FunctionSpaceType :: RangeFieldType    RangeFieldType;
  typedef typename FunctionSpaceType :: RangeType         RangeType;

//   typedef TestThermodynamics ThermodynamicsType;
  typedef BalancedThermodynamics ThermodynamicsType;

  SharpProblem() : 
    myName_( "Mixedtest Heatproblem" ),
    endTime_ ( Fem::Parameter::getValue<double>( "phasefield.endTime",1.0 )), 
    mu_( Fem::Parameter :: getValue< double >( "phasefield.mu1" )),
    delta_(Fem::Parameter::getValue<double>( "phasefield.delta" )),
    A_(Fem::Parameter::getValue<double>("phasefield.A")),
    rhofactor_( Fem::Parameter::getValue<double> ("phasefield.rhofactor")),
    rho1_( Fem::Parameter::getValue<double> ("phasefield.mwpliq")),
    rho2_( Fem::Parameter::getValue<double> ("phasefield.mwpvap")),
    phiscale_(Fem::Parameter::getValue<double> ("phasefield.phiscale")),
    radius_(Fem::Parameter::getValue<double>("phasefield.radius")),
    veloright_(Fem::Parameter::getValue<double>("phasefield.veloright")),
    veloleft_(Fem::Parameter::getValue<double>("phasefield.veloleft")),
    thermodyn_()
    {
    }



  // initialize A and B 
  double init(const bool returnA ) const ;

  // print info 
  void printInitInfo() const;

  // source implementations 
  inline bool hasStiffSource() { return true; }
  inline bool hasNonStiffSource() { return false; }
  // this is the initial data
  inline void evaluate( const DomainType& arg , RangeType& res ) const 
  {
    evaluate( 0.,arg ,res);
  }


  // evaluate function 
  inline void evaluate( const double t, const DomainType& x, RangeType& res ) const;

  // cloned method 
  inline void evaluate( const DomainType& x, const double t, RangeType& res ) const
  {
    evaluate( t, x, res );
  }

  inline double dxr( const DomainType& x) const
  {
    return x[0]/x.two_norm();
  }


  inline double dxdxr( const DomainType& x) const
  {
    return 0; 
  }

  template< class DiscreteFunctionType >
  void finalizeSimulation( DiscreteFunctionType& variablesToOutput,
                           const int eocloop) const
  {}


  const ThermodynamicsType& thermodynamics() const {return thermodyn_;}
  void printmyInfo( std::string filename ) const {}
  inline double endtime() const { return endTime_; }
  inline std::string myName() const { return myName_; }
  void paraview_conv2prim() const {}
  std::string description() const;
 
  inline double mu() const { abort(); return mu_; }
  inline double delta() const{return delta_;}
  protected:
  const std::string myName_;
  const double endTime_;
  const double mu_;
  const double delta_;
  double A_;
  double rhofactor_; 
  double rho1_;
  double rho2_;
  const double phiscale_;
  const double radius_;
  const double veloright_;
  const double veloleft_;
  const ThermodynamicsType thermodyn_;
  
};


template <class GridType,class RangeProvider>
inline double SharpProblem<GridType,RangeProvider>
:: init(const bool returnA ) const 
{
  return 0;
}



template <class GridType,class RangeProvider>
inline void SharpProblem<GridType,RangeProvider>
:: printInitInfo() const
{}

template <class GridType,class RangeProvider>
inline void SharpProblem<GridType,RangeProvider>
:: evaluate( const double t, const DomainType& arg, RangeType& res ) const 
{

  double rhodiff=rho2_-rho1_; 
  double rhomean=0.5*(rho2_+rho1_);
  double deltaInv=rhomean/(delta_*phiscale_);

  double phi=0;
  double rho=rho1_;
  double v=veloleft_;

  if( arg[ 0 ]>0.3 && arg[ 0 ]<0.7 &&  arg[ 1 ]>0.1 && arg[ 1 ]<0.9)
  {
    phi=1;
    rho=rho2_;
    v=veloright_;
  }
  else
  {
    phi=0;
    rho=rho1_;
    veloleft_;
  }
  //rho
  res[0]= rho;
 
  for(int i=1;i<=dimension;i++)
  {
    res[i]=v;
  }
  if(dimension==2)
    res[2]=0;

  
  res[dimension+1]=phi;
 
#if MIXED
  double sigma=0;

#if RHOMODEL     
  //mu
  res[dimension+2]=0.;
  //tau
  res[dimension+3]=0.;
  //sigma_x
  res[dimension+4]=0.;
#else
  //mu
  res[dimension+2]=0.;
  //tau
  res[dimension+3]=0.;
  //sigma_x
  res[dimension+4]=0.;
#endif
#if LAMBDASCHEME
  //lambda_x
  res[dimension+4+dimension]=0;
#endif
  if(dimension==2)
  {
    //sigma_y
    res[dimension+5]=0.;
#if LAMBDASCHEME
    //alpha_x
    res[dimension+7]=0.;
#endif
  }
    
#else
#if NONCONTRANS
#else
    res[dimension+1]*=res[0];
#endif
#endif
}





template <class GridType,class RangeProvider>
inline std::string SharpProblem<GridType,RangeProvider>
:: description() const
{
  std::ostringstream stream;
  std::string returnString = stream.str();
  return returnString;
}

} // end namespace Dune
#endif

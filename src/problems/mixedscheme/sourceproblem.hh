#ifndef SOURCERROBLEM_HH
#define SOURCEPROBLEM_HH
#include <dune/common/version.hh>

// dune-fem includes
#include <dune/fem/misc/linesegmentsampler.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/space/common/functionspace.hh>

// local includes

#include <dune/phasefield/modelling/thermodynamicsbalancedphases.hh>
//#include <dune/phasefield/modelling/thermodynamicsvanderWaals.hh>
//#include <dune/fem/probleminterfaces.hh>

#include <dune/fem-dg/models/defaultprobleminterfaces.hh>


namespace Dune {

template <class GridType, class RangeProvider>
class SourceProblem : public EvolutionProblemInterface<
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

  SourceProblem() : 
    myName_( "Mixedtest Heatproblem" ),
    endTime_ ( Fem::Parameter::getValue<double>( "phasefield.endTime",1.0 )), 
    mu_( Fem::Parameter :: getValue< double >( "phasefield.mu1" )),
    delta_(Fem::Parameter::getValue<double>( "phasefield.delta" )),
    deltaInv_(1./delta_),
    alpha_(Fem::Parameter::getValue<double>( "phasefield.alpha", 0 )),
    beta_(Fem::Parameter::getValue<double>( "phasefield.beta",1 )),
    rho_( Fem::Parameter::getValue<double> ("phasefield.rho0")),
    rho1_( Fem::Parameter::getValue<double> ("phasefield.mwp1")),
    rho2_( Fem::Parameter::getValue<double> ("phasefield.mwp2")),
    lambda_(1.),
    radius_(Fem::Parameter::getValue<double>("phasefield.radius")),
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
#if RHOMODEL
#include "sourceprobCODEGEN/balancedh.cc"
#else
#include "sourceprobCODEGEN/balanced.cc"
#endif

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
  const double deltaInv_;
  const double alpha_;
  const double beta_;
  double rho_; 
  double rho1_;
  double rho2_;
  const double lambda_; 
  const double radius_;
  const ThermodynamicsType thermodyn_;
  
};


template <class GridType,class RangeProvider>
inline double SourceProblem<GridType,RangeProvider>
:: init(const bool returnA ) const 
{
  return 0;
}



template <class GridType,class RangeProvider>
inline void SourceProblem<GridType,RangeProvider>
:: printInitInfo() const
{}

template <class GridType,class RangeProvider>
inline void SourceProblem<GridType,RangeProvider>
:: evaluate( const double t, const DomainType& arg, RangeType& res ) const 
{
 double  velo_=0.; 
  double x=arg[0];
  double ksi=x-velo_*t;
  
  double phi=0.5*tanh(ksi*lambda_*deltaInv_)+0.5;
  
  //rhoo
  res[0]=rhosol(ksi);
  //v
  res[1]=velo_;
  //phi
  res[2]=phi;
  //mu
  res[3]=musol(ksi)+0.5*velo_*velo_;
  //theta
  res[4]=thetasol(ksi);
  //sigma
  res[5]=gradphi(ksi);
#if LAMBDASCHEME
  res[6]=res[5]/res[0];
#endif
}



template <class GridType,class RangeProvider>
inline std::string SourceProblem<GridType,RangeProvider>
:: description() const
{
  std::ostringstream stream;
  std::string returnString = stream.str();
  return returnString;
}

} // end namespace Dune
#endif

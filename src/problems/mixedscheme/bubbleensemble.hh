#ifndef PHASEFIELD_BUBBLEENSEMBLE_HH
#define PHASEFIELD_BUBBLEENSEMBLE_HH
#include <dune/common/version.hh>

// dune-fem includes
#include <dune/fem/io/parameter.hh>
#include <dune/fem/space/common/functionspace.hh>

// local includes

#include <dune/phasefield/modelling/thermodsurface.hh>

#include <dune/fem-dg/models/defaultprobleminterfaces.hh>


namespace Dune {

template <class GridType, class RangeProvider>
class BubbleEnsemble : public EvolutionProblemInterface<
                       Dune::Fem::FunctionSpace< double, double, GridType::dimension,RangeProvider::rangeDim >, true >
{
 
public:
  constexpr static int dimension = GridType::dimensionworld;
  constexpr static int phasefieldId = dimension + 1 ;
  constexpr static int dimRange=RangeProvider::rangeDim;
  
  using FunctionSpaceType = typename Fem::FunctionSpace<typename GridType::ctype,
                                                        double,
                                                        GridType::dimensionworld,
                                                        dimRange >;

  using DomainFieldType = typename FunctionSpaceType :: DomainFieldType;
  using DomainType      = typename FunctionSpaceType :: DomainType;
  using RangeFieldType  = typename FunctionSpaceType :: RangeFieldType;
  using RangeType       = typename FunctionSpaceType :: RangeType;

  using ThermodynamicsType = BalancedThermodynamics;

 BubbleEnsemble() : 
    myName_( "BubbleProblem" ),
    endTime_ ( Fem::Parameter::getValue<double>( "phasefield.endTime",1.0 )), 
    delta_(Fem::Parameter::getValue<double>( "phasefield.delta" )),
    A_(Fem::Parameter::getValue<double>("phasefield.A")),
    rho_( Fem::Parameter::getValue<double> ("phasefield.rho0")),
    mwpliq_( Fem::Parameter::getValue<double> ("phasefield.mwpliq")),
    mwpvap_( Fem::Parameter::getValue<double> ("phasefield.mwpvap")),
    phiscale_(Fem::Parameter::getValue<double> ("phasefield.phiscale")),
    rhofactor_(Fem::Parameter::getValue<double>("phasefield.rhofactor")),
    bubblefilename_(Fem::Parameter::getValue<std::string>("phasefield.bubbles")),
    bubblevector_(0),
    thermodyn_()
    {
      distortion_[0]=1;
      distortion_[1]=1;

      readDataFile( bubblefilename_ );
    }      

  void readDataFile(std::string filename)
  {
    std::ifstream bubblefile;
    std::string tmp;
    bubblefile.open(filename.c_str());
    double a;
     
    if(! bubblefile )
    {
      std::cout<<"There was a problem opening the file\n";
      abort();
    }
    else
    {
      while(  bubblefile >> a )
      {
        std::cout<<a<<"\n";
        bubblevector_.push_back(a);
      }

      numBubbles_=bubblevector_.size()/(dimension+1);
      bubblefile.close();
    }
  }




  // initialize A and B 
  double init(const bool returnA ) const ;

  // print info 
  void printInitInfo() const;
#if THERMO == 1
#include <dune/phasefield/modelling/CoquelTaylorSources/coquelTaylorRho.cc>
#elif THERMO == 2
#include <dune/phasefield/modelling/CoquelTaylorSources/realRho.cc>
#elif THERMO == 3
#include <dune/phasefield/modelling/PhasefieldvanderWaalsSources/pfvdWaalRho.cc>
#endif
// source implementations 
  inline bool hasStiffSource() { return true; }
  inline bool hasNonStiffSource() { return false; }
  // this is the initial data
  inline void evaluate( const DomainType& arg , RangeType& res ) const 
  {
    evaluate( 0.,arg ,res);
  }

  
  inline double mynorm( const DomainType& x) const
  {
    double res(0.);
    for( int ii=0 ; ii < dimension ; ++ii)
      res+=distortion_[ii]*x[ii]*x[ii];
    return std::sqrt(res);
  }


  // evaluate function 
  inline void evaluate( const double t, const DomainType& x, RangeType& res ) const;

  // cloned method 
  inline void evaluate( const DomainType& x, const double t, RangeType& res ) const
  {
    evaluate( t, x, res );
  }

  inline double dxdi( const DomainType& x , const int i) const
  {
    return distortion_[i]*x[i]/mynorm(x);
  }

  inline double dxdidi( const DomainType& x , const int i) const
  {
    double norm=mynorm(x);
    double a=distortion_[i];
    double res=1-a*x[i]*x[i]/(norm*norm);
    res*=a;
    res/=norm;

    return res;
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
 
  inline double delta() const{return delta_;}
  protected:
  const std::string myName_;
  const double endTime_;
  const double delta_;
  const double A_;
  double rho_;
  double mwpliq_;
  double mwpvap_;
  const double phiscale_;
  const double rhofactor_;
  std::string bubblefilename_; 
  std::vector<double> bubblevector_;
  const ThermodynamicsType thermodyn_;
  DomainType distortion_;
  int numBubbles_;
};


template <class GridType,class RangeProvider>
inline double BubbleEnsemble<GridType,RangeProvider>
:: init(const bool returnA ) const 
{
  return 0;
}



template <class GridType,class RangeProvider>
inline void BubbleEnsemble<GridType,RangeProvider>
:: printInitInfo() const
{}

template <class GridType,class RangeProvider>
inline void BubbleEnsemble<GridType,RangeProvider>
:: evaluate( const double t, const DomainType& arg, RangeType& res ) const 
{
  double phi; 
  double width=1/delta_;
  width*=phiscale_;
  //Outside: phi=1
  res[dimension+1]=1;
  // double rho= rhofactor_*(tanh(-50*(arg[0]-0.19))+1)+mwpliq(0.);
  double rho=mwpliq(0.);

#if MIXED 
  for(int ii = 0; ii < dimension ; ++ii)
    res[dimension+4+ii]=0.;

  double dFdphi= 0;
  double dFdrho=thermodyn_.chemicalPotential(rho,phi, rho);
  //mu
  res[dimension+2]=dFdrho;
  //tau
  res[dimension+3]=dFdphi;
#endif 
  const int offset=dimension+1;
  double bubblephi;
  double laplace(0.);
  for(size_t i=0 ; i<numBubbles_ ;++i)
  { 

    DomainType center;
    double mirror;
    for(int ii = 0 ; ii <dimension ; ++ ii)
      center[ii]=bubblevector_[i*offset+ii];


    double radius=bubblevector_[dimension+(i*offset)];

    DomainType vector=center;
    vector-=arg;
    double r=mynorm(vector);
    
    double tanr=tanh((radius-r) * (width ));
    double tanhr=tanr;
    double dtanhr=1-tanhr*tanhr; 
    double ddtanhr=-2*tanhr*dtanhr;
    bubblephi=0.5*(tanr+1);   
    res[dimension+1]-=bubblephi;

    for( int ii=0 ; ii< dimension ;++ii)
    {
      res[dimension+4+ii]+=-0.5*dtanhr*(width)*dxdi(vector,ii);
      laplace+=dxdidi(vector,ii);
    }

    laplace*=thermodyn_.delta()*ddtanhr*(width*width);

    continue;
  }
  phi=res[dimension+1];
  double v=0;

  //rho
  rho=evalRho(phi);
  rho+=rhofactor_*(tanh(-50*(arg[0]-0.19))+1);

  res[0]=rho;
  res[dimension+2]=thermodyn_.chemicalPotential(rho,phi,rho); 
  res[dimension+3]=thermodyn_.reactionSource(rho,phi,phi)+laplace;

  for(int i=1;i<=dimension;i++)
  {
#if MIXED
     res[i]=v;
#else
     res[i]=v*rho;
#endif
   }

#if MIXED || NONCONTRANS  
#else
  res[dimension+1]*=rho;
#endif
 
}


template <class GridType,class RangeProvider>
inline std::string BubbleEnsemble<GridType,RangeProvider>
:: description() const
{
  std::ostringstream stream;
  std::string returnString = stream.str();
  return returnString;
}

////////////////////////////////////////////////////////////////////////////////////
//                      WRAPPERS FOR SPLIT OPERATOR                               //
////////////////////////////////////////////////////////////////////////////////////

template <class GridType, class RangeProvider>
class NvStBubbleEnsemble:
private BubbleEnsemble< GridType, RangeProvider>
{
  using BaseType = BubbleEnsemble< GridType, RangeProvider> ;
  using BaseType::fixedTimeFunction;
  constexpr static int dimension = GridType::dimensionworld;
  constexpr static int dimDomain = dimension;
  constexpr static int dimRange  = RangeProvider::rangeDim/2;
  using FunctionSpaceType = typename Fem::FunctionSpace<typename GridType::ctype,
                                                        double,
                                                        GridType::dimensionworld,
                                                        dimRange >;

  using DomainFieldType = typename FunctionSpaceType :: DomainFieldType;
  using DomainType      = typename FunctionSpaceType :: DomainType;
  using RangeFieldType  = typename FunctionSpaceType :: RangeFieldType;
  using RangeType       = typename FunctionSpaceType :: RangeType;

  public:
  NvStBubbleEnsemble():
  BaseType()
  {}
  template< class DomainType , class RangeType>
  void evaluate( const DomainType& arg, RangeType& res) const
  {
    evaluate(0,arg,res);
  }


  template< class DomainType , class RangeType>
  void evaluate( const double time, const DomainType& arg, RangeType& res ) const
  {
    typename BaseType::RangeType resfull;
    BaseType::evaluate( time , arg , resfull );
    res[0]=resfull[0];
    for( int ii = 0 ; ii<dimension ; ++ii)
      res[1+ii]=resfull[1+ii];
    res[dimension+1]=resfull[dimension+2];
  }

};

template <class GridType, class RangeProvider>
class AcBubbleEnsemble :
private BubbleEnsemble< GridType, RangeProvider>
{
  using BaseType = BubbleEnsemble< GridType, RangeProvider> ;
  using BaseType::fixedTimeFunction;
  constexpr static int  dimension = GridType::dimensionworld;
  constexpr static int dimDomain = dimension;
  constexpr static int dimRange=RangeProvider::rangeDim/2;
  using FunctionSpaceType = typename Fem::FunctionSpace<typename GridType::ctype,
                                                        double,
                                                        GridType::dimensionworld,
                                                        dimRange >;

  using DomainFieldType = typename FunctionSpaceType :: DomainFieldType;
  using DomainType      = typename FunctionSpaceType :: DomainType;
  using RangeFieldType  = typename FunctionSpaceType :: RangeFieldType;
  using RangeType       = typename FunctionSpaceType :: RangeType;

  public:
  AcBubbleEnsemble():
  BaseType()
  {}
  
  template< class DomainType , class RangeType>
  void evaluate( const DomainType& arg, RangeType& res) const
  {
    evaluate(0,arg,res);
  } 
  
  template< class DomainType , class RangeType>
  void evaluate( const double time, const DomainType& arg, RangeType& res) const
  {
    typename BaseType::RangeType resfull;
    BaseType::evaluate( time , arg , resfull );
    res[0]=resfull[dimension+1];
    res[1]=resfull[dimension+3];
    for( int ii = 0 ; ii < dimension ; ++ii)
     res[2+ii]=resfull[dimension+4+ii];
  }

};

} // end namespace Dune
#endif

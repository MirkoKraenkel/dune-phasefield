#ifndef PHASEFIELD_BUBBLEENSEMBLE_HH
#define PHASEFIELD_BUBBLEENSEMBLE_HH
#include <dune/common/version.hh>

// dune-fem includes
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
    rho1_( Fem::Parameter::getValue<double> ("phasefield.mwp1")),
    rho2_( Fem::Parameter::getValue<double> ("phasefield.mwp2")),
    phiscale_(Fem::Parameter::getValue<double> ("phasefield.phiscale")),
    bubblefilename_(Fem::Parameter::getValue<std::string>("phasefield.bubbles")),
    bubblevector_(0),
    thermodyn_()
    {
      distortion_[0]=1;
      distortion_[1]=1.5;
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
 //     return x[i]/x.two_norm();
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
  double rho1_;
  double rho2_;
  const double phiscale_;
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
  double width=6*delta_;
  double width2=delta_;
#if SURFACE
#else
  width/=sqrt(A_);
#endif
  phi=0.5-phiscale_*0.5;
  double rho=rho1_;
#if MIXED 
  for(int ii = 0; ii < dimension ; ++ii)
    res[dimension+4+ii]=0.;

  double dFdphi= 0;
  double dFdrho=0;
 //mu
 res[dimension+2]=dFdrho;
  //tau
 res[dimension+3]=dFdphi;
#endif 
  const int offset=dimension+1;

  for(size_t i=0 ; i<numBubbles_ ;++i)
  { 
    DomainType center;
    for(int ii = 0 ; ii <dimension ; ++ ii)
      center[ii]=bubblevector_[i*offset+ii];
    double radius=bubblevector_[dimension+(i*offset)];
  
    DomainType vector=center;
    vector-=arg;
    double r=mynorm(vector);
    //double r=vector.two_norm();
    
    double tanr=tanh((radius-r) * (1 /width2 ));
    double tanhr=tanr;
    double dtanr=1+tanr*tanr;
    double dtanhr=1-tanhr*tanhr; 
    double ddtanr=2*tanr*dtanr;
    double ddtanhr=-2*tanhr*dtanhr;
 
  
    if( r < radius+(0.5*width))
      {
        if( r < radius-(0.5*width))
          {//Inside bubble
            phi=0.5+0.5*phiscale_;
            rho=rho2_;
#if MIXED
            //mu
            res[dimension+2]=dFdrho;
            //tau
            res[dimension+3]=dFdphi;
#endif 
            continue;
          }
        else
          {
            phi=phiscale_*0.5*( tanhr )+0.5;
            double rhodiff=rho1_-rho2_;
            rho=(rhodiff)*(-0.5*tanhr+0.5)+rho2_;
#if MIXED
            for( int ii=0 ; ii< dimension ;++ii)
              res[dimension+4+ii]=0.5*dtanhr*(1/width2)*dxdi(vector,ii);
            //mu
            res[dimension+2]=0;
            //tau
            res[dimension+3]=0;
#endif
            continue;
          }
        }
      }
      
  double v=0;
  //rho
  res[0]= rho;
   
   for(int i=1;i<=dimension;i++)
   {
     res[i]=v;
   }
  res[dimension+1]=phi;

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

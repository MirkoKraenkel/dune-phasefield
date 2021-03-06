#ifndef PHASEFIELD_HANDLEPROBLEM_HH
#define PHASEFIELD_HANDLEPROBLEM_HH
#include <dune/common/version.hh>

// dune-fem includes
#include <dune/fem/io/parameter.hh>
#include <dune/fem/space/common/functionspace.hh>

// local includes
#include <dune/phasefield/modelling/thermodsurface.hh>
//#include <dune/phasefield/modelling/thermodscaled.hh>

#include <dune/fem-dg/models/defaultprobleminterfaces.hh>


namespace Dune {

template <class GridType, class RangeProvider>
class HandleProblem : public EvolutionProblemInterface<
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

  HandleProblem() : 
    myName_( "HandleProblem" ),
    endTime_ ( Fem::Parameter::getValue<double>( "phasefield.endTime",1.0 )), 
    delta_(Fem::Parameter::getValue<double>( "phasefield.delta" )),
    A_(Fem::Parameter::getValue<double>("phasefield.A")),
    rhofactor_( Fem::Parameter::getValue<double> ("phasefield.rhofactor")),
    rholiq_( Fem::Parameter::getValue<double> ("phasefield.mwpliq")),
    rhovap_( Fem::Parameter::getValue<double> ("phasefield.mwpvap")),
    phiscale_(Fem::Parameter::getValue<double> ("phasefield.phiscale")),
    pos1_(Fem::Parameter::getValue<double>("phasefield.p1")),
    pos2_(Fem::Parameter::getValue<double>("phasefield.p2")),
    length_(Fem::Parameter::getValue<double>("phasefield.k")),
    radius1_(Fem::Parameter::getValue<double>("phasefield.r1")),
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
#if THERMO == 1
#include <dune/phasefield/modelling/CoquelTaylorSources/coquelTaylorRho.cc>
#elif THERMO == 2
#include <dune/phasefield/modelling/CoquelTaylorSources/realRho.cc>
#elif THERMO == 3
#include <dune/phasefield/modelling/PhasefieldvanderWaalsSources/phasefieldvanderWaalsRho.cc>
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
 
  inline double delta() const{return delta_;}
  protected:
  const std::string myName_;
  const double endTime_;
  const double delta_;
  double A_;
  double rhofactor_; 
  double rholiq_;
  double rhovap_;
  const double phiscale_;
  const double pos1_;
  const double pos2_;
  const double length_;
  const double radius1_;
  const ThermodynamicsType thermodyn_;
  
};






template <class GridType,class RangeProvider>
inline double HandleProblem<GridType,RangeProvider>
:: init(const bool returnA ) const 
{
  return 0;
}



template <class GridType,class RangeProvider>
inline void HandleProblem<GridType,RangeProvider>
:: printInitInfo() const
{}

template <class GridType,class RangeProvider>
inline void HandleProblem<GridType,RangeProvider>
:: evaluate( const double t, const DomainType& arg, RangeType& res ) const 
{

  double rhodiff=mwpliq(0.)-mwpvap(0.);
  const double rhomean=1;//0.5*(rhovap_+rholiq_);

  double deltaInv=1./delta_;
  double k=length_;
  double r,r1,r2,r3;
  double z=sqrt(radius1_*radius1_ - std::pow( pos1_-(0.5+0.5*k),2));

  double radius2 = sqrt( z*z + std::pow( (0.5-0.5*k)-pos2_, 2));
  double myradius;
  // drop1
  DomainType vector{ 0.5,pos2_ };
  vector-=arg;
  r1=vector.two_norm();

  DomainType vector2{ 0.5,pos1_ };
  vector2-=arg;
  r2=vector2.two_norm();

  r3=std::abs( arg[0]-0.5 );
  if( std::abs(arg[1]-0.5)<0.2 )
    r=std::min(r1-radius1_,std::min( r2-radius2,r3-z));
  else
    r=std::min(r1-radius1_, r2-radius2);


  double tanhr=tanh( -( r  )*deltaInv );
  double tanhrho=-tanh( -( r )*rhofactor_*deltaInv );

  //(1-tanh^2)/delta
  double drtanhr=deltaInv*(1-tanhr*tanhr);
  double drdrtanhr=-tanhr*drtanhr;
  //-2/delta^2 *(tanh*(1-tanh^2)
  drdrtanhr*=-2*deltaInv;
  double phi=-0.5*phiscale_*tanhr+0.5;

  double rho=evalRho(phi);
  //double rho=(rhodiff)*(0.5*-tanhrho+0.5)+rhovap_;

  //rho
  res[0]= rho;

  

  res[dimension+1]=phi;

#if MIXED
  double dFdphi= thermodyn_.reactionSource(rho,phi,phi);
  double dFdrho= thermodyn_.chemicalPotential(rho,phi,rho);
  double sigma=0.5*phiscale_*drtanhr;
  double laplacePhi=0.5*drdrtanhr;
  if(arg[0]<0.)
    {
      sigma*=-1;
    }

#if RHOMODEL     
  double gradrho=-0.5*rhodiff*drtanhr;
  double rhoInv=thermodyn_.h2(rho);            //1./rho;
  double hprime=thermodyn_.h2prime(rho);//-1./(rho*rho);
  //mu
  res[dimension+2]=0;//0.5*v*v+dFdrho-delta_*rhoInv *0.5*sigma*sigma;
  //tau
  res[dimension+3]=0;//dFdphi-delta_*rhoInv*laplacePhi+hprime*gradrho*sigma;
  //sigma_x
  res[dimension+4]=sigma;
#else
  //mu
  res[dimension+2]=dFdrho;
  //tau
  res[dimension+3]=dFdphi-delta_*A_*laplacePhi;
  //sigma_x
  res[dimension+4]=sigma;
#endif
#if LAMBDASCHEME
  //lambda_x
  res[dimension+4+dimension]=res[dimension+4]*rhoInv;
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
  for(int i=1;i<=dimension;i++)
    {
      res[i]*=rho;
    }
#if NONCONTRANS
#else
    res[dimension+1]*=res[0];
#endif
#endif
}


template <class GridType,class RangeProvider>
inline std::string HandleProblem<GridType,RangeProvider>
:: description() const
{
  std::ostringstream stream;
  stream<<myName_;
  std::string returnString = stream.str();
  return returnString;
}

////////////////////////////////////////////////////////////////////////////////////
//                      WRAPPERS FOR SPLIT OPERATOR                               //
////////////////////////////////////////////////////////////////////////////////////

template <class GridType, class RangeProvider>
class NvStHandleProblem:
private HandleProblem< GridType, RangeProvider>
{
  using BaseType=HandleProblem< GridType, RangeProvider> ;
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
  NvStHandleProblem():
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
class AcHandleProblem:
private HandleProblem< GridType, RangeProvider>
{
  using BaseType=HandleProblem< GridType, RangeProvider> ;
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
  AcHandleProblem():
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

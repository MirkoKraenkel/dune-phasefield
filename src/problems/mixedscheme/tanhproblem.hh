#ifndef PHASEFIELD_TANHPROBLEM_HH
#define PHASEFIELD_TANHPROBLEM_HH
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
class TanhProblem : public EvolutionProblemInterface<
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

  TanhProblem() : 
    myName_( "TanhProblem" ),
    endTime_ ( Fem::Parameter::getValue<double>( "phasefield.endTime",1.0 )), 
    delta_(Fem::Parameter::getValue<double>( "phasefield.delta" )),
    A_(Fem::Parameter::getValue<double>("phasefield.A")),
    rhofactor_( Fem::Parameter::getValue<double> ("phasefield.rhofactor")),
    rho1_( Fem::Parameter::getValue<double> ("phasefield.mwp1")),
    rho2_( Fem::Parameter::getValue<double> ("phasefield.mwp2")),
    phiscale_(Fem::Parameter::getValue<double> ("phasefield.phiscale")),
    radius_(Fem::Parameter::getValue<double>("phasefield.radius")),
    shift_(Fem::Parameter::getValue<double>("phasefield.shift",0.)),
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
  double rho1_;
  double rho2_;
  const double phiscale_;
  const double radius_;
  const double shift_;
  const double veloright_;
  const double veloleft_;
  const ThermodynamicsType thermodyn_;
  
};






template <class GridType,class RangeProvider>
inline double TanhProblem<GridType,RangeProvider>
:: init(const bool returnA ) const 
{
  return 0;
}



template <class GridType,class RangeProvider>
inline void TanhProblem<GridType,RangeProvider>
:: printInitInfo() const
{}

template <class GridType,class RangeProvider>
inline void TanhProblem<GridType,RangeProvider>
:: evaluate( const double t, const DomainType& arg, RangeType& res ) const 
{

  double rhodiff=rho2_-rho1_;
#if UNBALMODEL
  const double rhomean=1;//0.5*(rho2_+rho1_);
#else
  const double rhomean=1;
#endif
#if SURFACE
  double deltaInv=rhomean/(delta_*phiscale_);
#else  
  double deltaInv=sqrt(A_)/(delta_*phiscale_);
#endif  
  double r=std::abs( arg[0]);
  double tanhr=-tanh( ( r-radius_ )*deltaInv ); 
  double tanhrho=-tanh( ( r-(radius_+shift_))*rhofactor_*deltaInv ); 

  double velopos=r-(radius_+(delta_/rhomean));
  double peak=exp(-velopos*velopos*deltaInv*deltaInv);

  //(1-tanh^2)/delta
  double drtanhr=-deltaInv*(1-tanhr*tanhr);
  double drdrtanhr=tanhr*drtanhr;
  //-2/delta^2 *(tanh*(1-tanh^2)
  drdrtanhr*=-2*deltaInv;
  double phi=-0.5*tanhr+0.5;

  

  //double rho=thermodyn_.evalRho(phi);
  double rho=(rhodiff)*(-0.5*tanhrho+0.5)+rho1_;
  
  double v=0;//(veloright_-veloleft_)*(peak)+veloleft_;
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
  double dFdphi= thermodyn_.reactionSource(rho,phi); 
  double dFdrho=thermodyn_.chemicalPotential(rho, phi);
  double sigma=-0.5*drtanhr;
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
  res[dimension+2]=0.5*v*v+dFdrho-delta_*rhoInv *0.5*sigma*sigma;
  //tau
  res[dimension+3]=dFdphi-delta_*rhoInv*laplacePhi+hprime*gradrho*sigma;
  //sigma_x
  res[dimension+4]=sigma;
#else
  //mu
  res[dimension+2]=0.5*v*v+dFdrho;
  //tau
  res[dimension+3]=dFdphi-delta_*laplacePhi;
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
#if NONCONTRANS
#else
    res[dimension+1]*=res[0];
#endif
#endif
}


template <class GridType,class RangeProvider>
inline std::string TanhProblem<GridType,RangeProvider>
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
class NvStTanhProblem:
private TanhProblem< GridType, RangeProvider>
{
  using BaseType=TanhProblem< GridType, RangeProvider> ;
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
  NvStTanhProblem():
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
class AcTanhProblem:
private TanhProblem< GridType, RangeProvider>
{
  using BaseType=TanhProblem< GridType, RangeProvider> ;
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
  AcTanhProblem():
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

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
class TanhProblem : public EvolutionProblemInterface<
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

  TanhProblem() : 
    myName_( "Mixedtest Heatproblem" ),
    endTime_ ( Fem::Parameter::getValue<double>( "phasefield.endTime",1.0 )), 
    mu_( Fem::Parameter :: getValue< double >( "phasefield.mu1" )),
    delta_(Fem::Parameter::getValue<double>( "phasefield.delta" )),
    A_(Fem::Parameter::getValue<double>("phasefield.A")),
    rhofactor_( Fem::Parameter::getValue<double> ("phasefield.rhofactor")),
    rho1_( Fem::Parameter::getValue<double> ("phasefield.mwp1")),
    rho2_( Fem::Parameter::getValue<double> ("phasefield.mwp2")),
    phiscale_(Fem::Parameter::getValue<double> ("phiscale")),
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
#if SURFACE
  double deltaInv=1/delta_;
#else  
  double deltaInv=sqrt(A_)/delta_;
#endif  
  double r=std::abs( arg[0]);
  double tanhr=-tanh( ( r-radius_ )*deltaInv ); 
  double tanhrho=-tanh( ( r-radius_ )*rhofactor_*deltaInv ); 

  //(1-tanh^2)/delta
  double drtanhr=-deltaInv*(1-tanhr*tanhr);
  double drdrtanhr=tanhr*drtanhr;
  //-2/delta^2 *(tanh*(1-tanh^2)
  drdrtanhr*=-2*deltaInv;
    double phi=0.5*tanhr+0.5;
 
  double rhodiff=rho2_-rho1_;
  //double rho=thermodyn_.evalRho(phi);
  double rho=(rhodiff)*(-0.5*tanhrho+0.5)+rho1_;
  
  double v=0;//0.1*sin(2*M_PI*arg[0]);
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
  double sigma=0.5*drtanhr;
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

} // end namespace Dune
#endif

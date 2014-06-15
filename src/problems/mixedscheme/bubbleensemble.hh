#ifndef HEATPROBLEM_HH
#define HEATPROBLEM_HH
#include <dune/common/version.hh>

// dune-fem includes
#include <dune/fem/misc/linesegmentsampler.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/space/common/functionspace.hh>

// local includes

#include <dune/phasefield/modelling/thermodynamicsbalancedphases.hh>

//#include <dune/fem/probleminterfaces.hh>

#include <dune/fem-dg/models/defaultprobleminterfaces.hh>


namespace Dune {

template <class GridType, class RangeProvider>
class HeatProblem : public EvolutionProblemInterface<
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

  HeatProblem() : 
    myName_( "Mixedtest Heatproblem" ),
    endTime_ ( Fem::Parameter::getValue<double>( "phasefield.endTime",1.0 )), 
    mu_( Fem::Parameter :: getValue< double >( "phasefield.mu1" )),
    delta_(Fem::Parameter::getValue<double>( "phasefield.delta" )),
    rho_( Fem::Parameter::getValue<double> ("phasefield.rho0")),
    phiscale_(Fem::Parameter::getValue<double> ("phiscale")),
    bubblefilename_(Fem::Parameter::getValue<std::string>("phasefield.bubbles")),
    bubblevector_(0),
    thermodyn_()
    {
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
     
      std::cout<<"size="<<bubblevector_.size()<<"\n";
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

  inline double dyr( const DomainType& x) const
  {
    return x[1]/x.two_norm();
  }

  inline double dxdxr( const DomainType& x) const
  {
    return pow( x[1], 2)*pow(x.two_norm(),-3./2.);
  }

 inline double dydyr( const DomainType& x) const
  {
    return pow( x[0], 2)*pow(x.two_norm(),-3./2.);
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
  double rho_;
  const double phiscale_;
  std::string bubblefilename_; 
  std::vector<double> bubblevector_;
  const ThermodynamicsType thermodyn_;
  
};


template <class GridType,class RangeProvider>
inline double HeatProblem<GridType,RangeProvider>
:: init(const bool returnA ) const 
{
  return 0;
}



template <class GridType,class RangeProvider>
inline void HeatProblem<GridType,RangeProvider>
:: printInitInfo() const
{}

template <class GridType,class RangeProvider>
inline void HeatProblem<GridType,RangeProvider>
:: evaluate( const double t, const DomainType& arg, RangeType& res ) const 
{
  double deltaInv=1./delta_;
  double phi; 
  double width=6*delta_;
  phi=0.;
  res[dimension+4]=0.;
  res[dimension+5]=0.;
 
  for( int i=0 ; i<(bubblevector_.size())/3 ;++i)
  { 

    
  DomainType center;
  center[0]=bubblevector_[i*3];
 
  center[1]=bubblevector_[1+(i*3)];
  double radius=bubblevector_[2+(i*3)];
  
  DomainType vector=center;
  vector-=arg;
  double r=vector.two_norm();
  double tanr=tan((radius-r) * ( M_PI / width ));
  double tanhr=tanh(tanr);
  double dtanr=1+tanr*tanr;
  double dtanhr=1-tanhr*tanhr; 
  if( r < radius+(0.5*width))
   {
      if( r < radius-(0.5*width))
      {
        phi=1.;
        continue;
      }
      else
      {
        phi=0.5*( tanhr )+0.5;
        res[dimension+4]=dtanhr*dtanr*(M_PI/width)*dxr(arg);
        res[dimension+5]=dtanhr*dtanr*(M_PI/width)*dyr(arg);
         continue;
      }
    }
 
  }
      
   

  double rho=1.;

  double v=0;
  //rho
  res[0]= rho;
 
   for(int i=1;i<=dimension;i++)
   {
     res[i]=v;
   }

   
  res[dimension+1]=phi;
  
  double laplacePhi=0;

  double dFdphi= thermodyn_.reactionSource(rho,phi); 
  double dFdrho=thermodyn_.chemicalPotential(rho, phi);
     
   //mu
  res[dimension+2]=0;
  //tau
  res[dimension+3]=0;
  
    //sigma
 }





template <class GridType,class RangeProvider>
inline std::string HeatProblem<GridType,RangeProvider>
:: description() const
{
  std::ostringstream stream;
  std::string returnString = stream.str();
  return returnString;
}

} // end namespace Dune
#endif

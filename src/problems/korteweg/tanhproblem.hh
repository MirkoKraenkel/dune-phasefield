#ifndef TANHBUBBLEPROBLEM_HH
#define TANHBUBBLEPROBLEM_HH
#include <dune/common/version.hh>

// dune-fem includes
#include <dune/fem/misc/linesegmentsampler.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/space/common/functionspace.hh>

// local includes

#include <dune/phasefield/modelling/thermodynamicsKorteweg.hh>


#include <dune/fem-dg/models/defaultprobleminterfaces.hh>


namespace Dune {

template <class GridType,class RangeProvider >
class TanhBubbleProblem : public EvolutionProblemInterface<
                        Dune::Fem::FunctionSpace< double, double, GridType::dimension, RangeProvider::rangeDim >, true >
{
 
public:
  enum{ dimension = GridType::dimensionworld };
  enum{ dimDomain = dimension };
  enum{ dimRange = RangeProvider::rangeDim };
  typedef Fem::FunctionSpace<typename GridType::ctype, double, GridType::dimensionworld,dimRange > FunctionSpaceType ;
  


  typedef typename FunctionSpaceType :: DomainFieldType   DomainFieldType;
  typedef typename FunctionSpaceType :: DomainType        DomainType;
  typedef typename FunctionSpaceType :: RangeFieldType    RangeFieldType;
  typedef typename FunctionSpaceType :: RangeType         RangeType;

//   typedef TestThermodynamics ThermodynamicsType;
  typedef BalancedThermodynamics ThermodynamicsType; 

  TanhBubbleProblem() : 
    myName_( "Mixedtest Heatproblem" ),
    endTime_ ( Fem::Parameter::getValue<double>( "phasefield.endTime",1.0 ) ), 
    mu_( Fem::Parameter :: getValue< double >( "phasefield.mu1" ) ),
    delta_(Fem::Parameter::getValue<double>( "phasefield.delta" ) ) ,
    rho1_( Fem::Parameter::getValue<double> ("korteweg.mwp1" ) ),
    rho2_( Fem::Parameter::getValue<double> ("korteweg.mwp2") ),
    phiscale_(Fem::Parameter::getValue<double> ("phiscale") ),
    radius_(Fem::Parameter::getValue<double>("phasefield.radius") ),
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
  double rho1_;
  double rho2_;
  const double phiscale_;
  const double radius_;
  const ThermodynamicsType thermodyn_;
  
};


template < class GridType , class RangeProvider >
inline double TanhBubbleProblem<GridType, RangeProvider>
:: init(const bool returnA ) const 
{
  return 0;
}



template <class GridType ,class RangeProvider>
inline void TanhBubbleProblem<GridType , RangeProvider >
:: printInitInfo() const
{}

template <class GridType , class RangeProvider >
inline void TanhBubbleProblem<GridType, RangeProvider >
:: evaluate( const double t, const DomainType& arg, RangeType& res ) const 
{
  double deltaInv=sqrt(1./delta_);
  double r= std::abs(arg[0]);
  double tanhr=tanh( ( r-0.5 )*deltaInv ); 
  //(1-tanh^2)/delta
  double drtanhr=deltaInv*(1-tanhr*tanhr);
  double drdrtanhr=tanhr*drtanhr;
  //-2/delta^2 *(tanh*(1-tanh^2)
  drdrtanhr*=-2*deltaInv;
 
  // double rho=1.85*(-0.5*tanhr+0.5)+
  double rho=(rho2_-rho1_)*(-0.5*tanhr+0.5)+rho1_;
  
  double v=0;//0.1*sin(2*M_PI*arg[0]);
  //rho
  res[0]= rho;
 
   for(int i=1;i<=dimension;i++)
   {
     res[i]=v;
   }
   if(dimension==2)
     res[2]=0;
   
 
#if MIXED
  double dFdrho=thermodyn_.chemicalPotential(rho);
  double sigma=0.5*drtanhr;
  double gradrho=-0.5*(rho2_-rho1_)*drtanhr;
  double laplacePhi=0.5*drdrtanhr;
  if(arg[0]>0.5)
     laplacePhi*=-1;



  //mu
  res[dimension+2]=0.5*v*v+dFdrho;
  //sigma_x
  res[dimension+3]=sigma;
  if(dimension==2)
    {
      //sigma_y
      res[dimension+4]=0.;
    }
    
#endif
}





template <class GridType , class RangeProvider>
inline std::string TanhBubbleProblem<GridType , RangeProvider>
:: description() const
{
  std::ostringstream stream;
  std::string returnString = stream.str();
  return returnString;
}

} // end namespace Dune
#endif

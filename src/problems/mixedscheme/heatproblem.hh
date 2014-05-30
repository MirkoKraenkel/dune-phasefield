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
 
  inline double mu() const { abort(); return mu_; }
  inline double delta() const{return delta_;}
  protected:
  const std::string myName_;
  const double endTime_;
  const double mu_;
  const double delta_;
  double rho_;
  const double phiscale_;
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
  double x=arg[0];
  double cost=cos(M_PI*t);
  double cosx=cos(2*M_PI*x);
  double sinx=sin(2*M_PI*x);
    
  double rho=rho_;
   //double v=sinx*cost;
   double v=0;
   //rho
   res[0]= rho;
   //v
   for(int i=1;i<=dimension;i++)
   {
     res[i]=v;
   }
   if(dimension==2)
     res[2]=0;
   double  phi=0.05*cosx+0.5;
   res[dimension+1]=phi;
   
  double dFdphi= thermodyn_.reactionSource(rho,phi); 
  double dFdrho=thermodyn_.chemicalPotential(rho, phi);
     
   if( dimRange > dimDomain+2)
    {
      //mu
      res[dimension+2]=0.5*v*v+dFdrho;
      //tau
      res[dimension+3]=thermodyn_.delta()*4*0.05*M_PI*M_PI*cosx+dFdphi;


#if SCHEME==DG
    //sigma
    if(dimension==1)
      {
        res[dimension+4]=-2*M_PI*sinx*cost*0.05;
      #if RHOMODEL
        res[dimension+5]=res[dimension+4]*rho;
      #endif
      }
     else
      { 
        res[dimension+4]=-2*M_PI*sinx*cost*0.5;
        res[dimension+5]=0; 
      }
    }
   else
    {
#if NONCONTRANS
#else
   //   res[1]*=res[0];
      res[2]*=res[0];
#endif
    }
#elif SCHEME==FEM
#endif
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

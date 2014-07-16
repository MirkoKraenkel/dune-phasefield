#ifndef STATIONARYPROBLEM_HH
#define STATIONARYPROBLEM_HH
#include <dune/common/version.hh>

// dune-fem includes
#include <dune/fem/misc/linesegmentsampler.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/space/common/functionspace.hh>

// local includes

//#include <dune/phasefield/modelling/thermodynamicsexactequi.hh>
#warning "USE thermodynamicsexactequi.hh for EOC "
#include <dune/phasefield/modelling/thermodynamicsbalancedphases.hh>


//#include <dune/fem/probleminterfaces.hh>

#include <dune/fem-dg/models/defaultprobleminterfaces.hh>


namespace Dune {

template <class GridType, class RangeProvider>
class StationaryProblem : public EvolutionProblemInterface<
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

  StationaryProblem() : 
    myName_( "Mixedtest Heatproblem" ),
    endTime_ ( Fem::Parameter::getValue<double>( "phasefield.endTime",1.0 )), 
    mu_( Fem::Parameter :: getValue< double >( "phasefield.mu1" )),
    delta_(Fem::Parameter::getValue<double>( "phasefield.delta" )),
    rho_( Fem::Parameter::getValue<double> ("phasefield.rho0")),
    phiscale_(Fem::Parameter::getValue<double> ("phiscale")),
    radius_(Fem::Parameter::getValue<double>("phasefield.radius")),
    thermodyn_()
    {
    }



  // initialize A and B 
  double init (const bool returnA ) const ;

  // print info 
  void printInitInfo () const;

  // source implementations 
  inline bool hasStiffSource () { return true; }
  inline bool hasNonStiffSource () { return false; }
  // this is the initial data
  inline void evaluate ( const DomainType& arg , RangeType& res ) const 
  {
    evaluate( 0.,arg ,res);
  }


  // evaluate function 
  inline void evaluate ( const double t, const DomainType& x, RangeType& res ) const;

  // cloned method 
  inline void evaluate ( const DomainType& x, const double t, RangeType& res ) const
  {
    evaluate( t, x, res );
  }

  inline double evalrho ( const double phi ) const;

  template< class DiscreteFunctionType >
  void finalizeSimulation ( DiscreteFunctionType& variablesToOutput,
                             const int eocloop) const
  {}


  const ThermodynamicsType& thermodynamics () const {return thermodyn_;}
  void printmyInfo ( std::string filename ) const {}
  inline double endtime () const { return endTime_; }
  inline std::string myName () const { return myName_; }
  void paraview_conv2prim () const {}
  std::string description () const;
 
  inline double mu () const { abort(); return mu_; }
  inline double delta () const{ return delta_; }
  protected:
  const std::string myName_;
  const double endTime_;
  const double mu_;
  const double delta_;
  double rho_;
  const double phiscale_;
  const double radius_;
  const ThermodynamicsType thermodyn_;
  
};


template <class GridType,class RangeProvider>
inline double StationaryProblem<GridType,RangeProvider>
:: init (const bool returnA ) const 
{
  return 0;
}



template <class GridType,class RangeProvider>
inline void StationaryProblem<GridType,RangeProvider>
:: printInitInfo () const
{}

template <class GridType,class RangeProvider>
inline void StationaryProblem<GridType,RangeProvider>
:: evaluate ( const double t, const DomainType& arg, RangeType& res ) const 
{
  double deltaInv=1./delta_;
  double phi=tanh( arg[0]*deltaInv ); 

  //(1-tanh^2)/delta
  double dxphi=0.5*deltaInv*(1-phi*phi);
  double dxdxphi=phi*dxphi;
  //-2/delta^2 *(tanh*(1-tanh^2)
  dxdxphi*=-2*deltaInv;

  phi*=0.5;
  phi+=0.5;
  double rho=evalrho(phi);
  res[0] = rho;
  for( int ii=0 ; ii<dimension ; ++ii)
    {
      res[1+ii] = 0;
    }
  res[2] = phi;
   
  double laplacePhi=delta_*dxdxphi;
  double dFdphi = thermodyn_.reactionSource(rho,phi); 
  double dFdrho = thermodyn_.chemicalPotential(rho, phi);
  
  
  if( dimRange > dimDomain+2)
    {
      //mu
      res[dimension+2]=dFdrho;
      //tau
      res[dimension+3]=0;
      //sigma
      for( int ii = 0 ; ii<dimension; ++ii )
        res[dimension+4+ii]=0;
        
      res[dimension+4]=dxphi;
      
    }
   else
    {
#if NONCONTRANS

#else
      res[2]*=res[0];
#endif
    }
}

template <class GridType,class RangeProvider>
inline double StationaryProblem<GridType,RangeProvider>
::evalrho ( const double phi ) const
{
  double t1;
  double t2;
  double t3;
  double t6;
  t1 = phi*phi;
  t2 = t1*t1;
  t3 = t2*phi;
  t6 = t1*phi;
  return(exp((4.0-18.0*t3+45.0*t2-30.0*t6)/(-0.9E1*t3+0.225E2*t2-0.15E2*t6+3.0)));
}



template <class GridType,class RangeProvider>
inline std::string StationaryProblem<GridType,RangeProvider>
:: description() const
{
  std::ostringstream stream;
  std::string returnString = stream.str();
  return returnString;
}

} // end namespace Dune
#endif

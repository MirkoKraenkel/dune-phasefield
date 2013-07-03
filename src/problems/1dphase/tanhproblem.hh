#ifndef PROBLEM_HH
#define PROBLEM_HH

#include <dune/common/version.hh>

// dune-fem includes
#include <dune/fem/misc/linesegmentsampler.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/space/common/functionspace.hh>

// local includes
// #include "idealthermodynamics2interpol.hh"
#include <dune/phasefield/modelling/thermodynamicsTest.hh>
#include <dune/fem/probleminterfaces.hh>




namespace Dune {

template <class GridType>
class TanhProblem : public EvolutionProblemInterface<
                    Dune::Fem::FunctionSpace< double, double, GridType::dimension, GridType::dimension + 2 >,
                    true >
                  
{
 
  public:
  typedef Fem::FunctionSpace<typename GridType::ctype,
                        double, GridType::dimensionworld,
                        GridType::dimensionworld + 2 > FunctionSpaceType ;

  enum{ dimension = GridType::dimensionworld };
  enum { energyId = dimension + 1 };
  typedef typename FunctionSpaceType :: DomainFieldType   DomainFieldType;
  typedef typename FunctionSpaceType :: DomainType        DomainType;
  typedef typename FunctionSpaceType :: RangeFieldType    RangeFieldType;
  typedef typename FunctionSpaceType :: RangeType         RangeType;

  typedef PGThermodynamics  ThermodynamicsType;
  
  TanhProblem() : 
		 myName_( "Constant Problem" ),
		 endTime_ ( Fem::Parameter::getValue<double>( "phasefield.endTime", 1.)), 
		 mu_( Fem::Parameter :: getValue< double >( "phasefield.mu1" )),
		 delta_(Fem::Parameter::getValue<double>( "phasefield.delta" )),
     rho1_(Fem::Parameter :: getValue< double >( "phasefield.rho1")),
		 rho2_(Fem::Parameter ::getValue< double> ("phasefield.rho2")),
		 smear_( Fem::Parameter::getValue<double> ("phasefield.smear")),
		 persistentsmear_(smear_),
     phiscale_(Fem::Parameter::getValue<double> ("phiscale")),
		 gamma_(Fem::Parameter::getValue<double> ("gamma")),
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
    evaluate(0.,arg,res);
  }


  // evaluate function 
  inline void evaluate( const double t, const DomainType& x, RangeType& res ) const;

  // cloned method 
  inline void evaluate( const DomainType& x, const double t, RangeType& res ) const
  {
    evaluate( t, x, res );
  }

  inline double rhoval(const double phi) const;

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
 
  inline double mu() const { return mu_; }
  inline double delta() const{return delta_;}
  inline double rho1() const {return rho1_;}
  inline double rho2() const {return rho2_;}
  inline double gamma() const {return gamma_;}
  inline void smearone()  {smear_=1.;}
  inline void resetsmear()  {smear_=persistentsmear_;}
  protected:
  const std::string myName_;
  const double endTime_;
  const double mu_;
  const double delta_;
  const double rho1_,rho2_;
  double smear_;
  double persistentsmear_;
  const double phiscale_;
  const double gamma_;
  const ThermodynamicsType thermodyn_;
  
};


template <class GridType>
inline double TanhProblem<GridType>
:: init(const bool returnA ) const 
{

  return 0;
}



template <class GridType>
inline void TanhProblem<GridType>
:: printInitInfo() const
{}

template <class GridType>
inline void TanhProblem<GridType>
:: evaluate( const double t, const DomainType& arg, RangeType& res ) const 
{
  double x=arg[0];
 
  for(int i=1;i<=dimension;i++)
      res[i]=0;

   double tanx;
  
  if(t==0.)
    tanx=0.5*tanh(x/(smear_*delta_))+0.5;
  else
    tanx=0.5*tanh(x/delta_)+0.5;
  
   res[0]=1.;
   
   res[dimension+1]=tanx;
 
}

template <class GridType>
inline double TanhProblem<GridType>
::rhoval( const double x ) const 
{
abort();

}



template <class GridType>
inline std::string TanhProblem<GridType>
:: description() const
{
  std::ostringstream stream;
  std::string returnString = stream.str();
  return returnString;
}

} // end namespace Dune
#endif

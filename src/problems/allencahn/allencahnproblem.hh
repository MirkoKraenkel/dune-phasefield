#ifndef PROBLEM_HH
#define PROBLEM_HH

#include <dune/common/version.hh>

// dune-fem includes
#include <dune/fem/misc/linesegmentsampler.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/space/common/functionspace.hh>




#include <dune/fem/probleminterfaces.hh>




namespace Dune {

template <class GridType>
class AllenCahnProblem : public EvolutionProblemInterface<
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

  typedef BalancedThermodynamics  ThermodynamicsType;
  
  AllenCahnProblem() : 
    myName_( "AllenCahnBalanced Problem" ),
    endTime_ ( Fem::Parameter::getValue<double>( "phasefield.endTime",1.0 )), 
    delta_(Fem::Parameter::getValue<double>( "phasefield.delta" ))
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
 
  inline double mu() const { abort(); return mu_; }
  inline double delta() const{return delta_;}
  inline double rho1() const {return rho1_;}
  inline double rho2() const {return rho2_;}
  inline double gamma() const {return gamma_;}
  inline void smearone()   {smear_=1.;}
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
inline double AllenCahnProblem<GridType>
:: init(const bool returnA ) const 
{

  return 0;
}



template <class GridType>
inline void AllenCahnProblem<GridType>
:: printInitInfo() const
{}

template <class GridType>
inline void AllenCahnProblem<GridType>
:: evaluate( const double t, const DomainType& arg, RangeType& res ) const 
{

}





template <class GridType>
inline std::string AllenCahnProblem<GridType>
:: description() const
{
  std::ostringstream stream;
  std::string returnString = stream.str();
  return returnString;
}

} // end namespace Dune
#endif

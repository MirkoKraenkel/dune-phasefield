#ifndef PROBLEM_HH
#define PROBLEM_HH
#warning "ACPROBLEM"
#include <dune/common/version.hh>

// dune-fem includes
#include <dune/fem/io/parameter.hh>
#include <dune/fem/space/common/functionspace.hh>

#include <dune/phasefield/modelling/thermodynamicsTest.hh>
#include <dune/fem/probleminterfaces.hh>

namespace Dune {
#if 0  
typedef GridSelector::GridType GridType;
typedef typename Dune::Fem::FunctionSpace< double, double, GridType::dimension, 1> FunctionSpaceType;
typedef typename Dune::Fem::FunctionSpace< double, double, GridType::dimension, 1> VelocityFunctionSpaceType;
#endif


template<class FuncSpace>
class Velocity
: public Fem::Function<FuncSpace,Velocity<FuncSpace>>
{
  public:
    typedef FuncSpace VelocityFunctionSpaceType;
    typedef typename VelocityFunctionSpaceType::RangeType RangeType;
    typedef typename VelocityFunctionSpaceType::RangeFieldType RangeFieldType;
    typedef typename VelocityFunctionSpaceType::DomainType DomainType;

    Velocity():gamma_(Fem::Parameter::getValue<double>("phasefield.gamma")){}
  

    void evaluate(const DomainType& x, RangeType &ret) const
    {
      ret[0]=gamma_*1.;
    }

    void evaluate(const DomainType& x, double time, RangeType& ret) const
    { 
      evaluate(x,ret);
    }
 
  private:
    double gamma_;
};

template <class GridType>
class AllenCahnProblem : public EvolutionProblemInterface<
                         Fem::FunctionSpace< double, double, GridType::dimension, 1>,
                         true >               
{
 
  public:
 typedef Fem::FunctionSpace<typename GridType::ctype,
                       double, GridType::dimensionworld,
                        1  > FunctionSpaceType ;

  enum{ dimension = FunctionSpaceType::dimDomain };
  typedef typename FunctionSpaceType :: DomainFieldType   DomainFieldType;
  typedef typename FunctionSpaceType :: DomainType        DomainType;
  typedef typename FunctionSpaceType :: RangeFieldType    RangeFieldType;
  typedef typename FunctionSpaceType :: RangeType         RangeType;

  
  AllenCahnProblem() : 
    myName_( "AllenCahn Problem" ),
    endTime_ ( Fem::Parameter::getValue<double>( "phasefield.endTime",1.0 )), 
    delta_(Fem::Parameter::getValue<double>( "phasefield.delta" )),
    smear_(Fem::Parameter::getValue<double>( "phasefield.smear" )),
    gamma_(Fem::Parameter::getValue<double>( "phasefield.gamma" ))
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
  inline void evaluate( const DomainType& arg, const double t, RangeType& res ) const
  {
   // res=arg[0];
   evaluate( t, arg, res );
  }


  template< class DiscreteFunctionType >
  void finalizeSimulation( DiscreteFunctionType& variablesToOutput,
                           const int eocloop) const
  {}


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
  double smear_;
  const double gamma_; 
};


template <class GridType>
inline double AllenCahnProblem<GridType>
:: init(const bool returnA ) const 
{
  abort();
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
  double x=arg[0];
  res=0.5*tanh((x-gamma_*t)/delta_)+0.5;
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

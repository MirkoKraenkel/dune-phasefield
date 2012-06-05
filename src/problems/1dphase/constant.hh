#ifndef DUNE_CONSTANT_HH
#define DUNE_CONSTANT_HH

// dune-fem include
#include <dune/fem/io/parameter.hh>
#include <dune/fem/space/common/functionspace.hh>

// local includes
#include <dune/fem/probleminterfaces.hh>
#include "../coldbubble/thermodynamics.hh"


/*****************************************************
 *                                                   *
 * 2d problem                                        *
 * NSConstant solution problem                       *
 *   with exact analytical solution for Euler and    *
 *   NS eqns with source term equal 0                *
 *                                                   *
 ****************************************************/

namespace Dune {

template <class GridType>
class NSConstant : public EvolutionProblemInterface<
                  Dune::FunctionSpace< double, double, GridType::dimension, GridType::dimension + 2 >,
                  true >
{
  NSConstant( const NSConstant& );
 public:
  typedef FunctionSpace<typename GridType::ctype, 
                        double, GridType::dimensionworld,
                        GridType::dimensionworld + 2 > FunctionSpaceType ;

  typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType :: DomainType DomainType;
  typedef typename FunctionSpaceType :: RangeFieldType  RangeFieldType;
  typedef typename FunctionSpaceType :: RangeType  RangeType;


  NSConstant() 
    : myName_( "NSConstant" )

    , Re_( Parameter::getValue< double >( "reynold" ) )
    , Re_inv_( 1. / Re_ )
    , Pr_( Parameter::getValue< double >( "prandtl" ) )
    , Pr_inv_( 1. / Pr_ )
  {
    init();
  }


  void init();


  void printInitInfo() const;

  // source implementations 
  inline bool hasStiffSource() { return false; }
  inline bool hasNonStiffSource() { return false; }
  inline double stiffSource( const double t, const DomainType& x, RangeType& res ) const;
  inline double nonStiffSource( const double t, const DomainType& x, RangeType& res ) const;

  // this is U0
  inline void evaluate( const DomainType& arg, RangeType& res ) const
  {
    evaluate( 0, arg, res );
  }


  inline void evaluate( const DomainType& arg, const double t, RangeType& res ) const
  {
    evaluate( t, arg, res );
  }


  inline void evaluate( const double t, const DomainType& arg, RangeType& res ) const
  {
    res[0] = 1.;
    res[1] = 0.;
    res[2] = 0.;
#if GRIDDIM==3    
    res[3] = 0.;
    res[4] = 3.;
#elif GRIDDIM==2
    res[3] = 3.;
#endif
  }


  template <class RangeImp>
  inline void pressAndTemp( const RangeImp& u, double& p, double& T ) const;


  Thermodynamics thermodynamics() const { return thermodynamics_; }
  inline double Re() const { return Re_; }
  inline double Re_inv() const { return Re_inv_; }
  inline double Pr() const { return Pr_; }
  inline double Pr_inv() const { return Pr_inv_; }
  inline double c_p() const { return c_p_; }
  inline double c_v_inv() const { return c_v_inv_; }
  inline double R_d() const { return R_d_; }
  inline double R_v() const { return R_v_; }
  inline double endTime() const { return endTime_; }
  inline double gamma() const { return gamma_; }
  inline double g() const { return g_; }
  inline std::string myName() const { return myName_; }
  inline double lambda() const { return -2./3.*mu_; }
  inline double lambda( const double T ) const { return -2./3.*mu_; }
  inline double mu() const { return mu_; }
  inline double mu( const double T ) const { return mu_; }
  inline double k() const { return c_p_*mu()*Pr_inv_; }
  inline double k( const double T ) const { return c_p_*mu(T)*Pr_inv_; }
  void printmyInfo( const std::string& filename )const {}
  std::string description() const { return ""; }
  void paraview_conv2prim() const {}

private:
  const Thermodynamics thermodynamics_;
  const std::string myName_;
  const double Re_;
  const double Re_inv_;
  const double Pr_;
  const double Pr_inv_;
  double gamma_;
  double g_;
  double endTime_;
  double c_v_;
  double c_v_inv_;
  double c_p_;
  double c_p_inv_;
  double R_d_;
  double R_d_inv_;
  double R_v_;
  double R_v_inv_;
  double lambda_;
  double mu_;
};




template <class GridType> 
inline void NSConstant<GridType>
:: init()
{
  g_ = Dune::Parameter::getValue< double >( "g" );
  endTime_ = Dune::Parameter::getValue< double >( "femhowto.endTime" );
  mu_ = Dune::Parameter :: getValue< double >( "mu" ); 
  c_v_ = Dune::Parameter :: getValue< double >( "c_v" );
  c_p_ = Dune::Parameter :: getValue< double >( "c_p" );

  // additional coefficients
  assert( g_ == 0. );
  c_v_inv_ = 1. / c_v_;
  c_p_inv_ = 1. / c_p_;
  R_d_ = c_p_ - c_v_;
  R_d_inv_ = 1. / R_d_;
  gamma_ = c_p_ * c_v_inv_;
}



template <class GridType> 
inline void NSConstant<GridType>
:: printInitInfo() const
{
  std::cout <<R_d_ <<" = R_d.\n";
  std::cout <<R_v_ <<" = R_v.\n";
  std::cout <<c_v_inv_ <<" = c_v_inv.\n";
  std::cout <<mu_ <<" = mu\n";
  std::cout <<c_p_ <<" = c_p_.\n";
  std::cout <<gamma_ <<" = gamma.\n";
  std::cout <<"DONE.\n";
}


template <class GridType>
inline double NSConstant<GridType>
:: stiffSource( const double t, const DomainType& x, RangeType& res ) const
{
  res = 0.;
  return 0.;
}


template <class GridType>
inline double NSConstant<GridType>
:: nonStiffSource( const double t, const DomainType& x, RangeType& res ) const
{
  res = 0.;
  return 0.;
}


template <class GridType> 
template <class RangeImp>
inline void NSConstant<GridType>
:: pressAndTemp( const RangeImp& u, double& p, double& T ) const
{
  const double rho_inv = 1./u[0];
#if GRIDDIM==2
  const double rhoeps = u[3] - 0.5*rho_inv*(u[1]*u[1] + u[2]*u[2]);
#elif GRIDDIM==3
  const double rhoeps = u[4] - 0.5*rho_inv*(u[1]*u[1] + u[2]*u[2]+u[3]*u[3]);
#else
#error "Grid dimension must be 2 or 3"
#endif
  assert( rhoeps > 1e-15 );
  p = (gamma_ - 1.)*rhoeps;
  T = p * rho_inv * R_d_inv_;
}

}
#endif

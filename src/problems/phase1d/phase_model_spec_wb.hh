#ifndef Phase_MODEL_SPEC_HH
#define Phase_MODEL_SPEC_HH

#include <dune/common/fvector.hh>

#include <dune/fem/misc/fmatrixconverter.hh>

#include "problemtype.hh"

namespace Dune {
  
template< int dimDomain >
class PhaseFlux
{
  enum { e = dimDomain + 1 };
  enum { dimRange = dimDomain + 2 };
  enum { dimGradRange = dimRange * dimDomain };
  enum { dimension = dimDomain};
  public:
  typedef Dune::FieldVector< double, dimDomain > DomainType;
  typedef Dune::FieldVector< double, dimDomain + 2 > RangeType;
  typedef Dune::FieldVector< double, dimGradRange > GradientRangeType;
  typedef Dune::FieldMatrix< double, dimRange, dimDomain > JacobianRangeType;
  typedef Dune::FieldMatrix< double, dimGradRange, dimDomain > JacobianFluxRangeType;
  typedef Dune::FieldMatrixConverter< GradientRangeType, JacobianRangeType > 
    ConvertedJacobianRangeType;

  PhaseFlux( const PhaseProblemType& problem )
    : problem_( problem )
    , delta_( problem.thermodynamics().delta()),
      gamma_( problem.gamma())
  {}

  inline void analyticalFlux( const RangeType& u, JacobianRangeType& f ) const;
  inline void jacobian( const RangeType& u, JacobianFluxRangeType& du ) const;
  inline double maxSpeed( const DomainType& n, const RangeType& u ) const;

  inline void diffusion( const RangeType& u,
                         const GradientRangeType& du,
												 JacobianRangeType& f ) const;

  template <class JacobianRangeImp>
  inline void diffusion( const RangeType& u,
												 const JacobianRangeImp& du,
												 JacobianRangeType& f ) const;
  
  template< class JacobianRangeImp >  
  inline void tension ( const RangeType& u,
												const JacobianRangeImp& du,
												GradientRangeType& f ) const;

 

  inline double mu(  ) const { return problem_.mu(); }
  inline double lambda( const double T ) const {abort(); }
  inline double k( const double T ) const { abort(); }
  inline double delta() const {return delta_;}

private:
  const PhaseProblemType& problem_;
  const  double delta_;
  const double gamma_;
};




/********************
 * 1D
 *******************/
 
template <>
  inline  void PhaseFlux<1>  
  :: analyticalFlux( const FieldVector<double,1+2>& u
		     , FieldMatrix<double,1+2,1>& f ) const 
  {
 
    assert( u[0] > 1e-10 );
    double rho_inv = 1. / u[0];
    const double v = u[1]*rho_inv;
    double p;
    double W;
    
    problem_.pressAndTemp( u, p, W );
    
    JacobianRangeType phiT(0.);
    f[0][0] = u[1];
    f[1][0] = v*u[1]+p;
    f[e][0] = u[e]*v;
      
  }


/*e
 * @brief divergence of this method is grad(u)
 */
template< >
inline void PhaseFlux<1>
:: jacobian( const RangeType& u,
             JacobianFluxRangeType& a ) const
{
  assert(u[0] > 1e-10);

  a[0][0] = u[0];
  a[1][0] = u[1];
  a[2][0] = u[2];
 
}


template< >
template< class JacobianRangeImp >
void PhaseFlux<1> 
:: diffusion( const RangeType& u,
              const JacobianRangeImp& du,
	    //   const JacobianRangeImp& tension,
	      JacobianRangeType& diff ) const 
{
  assert( u[0] > 1e-10 );
  double rho_inv = 1. / u[0];
  double p;
  double T;
  problem_.pressAndTemp( u, p, T );
  const double muLoc = mu(  );

  const double v   =  u[1]*rho_inv;
  const double phi =  u[2]*rho_inv;
  const double dxrho     = du[0][0]; //drho/dx
  const double dxrhou    = du[1][0]; //d(rho*v)/dx
  const double dxrhophi  = du[2][0]; //d(rho*phi)/dx
  
  const double dxv   = rho_inv*(dxrhou - v*dxrho);
  const double dxphi = rho_inv*(dxrhophi - phi*dxrho);
 
  diff[0][0]=0.;
  diff[1][0]=muLoc*dxv;
  diff[2][0]=dxphi*gamma_;
  

}

template< >
template< class JacobianRangeImp >
void PhaseFlux<1> 
:: tension(const RangeType& u,
	   const JacobianRangeImp& du,
	   GradientRangeType& diff ) const 
{
  assert( u[0] > 1e-10 );
  double rho_inv = 1. / u[0];
  
  const double phi =  u[2]*rho_inv;
  const double dxrho     = du[0][0]; //drho/dx
  const double dxrhophi  = du[2][0]; //d(rho*phi)/dx
  
  const double rhodxphi = (dxrhophi - phi*dxrho);
  const double dxphi = rho_inv*(dxrhophi - phi*dxrho);
  double tension = gamma_*(dxphi*dxphi-dxphi*dxphi*0.5);
  
  diff[0]=0.;
  diff[1]=tension;
  diff[2]=0;
  
}

template <>
inline
double PhaseFlux<1> :: maxSpeed( const FieldVector<double,1>& n
                              , const FieldVector<double,1+2>& u) const
{

  assert( u[0] > 1e-10 );
  double u_normal = u[1]*n[0]/u[0];
 
  
  return fabs( u_normal ) + sqrt( 1.4 );
  
}







} // end namespace Dune
#endif

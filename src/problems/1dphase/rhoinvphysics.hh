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
    , delta_( problem.thermodynamics().delta())
    , gamma_(problem.gamma())
	   //     , Re_inv_( problem.Re_inv() )
//     , c_v_inv_( problem.thermodynamics().c_vd_inv() )
//     , c_v_( problem.thermodynamics().c_vd() )
  {}
  inline void analyticalFlux( const RangeType& u, JacobianRangeType& f ) const;
  inline void jacobian( const RangeType& u, JacobianFluxRangeType& du ) const;
  inline double maxSpeed( const DomainType& n, const RangeType& u ) const;

  inline void diffusion( const RangeType& u,
                         const GradientRangeType& du,
		// 	 const GradientRangeType& tension
			 JacobianRangeType& f ) const;

  template <class JacobianRangeImp>
  inline void diffusion( const RangeType& u,
			 const JacobianRangeImp& du,
		// 	 const JacobianRangeImp& tension,
			 JacobianRangeType& f ) const;
  
  template< class JacobianRangeImp >  
  inline void tension ( const RangeType& u,
			const JacobianRangeImp& du,
			GradientRangeType& f ) const;

 

  inline double mu(  ) const { return problem_.mu(); }
  inline double lambda( const double T ) const {abort(); }
  inline double k( const double T ) const { abort(); }
  inline double delta() const {return delta_;}
  inline double gamma() const {return gamma_;}

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
  //   double phi=u[2]*rho_inv;
    
    
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
//   const double v = u[1]*rho_inv;
//   const double v2 = v*v;
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
   const double rhodxphi = (dxrhophi - phi*dxrho)*rho_inv*rho_inv;
   

  diff[0][0]=0.;
  diff[1][0]=muLoc*dxv;
  diff[2][0]=rhodxphi*gamma_;
  
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
  double tension =gamma_*(dxphi*dxphi*rho_inv-0.5*dxphi*dxphi*rho_inv);
  
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
 

 

 // n is unitOuterNormal
 // double c2 = gamma_ * p* rho_inv;// * n.two_norm2();
  //  assert( c2 > 1e-10 );

  return fabs( u_normal ) + sqrt( 1.4 );
  
}



/********************
 * 2D
 *******************/
template <>
template< class JacobianRangeImp >
void PhaseFlux<2> 
:: diffusion( const RangeType& u,
              const JacobianRangeImp& du,
	      JacobianRangeType& diff ) const 
{
  // du is grad(u) which is 4x2 matrix (for 2d case)
  assert( u[0] > 1e-10 );
  double rho_inv = 1. / u[0];
  const double v[2] = { u[1]*rho_inv, u[2]*rho_inv };
  const double phi = u[3]*rho_inv;
  
  //double p;
  //double T;
  // problem_.pressAndTemp( u, p, T );
  const double muLoc = mu();
  const double lambdaLoc =mu();

  // get dx_u, dz_u, dx_w, dz_w, dx_T, dz_T (in 2d case) for du
  const double du00=du[0][0];
  const double du01=du[0][1];
  const double du10=du[1][0];//dx rho*U_1
  const double du11=du[1][1];//dy rho*U_1
  const double du20=du[2][0];//dx rho*U_2
  const double du21=du[2][1];//dy rho*U_2
  const double du30=du[3][0];//dx rho*phi
  const double du31=du[3][1];;//dz rho*phi
  
  const double dxu = rho_inv*(du10 - v[0]*du00);//=1/rho(dx(rho*v1)-v1*dx(rho))=dx(v1);
  const double dzu = rho_inv*(du11 - v[0]*du01);
  const double dxw = rho_inv*(du20 - v[1]*du00);
  const double dzw = rho_inv*(du21 - v[1]*du01);
  
  const double rhodxphi = du30 - phi*du00; //dx(rho*phi)-phi*dx(rho)=rho*dx(phi);
  const double rhodzphi = du31 - phi*du01;
  
  
  const double tau00 = (2.*muLoc+lambdaLoc)*dxu + lambdaLoc*dzw;
  const double tau01 = muLoc*(dxw + dzu);
  const double tau10 = tau01;
  const double tau11 = lambdaLoc*dxu + (2.*muLoc+lambdaLoc)*dzw;

  // 1st row
  diff[0][0] = 0.;                   diff[0][1] = 0.;

  // 2nd row
  diff[1][0] = tau00;                diff[1][1] = tau01;

  // 3rd row
  diff[2][0] = tau10;                diff[2][1] = tau11;

  // 4th row
  diff[3][0] = gamma_*rhodxphi;           diff[3][1] =gamma_*rhodzphi;
  
}


template<>
inline void PhaseFlux<2> 
:: analyticalFlux( const RangeType& u,
                   JacobianRangeType& f ) const 
{
  assert( u[0] > 1e-10 );
  double rho_inv = 1. / u[0];
  double v[2] = { u[1]*rho_inv, u[2]*rho_inv};
  double p;
  double T;
 //  double phi=u[3]*rho_inv;

  
  
  problem_.pressAndTemp( u, p, T );
  f=0; 
  
  f[0][0] = u[1];                              f[0][1] = u[2];
  f[1][0] = v[0]*u[1]+p;                       f[1][1] = v[1]*u[1];
  f[2][0] = f[1][1];                           f[2][1] = v[1]*u[2] + p;
  f[3][0] = v[0]*u[3];                         f[3][1] = v[1]*u[3] ;
  
  
  
  
}



//converts Range in JacobianFluxRange allos grad u <-> div A(u);
template <>
inline void PhaseFlux<2>
:: jacobian( const PhaseFlux<2>::RangeType& u, 
             PhaseFlux<2>::JacobianFluxRangeType& a ) const 
{ 
  assert( u[0] > 1e-10 );
  a[0][0] = u[0];       a[0][1] = 0.;
  a[1][0] = 0.;         a[1][1] = u[0];
  a[2][0] = u[1];       a[2][1] = 0.;
  a[3][0] = 0.;         a[3][1] = u[1]; 
  a[4][0] = u[2];       a[4][1] = 0.;
  a[5][0] = 0.;         a[5][1] = u[2];
  a[6][0] = u[3];       a[6][1] = 0.;
  a[7][0] = 0.;         a[7][1] = u[3];
}

template< >
template< class JacobianRangeImp >
void PhaseFlux<2> 
:: tension(const RangeType& u,
	   const JacobianRangeImp& du,
	   GradientRangeType& diff ) const 
{
  
  
  assert( u[0] > 1e-10 );
  double rho_inv = 1. / u[0];
  
  const double phi =  u[2]*rho_inv;
  const double dxrho     = du[0][0]; //drho/dx
  const double dyrho     = du[0][1]; //drho/dy
  const double dxrhophi  = du[3][0]; //d(rho*phi)/dx
  const double dyrhophi  = du[3][1]; //d(rho*phi)/dy
  
  const double rhodxphi = (dxrhophi - phi*dxrho);
  const double rhodyphi = (dyrhophi - phi*dyrho);
  
  const double dxphi = rho_inv*(dxrhophi - phi*dxrho);
  const double dyphi = rho_inv*(dyrhophi - phi*dyrho);

  // double tension = delta_*rhodxphi*dxphi;
  
  diff[0]=0; 
  diff[1]=0;
  diff[2]=rhodxphi*dxphi;
  diff[3]=rhodyphi*dxphi;
  diff[4]=rhodyphi*dxphi;
  diff[5]=rhodyphi*dyphi;
  diff[6]=0;
  diff[7]=0;
  diff*=gamma_;
}


template <>
inline double PhaseFlux<2> 
:: maxSpeed( const FieldVector<double,2>& n, 
             const FieldVector<double,2+2>& u ) const 
{
  assert( u[0] > 1e-10 );
  double rho_inv = 1. / u[0];
  double u_normal = ( u[1]*n[0] + u[2]*n[1] ) * rho_inv;
  double p;
  double T;
  problem_.pressAndTemp( u, p, T );

  // n is unitOuterNormal
 

  
  return fabs( u_normal );
}




} // end namespace Dune
#endif

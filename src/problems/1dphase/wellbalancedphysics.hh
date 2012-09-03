#ifndef WELLBALANCEDPHYSICS_HH
#define WELLBALANCEDPHYSICS_HH

#include <dune/common/fvector.hh>

#include <dune/fem/misc/fmatrixconverter.hh>

#include "problemtype.hh"

namespace Dune {
  
template< int dimDomain >
class PhasePhysics
{
  enum { e = dimDomain + 1 };
  enum { dimRange = dimDomain + 2 };
 	enum { dimThetaRange = 2 };
	enum { dimGradRange = dimRange * dimDomain };
	enum { dimThetaGradRange = dimThetaRange * dimDomain };
	enum { dimension = dimDomain};
  public:
  typedef Dune::FieldVector< double, dimDomain > DomainType;
 
	typedef Dune::FieldVector< double, dimDomain + 2 > RangeType;
	typedef Dune::FieldVector< double,  2 > ThetaRangeType;
  
	typedef Dune::FieldVector< double, dimGradRange > GradientRangeType;
  typedef Dune::FieldVector< double, dimThetaGradRange > ThetaGradientRangeType;
  
	typedef Dune::FieldMatrix< double, dimRange, dimDomain > JacobianRangeType;
	typedef Dune::FieldMatrix< double, dimThetaRange, dimDomain > ThetaJacobianRangeType;

  typedef Dune::FieldMatrix< double, dimGradRange, dimDomain > JacobianFluxRangeType;
  typedef Fem::FieldMatrixConverter< GradientRangeType, JacobianRangeType > 
    ConvertedJacobianRangeType;

  PhasePhysics( const PhaseProblemType& problem )
    : problem_( problem )
    , delta_( problem.thermodynamics().delta()),
      gamma_( problem.gamma())
  {}

  inline void analyticalFlux( const RangeType& u, JacobianRangeType& f ) const;
  inline void jacobian( const RangeType& u, JacobianFluxRangeType& du ) const;
  inline double maxSpeed( const DomainType& n, const RangeType& u ) const;

	inline void nonConProduct(const RangeType & uL, 
														const RangeType & uR,
														const ThetaRangeType& thetaL,
														const ThetaRangeType& thetaR,
														RangeType& ret) const;


	inline double stiffSource(const RangeType& u,
														const GradientRangeType& du,
														const ThetaRangeType& theta,
														const ThetaJacobianRangeType& dtheta,
														RangeType& f) const;
//   template< class JacobianRangeImp >  
// 	inline double stiffSource( const RangeType& u,
// 														 const JacobianRangeType& du,
// 														 RangeType& f ) const;



//   inline void diffusion( const RangeType& u,
//                          const GradientRangeType& du,
// 												 JacobianRangeType& f ) const;

  template <class JacobianRangeImp>
  inline void diffusion( const RangeType& u,
												 const JacobianRangeImp& du,
												 JacobianRangeType& f ) const;
	
	
  inline void allenCahn( const RangeType& u,
												const GradientRangeType& du,
												ThetaJacobianRangeType& f ) const;

  
  template< class JacobianRangeImp >  
  inline void allenCahn(const RangeType& u,
												const JacobianRangeImp& du,
												ThetaJacobianRangeType& f ) const;

 




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
	template<>
 	inline void   PhasePhysics<1>
	::nonConProduct(const RangeType & uL, 
									const RangeType & uR,
									const ThetaRangeType& thetaL,
									const ThetaRangeType& thetaR,
									RangeType& ret) const
	{
		ret[1]=0.5*(uL[0]+uR[0])*(thetaL[0]-thetaR[0]);
 	}
	
	template<>
	inline double PhasePhysics<1>
	::stiffSource(const RangeType& u,
								const GradientRangeType& du,
								const ThetaRangeType& theta,
								const ThetaJacobianRangeType& dtheta,
								RangeType& f) const
	{
			f[0]=0;
			f[1]=dtheta[0]*u[0]+du[2]*theta[1];
			f[2]=theta[1];

	
	}

	
	template <>
	inline  void PhasePhysics<1>  
	:: analyticalFlux( const FieldVector<double,1+2>& u
										 , FieldMatrix<double,1+2,1>& f ) const 
	{
		assert( u[0] > 1e-10 );
		double rho_inv = 1. / u[0];
		const double v = u[1]*rho_inv;
		double p;
		double W;
   
		JacobianRangeType phiT(0.);
		f[0][0] = u[1];
		f[1][0] = v*u[1];
		f[e][0] = u[e]*v;
	}


	/*e
	 * @brief divergence of this method is grad(u)
	 */
	template< >
	inline void PhasePhysics<1>
	:: jacobian( const RangeType& u,
							 JacobianFluxRangeType& a ) const
	{
		assert(u[0] > 1e-10);
		
		a[0][0] = u[0]; //rho
		a[1][0] = u[1]/u[0];//(rho v)/rho
		a[2][0] = u[2]/u[0];//(rho phi)/rho
	}
	
	//Navier Stokes Tensor
	template< >
	template< class JacobianRangeImp >
	void PhasePhysics<1> 
	:: diffusion( const RangeType& u,
								const JacobianRangeImp& du,
								JacobianRangeType& diff ) const 
	{
		assert( u[0] > 1e-10 );
		double rho_inv = 1. / u[0];
		double p;
		double T;
		problem_.pressAndTemp( u, p, T );
		const double muLoc = mu(  );
		
		const double dxv   = du[1][0]; //dv/dx
	
		diff[0][0]=0.;
		diff[2][0]=0.;
		diff[1][0]=muLoc*dxv;
	}


template< >
template< class JacobianRangeImp >
void PhasePhysics<1> 
:: allenCahn( const RangeType& u,
              const JacobianRangeImp& du,
							ThetaJacobianRangeType& diff ) const 
{
   assert( u[0] > 1e-10 );
	 
	 diff[0][0]=0.;
	 diff[1][0]=delta_*du[2][0];

}

template <>
inline
double PhasePhysics<1> :: maxSpeed( const FieldVector<double,1>& n
                              , const FieldVector<double,1+2>& u) const
{

  assert( u[0] > 1e-10 );
  double u_normal = u[1]*n[0]/u[0];
 
  
  return fabs( u_normal ) + sqrt( 1.4 );
  
}


/********************
 * 2D
 *******************/
	template<>
 	inline void   PhasePhysics<2>
	::nonConProduct(const RangeType & uL, 
									const RangeType & uR,
									const ThetaRangeType& thetaL,
									const ThetaRangeType& thetaR,
									RangeType& ret) const
	{
		ret[1]=0.5*(uL[0]+uR[0])*(thetaL[0]-thetaR[0]);
 	}
	
	template<>
	inline double PhasePhysics<2>
	::stiffSource(const RangeType& u,
								const GradientRangeType& du,
								const ThetaRangeType& theta,
								const ThetaJacobianRangeType& dtheta,
								RangeType& f) const
	{
			f[0]=0;
			const double dxmu = dtheta[0][0];
			const double dymu = dtheta[0][1];
			const double dxphi= du[6];
			const double dyphi= du[7];


			//dxmu*rho+dxphi*theta2
			f[1]=dxmu*u[0]+dxphi*theta[1];
			//dymu*rho+dyphi*theta2
			f[2]=dymu*u[0]+dyphi*theta[1];
			
			f[3]=theta[1];

	
	}

	
	template <>
	inline  void PhasePhysics<2>  
	:: analyticalFlux(const RangeType& u,
										JacobianRangeType& f ) const 
	{
		assert( u[0] > 1e-10 );
		double rho_inv = 1. / u[0];
		double v[2] = { u[1]*rho_inv, u[2]*rho_inv};
		
		f[0][0] = u[1];                              f[0][1] = u[2];
		f[1][0] = v[0]*u[1];                         f[1][1] = v[1]*u[1];
		f[2][0] = f[1][1];                           f[2][1] = v[1]*u[2];
		f[3][0] = v[0]*u[3];                         f[3][1] = v[1]*u[3] ;
  
		
	}


	/*e
	 * @brief divergence of this method is grad(u)
	 */
	template< >
	inline void PhasePhysics<2>
	:: jacobian( const RangeType& u,
							 JacobianFluxRangeType& a ) const
	{
		assert(u[0] > 1e-10);
		const double rhoInv=1./u[0];
		a[0][0] = u[0];              a[0][1] = 0.;
		a[1][0] = 0.;                a[1][1] = u[0];
		a[2][0] = u[1]*rhoInv;       a[2][1] = 0.;
		a[3][0] = 0.;                a[3][1] = u[1]*rhoInv; 
		a[4][0] = u[2]*rhoInv;       a[4][1] = 0.;
		a[5][0] = 0.;                a[5][1] = u[2]*rhoInv;
		a[6][0] = u[3]*rhoInv;       a[6][1] = 0.;
		a[7][0] = 0.;                a[7][1] = u[3]*rhoInv;
		
	}
	
	//Navier Stokes Tensor
	template< >
	template< class JacobianRangeImp >
	void PhasePhysics<2> 
	:: diffusion( const RangeType& u,
								const JacobianRangeImp& du,
								JacobianRangeType& diff ) const 
	{
		// du is grad(rho,v,phi) which is 4x2 matrix (for 2d case)
		assert( u[0] > 1e-10 );
		double rho_inv = 1. / u[0];
		double p;
		double T;
		problem_.pressAndTemp( u, p, T );
		const double muLoc = 1.;
		const double lambdaLoc = 1.;
		
			const double dxu1=du[1][0];//dx U_1
		const double dyu1=du[1][1];//dy U_1
		const double dxu2=du[2][0];//dx U_2
		const double dyu2=du[2][1];//dy U_2
	 
		const double tau00 = (2.*muLoc+lambdaLoc)*dxu1 + lambdaLoc*dyu2;
		const double tau01 = muLoc*(dxu2 + dyu1);
		const double tau10 = tau01;
		const double tau11 = lambdaLoc*dxu1 + (2.*muLoc+lambdaLoc)*dyu2;
	  // 1st row
		diff[0][0] = 0.;                   diff[0][1] = 0.;
		// 2nd row
		diff[1][0] = tau00;                diff[1][1] = tau01;
		// 3rd row
		diff[2][0] = tau10;                diff[2][1] = tau11;
		// 4th row
		diff[3][0] = 0.;                   diff[3][1] = 0.;
	}


template< >
template< class JacobianRangeImp >
void PhasePhysics<2> 
:: allenCahn( const RangeType& u,
              const JacobianRangeImp& du,
							ThetaJacobianRangeType& diff ) const 
{
   assert( u[0] > 1e-10 );
	 
	 diff[0][0]=0.; 
	 diff[0][1]=0.;
	 diff[1][0]=du[2][0];
	 diff[1][1]=du[2][1];
	 diff*=delta_;
}

template <>
inline
double PhasePhysics<2> :: maxSpeed( const DomainType& n
                              , const RangeType& u) const
{

  assert( u[0] > 1e-10 );
  double u_normal = u[1]*n[0]/u[0];
 
  
  return fabs( u_normal ) + sqrt( 1.4 );
  
}









} // end namespace Dune
#endif

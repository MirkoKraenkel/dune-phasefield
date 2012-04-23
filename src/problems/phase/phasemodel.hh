#ifndef NS_MODEL_HH
#define NS_MODEL_HH

// DUNE includes
#include <dune/common/version.hh>
#include <dune/fem/misc/fmatrixconverter.hh>

// local includes
#include "phase_model_spec.hh"


namespace Dune {

// forward declaration of NSFlux
template< int >
struct PhaseFlux;

//////////////////////////////////////////////////////
//
//                 NAVIER-STOKES EQNS 
//
//////////////////////////////////////////////////////
template< class GridPart , int dimRangeImp ,
          int dimGradRangeImp = dimRangeImp * GridPart::GridType::dimensionworld >
class PhaseModelTraits 
{
  public:
  typedef GridPart GridPartType;
  typedef typename GridPart :: GridType GridType;
  enum { dimDomain = GridType :: dimensionworld };
  enum { dimRange = dimRangeImp };
  enum { dimGradRange = dimGradRangeImp };
  enum { dimGradient = dimDomain + 1 };
  
  typedef FieldVector< double, dimDomain > DomainType;
  typedef FieldVector< double, dimDomain - 1 > FaceDomainType;
  typedef FieldVector< double, dimRange > RangeType;
  typedef FieldVector< double, dimGradRange > GradientType;
  typedef FieldMatrix< double, dimRange, dimDomain > FluxRangeType;
  typedef FieldVector< double, dimGradRange > GradientRangeType;
  typedef FieldMatrix< double, dimRange, dimDomain > JacobianRangeType;
  typedef FieldMatrix< double, dimGradRange, dimDomain > JacobianFluxRangeType;

  typedef typename GridPart :: IntersectionIteratorType IntersectionIterator;
  typedef typename IntersectionIterator::Intersection Intersection;
  typedef typename GridType :: template Codim<0> :: Entity EntityType;
  typedef typename EntityType :: EntityPointer EntityPointerType;
};



template< class GridPartType , class ProblemImp >
class PhaseModel 
{
  public:
  typedef ProblemImp ProblemType;
  enum { dimDomain = GridPartType :: GridType :: dimensionworld };
  enum { dimRange = dimDomain + 2 };
  enum { dimGradRange = dimRange * dimDomain };
  typedef PhaseModelTraits< GridPartType, dimRange > Traits;

  typedef typename Traits :: GridType GridType;
  typedef typename Traits :: EntityType EntityType;
  typedef typename Traits :: EntityPointerType EntityPointerType;
  typedef typename Traits :: IntersectionIterator IntersectionIterator;
  typedef typename Traits :: Intersection IntersectionType;
  typedef typename Traits :: FaceDomainType FaceDomainType;

  typedef typename Traits :: DomainType DomainType;
  typedef typename Traits :: RangeType RangeType;
  typedef typename Traits :: GradientRangeType GradientRangeType;
  typedef typename Traits :: JacobianRangeType JacobianRangeType;
  typedef typename Traits :: JacobianFluxRangeType JacobianFluxRangeType;

 public:
  PhaseModel( const ProblemType& problem ) 
    : mu_( problem.mu() )
    , problem_( problem )
    , nsFlux_( problem )
    , c1_(problem_.c1())
    , c2_(problem_.c2())
      {}

  inline bool hasStiffSource() const { return false; }
  inline bool hasNonStiffSource() const { return true; }
  inline bool hasFlux()const{return true;}

  inline double stiffSource( const EntityType& en,
                        const double time,
                        const DomainType& x,
                        const RangeType& u,
                        const GradientRangeType& du,
                        RangeType & s) const
  {

    return stiffSource( en, time, x, u, s ); 

  }



  inline double stiffSource( const EntityType& en
                      , const double time
                      , const DomainType& x
                      , const RangeType& u
                      , RangeType& s ) const 
  {abort();
    
      s=0.;
      //double tmp=0;
      double delta=problem_.delta();
      double delta_inv=1./delta;
      double phi=u[dimDomain+1];
      double rho=u[0];
      //double c1=problem_.c1();
      double c1=0.1;
      double c2=0.9;
      //double c2=problem_.c2();

      phi/=rho;

      s[dimDomain+1]= 0.4e1 * pow(phi, 0.3e1)
        - 0.6e1 * phi * phi
        + 0.2e1 * phi
        + 0.4e1 * phi * rho * rho
        - 0.4e1 * phi * rho * c2
        + 0.2e1 * phi * c2 * c2
        - 0.2e1 * rho * rho
        - 0.4e1 * rho * c1 * phi
        + 0.4e1 * rho * c1
        + 0.2e1 * c1 * c1 * phi
        - 0.2e1 * c1 * c1;
    s[dimDomain+1]*=-delta_inv*delta_inv;
#if MYDEBUG
    std::cout<<"SOURCE=("<<s<<")\n******\n";
#endif
    return 1e+10;

 }


  inline double nonStiffSource( const EntityType& en,
                        const double time,
                        const DomainType& x,
                        const RangeType& u,
                        const GradientRangeType& du,
                        RangeType & s) const
  {
    //FieldMatrixConverter< GradientType, JacobianRangeType> jac( du );
    return nonStiffSource( en, time, x, u, s );
  }



  inline double nonStiffSource( const EntityType& en,
                        const double time,
                        const DomainType& x,
                        const RangeType& u,
                        const JacobianRangeType& jac,
                        RangeType & s) const
  {
    return nonStiffSource( en, time, x, u, s ); 
  }

xe
  inline double nonStiffSource( const EntityType& en
                      , const double time
                      , const DomainType& x
                      , const RangeType& u
                      , RangeType& s ) const 
  {
    const DomainType& xgl = en.geometry().global(x);
    return problem_.nonStiffSource( time, xgl, s );
  }
  

  inline void advection( const EntityType& en,
                         const double time,
                         const DomainType& x,
                         const RangeType& u,
                         JacobianRangeType& f ) const 
  {
    abort();
    nsFlux_.analyticalFlux( u, f );
  }

  inline double diffusionTimeStep( const IntersectionType& it,
                                   const double enVolume,
                                   const double circumEstimate,
                                   const double time,
                                   const FaceDomainType& x,
                                   const RangeType& u ) const
  {
    // look at Ch. Merkle Diplom thesis, pg. 38

    // get value of mu at temperature T
    double p;
    double T;
    problem_.pressAndTemp( u, p, T );
    // get mu 
    const double mu = problem_.mu( T );

    // ksi = 0.25 
    return mu * circumEstimate * alpha_ / (0.25 * u[0] * enVolume);
  }


  inline void jacobian( const EntityType& en
                        , const double time
                        , const DomainType& x
                        , const RangeType& u 
                        , JacobianFluxRangeType& a ) const 
  {
    nsFlux_.jacobian( u, a );
  }
  

  inline bool hasBoundaryValue( const IntersectionType& it
                                , const double time
                                , const FaceDomainType& x ) const 
  { 
    return true;
  }
 

  inline double boundaryFlux( const IntersectionType& it
                              , const double time
                              , const FaceDomainType& x
                              , const RangeType& uLeft
                              , const GradientRangeType& duLeft
                              , RangeType& gLeft ) const; 


  inline void boundaryValue( const IntersectionType& it,
                             const double time,
                             const FaceDomainType& x,
                             const RangeType& uLeft,
                             RangeType& uRight ) const; 
 
 

inline double diffusionBoundaryFlux( const IntersectionType& it,
                                     const double time,
                                     const FaceDomainType& x,
                                     const RangeType& uLeft,
                                     const GradientRangeType& gradLeft,
                                     RangeType& gLeft ) const;


template<class JacobianRangeImp>
double diffusionBoundaryFlux( const IntersectionType& it,
                                     const double time,
                                     const FaceDomainType& x,
                                     const RangeType& uLeft,
                                     const JacobianRangeImp& jacLeft,
                                     RangeType& gLeft ) const
    {
         std::cerr <<"diffusionBoundaryFlux shouldn't be used for this testcase" <<std::endl;
             abort();
     }
inline void eigenValues(const EntityType& en,
                        const double time,
                        const DomainType& x,
                        const RangeType& u,
                          RangeType& maxValue) const
    {
      abort();
          // TODO: calculate eigen values here
               }
          
 

  // here x is in global coordinates
  inline void maxSpeed( const DomainType& normal,
                        const double time,
                        const DomainType& x,
                        const RangeType& u,
                        double& advspeed,
                        double& totalspeed ) const 
  {
    advspeed = nsFlux_.maxSpeed( normal , u );
    totalspeed=advspeed;
  }


  void diffusion( const EntityType& en,
                  const double time,
                  const DomainType& x,
                  const RangeType& u,
                  const GradientRangeType& v,
                  JacobianRangeType& diff ) const
  {
    FieldMatrixConverter< GradientRangeType, JacobianRangeType> jac( v );
    diffusion( en, time, x, u, jac, diff );
  }


  template <class JacobianRangeImp>
  void diffusion( const EntityType& en,
                  const double time,
                  const DomainType& x,
                  const RangeType& u,
                  const JacobianRangeImp& jac,
                  JacobianRangeType& diff ) const
  {
    nsFlux_.diffusion( u, jac, diff );
  }

  inline double boundaryFlux( const IntersectionType& it
                   , const double time
                   , const FaceDomainType& x
                   , const RangeType& uLeft
                   //, const JacobianRangeType& duLeft
                   , RangeType& gLeft ) const  
  {
    gLeft = 0.;
    return 0.;
  }



  inline void conservativeToPrimitive( const DomainType& xgl,
                                       const RangeType& cons, 
                                       RangeType& prim
                                     ) const
  {
    thermodynamics_.conservativeToPrimitiveEnergyForm( cons, prim );
  }







  inline const ProblemType& problem() const { return problem_; }
  inline const double mu( const double T ) const { return problem_.mu(T); }
  inline const double lambda( const double T ) const { return problem_.lambda(T); }
  inline const double k( const double T ) const { return problem_.k(T); }
  inline const double c_v_inv() const { return c_v_inv_; }
  inline const double c_p() const { return problem_.c_p(); }

  

 protected:
  const double mu_; 
  const ProblemType& problem_;
  const PhaseFlux< dimDomain > nsFlux_;
  const double c1_;
  const double c2_;
};




template< class GridPartType, class ProblemImp >
inline double PhaseModel< GridPartType, ProblemImp >
:: boundaryFlux( const IntersectionType& it
                 , const double time
                 , const FaceDomainType& x
                 , const RangeType& uLeft
                 , const GradientRangeType& duLeft
                 , RangeType& gLeft ) const  
{
  DomainType xgl=it.intersectionGlobal().global(x);
  const typename Traits :: DomainType normal = it.integrationOuterNormal(x); 
  double p;
  double T;
  problem_.pressAndTemp( uLeft, p, T );
  gLeft = 0;

  // bnd. cond. from euler part
  for (int i=0;i<dimDomain; ++i)
    gLeft[i+1] = normal[i]*p;

#if 0
  double div_uw = duLeft[2] + duLeft[5];
  double tau11 = 2.*mu()*duLeft[2] + lambda_*div_uw;
  double tau12 = mu_*div_uw;
  double tau22 = 2.*mu_*duLeft[5] + lambda_*div_uw;
  // bnd. cond. from visc. part
  gLeft[1] += tau11*normal[0] + tau12*normal[1];
  gLeft[2] += tau12*normal[0] + tau22*normal[1];
  gLeft[3] += -2.*mu_*(duLeft[2]*normal[0]+duLeft[5]*normal[1]);
#endif

  return 0.;
}
template< class GridPartType, class ProblemImp >
inline double PhaseModel< GridPartType, ProblemImp >
::diffusionBoundaryFlux( const IntersectionType& it,
                                     const double time,
                                     const FaceDomainType& x,
                                     const RangeType& uLeft,
                                     const GradientRangeType& gradLeft,
                                     RangeType& gLeft ) const
   {
      FieldMatrixConverter< GradientRangeType, JacobianRangeType> jacLeft( gradLeft );
      return diffusionBoundaryFlux( it, time, x, uLeft, jacLeft, gLeft );
   }

  /** \brief boundary flux for the diffusion part
   *    */
/*
template< class GridPartType, class ProblemImp,class JacobianRangeImp>
inline double PhaseModel< GridPartType, ProblemImp >
::diffusionBoundaryFlux( const IntersectionType& it,
                                      const double time,
                                        const FaceDomainType& x,
                                        const RangeType& uLeft,
                                        const JacobianRangeImp& jacLeft,
                                        RangeType& gLeft ) const
      {
         std::cerr <<"diffusionBoundaryFlux shouldn't be used for this testcase" <<std::endl;
             abort();
     }
*/

template< class GridPartType, class ProblemImp >
inline void PhaseModel< GridPartType, ProblemImp >
:: boundaryValue( const IntersectionType& it
                  , const double time
                  , const FaceDomainType& x
                  , const RangeType& uLeft
                  , RangeType& uRight ) const 
{
  abort();
  DomainType xgl = it.geometry().global( x );
  problem_.evaluate( time , xgl , uRight );
}

} // end namespace Dune

#endif // file definition

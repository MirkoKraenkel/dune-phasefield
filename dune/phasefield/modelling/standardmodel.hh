#ifndef NS_MODEL_HH
#define NS_MODEL_HH

// DUNE includes
#include <dune/common/version.hh>
#include <dune/fem/misc/fmatrixconverter.hh>

// local includes
#include <dune/fem-dg/models/defaultmodel.hh>
#include "physics.hh"

namespace Dune {

// forward declaration of PhaseFlux
template< int >
struct PhaseFlux;

//////////////////////////////////////////////////////
//
//                 NAVIER-STOKES EQNS 
//
//////////////////////////////////////////////////////
template< class GridPart > 
class PhaseModelTraits 
{
  public:
  typedef GridPart GridPartType;
  typedef typename GridPart :: GridType                     GridType;
  enum { dimDomain = GridType :: dimensionworld };
  enum { dimRange = dimDomain + 2 };
  enum { dimScalGradRange = dimDomain };
	enum { dimGradRange = dimRange * dimDomain };
	enum { dimGradient = dimDomain + 1 };
  
  typedef FieldVector< double, dimDomain >                  DomainType;
  typedef FieldVector< double, dimDomain - 1 >              FaceDomainType;
  typedef FieldVector< double, dimRange >                   RangeType;
  typedef FieldVector< double, dimGradRange >               GradientType;
  typedef FieldMatrix< double, dimRange, dimDomain >        FluxRangeType;
  typedef FieldVector< double, dimGradRange >               GradientRangeType;
  typedef FieldVector< double, dimScalGradRange >           ScalarGradientRangeType;
	typedef FieldMatrix< double, dimRange, dimDomain >        JacobianRangeType;
  typedef FieldMatrix< double, dimGradRange, dimDomain >    JacobianFluxRangeType;

  typedef typename GridPart :: IntersectionIteratorType     IntersectionIterator;
  typedef typename IntersectionIterator::Intersection       Intersection;
  typedef Intersection       IntersectionType;
  typedef typename GridType :: template Codim<0> :: Entity  EntityType;
  typedef typename EntityType :: EntityPointer              EntityPointerType;

  typedef Thermodynamics< dimDomain >                       ThermodynamicsType;
};



template< class GridPartType , class ProblemImp >
class PhaseModel : public DefaultModel < PhaseModelTraits< GridPartType > >
{
  public:
  typedef ProblemImp                                        ProblemType;
  typedef typename GridPartType::GridType                   GridType;
  typedef PhaseModelTraits< GridPartType >                     Traits;

  enum { dimDomain = Traits :: dimDomain };
  enum { dimRange  = Traits :: dimRange  };
  enum { dimGradRange = Traits::dimGradRange } ;

  typedef typename Traits :: EntityType                     EntityType;
  typedef typename Traits :: EntityPointerType              EntityPointerType;
  typedef typename Traits :: IntersectionIterator           IntersectionIterator;
  typedef typename Traits :: Intersection                   IntersectionType;
  typedef typename Traits :: FaceDomainType                 FaceDomainType;

  typedef typename Traits :: DomainType                     DomainType;
  typedef typename Traits :: RangeType                      RangeType;
  typedef typename Traits :: GradientRangeType              GradientRangeType;
  typedef typename Traits :: JacobianRangeType              JacobianRangeType;
  typedef typename Traits :: JacobianFluxRangeType          JacobianFluxRangeType;
  typedef typename Traits :: ThermodynamicsType             ThermodynamicsType;

 public:
  PhaseModel( const ProblemType& problem ) 
    : thermodynamics_( problem.thermodynamics() )
    , problem_( problem )
    , nsFlux_( problem )
    , visc_(Dune::Fem::Parameter::getValue<double>("visc"))
  {
    
  }

  inline bool hasStiffSource() const { return true; }
  inline bool hasNonStiffSource() const { return false; }
  inline bool hasFlux() const { return true ; }

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
  {
   
		double reaction,pressure;
		thermodynamics_.pressureAndReaction(u,pressure,reaction);
	  s[0]=0.;
		s[1]=0.;
		s[2]=reaction;
		return 1;
  }

 


  inline double nonStiffSource( const EntityType& en,
                        const double time,
                        const DomainType& x,
                        const RangeType& u,
                        const GradientRangeType& du,
                        RangeType & s) const
  {
    Dune::Fem::FieldMatrixConverter< GradientRangeType, JacobianRangeType > jac( du );
    return nonStiffSource( en, time, x, u, jac, s );
  }



  template< class JacobianRangeTypeImp >
  inline double nonStiffSource( const EntityType& en,
                        const double time,
                        const DomainType& x,
                        const RangeType& u,
                        const JacobianRangeTypeImp& jac,
                        RangeType& s) const
  {

    abort();
    return 0.;
  }

  

  inline void advection( const EntityType& en,
                         const double time,
                         const DomainType& x,
                         const RangeType& u,
			 JacobianRangeType& f ) const 
  {
  
    nsFlux_.analyticalFlux(u, f);
    
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
    pressAndTemp( u, p, T );
    // get mu 
    const double mu = problem_.mu(  );

    // ksi = 0.25 
    
    return mu * circumEstimate * 1 / (0.25 * u[0] * enVolume);
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
  

  inline double diffusionBoundaryFlux( const IntersectionType& it,
                                       const double time,
                                       const FaceDomainType& x,
                                       const RangeType& uLeft,
                                       const GradientRangeType& gradLeft,
                                       RangeType& gLeft ) const  
  {
    Dune::Fem::FieldMatrixConverter< GradientRangeType, JacobianRangeType> jacLeft( gradLeft );
    return diffusionBoundaryFlux( it, time, x, uLeft, jacLeft, gLeft );
  }

  /** \brief boundary flux for the diffusion part
   */
  template <class JacobianRangeImp>
  inline double diffusionBoundaryFlux( const IntersectionType& it,
                                       const double time,
                                       const FaceDomainType& x,
                                       const RangeType& uLeft,
                                       const JacobianRangeImp& jacLeft,
                                       RangeType& gLeft ) const  
  {
    return 1.;
  }


  inline void boundaryValue( const IntersectionType& it,
                             const double time,
                             const FaceDomainType& x,
                             const RangeType& uLeft,
                             RangeType& uRight ) const; 
 
  
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
    Dune::Fem:: FieldMatrixConverter< GradientRangeType, JacobianRangeType> jac( v );
  //   FieldMatrixConverter< GradientRangeType, JacobianRangeType> tension( tens );
    
     diffusion( en, time, x, u, jac,diff );
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
   void tension( const EntityType& en,
		 const double time,
		 const DomainType& x,
		 const RangeType& u,
		 const GradientRangeType& v,
		 JacobianRangeType& diff )const
  {    
    Dune::Fem::FieldMatrixConverter< GradientRangeType, JacobianRangeType> jac( v );
    diff=jac;
  }
  
   void tension( const EntityType& en,
		 const double time,
		 const DomainType& x,
		 const RangeType& u,
		 const GradientRangeType& v,
		 GradientRangeType& diff )const
  {
    Dune::Fem::FieldMatrixConverter< GradientRangeType, JacobianRangeType> jac( v );
     tension( u, jac, diff );
  }
  
  template <class JacobianRangeImp>
   void tension( const RangeType& u,
		 const JacobianRangeImp& jac,
		 GradientRangeType& diff ) const
  { 
    nsFlux_.tension( u, jac, diff );
  }


  inline double boundaryFlux( const IntersectionType& it
                   , const double time
                   , const FaceDomainType& x
                   , const RangeType& uLeft
		   , RangeType& gLeft ) const  
  {
    gLeft = 0.;
    return 0.;
  }

  inline void pressAndTemp( const RangeType& u, double& p, double& T ) const
  {
    thermodynamics_.pressureAndReaction( u, p, T );
  }

  inline void conservativeToPrimitive( const DomainType& xgl,
                                       const RangeType& cons, 
                                       RangeType& prim ) const
  {
    thermodynamics_.conservativeToPrimitive( cons, prim );
  }
  inline void totalEnergy( const DomainType& xgl,
			   const RangeType& cons, 
			   const GradientRangeType& grad,
			   double& res ) const
  {
    
    thermodynamics_.totalEnergy(cons, grad,res );
  }


  inline const ProblemType& problem() const { return problem_; }
  inline const double mu( const double T ) const { return problem_.mu(T); }
  inline const double delta( ) const {return problem_.delta();}
	    inline const double visc( ) const {return visc_;}
 

	    protected:
	    const ThermodynamicsType& thermodynamics_;
	    const ProblemType& problem_;
	    const PhaseFlux< dimDomain > nsFlux_;
	    const double visc_;
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
  abort();
  DomainType xgl=it.intersectionGlobal().global(x);
  const typename Traits :: DomainType normal = it.integrationOuterNormal(x); 
  double p;
  double T;
  pressAndTemp( uLeft, p, T );
  gLeft = 0;

  // bnd. cond. from euler part
  for (int i=0;i<dimDomain; ++i)
    gLeft[i+1] = 0;


  return 0.;
}



template< class GridPartType, class ProblemImp >
inline void PhaseModel< GridPartType, ProblemImp >
:: boundaryValue( const IntersectionType& it
                  , const double time
                  , const FaceDomainType& x
                  , const RangeType& uLeft
                  , RangeType& uRight ) const 
{
  
  DomainType xgl = it.geometry().global( x );
  // uRight[0]=uLeft[0];     
  
  for(int i=0;i<dimDomain+1;i++)
    uRight[i]=0.;      
      
  uRight[0]=uLeft[0];
  uRight[dimDomain+1]=uLeft[dimDomain+1];
  // uRight[dimDomain+1]=uLeft[dimDomain+1];     
 
}

} // end namespace Dune

#endif // file definition
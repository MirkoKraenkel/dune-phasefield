#ifndef ALLENCAHNMODEL_HH
#define ALLENCAHNMODEL_HH

// DUNE includes
#include <dune/common/version.hh>
#include <dune/fem/misc/fmatrixconverter.hh>

// local includes
#include <dune/fem-dg/models/defaultmodel.hh>

namespace Dune {

// forward declaration of PhaseFlux
//template< int >
//struct PhaseFlux;

//////////////////////////////////////////////////////
//
//                 NAVIER-STOKES EQNS 
//
//////////////////////////////////////////////////////
template< class GridPart > 
class AllenCahnModelTraits 
{
  public:
  typedef GridPart GridPartType;
  typedef typename GridPart :: GridType                     GridType;
  enum { dimDomain = GridType :: dimensionworld };
  enum { dimRange = 1 };
  //enum { dimScalGradRange = dimDomain };
	enum { dimGradRange = dimRange * dimDomain };
	enum { dimGradient = dimDomain };
  
  typedef FieldVector< double, dimDomain >                  DomainType;
  typedef FieldVector< double, dimDomain - 1 >              FaceDomainType;
  typedef FieldVector< double, dimRange >                   RangeType;
  typedef FieldVector< double, dimGradRange >               GradientType;
  typedef FieldMatrix< double, dimRange, dimDomain >        FluxRangeType;
  typedef FieldVector< double, dimGradRange >               GradientRangeType;
	typedef FieldMatrix< double, dimRange, dimDomain >        JacobianRangeType;
  typedef FieldMatrix< double, dimGradRange, dimDomain >    JacobianFluxRangeType;

  typedef typename GridPart :: IntersectionIteratorType     IntersectionIterator;
  typedef typename IntersectionIterator::Intersection       Intersection;
  typedef Intersection       IntersectionType;
  typedef typename GridType :: template Codim<0> :: Entity  EntityType;
  typedef typename EntityType :: EntityPointer              EntityPointerType;

};



template< class GridPartType , class Problem >
class AllenCahnModel : public DefaultModel < AllenCahnModelTraits< GridPartType > >
{
  public:
  typedef Problem                                           ProblemType;
  typedef typename GridPartType::GridType                   GridType;
  typedef AllenCahnModelTraits< GridPartType >                  Traits;

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

 public:
  AllenCahnModel( const ProblemType& problem ) 
    : problem_( problem ),
      delta_(Dune::Fem::Parameter::getValue<double>("phasefield.delta")) 
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

    DomainType xgl = en.geometry().global( x );
 
    double deltainv=1/delta_;
 


    s[0]=u[0]*u[0]*u[0]-u[0];
  	
    return deltainv;
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
     abort(); 
  }

  inline double diffusionTimeStep( const IntersectionType& it,
                                   const double enVolume,
                                   const double circumEstimate,
                                   const double time,
                                   const FaceDomainType& x,
                                   const RangeType& u ) const
  {
    // look at Ch. Merkle Diplom thesis, pg. 38
    return 1/delta_; 
  }


  inline void jacobian( const EntityType& en
                        , const double time
                        , const DomainType& x
                        , const RangeType& u 
                        , JacobianFluxRangeType& a ) const 
  {
    a=0;
    for(int i=0; i<dimDomain; ++i)
      a[i][i]=u[0];
    
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

  template <class JacobianRangeImp>
  inline double diffusionBoundaryFlux( const IntersectionType& it,
                                       const double time,
                                       const FaceDomainType& x,
                                       const RangeType& uLeft,
                                       const JacobianRangeImp& jacLeft,
                                       RangeType& gLeft ) const  
  {
    abort();
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
    abort();
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

    for(int i = 0; i<dimDomain; ++i)
    diff[i][i]=delta_*jac[i][i];
  
  
  }
   void boundarydiffusion( const EntityType& en,
                  const double time,
                  const DomainType& x,
                  const RangeType& u,
                  const GradientRangeType& v,
		              JacobianRangeType& diff ) const
  {
    Dune::Fem:: FieldMatrixConverter< GradientRangeType, JacobianRangeType> jac( v );
    boundarydiffusion( en, time, x, u, jac,diff );
  }


  template <class JacobianRangeImp>
  void boundarydiffusion( const EntityType& en,
                  const double time,
                  const DomainType& x,
                  const RangeType& u,
	            	  const JacobianRangeImp& jac,
		              JacobianRangeType& diff ) const
  { 
  diff=0.;
  }
  
  
 
 

  inline double boundaryFlux( const IntersectionType& it
                   , const double time
                   , const FaceDomainType& x
                   , const RangeType& uLeft
		   , RangeType& gLeft ) const  
  {
   abort();
    gLeft = 0.;
    return 0.;
  }



 
  
  inline void totalEnergy( const DomainType& xgl,
		                  	   const RangeType& cons, 
			                     const GradientRangeType& grad, 
                           FieldVector<double,1>& kin,
			                     FieldVector<double,1>& total ) const
  {
    Fem::FieldMatrixConverter< GradientRangeType, JacobianRangeType> jac( grad);
    totalEnergy(cons, jac,kin,total );
  }
  template<class JacobianRangeImp>
  inline void totalEnergy( const RangeType& cons, 
			                     const JacobianRangeImp& grad, 
                           FieldVector<double,1>& kin,
			                     FieldVector<double,1>& total ) const
  {
   abort();
  }

  inline double delta()const
  {
    return delta_;
  }

  
///Data Members 
 protected:
  const ProblemType problem_;
  const double delta_;

};

/////////////////////////////////////////
//Implementations
/////////////////////////////////////////
template< class GridPartType, class ProblemImp >
inline double AllenCahnModel< GridPartType, ProblemImp >
:: boundaryFlux( const IntersectionType& it
                 , const double time
                 , const FaceDomainType& x
                 , const RangeType& uLeft
                 , const GradientRangeType& duLeft
                 , RangeType& gLeft ) const  
{
abort();


  return 0.;
}



template< class GridPartType, class ProblemImp >
inline void AllenCahnModel< GridPartType, ProblemImp >
:: boundaryValue( const IntersectionType& it
                  , const double time
                  , const FaceDomainType& x
                  , const RangeType& uLeft
                  , RangeType& uRight ) const 
{
  uRight[0]=uLeft[0];
}

} // end namespace Dune

#endif // file definition

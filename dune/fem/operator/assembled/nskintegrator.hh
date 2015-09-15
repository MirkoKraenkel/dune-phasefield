#ifndef PHASEFIELD_Nsk_INTEGRATOR_HH
#define PHASEFIELD_Nsk_INTEGRATOR_HH
#include "nskfilter.hh"
template< class DiscreteFunction, class Model, class Flux>
class NskIntegrator
{ 
  public:
  using DiscreteFunctionType=DiscreteFunction;
  using ModelType=Model;
  using NumericalFluxType=Flux;
  using DiscreteFunctionSpaceType=typename DiscreteFunctionType::DiscreteFunctionSpaceType;
  
  typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType;
  typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;
  typedef typename LocalFunctionType::RangeType RangeType;
  typedef typename LocalFunctionType::JacobianRangeType JacobianRangeType;

  typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
  typedef typename IteratorType::Entity       EntityType;
  typedef typename EntityType::EntityPointer  EntityPointerType;

  typedef typename EntityType::Geometry       GeometryType;

  typedef typename DiscreteFunctionSpaceType::DomainType DomainType; 

  typedef typename DiscreteFunctionSpaceType::GridPartType  GridPartType;
  typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
  typedef typename IntersectionIteratorType::Intersection IntersectionType;
  typedef typename IntersectionType::Geometry  IntersectionGeometryType;

  typedef Dune::Fem::CachingQuadrature< GridPartType, 1 > FaceQuadratureType;

  typedef Dune::Fem::CachingQuadrature< GridPartType, 0 > QuadratureType;

  typedef Dune::Fem::TemporaryLocalFunction<DiscreteFunctionSpaceType> TemporaryLocalType;

  static const int dimDomain = LocalFunctionType::dimDomain;
  static const int dimRange = LocalFunctionType::dimRange;

  typedef  PhasefieldFilter<RangeType> Filter; 



  public:
    NskIntegrator( const ModelType& model,
                              const DiscreteFunctionSpaceType& space):
                              model_( model),
                              flux_( model ),
                              theta_(Dune::Fem::Parameter::getValue<double>("phasefield.mixed.theta")),
                              time_(0.),
                              deltaT_(0.),
                              deltaTInv_(0.),
                              maxSpeed_(1.),
                              lastSpeed_(1.),
                              uOld_("uOld" , space ),
                              uOldLocal_(space),
                              uOldNeighbor_(space),
                              outflow_(Dune::Fem::Parameter::getValue<bool>("phasefield.outflow")),
                              minArea_( std::numeric_limits<double>::max() ),
                              rescale_(Dune::Fem::Parameter::getValue<bool>("phasefield.rescale")),
                              precFactor_(1.)
      {
        uOld_.clear();
        factorImp_=0.5*(1+theta_);
        factorExp_=0.5*(1-theta_);
      }
  
  
  //prepares data on a given entity 
  void setEntity ( const EntityType& entity ) 
  {
    uOldLocal_.init(entity);
    uOldLocal_.assign(uOld_.localFunction( entity ));
    areaEn_=entity.geometry().volume();
    minArea_=std::min( areaEn_ , minArea_ );
  }
  
  //prepares data on a given entity
  void setNeighbor ( const EntityType& entity ) 
  {
    uOldNeighbor_.init(entity);
    uOldNeighbor_.assign(uOld_.localFunction( entity ));
    areaNb_=entity.geometry().volume();
    minArea_=std::min( areaNb_ , minArea_ ); 
  }
  
  void setTime ( const double time) { time_=time; }
  void setDeltaT ( const double deltat)
  { 
    deltaT_=deltat;

    deltaTInv_=1./deltaT_;
    if( rescale_ )
      precFactor_*=deltaTInv_;
  }

  double timeStepEstimate () { return std::min( 1./maxSpeed_,lipschitzC()); }

  double maxSpeed() { return maxSpeed_; }
  double lastSpeed() { return lastSpeed_;}
  double lipschitzC() { return 1; }

  void setSpeed() { lastSpeed_=maxSpeed_; maxSpeed_=0;}

  void setPreviousTimeStep( DiscreteFunctionType& uOld)  { uOld_.assign(uOld);} 
 
  void localIntegral( size_t  pt,
                      const GeometryType& geometry,
                      const QuadratureType& quadrature,
                      RangeType& vu,
                      JacobianRangeType& du,
                      RangeType& avu, // to be added to the result local function
                      JacobianRangeType& advu) const;

  template< class IntersectionQuad>
  void intersectionIntegral ( const IntersectionType& intersection,
                              const size_t pt,
                              const IntersectionQuad& quadInside,
                              const IntersectionQuad& quadOutside,
                              const RangeType& vuEn,
                              const RangeType& vuNb,
                              const JacobianRangeType& duEn,
                              const JacobianRangeType& duNb,
                              RangeType& avuLeft,
                              RangeType& avuRight,
                              JacobianRangeType& aduLeft,
                              JacobianRangeType& aduRight) const;

  void boundaryIntegral( const IntersectionType& intersection,
                        const size_t pt,  
                        const FaceQuadratureType& quadInside,
                        const RangeType& vuEn,
                        const JacobianRangeType& duEn,
                        RangeType& avuLeft,
                        JacobianRangeType& aduLeft) const;

  protected:
  ModelType model_;    
  NumericalFluxType flux_;
  const double  theta_;
  double time_;
  double deltaT_;
  double deltaTInv_;
  mutable double maxSpeed_;
  mutable double lastSpeed_;
  DiscreteFunctionType uOld_;
  TemporaryLocalType uOldLocal_;
  TemporaryLocalType uOldNeighbor_;
  const bool outflow_;
  mutable double minArea_;
  mutable double areaEn_;
  mutable double areaNb_;
  double factorImp_;
  double factorExp_;
  bool rescale_;
  double precFactor_;

};






template<class DiscreteFunction, class Model, class Flux >
void NskIntegrator<DiscreteFunction, Model,Flux>
::localIntegral( size_t  pt,
    const GeometryType& geometry,
    const QuadratureType& quadrature,
    RangeType& vu,
    JacobianRangeType& du,
    RangeType& avu, // to be added to the result local function
    JacobianRangeType& adu) const
{

  const typename QuadratureType::CoordinateType &x = quadrature.point( pt );
  const double weight = quadrature.weight( pt )* geometry.integrationElement( x );
  const DomainType xgl = geometry.global(x);
  RangeType vuOld, vuMid(0.);

  //this should stay instide local Integral as it is operator specific
  uOldLocal_.evaluate( quadrature[ pt ], vuOld); 

  double deltaInv=1./deltaT_;
  //(1+theta)/2*U^n+(1-theta)/2*U^(n-1)
  vuMid.axpy(factorImp_,vu);
  vuMid.axpy(factorExp_,vuOld);

  JacobianRangeType duOld,duMid(0.);
  uOldLocal_.jacobian( quadrature[ pt ], duOld);

  //(1+theta)/2*DU^n+(1-theta)/2*DU^(n-1)
  duMid.axpy(factorImp_,du);
  duMid.axpy(factorExp_,duOld);
  double precInv=1.;
  
  RangeType  source(0.);
  //model_.systemSource(time_,deltaT_,vuOld, xgl, source);
  //rho------------------------------------------------------------- 
  //d_t rho=(rho^n-rho^n-1)/delta t
  Filter::rho(avu)=Filter::rho(vu);
  Filter::rho(avu)-=Filter::rho(vuOld);
  Filter::rho(avu)*=deltaInv;

  RangeFieldType div(0.),gradrhodotv(0.);
  //div(rho v)=rho*div v+gradrho v
  for(int ii = 0; ii <dimDomain ; ++ii )
  { 
    //sum d_i v_i 
    div+=Filter::dvelocity( duMid, ii , ii );
    // sum d_i rho*v_i
    gradrhodotv+=Filter::drho( duMid , ii )*Filter::velocity( vuMid, ii );
  }

  Filter::rho(avu)+=div*Filter::rho(vuMid)+gradrhodotv;
  //---------------------------------------------------------------

  //v--------------------------------------------------------------

  for( int ii = 0; ii < dimDomain ; ++ii )
  {
    //dtv
    Filter::velocity(avu,ii)=Filter::velocity(vu,ii);
    Filter::velocity(avu,ii)-=Filter::velocity(vuOld,ii);
    Filter::velocity(avu,ii)*=deltaInv;
    RangeFieldType sgradv(0);

    //sum_j v_j( d_j v_i - d_i v_j)
    for(int jj = 0;jj < dimDomain ; ++jj)
      sgradv+=Filter::velocity( vuMid , jj )
              *(Filter::dvelocity(duMid, ii , jj )-Filter::dvelocity( duMid, jj , ii));

    Filter::velocity( avu , ii )+=sgradv;
    Filter::velocity( avu , ii )+=precInv*Filter::dmu( duMid, ii);
    Filter::velocity( avu , ii )*=Filter::rho( vuMid);

  }
  // A(dv) 
  model_.diffusion( vuMid, duMid , adu );
  //------------------------------------------------------------------

 //mu-----------------------------------------------------------------
  
  double dFdrho;
  model_.muSource(Filter::rho(vu),Filter::rho(vuOld),dFdrho);

  Filter::mu(avu)=precFactor_*Filter::mu( vuMid );
  Filter::mu(avu)-=dFdrho;
  RangeFieldType usqr(0.) , uOldsqr(0.) , sigmasqr(0.) , sigmaOldsqr(0.);

  for( int ii = 0; ii < dimDomain ; ++ii) 
  {
    // |v^n|^2
    usqr+=Filter::velocity( vu , ii )*Filter::velocity( vu , ii );

    // |v^{n-1}|^2
    uOldsqr+=Filter::velocity( vuOld , ii )*Filter::velocity( vuOld , ii );
  }

  Filter::mu(avu)-=0.25*(usqr+uOldsqr);
  RangeFieldType divsigma(0.), gradrhosigma(0.);

  for( int ii = 0 ; ii < dimDomain ; ++ii) 
  {
    divsigma+=precInv*Filter::dsigma( duMid, ii , ii );
  }
  
  Filter::mu( avu )+=model_.delta()*divsigma;
  //-------------------------------------------------------------------

  //sigma--------------------------------------------------------------
  //\sigma-\nabla\phi
  for( int ii = 0 ; ii < dimDomain ; ++ii) 
  {
    //sigma^n
    Filter::sigma( avu , ii )=precFactor_*Filter::sigma( vu , ii );
    //\nabla\rho^n
    Filter::sigma( avu , ii )-=Filter::drho( du , ii );
 }
  
  //------------------------------------------------------------------        
  
  for(int ii = 0; ii < dimRange ; ii++)
  {
    assert( avu[ii]==avu[ii]) ;
  }
  avu-=source;
  avu*=weight;
  adu*=weight;

}



template<class DiscreteFunction, class Model, class Flux>
template< class IntersectionQuad>
void NskIntegrator<DiscreteFunction, Model,Flux>
::intersectionIntegral( const IntersectionType& intersection,
                        const size_t pt,  
                        const IntersectionQuad& quadInside,
                        const IntersectionQuad& quadOutside,
                        const RangeType& vuEn,
                        const RangeType& vuNb,
                        const JacobianRangeType& duEn,
                        const JacobianRangeType& duNb,
                        RangeType& avuLeft,
                        RangeType& avuRight,
                        JacobianRangeType& aduLeft,
                        JacobianRangeType& aduRight) const
{
  RangeType vuOldEn(0.),vuMidEn(0.),vuOldNb(0.),vuMidNb(0.);
  RangeType vuMidEn2(0.),vuMidNb2(0.);
  JacobianRangeType duOldEn(0.),duOldNb(0.),duMidEn(0.), duMidNb(0.);
  JacobianRangeType duMidEn2(0.),duMidNb2(0.);
  //calc vuOldEn....
  uOldLocal_.evaluate( quadInside[ pt ] , vuOldEn );
  uOldLocal_.jacobian( quadInside[ pt ] , duOldEn );
  uOldNeighbor_.evaluate( quadOutside[ pt ] , vuOldNb );
  uOldNeighbor_.jacobian( quadOutside[ pt ] , duOldNb );

  vuMidEn.axpy( factorImp_ , vuEn );
  vuMidEn.axpy( factorExp_ , vuOldEn );

  vuMidNb.axpy( factorImp_ , vuNb );
  vuMidNb.axpy( factorExp_ , vuOldNb);

  duMidEn.axpy( factorImp_ , duEn );
  duMidEn.axpy( factorExp_ , duOldEn );

  duMidNb.axpy( factorImp_ , duNb );
  duMidNb.axpy( factorExp_ , duOldNb );

  const typename FaceQuadratureType::LocalCoordinateType &x = quadInside.localPoint( pt );
  const DomainType normal = intersection.integrationOuterNormal( x );

  // compute penalty factor
  const double intersectionArea = normal.two_norm();
  DomainType unitnormal=normal;
  unitnormal/=intersectionArea;
  const double localMinArea=std::min(areaNb_,areaEn_);
  const double penaltyFactor = intersectionArea / localMinArea;
  const double localMaxSpeed =  std::max( model_.maxSpeed( unitnormal , vuOldEn ), model_.maxSpeed( unitnormal , vuOldNb ) );

  double viscFactor=0.5*localMinArea*lastSpeed_;
  maxSpeed_=std::max(localMaxSpeed*intersectionArea/localMinArea,maxSpeed_);

  JacobianRangeType dvalue(0.),advalue(0.);
  double fluxRet;
  fluxRet=flux_.numericalFlux( normal,
                              viscFactor,
                              penaltyFactor,
                              vuEn,
                              vuNb,
                              vuMidEn,  
                              vuMidNb, 
                              avuLeft,
                              avuRight); 
 
  RangeType value(0);
  
  fluxRet+=flux_.diffusionFlux( normal,
                                penaltyFactor,
                                vuMidEn,
                                vuMidNb,
                                duMidEn,
                                duMidNb,
                                value,
                                aduLeft);

  avuLeft+=value;
  avuRight-=value; 
  aduRight=aduLeft;
  aduRight*=-1.;  
}

//Boundary Intgral
template<class DiscreteFunction, class Model, class Flux>
void NskIntegrator<DiscreteFunction, Model,Flux>
::boundaryIntegral( const IntersectionType& intersection,
                    const size_t pt,  
                    const FaceQuadratureType& quadInside,
                    const RangeType& vuEn,
                    const JacobianRangeType& duEn,
                    RangeType& avuLeft,
                    JacobianRangeType& aduLeft) const
{

  size_t boundaryIndex=intersection.boundaryId();
  typedef typename IntersectionType::Geometry  IntersectionGeometryType;
  const IntersectionGeometryType &intersectionGeometry = intersection.geometry();

  RangeType vuOldEn(0.),vuMidEn(0.), bndValue(0.);
  JacobianRangeType duOldEn(0.),duMidEn(0.);

  //calc vuOldEn....
  uOldLocal_.evaluate( quadInside[ pt ] , vuOldEn );
  uOldLocal_.jacobian( quadInside[ pt ] , duOldEn );

  vuMidEn.axpy( factorImp_ , vuEn );
  vuMidEn.axpy( factorExp_ , vuOldEn );


  duMidEn.axpy( factorImp_ , duEn );
  duMidEn.axpy( factorExp_ , duOldEn );


  const typename FaceQuadratureType::LocalCoordinateType &x = quadInside.localPoint( pt );
  const DomainType normal = intersection.integrationOuterNormal( x );

  // compute penalty factor
  const double intersectionArea = normal.two_norm();

  const double penaltyFactor = intersectionArea /  areaEn_; 
  const double localMaxSpeed = model_.maxSpeed(normal,vuOldEn);

  double viscFactor =  areaEn_*lastSpeed_;
  const DomainType xgl = intersectionGeometry.global(x);
  model_.dirichletValue( time_,xgl, bndValue);

  JacobianRangeType dvalue(0.),advalue(0.);
  double fluxRet;
  RangeType gLeft(0.),dummy(0.);
#if 1 
  if( boundaryIndex==1 || !outflow_)
    {
      fluxRet=flux_.boundaryFlux( normal,
                                  viscFactor,
                                  vuEn,
                                  vuMidEn,gLeft);
    }
  else
    { 
    std::cout<<"implement outflowflux\n";
    abort();
#if 0
flux_.outFlowFlux( normal,
                         viscFactor,
                         vuEn,
                         vuMidEn,gLeft);
#endif
    }
#else
   fluxRet=flux_.numericalFlux( normal,
                                area,
                                vuEn,
                                bndValue,
                                vuMidEn,
                                bndValue,
                                gLeft,
                                dummy);
#endif
  avuLeft+=gLeft;
  RangeType value(0.);

  if( boundaryIndex==1 || !outflow_)
    {
      fluxRet+=flux_.diffusionBoundaryFlux( normal,
                                            penaltyFactor,
                                            vuMidEn,
                                            duMidEn,
                                            value,
                                            aduLeft);
    }
  else
    {
      std::cout<<"implement outflowflux\n";
      abort();
#if 0
      flux_.diffusionOutFlowFlux( value , aduLeft);
#endif
    }

  avuLeft+=value;



}




#endif





#ifndef PHASEFIELD_NAVIERSTOKESINTEGRATOR_HH
#define PHASEFIELD_NAVIERSTOKESINTEGRATOR_HH
#include "splitfilter.hh"
template< class DiscreteFunction, class AddFunction, class Model, class Flux>
class PhasefieldNavierStokesIntegrator
{ 
  public:
  using DiscreteFunctionType=DiscreteFunction;
  using AddFunctionType=AddFunction;
  using ModelType=Model;
  using NumericalFluxType=Flux;
  using DiscreteFunctionSpaceType=typename DiscreteFunctionType::DiscreteFunctionSpaceType;
  
  typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType;
  typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;
  typedef typename LocalFunctionType::RangeType RangeType;
  typedef typename LocalFunctionType::JacobianRangeType JacobianRangeType;

  using AddLocalFunctionType=typename AddFunctionType::LocalFunctionType;
  using AddRangeType=typename AddLocalFunctionType::RangeType;
  using AddJacobianRangeType=typename AddLocalFunctionType::JacobianRangeType;
  using AddDiscreteFunctionSpaceType=typename AddFunctionType::DiscreteFunctionSpaceType;


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
  typedef Dune::Fem::TemporaryLocalFunction<AddDiscreteFunctionSpaceType> TemporaryAddLocalType;


  static const int dimDomain = LocalFunctionType::dimDomain;
  static const int dimRange = LocalFunctionType::dimRange;
  static const int combinedRange= dimRange+dimRange;

  typedef Dune::FieldVector<RangeFieldType, combinedRange> CombinedRangeType;

  typedef  NvStFilter<RangeType> Filter; 
  typedef  AcFilter< AddRangeType> AddFilter;


  public:
    PhasefieldNavierStokesIntegrator( const ModelType& model,
                              const DiscreteFunctionSpaceType& space):
                              model_( model),
                              flux_( model ),
                              theta_(Dune::Fem::Parameter::getValue<double>("phasefield.mixed.theta")),
                              time_(0.),
                              deltaT_(0.),
                              deltaTInv_(0.),
                              maxSpeed_(0.),
                              lastSpeed_(1.),
                              uOld_("uOld" , space ),
                              uOldLocal_(space),
                              uOldNeighbor_(space),
                              addFunction_("uAdd",space),
                              addLocal_(space),
                              addNeighbor_(space),
                              outflow_(Dune::Fem::Parameter::getValue<bool>("phasefield.outflow")),
                              minArea_( std::numeric_limits<double>::max() )
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
    addLocal_.init(entity);
    addLocal_.assign(addFunction_.localFunction( entity ));
    areaEn_=entity.geometry().volume();
    minArea_=std::min( areaEn_ , minArea_ );
    
  }
  
  //prepares data on a given entity
  void setNeighbor ( const EntityType& entity ) 
  {
    uOldNeighbor_.init(entity);
    uOldNeighbor_.assign(uOld_.localFunction( entity ));
    addNeighbor_.init(entity);
    addNeighbor_.assign(addFunction_.localFunction( entity ));
    areaNb_=entity.geometry().volume();
    minArea_=std::min( areaNb_ , minArea_ ); 
  }
  
  void setTime( const double time) { time_=time;}
  void setDeltaT( const double deltat) {deltaT_=deltat;}

  double timeStepEstimate() { return std::min( minArea_/maxSpeed_,lipschitzC());}

  double maxSpeed() { return maxSpeed_; }
  
  double lipschitzC() { return model_.lipschitzC(); }

  void setSpeed() { lastSpeed_=maxSpeed_; maxSpeed_=0;}

  void setPreviousTimeStep( DiscreteFunctionType& uOld)  { uOld_.assign(uOld); }
  void setAddVariables( AddFunctionType& addF) {addFunction_.assign(addF);} 
  
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
  double factorImp_;
  double factorExp_;
  DiscreteFunctionType uOld_;
  TemporaryLocalType uOldLocal_;
  TemporaryLocalType uOldNeighbor_;
  AddFunctionType addFunction_;
  TemporaryAddLocalType addLocal_;
  TemporaryAddLocalType addNeighbor_;
  const bool outflow_;
  mutable double minArea_;
  mutable double areaEn_;
  mutable double areaNb_;
  };






template<class DiscreteFunction, class AddFunction,class Model, class Flux >
void PhasefieldNavierStokesIntegrator<DiscreteFunction,AddFunction, Model,Flux>
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
  RangeType vuOld,vuMid(0);
  AddRangeType vuAdd; 
  
  //this should stay instide local Integral as it is operator specific
  uOldLocal_.evaluate( quadrature[ pt ], vuOld); 
  addLocal_.evaluate( quadrature[ pt ], vuAdd);
  double deltaInv=1./deltaT_;
  //(1+theta)/2*U^n+(1-theta)/2*U^(n-1)
  vuMid.axpy(factorImp_,vu);
  vuMid.axpy(factorExp_,vuOld);
  JacobianRangeType duOld,duMid ;
  AddJacobianRangeType duAdd;
  uOldLocal_.jacobian( quadrature[ pt ], duOld);
  addLocal_.jacobian( quadrature[ pt ], duAdd);
  //(1+theta)/2*DU^n+(1-theta)/2*DU^(n-1)
  duMid.axpy(factorImp_,du);
  duMid.axpy(factorExp_,duOld);



  CombinedRangeType  source(0.);
  model_.systemSource(time_, xgl, source);        
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
    div+=Filter::dvelocity(duMid, ii , ii );
    // sum d_i rho*v_i
    gradrhodotv+=Filter::drho(duMid , ii )*Filter::velocity(vuMid, ii );
  }

  Filter::rho(avu)+=div*Filter::rho(vuMid)+gradrhodotv;
  //---------------------------------------------------------------

  //v--------------------------------------------------------------

  for( int ii = 0; ii < dimDomain ; ++ii )
  {
    Filter::velocity(avu,ii)=Filter::velocity(vu,ii);
    Filter::velocity(avu,ii)-=Filter::velocity(vuOld,ii);
    Filter::velocity(avu,ii)*=deltaInv;
    RangeFieldType sgradv(0);

    //sum_j v_j( d_j v_i - d_i v_j)
    for(int jj = 0;jj < dimDomain ; ++jj)
      sgradv+=Filter::velocity( vuMid , jj )
              *(Filter::dvelocity(duMid, ii , jj )-Filter::dvelocity( duMid, jj , ii));

    Filter::velocity( avu , ii )+=sgradv;
    Filter::velocity( avu , ii )+=Filter::dmu( du, ii);
    Filter::velocity( avu , ii )*=Filter::rho( vuMid);

    //-tau\nabla phi
    Filter::velocity( avu , ii )-=AddFilter::tau( vuAdd )*AddFilter::dphi( duAdd , ii );
  }
  // A(dv)
  
  model_.diffusion( duMid , adu );
  //------------------------------------------------------------------

 //mu-----------------------------------------------------------------
  
  double dFdrho;
  //model_.muSource(Filter::rho(vu),Filter::rho(vu),Filter::phi(vu),dFdrho);
  //old version like Paris talk
  //  model_.muSource(Filter::rho(vuOld),Filter::rho(vu),Filter::phi(vu),dFdrho);
  model_.muSource(Filter::rho(vu),Filter::rho(vuOld),AddFilter::phi(vuAdd),dFdrho);
 
  
  Filter::mu(avu)=Filter::mu( vu );
  Filter::mu(avu)-=dFdrho;
  RangeFieldType usqr(0.) , uOldsqr(0.) , sigmasqr(0.) , sigmaOldsqr(0.);

  for( int ii = 0; ii < dimDomain ; ++ii) 
  {
    // |v^n|^2
    usqr+=Filter::velocity( vu , ii )*Filter::velocity( vu , ii );

    // |v^{n-1}|^2
    uOldsqr+=Filter::velocity( vuOld , ii )*Filter::velocity( vuOld , ii );
#if RHOMODEL
    //\sigma^n*\sigma^{n-1}
    sigmasqr+=AddFilter::sigma( vuAdd , ii )*AddFilter::sigma( vuAdd , ii );
#endif
  }

  Filter::mu(avu)-=0.25*(usqr+uOldsqr);
#if RHOMODEL
#if DIFFQUOT
  double rhodiff=(Filter::rho(vu)-Filter::rho(vuOld));
  
  if( std::abs(rhodiff)<1e-9)
    Filter::mu(avu)-=0.5*model_.delta()*model_.h2prime(Filter::rho(vuOld))*sigmasqr;
  else
    Filter::mu(avu)+=0.5*model_.delta()*(1/rhodiff)*(model_.h2(Filter::rho(vu))-model_.h2(Filter::rho(vuOld)))*sigmasqr;
#else
    Filter::mu(avu)-=0.5*model_.delta()*model_.h2prime(Filter::rho(vuOld))*sigmasqr;
#endif
#endif
  //------------------------------------------------------------------

  
 
  for(int ii = 0; ii < dimRange ; ii++)
  {
    assert( avu[ii]==avu[ii]) ;
  }
  //avu-=source;
  avu*=weight;
  adu*=weight;

}



template<class DiscreteFunction,class AddFunction, class Model, class Flux>
template< class IntersectionQuad>
void PhasefieldNavierStokesIntegrator<DiscreteFunction,AddFunction, Model,Flux>
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
  typedef typename IntersectionType::Geometry  IntersectionGeometryType;

  RangeType vuOldEn(0.),vuMidEn(0.),vuOldNb(0.),vuMidNb(0.),vuAddEn(0.),vuAddNb(0.);
  JacobianRangeType duOldEn(0.),duOldNb(0.),duMidEn(0.), duMidNb(0.), duAddEn(0.),duAddNb(0.);

  //calc vuOldEn....
  uOldLocal_.evaluate( quadInside[ pt ] , vuOldEn );
  uOldLocal_.jacobian( quadInside[ pt ] , duOldEn );
  uOldNeighbor_.evaluate( quadOutside[ pt ] , vuOldNb );
  uOldNeighbor_.jacobian( quadOutside[ pt ] , duOldNb );
  
  addLocal_.evaluate( quadInside[ pt ] , vuAddEn );
  addNeighbor_.evaluate( quadInside[ pt ] , vuAddNb );


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
  const double penaltyFactor = intersectionArea / std::min( areaEn_, areaNb_ ); 
  const double area=lastSpeed_*std::min(areaEn_,areaNb_)/intersectionArea; 

  maxSpeed_ = 0;//std::max( std::max( model_.maxSpeed(normal,vuOldEn), model_.maxSpeed(normal,vuOldNb)),maxSpeed_);

  JacobianRangeType dvalue(0.),advalue(0.);
  double fluxRet;
  fluxRet=flux_.numericalFlux( normal,
                              area,
                              penaltyFactor,
                              vuEn,
                              vuNb,
                              vuMidEn,  
                              vuMidNb, 
                              vuAddEn,
                              vuAddNb,
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
template<class DiscreteFunction, class AddFunction,class Model, class Flux>
void PhasefieldNavierStokesIntegrator<DiscreteFunction,AddFunction, Model,Flux>
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
  AddRangeType vuAddEn(0.);
  //calc vuOldEn....
  uOldLocal_.evaluate( quadInside[ pt ] , vuOldEn );
  uOldLocal_.jacobian( quadInside[ pt ] , duOldEn );

  addLocal_.evaluate( quadInside[ pt ] , vuAddEn );
 
  vuMidEn.axpy( factorImp_ , vuEn );
  vuMidEn.axpy( factorExp_ , vuOldEn );


  duMidEn.axpy( factorImp_ , duEn );
  duMidEn.axpy( factorExp_ , duOldEn );



  // compute penalty factor
  const double intersectionArea = intersectionGeometry.volume();
  const double penaltyFactor = intersectionArea /  areaEn_; 
  const double area=lastSpeed_*areaEn_/intersectionArea; 
  const typename FaceQuadratureType::LocalCoordinateType &x = quadInside.localPoint( pt );
  const DomainType normal = intersection.integrationOuterNormal( x );
  const DomainType xgl = intersectionGeometry.global(x);
//  model_.dirichletValue( time_,xgl, bndValue);

  JacobianRangeType dvalue(0.),advalue(0.);
  double fluxRet;
  RangeType gLeft(0.),dummy(0.);
#if 1 
  if( boundaryIndex==1 || !outflow_)
    {
      fluxRet=flux_.boundaryFlux( normal,
                                  area,
                                  vuEn,
                                  vuMidEn,
                                  gLeft);
    }
  else
    {
      flux_.outFlowFlux( normal,
                         area,
                         vuEn,
                         vuMidEn,
                         gLeft);

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
      flux_.diffusionOutFlowFlux( value , aduLeft);
    }

  avuLeft+=value;



}




#endif





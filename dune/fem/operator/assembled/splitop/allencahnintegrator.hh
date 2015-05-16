#ifndef PHASEFIELD_ALLENCAHNINTEGRATOR_HH
#define PHASEFIELD_ALLENCAHNINTEGRATOR_HH
#include "splitfilter.hh"
#include "../matrixhelper.hh"
template< class DiscreteFunction, class AddFunction ,class Model, class Flux>
class PhasefieldAllenCahnIntegrator
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

  typedef typename EntityType::Geometry       GeometryType;

  typedef typename DiscreteFunctionSpaceType::DomainType DomainType; 

  typedef typename DiscreteFunctionSpaceType::GridPartType  GridPartType;
  typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
  typedef typename IntersectionIteratorType::Intersection IntersectionType;

  typedef Dune::Fem::CachingQuadrature< GridPartType, 1 > FaceQuadratureType;

  typedef Dune::Fem::CachingQuadrature< GridPartType, 0 > QuadratureType;

  typedef Dune::Fem::TemporaryLocalFunction<DiscreteFunctionSpaceType> TemporaryLocalType;
  typedef Dune::Fem::TemporaryLocalFunction<AddDiscreteFunctionSpaceType> TemporaryAddLocalType;


  static const int dimDomain = LocalFunctionType::dimDomain;
  static const int dimRange = LocalFunctionType::dimRange;
  static const int combinedRange= dimRange+dimRange;
  typedef Dune::FieldVector<RangeFieldType, combinedRange> CombinedRangeType;


  typedef  AcFilter<RangeType> Filter; 
  typedef  NvStFilter< AddRangeType > AddFilter;


  public:
    PhasefieldAllenCahnIntegrator( const ModelType& model,
                              const DiscreteFunctionSpaceType& space):
                              model_( model),
                              flux_( model ),
                              theta_(Dune::Fem::Parameter::getValue<double>("phasefield.mixed.theta")),
                              time_(0.),
                              deltaT_(0.),
                              deltaTInv_(0.),
                              maxSpeed_(1.),
                              lastSpeed_(0.),
                              uOld_("uOld" , space ),
                              uOldLocal_(space),
                              uOldNeighbor_(space),
                              addFunction_("uuAdd",space),
                              addLocal_( space ),
                              addNeighbor_( space),
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
    addNeighbor_.assign(addFunction_.localFunction( entity));
    areaNb_=entity.geometry().volume();
    minArea_=std::min( areaNb_ , minArea_ ); 
  }
  
  void setTime( const double time) { time_=time;}
  void setDeltaT( const double deltat) {deltaT_=deltat, deltaTInv_=1./deltaT_;}

  double timeStepEstimate() { return std::min( minArea_/maxSpeed_,lipschitzC());}

  double maxSpeed() { return maxSpeed_; }
  
  double lipschitzC() { return model_.lipschitzC(); }

  void setSpeed() { lastSpeed_=maxSpeed_; maxSpeed_=0;}

  void setPreviousTimeStep( DiscreteFunctionType& uOld)  { uOld_.assign(uOld);} 
  void setAddVariables( AddFunctionType& addF) { addFunction_.assign(addF);} 
 
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
  double theta_;
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






template<class DiscreteFunction,class AddFunction, class Model, class Flux >
void PhasefieldAllenCahnIntegrator<DiscreteFunction,AddFunction,Model,Flux>
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
  RangeType vuOld,vuMid(0.);
  AddRangeType vuAdd;

  //this should stay instide local Integral as it is operator specific
  uOldLocal_.evaluate( quadrature[ pt ], vuOld); 
  addLocal_.evaluate( quadrature[ pt ], vuAdd);
  double deltaInv=1./deltaT_;
  //(1+theta)/2*U^n+(1-theta)/2*U^(n-1)
  vuMid.axpy(factorImp_,vu);
  vuMid.axpy(factorExp_,vuOld);


  JacobianRangeType duOld,duMid(0);
  AddJacobianRangeType duAdd;
  uOldLocal_.jacobian( quadrature[ pt ], duOld);
  addLocal_.jacobian( quadrature[pt],duAdd);
  duMid.axpy(factorImp_,du);
  duMid.axpy(factorExp_,duOld);


  CombinedRangeType  source(0.);
  model_.systemSource(time_, xgl, source);        
  //phi---------------------------------------------------------------

  Filter::phi( avu )=Filter::phi( vu )-Filter::phi( vuOld );
  Filter::phi( avu )*=deltaInv;

  RangeFieldType transport(0.);

  // \nabla phi\cdot v
  for( int ii = 0; ii < dimDomain ; ++ii ) 
  { 
    transport+=AddFilter::velocity( vuAdd , ii )*Filter::dphi(duMid, ii );
  }
  Filter::phi( avu )+=transport+model_.reactionFactor()*Filter::tau( vu )/AddFilter::rho(vuAdd);

  //tau---------------------------------------------------------------
  // dF/dphi
  double dFdphi;
  /*
    model_.tauSource( Filter::phi(vuOld),
                      Filter::phi(vu),
                      Filter::rho(vuOld),
                      dFdphi);
   */
  model_.tauSource( Filter::phi(vu),
                    Filter::phi(vuOld),
                    AddFilter::rho(vuAdd),
                    dFdphi);


  Filter::tau( avu )=Filter::tau( vu );
  Filter::tau( avu )-=dFdphi;

  RangeFieldType divsigma(0.), gradrhosigma(0.);

  for( int ii = 0 ; ii < dimDomain ; ++ii) 
    {
#if RHOMODEL
#if LAMBDASCHEME 
      divsigma+=Filter::dalpha( du, ii , ii );
#else
    #warning "USE RHOMODEL ONLY WITH LAMBDASCHEME"
    #error
    abort();
#endif
#else
      divsigma+=Filter::dsigma( duMid, ii , ii );
#endif
    }
#if RHOMODEL && !LAMBDASCHEME
  Filter::tau( avu )+=model_.delta()*(divsigma+gradrhosigma);
#else
  Filter::tau( avu )+=model_.delta()*divsigma;
#endif
  //-------------------------------------------------------------------

  //sigma--------------------------------------------------------------
  //\sigma-\nabla\phi
  for( int ii = 0 ; ii < dimDomain ; ++ii) 
  {
    //sigma^n
    Filter::sigma( avu , ii )=Filter::sigma( vu , ii );
    //\nabla\phi^n
    Filter::sigma( avu , ii )-=Filter::dphi( du , ii );
#if RHOMODEL && LAMBDASCHEME
    Filter::alpha( avu, ii )=Filter::alpha(vu,ii)-model_.h2(AddFilter::rho(vuAdd))*Filter::sigma(vu, ii);
#endif
 }
  
  //------------------------------------------------------------------        
  
  for(int ii = 0; ii < dimRange ; ii++)
  {
    assert( avu[ii]==avu[ii]) ;
  }
 // avu-=source;
  avu*=weight;
  adu*=weight;

}



template<class DiscreteFunction,class AddFunction, class Model, class Flux>
template< class IntersectionQuad>
void PhasefieldAllenCahnIntegrator<DiscreteFunction,AddFunction, Model,Flux>
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
  RangeType vuOldEn(0.),vuOldNb(0.), vuAddEn(0.),vuAddNb(0.), vuMidEn(0.),vuMidNb(0.);
  JacobianRangeType duOldEn(0.),duOldNb(0.),duAddEn(0.),duAddNb(0.);

  //calc vuOldEn....
  uOldLocal_.evaluate( quadInside[ pt ] , vuOldEn );
  uOldNeighbor_.evaluate( quadOutside[ pt ] , vuOldNb );
  //calc vuAddEn..
  addLocal_.evaluate( quadInside[ pt ] , vuAddEn );
  addNeighbor_.evaluate( quadOutside[ pt ] , vuAddNb );

  vuMidEn.axpy( factorImp_ , vuEn );
  vuMidEn.axpy( factorExp_ , vuOldEn );

  vuMidNb.axpy( factorImp_ , vuNb );
  vuMidNb.axpy( factorExp_ , vuOldNb);


  const typename FaceQuadratureType::LocalCoordinateType &x = quadInside.localPoint( pt );

  const DomainType normal = intersection.integrationOuterNormal( x );
  // compute penalty factor
  const double intersectionArea = normal.two_norm();
  const double penaltyFactor = intersectionArea / std::min( areaEn_, areaNb_ ); 
  const double area=lastSpeed_*std::min(areaEn_,areaNb_)/intersectionArea; 
  CombinedRangeType combinedEn,combinedNb;
  MatrixHelper::concat(vuAddEn,vuOldEn,combinedEn);
  MatrixHelper::concat(vuAddNb,vuOldNb,combinedNb);
  
  maxSpeed_ = std::max( std::max( model_.maxSpeed(normal,combinedEn), model_.maxSpeed(normal,combinedNb)),maxSpeed_);

 

  JacobianRangeType dvalue(0.),advalue(0.);
  flux_.numericalFlux( normal,
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
 
  
}

//Boundary Intgral
template<class DiscreteFunction, class AddFunction, class Model, class Flux>
void PhasefieldAllenCahnIntegrator<DiscreteFunction,AddFunction, Model,Flux>
::boundaryIntegral( const IntersectionType& intersection,
                    const size_t pt,  
                    const FaceQuadratureType& quadInside,
                    const RangeType& vuEn,
                    const JacobianRangeType& duEn,
                    RangeType& avuLeft,
                    JacobianRangeType& aduLeft) const
{

  size_t boundaryIndex=intersection.boundaryId();

  RangeType vuOldEn(0.),vuAddEn(0), bndValue(0.);
  
  //calc vuOldEn....
  uOldLocal_.evaluate( quadInside[ pt ] , vuOldEn );
  addLocal_.evaluate( quadInside[ pt ],vuAddEn);
  





  // compute penalty factor
  const typename FaceQuadratureType::LocalCoordinateType &x = quadInside.localPoint( pt );
  const DomainType normal = intersection.integrationOuterNormal( x );
  //const DomainType xgl = intersectionGeometry.global(x);
  const double intersectionArea = normal.two_norm();
  const double area=lastSpeed_*areaEn_/intersectionArea; 

  JacobianRangeType dvalue(0.),advalue(0.);
  RangeType gLeft(0.),dummy(0.);
  if( boundaryIndex==1 || !outflow_)
    {
      flux_.boundaryFlux( normal,
                          area,
                          vuEn,
                          vuEn,
                          vuAddEn,
                          gLeft);
    }
  else
    {
      flux_.outFlowFlux( normal,
                         area,
                         vuEn,
                         vuEn,
                         vuAddEn,
                         gLeft);

    }
  avuLeft+=gLeft;




}




#endif





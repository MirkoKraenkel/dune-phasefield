#ifndef DUNE_PHASEFIELD_MIXEDOPERATOR_HH
#define DUNE_PHASEFIELD_MIXEDOPERATOR_HH


#warning  "HEATOP2"
//globlas includes

//DUNE includes
#include <dune/common/fmatrix.hh>
#define HEATCHECK 1
#define OPCHECK  0 
#if OPCHECK
#warning "DEBUGGING VERSION"
#endif
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/common/operator.hh>

#include <dune/fem/function/localfunction/temporary.hh>

#include <dune/fem/operator/common/differentiableoperator.hh>
#include <dune/fem/operator/common/automaticdifferenceoperator.hh>
#include <dune/fem/operator/common/stencil.hh>

#include "phasefieldfilter.hh"
#include "flux.hh"
template<class DiscreteFunction, class Model, class Flux>
class DGPhasefieldOperator
: public virtual Dune::Fem::Operator<DiscreteFunction,DiscreteFunction>
{
  typedef Dune::Fem::Operator<DiscreteFunction,DiscreteFunction> BaseType;

protected:

  typedef DiscreteFunction DiscreteFunctionType;
  typedef Model            ModelType;
  typedef Flux             NumericalFluxType;
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
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

  typedef Dune::Fem::ElementQuadrature< GridPartType, 1 > FaceQuadratureType;
  typedef Dune::Fem::CachingQuadrature< GridPartType, 0 > QuadratureType;

  typedef Dune::Fem::TemporaryLocalFunction<DiscreteFunctionSpaceType> TemporaryLocalType;

  static const int dimDomain = LocalFunctionType::dimDomain;
  static const int dimRange = LocalFunctionType::dimRange;

  typedef  PhasefieldFilter<RangeType> Filter; 

public:
   //! constructor
   DGPhasefieldOperator(const ModelType &model,
                        const DiscreteFunctionSpaceType &space,
                        const NumericalFluxType &flux)
   : 
     model_(model),
     space_(space),
     flux_(flux),
     theta_(Dune::Fem::Parameter::getValue<double>("phasefield.mixed.theta")),
     time_(0.),
     deltaT_(0.),
     uOld_("uOld" , space )
    {
      assert(theta_>=0 && theta_<=1);
#if OPCHECK 
      factorImp_=0;
      factorExp_=1;
#else
      factorImp_=0.5*(1+theta_);
      factorExp_=0.5*(1-theta_);
#endif
    }

  // prepare the solution vector 
  template <class Function>
  void prepare( const Function &func, DiscreteFunctionType &u ) 
  { 
  }

  //! application operator 
  void operator() ( const DiscreteFunctionType &u, DiscreteFunctionType &w ) const;

  
  void setTime(const double time)  {time_=time;}
 
  double getTime(){return time_;}
  void setDeltaT( const double deltat) { deltaT_=deltat;}
  
  void setPreviousTimeStep( DiscreteFunctionType& uOld)  { uOld_.assign(uOld);} 
  DiscreteFunctionType& getPreviousTimeStep() { return uOld_;}
 
 protected:
  template< class ArgType, class LocalArgType,class LFDestType >
  void localOp( const EntityType& entity,
                const ArgType& u, 
                const LocalArgType& uLocal, 
                LFDestType& wlocal) const;

  template<class LocalArgType, class NeighborArgType, class LFDestType>
  void computeIntersection( const IntersectionType& intersection,
                            const EntityType& entity,
                            const EntityType& neighbor,
                            const double area,
                            const LocalArgType& uEn, 
                            const NeighborArgType& uNb, 
                            LFDestType& wlocal) const;
  
  template<class LocalArgType, class LFDestType>
  void computeBoundary( const IntersectionType& intersection,
                        const EntityType& entity,
                        const double area,
                        const LocalArgType& uEn, 
                        LFDestType& wlocal) const;

  const DiscreteFunctionSpaceType& space() const {return space_;}
  
  const ModelType& model() const{ return model_;}

  double penalty() const { return model_.penalty();}
  
protected:
  ModelType model_;
  const DiscreteFunctionSpaceType &space_;
  const NumericalFluxType flux_;
  const double  theta_;
  double time_;
  double deltaT_;
  double factorImp_;
  double factorExp_;
  mutable DiscreteFunctionType uOld_;
};





template<class DiscreteFunction, class Model, class Flux, class Jacobian>
class FDJacobianDGPhasefieldOperator
: public DGPhasefieldOperator<DiscreteFunction,Model,Flux>,
  public virtual Dune::Fem::AutomaticDifferenceOperator<DiscreteFunction,DiscreteFunction, Jacobian>
  {

    typedef DGPhasefieldOperator<DiscreteFunction,Model,Flux> MyOperatorType;
    typedef Dune::Fem::AutomaticDifferenceOperator<DiscreteFunction,DiscreteFunction,Jacobian> BaseType;
  
    typedef typename MyOperatorType::DiscreteFunctionType DiscreteFunctionType;
    typedef typename MyOperatorType::ModelType ModelType;
    typedef typename MyOperatorType::NumericalFluxType NumericalFluxType;
protected:
    typedef typename BaseType::RangeFunctionType RangeFunctionType;
    typedef typename BaseType::DomainFunctionType DomainFunctionType;
    typedef typename BaseType::RangeFieldType RangeFieldType;
    typedef typename BaseType::DomainFieldType DomainFieldType;
    typedef typename BaseType::DomainSpaceType DomainSpaceType;
    typedef typename BaseType::RangeSpaceType RangeSpaceType;
    typedef DomainSpaceType DiscreteFunctionSpaceType;

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

    typedef Dune::Fem::ElementQuadrature< GridPartType, 1 > FaceQuadratureType;
    typedef Dune::Fem::CachingQuadrature< GridPartType, 0 > QuadratureType;

    typedef Dune::Fem::TemporaryLocalFunction<DiscreteFunctionSpaceType> TemporaryLocalType;

    static const int dimDomain = LocalFunctionType::dimDomain;
    static const int dimRange = LocalFunctionType::dimRange;

    typedef  PhasefieldFilter<RangeType> Filter; 

public:
   //! constructor
   FDJacobianDGPhasefieldOperator(const ModelType &model,
                                  const DiscreteFunctionSpaceType &space,
                                  const NumericalFluxType &flux)
   :MyOperatorType(model,space,flux){} 
  
 using MyOperatorType::prepare;
  //! application operator 
  using MyOperatorType::operator();
  using MyOperatorType::setTime;
  using MyOperatorType::setDeltaT;
  using MyOperatorType::setPreviousTimeStep;
  using MyOperatorType::getPreviousTimeStep;
 
 protected:

  using MyOperatorType::localOp;
  using MyOperatorType::computeBoundary;
  using MyOperatorType::computeIntersection;
  using MyOperatorType::space;
  using MyOperatorType::model;
  using MyOperatorType::penalty;
  
protected:
  using MyOperatorType::model_;
  using MyOperatorType::space_;
  using MyOperatorType::flux_;
  using MyOperatorType::time_;
  using MyOperatorType::deltaT_;
  using MyOperatorType::theta_;
  using MyOperatorType::factorImp_;
  using MyOperatorType::factorExp_;
  using MyOperatorType::uOld_;
};


template<class DiscreteFunction, class Model, class Flux>
void DGPhasefieldOperator<DiscreteFunction, Model,Flux>
  ::operator() ( const DiscreteFunctionType &u, DiscreteFunctionType &w ) const 
{

  // clear destination 
  w.clear();
  assert(deltaT_>0);
  // iterate over grid 
  const IteratorType end = space().end();
  for( IteratorType it = space().begin(); it != end; ++it )
  {
    // get entity (here element) 
    const EntityType &entity = *it;
   
    // get local representation of the discrete functions 
    const LocalFunctionType uLocal = u.localFunction( entity);
    LocalFunctionType wLocal = w.localFunction( entity );

    localOp(entity,u, uLocal,wLocal);
    
  }
  // communicate data (in parallel runs)
  w.communicate();

}

template<class DiscreteFunction, class Model, class Flux >
template< class ArgType, class LocalArgType,class LFDestType >
void DGPhasefieldOperator<DiscreteFunction, Model,Flux>
::localOp(const EntityType& entity,
          const ArgType& uGlobal, 
          const LocalArgType& uLocal,
          LFDestType& wLocal) const
  {
    // get elements geometry 
    const GeometryType &geometry = entity.geometry();
    
    // obtain quadrature order 
    const int quadOrder = uLocal.order() + wLocal.order();
    LocalFunctionType uOldLocal=uOld_.localFunction(entity); 
    
    { // element integral
      QuadratureType quadrature( entity, quadOrder );
      const size_t numQuadraturePoints = quadrature.nop();
      for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
      {
        const typename QuadratureType::CoordinateType &x = quadrature.point( pt );
        const double weight = quadrature.weight( pt )* geometry.integrationElement( x );
        const DomainType xgl = geometry.global(x);
        RangeType vu,vuOld, vuMid{0};
        uLocal.evaluate( quadrature[ pt ], vu );
        uOldLocal.evaluate( quadrature[ pt ], vuOld); 
        
          double deltaInv=1./deltaT_;
        //(1+theta)/2*U^n+(1-theta)/2*U^(n-1)
        // #if OPCHECK vuMid=vuOld;      
        vuMid.axpy(factorImp_,vu);
        vuMid.axpy(factorExp_,vuOld);
  
        JacobianRangeType du,duOld,duMid, diffusion{0.};
        uLocal.jacobian( quadrature[ pt ], du );
        uOldLocal.jacobian( quadrature[ pt ], duOld);
       
        //(1+theta)/2*DU^n+(1-theta)/2*DU^(n-1)
        // #if OPCHECK vuMid=vuOld
        duMid.axpy(factorImp_,du);
        duMid.axpy(factorExp_,duOld);
       
        RangeType avu(0.), source(0.);
        model_.systemSource(time_, xgl, source);        
//rho------------------------------------------------------------- 

#if OPCHECK 
        Filter::rho(avu)=Filter::rho(vu);
        Filter::rho(avu)-=Filter::rho(vuOld);
#else
        //d_t rho=(rho^n-rho^n-1)/delta t
        Filter::rho(avu)+=Filter::rho(vu);
        Filter::rho(avu)-=Filter::rho(vuOld);
        Filter::rho(avu)*=deltaInv;
#endif
        RangeFieldType div(0.),gradrhodotv(0.);
        
//---------------------------------------------------------------

//v--------------------------------------------------------------

        for( int ii = 0; ii < dimDomain ; ++ii )
          {
#if  OPCHECK 
           Filter::velocity(avu,ii)=Filter::velocity(vu,ii);
           Filter::velocity(avu,ii)-=Filter::velocity(vuOld,ii);
               
#else 
            Filter::velocity(avu,ii)=Filter::velocity(vu,ii);
            Filter::velocity(avu,ii)-=Filter::velocity(vuOld,ii);
            Filter::velocity(avu,ii)/=deltaT_;
#endif            
            RangeFieldType sgradv(0);
            
         
         }
          // A(dv) 
          model_.diffusion( duMid , diffusion );
          //------------------------------------------------------------------

//phi---------------------------------------------------------------
#if  OPCHECK
        Filter::phi( avu )=Filter::phi( vu );
        Filter::phi( avu )-=Filter::phi( vuOld ); 
#else
        Filter::phi( avu )=Filter::phi( vu )-Filter::phi( vuOld );
        Filter::phi( avu )*=deltaInv;
#warning "TIMEDERIVATIVE"      
#endif
        RangeFieldType transport(0.);
 
        RangeFieldType laplace(0.);

         // \nabla phi\cdot v
        for( int ii = 0; ii < dimDomain ; ++ii ) 
        { 
          transport+=Filter::velocity( vuMid , ii )*Filter::dphi(duMid, ii );
          //remove after debug
          laplace+=Filter::dsigma( duMid, ii , ii );
        }
#if HEATCHECK
      //Filter::phi( avu )-=laplace;
       Filter::phi( avu) -=Filter::tau( vuMid );
#else       
   //     Filter::phi( avu )+=transport+Filter::tau( vuMid )/Filter::rho( vuMid );
#endif     
//------------------------------------------------------------------        
       
//tau---------------------------------------------------------------
        // dF/dphi
        double dFdphi;
        model_.tauSource( Filter::phi(vuOld),
                          Filter::phi(vu),
                          Filter::rho(vuOld),
                          dFdphi);
#if OPCHECK         
       Filter::tau(avu)=Filter::tau(vu);
       Filter::tau(avu)-=Filter::tau(vuOld);
#else
        Filter::tau(avu)=Filter::tau(vu);
#endif
        // Filter::tau(avu)-=dFdphi;
        RangeFieldType divsigma(0.);

        for( int ii = 0 ; ii < dimDomain ; ++ii) 
          divsigma+=Filter::dsigma( duMid, ii , ii );
         
        Filter::tau( avu )-=model_.delta()*divsigma;
//-------------------------------------------------------------------

//mu-----------------------------------------------------------------

        Filter::mu(avu)=Filter::mu(vu)-Filter::mu(vuOld);
//------------------------------------------------------------------

//sigma--------------------------------------------------------------
        //\sigma-\nabla\phi
        for( int ii = 0 ; ii < dimDomain ; ++ii) 
          {
            //sigma^n
            Filter::sigma( avu , ii )=Filter::sigma( vu , ii );
      //      Filter::sigma( avu , ii )-=Filter::sigma( vuOld, ii);
#if OPCHECK
            //\nabla\phi^n-1
           Filter::sigma( avu , ii )-=Filter::dphi( duOld , ii );
#else
            //\nabla\phi^n
           Filter::sigma( avu , ii )-=Filter::dphi( du , ii );
#endif 
          }
          //------------------------------------------------------------------        
          for(int ii = 0; ii < dimRange ; ii++)
          {
            assert( avu[ii]==avu[ii]) ;
          }
          
          avu-=source;
          avu*=weight;
          diffusion*=weight;
          //add to result
          //wlocal+=avu*phi+diffusion*dphi
          wLocal.axpy( quadrature[ pt ], avu,diffusion );
      }
    }
    
    if ( !space().continuous() )
    {
      const double area = entity.geometry().volume();
      const IntersectionIteratorType iitend = space().gridPart().iend( entity ); 
      for( IntersectionIteratorType iit = space().gridPart().ibegin( entity ); iit != iitend; ++iit ) // looping over intersections
      {
        const IntersectionType &intersection = *iit;
        if ( intersection.neighbor() ) 
        {
         const EntityPointerType pOutside = intersection.outside(); // pointer to outside element.
         const EntityType &outside = *pOutside;
         LocalFunctionType uOutLocal = uGlobal.localFunction( outside ); // local u on outisde element
 
         computeIntersection(intersection,entity,outside,area,uLocal,uOutLocal,wLocal);
        
        }  
        else if( intersection.boundary() )
        {
          //   if ( ! model().isDirichletIntersection( intersection ) )
          //     continue;
          
          computeBoundary(intersection, entity, area, uLocal, wLocal);
        }
      }
    }
  } 


template<class DiscreteFunction, class Model, class Flux>
template<class LocalArgType, class NeighborArgType, class LFDestType>
void DGPhasefieldOperator<DiscreteFunction, Model,Flux>
::computeIntersection(  const IntersectionType& intersection,
                        const EntityType& entity,
                        const EntityType& neighbor,
                        const double area,
                        const LocalArgType& uEn, 
                        const NeighborArgType& uNb, 
                        LFDestType& wLocal) const
  {
    typedef typename IntersectionType::Geometry  IntersectionGeometryType;
    const IntersectionGeometryType &intersectionGeometry = intersection.geometry();
 
    LocalFunctionType uOldEn=uOld_.localFunction(entity); 
    LocalFunctionType uOldNb=uOld_.localFunction(neighbor); 
  
    TemporaryLocalType ufMidEn(space());
    TemporaryLocalType ufMidNb(space());

    ufMidEn.init( entity );
    ufMidEn.clear();
    ufMidNb.init( neighbor );
    ufMidNb.clear();

    ufMidEn.axpy( factorImp_, uEn );
    ufMidEn.axpy( factorExp_, uOldEn );

    ufMidNb.axpy( factorImp_, uNb );
    ufMidNb.axpy( factorExp_, uOldNb);

    // compute penalty factor
    const double intersectionArea = intersectionGeometry.volume();
    const double penaltyFactor = penalty()*intersectionArea / std::min( area, neighbor.geometry().volume() ); 
   
    const int quadOrderEn = uEn.order() + wLocal.order();
    const int quadOrderNb = uNb.order() + wLocal.order();
    
    FaceQuadratureType quadInside( space().gridPart(), intersection, quadOrderEn, FaceQuadratureType::INSIDE );
    FaceQuadratureType quadOutside( space().gridPart(), intersection, quadOrderNb, FaceQuadratureType::OUTSIDE );

    const size_t numQuadraturePoints = quadInside.nop();
   
    std::vector<RangeType> valuesEn( numQuadraturePoints );
    std::vector<RangeType> valuesNb( numQuadraturePoints );
    
    std::vector<RangeType> midValuesEn( numQuadraturePoints );
    std::vector<JacobianRangeType> midJacobiansEn( numQuadraturePoints );
    std::vector<RangeType> midValuesNb( numQuadraturePoints );
    std::vector<JacobianRangeType> midJacobiansNb( numQuadraturePoints );
 
#if OPCHECK
    uOldEn.evaluateQuadrature( quadInside, valuesEn );
    uOldNb.evaluateQuadrature( quadOutside,valuesNb );
#else
    uEn.evaluateQuadrature( quadInside, valuesEn );
    uNb.evaluateQuadrature( quadOutside, valuesNb );
#endif    
    ufMidEn.evaluateQuadrature( quadInside,midValuesEn );
    ufMidNb.evaluateQuadrature( quadOutside,midValuesNb );
    ufMidEn.evaluateQuadrature( quadInside,midJacobiansEn );
    ufMidNb.evaluateQuadrature( quadOutside,midJacobiansNb );


    for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
      {
        const typename FaceQuadratureType::LocalCoordinateType &x = quadInside.localPoint( pt );
        const DomainType normal = intersection.integrationOuterNormal( x );

        DomainType xgl=intersectionGeometry.global(x); 
        const double weight = quadInside.weight( pt );
             
        RangeType valueEn, valuNb;
        JacobianRangeType dvalue{0.},advalue{0.};
        double fluxRet;
        RangeType gLeft{0.},gRight{0.};
        fluxRet=flux_.numericalFlux(normal,
                                    penaltyFactor,
                                    valuesEn[ pt ],
                                    valuesNb[ pt ],
                                    midValuesEn[ pt ],  
                                    midValuesNb[ pt ], 
                                    gLeft,
                                    gRight); 
   
        RangeType value{0.};
#if 1        
        fluxRet+=flux_.diffusionFlux(normal,
                                    penaltyFactor,
                                    midValuesEn[ pt ],
                                    midValuesNb[ pt ],
                                    midJacobiansEn[ pt ],
                                    midJacobiansNb[ pt ],
                                    value,
                                    advalue);
#endif
      
       gLeft+=value;
       gLeft*=weight;
     
        advalue*=weight;
        
     //   std::cout<<"DifffluxValue="<<value<<" DifffluxAValue="<<advalue<<"\n";
        wLocal.axpy(quadInside[pt],gLeft,advalue);
      }//end quad loop
  }

template< class DiscreteFunction,class Model, class Flux>
template< class LocalArgType, class LFDestType>
void DGPhasefieldOperator<DiscreteFunction, Model, Flux>
::computeBoundary( const IntersectionType& intersection,
                    const EntityType& entity,
                    const double area,
                    const LocalArgType& uEn,
                    LFDestType& wLocal) const
  {
    abort();
    typedef typename IntersectionType::Geometry  IntersectionGeometryType;

    const IntersectionGeometryType &intersectionGeometry = intersection.geometry();
    LocalFunctionType uOldEn=uOld_.localFunction(entity); 
    
    TemporaryLocalType ufMidEn(space());

    ufMidEn.init( entity );
    ufMidEn.clear();

    ufMidEn.axpy( factorImp_, uEn );
    ufMidEn.axpy( factorExp_, uOldEn );



    DomainType xglobal{0};

    // compute penalty factor
    const double intersectionArea = intersectionGeometry.volume();
    const double penaltyFactor=penalty()*intersectionArea / area;
    const int quadOrder = uEn.order() + wLocal.order();

    FaceQuadratureType quadInside( space().gridPart(), intersection, quadOrder, FaceQuadratureType::INSIDE );
    const size_t numQuadraturePoints = quadInside.nop();
    double fluxRet=0.;
    std::vector<RangeType> valuesEn( numQuadraturePoints );
    std::vector<RangeType> midValuesEn( numQuadraturePoints );
    std::vector<JacobianRangeType> midJacobiansEn( numQuadraturePoints );
    
    uEn.evaluateQuadrature( quadInside, valuesEn );
    ufMidEn.evaluateQuadrature( quadInside,midValuesEn );
    ufMidEn.evaluateQuadrature( quadInside,midJacobiansEn );
 
    for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
      {
        const typename FaceQuadratureType::LocalCoordinateType &x = quadInside.localPoint( pt );
        xglobal = intersectionGeometry.global(x);
      

        const DomainType normal = intersection.integrationOuterNormal( x );
        const double weight = quadInside.weight( pt );
      
        RangeType value, valueBnd,valueBndOld;
        JacobianRangeType dvalue{0.},advalue{0.};

        RangeType vuIn,jump, vuOld, vuMid,gLeft;
        JacobianRangeType duIn, aduIn, duMidNb{0.};

        model_.dirichletValue(time_,xglobal,valueBnd); 
        model_.dirichletValue(time_-deltaT_, xglobal,valueBndOld);
      
        valueBnd*=factorImp_;
        valueBnd.axpy( factorExp_,valueBndOld );


        
        flux_.boundaryFlux(normal,midValuesEn[ pt ],gLeft);

#warning "DIFFUSIONFLUX MISSING AT BOUNDARY"

#if 1     
        fluxRet+=flux_.diffusionFlux( normal,
                                      penaltyFactor,
                                      midValuesEn[ pt ],
                                      valueBnd,
                                      midJacobiansEn[ pt ],
                                      duMidNb,//=0.
                                      value,
                                      advalue);
        gLeft+=value;
        gLeft*=weight;

        advalue*=weight;
#endif 
        wLocal.axpy(quadInside[pt],gLeft,advalue);
     }
  }










#endif //DUNE_PHASEFIELD_MIXEDOPERATOR.HH

#ifndef DUNE_PHASEFIELD_MIXEDOPERATOR_HH
#define DUNE_PHASEFIELD_MIXEDOPERATOR_HH

//globlas includes

//DUNE includes
#include <dune/common/fmatrix.hh>



//DUNE-FEM include
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/common/operator.hh>

#include <dune/fem/function/localfunction/temporary.hh>

#include <dune/fem/operator/common/differentiableoperator.hh>
#include <dune/fem/operator/common/automaticdifferenceoperator.hh>
#include <dune/fem/operator/common/stencil.hh>

#include "phasefieldfilter.hh"
#include  "flux.hh"
template<class DiscreteFunction, class Model, class Flux>
class DGPhasefieldOperator
: public virtual Dune::Fem::AutomaticDifferenceOperator<DiscreteFunction>
{
  typedef DiscreteFunction DiscreteFunctionType;
  typedef Model            ModelType;
  typedef Flux             NumericalFluxType;
protected:
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
                        const NumericalFluxType &flux,
                        double theta=0.5)
   : model_(model),
     space_(space),
     flux_(flux),
     theta_(theta),
     time_(0.),
     deltaT_(0.),
     uOld_("uOld" , space )
    {
      factorImp_=(1+theta_)*0.5;
      factorExp_=(1-theta_)*0.5;
    }

  // prepare the solution vector 
  template <class Function>
  void prepare( const Function &func, DiscreteFunctionType &u ) 
  { 
  }

  //! application operator 
  virtual void
  operator() ( const DiscreteFunctionType &u, DiscreteFunctionType &w ) const;

  
  void setTime(const double time)  { time_=time;}
  
  void setDeltaT( const double deltat) { deltaT_=deltat;}
  
  void setPreviousTimeStep( DiscreteFunctionType uOld)  { uOld_.assign(uOld);} 
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
  void computeBoundary( const IntersectionType& intersectintersection,
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
  double time_;
  double deltaT_;
  const double theta_;
  double factorImp_;
  double factorExp_;
  DiscreteFunctionType uOld_;
};

template<class DiscreteFunction, class Model, class Flux>
void DGPhasefieldOperator<DiscreteFunction, Model,Flux>
  ::operator() ( const DiscreteFunctionType &u, DiscreteFunctionType &w ) const 
{
  // clear destination 
  w.clear();

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
          const ArgType& u, 
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
        const double weight = quadrature.weight( pt ) * geometry.integrationElement( x );

        RangeType vu,vuOld, vuMid{0};
        uLocal.evaluate( quadrature[ pt ], vu );
        uOldLocal.evaluate( quadrature[ pt ], vuOld); 
       
        vuMid.axpy(factorImp_,vu);
        vuMid.axpy(factorExp_,vuOld);
      
        JacobianRangeType du,duOld,duMid, diffusion;
        uLocal.jacobian( quadrature[ pt ], du );
        uOldLocal.jacobian( quadrature[ pt ], duOld);
        
        duMid.axpy(factorExp_,du);
        duMid.axpy(factorExp_,duOld);
        RangeType avu(0.);
//rho------------------------------------------------------------- 
        Filter::rho(avu)+=Filter::rho(vu);
        Filter::rho(avu)-=Filter::rho(vuOld);
        Filter::rho(avu)/=deltaT_;
        RangeFieldType div(0.),gradrhodotv(0.);
        for(int i= 0;i<dimDomain;++i)
          { 
            div+=Filter::dvelocity(duMid,i,i);
            gradrhodotv+=Filter::drho(duMid,i)*Filter::velocity(vuMid,i);
          }
        
        Filter::rho(avu)+=Filter::rho(vuMid)*div+gradrhodotv;
//---------------------------------------------------------------

//v--------------------------------------------------------------
        for(int i= 0;i<dimDomain;++i)
          {
            Filter::velocity(avu,i)=Filter::velocity(vu,i)-Filter::velocity(vuOld,i);
            Filter::velocity(avu,i)/=deltaT_;
            Filter::velocity(avu,i)+=Filter::dmu(du,i);
   
            RangeFieldType sgradv(0);
            
            for(int j=0;j<dimDomain;++j)
              sgradv+=Filter::velocity(vuMid,j)*(Filter::dvelocity(duMid,j,i)-Filter::dvelocity(duMid,i,j  ));
            
            Filter::velocity(avu,i)+=sgradv;
            Filter::velocity(avu,i)*=Filter::rho(vuMid);
            Filter::velocity(avu,i)-=Filter::tau(vuMid)*Filter::dphi(duMid,i); 
          }
    
          model_.diffusion(duMid,diffusion);

//------------------------------------------------------------------

//phi---------------------------------------------------------------
        Filter::phi(avu)=Filter::phi(vu)-Filter::phi(vuOld);
        Filter::phi(avu)/=deltaT_;
        
        RangeFieldType transport(0.);
        
        for( int i=0; i<dimDomain;++i) 
        { 
          transport+=Filter::velocity(vuMid,i)*Filter::dphi(duMid,i);
        }
   
        Filter::phi(avu)+=transport-Filter::tau(vuMid)/Filter::rho(vuMid);
//------------------------------------------------------------------        
       
//tau---------------------------------------------------------------
        
        // dF/dphi
        double dFdphi;
        model_.tauSource(Filter::phi(vuOld),Filter::phi(vu),Filter::rho(vuOld),dFdphi);
        
        Filter::tau(avu)=Filter::tau(vuMid)-dFdphi;
        RangeFieldType divsigma(0.);

        for( int i=0; i<dimDomain;++i) 
         divsigma+=Filter::dsigma(duMid,i,i);
        
        Filter::tau(avu)+=model_.delta()*divsigma;
//-------------------------------------------------------------------

//mu-----------------------------------------------------------------
        //dF/drho
        double dFdrho;
        model_.muSource(Filter::rho(vuOld),Filter::rho(vu),Filter::phi(vu),dFdrho);

        Filter::mu(avu)+=Filter::mu(vuMid)-dFdrho;


        RangeFieldType   usqr(0.),uOldsqr(0.);
        for( int i=0; i<dimDomain;++i) 
        {
          usqr+=Filter::velocity(vu,i)*Filter::velocity(vu,i);
          uOldsqr+=Filter::velocity(vuOld,i)*Filter::velocity(vuOld,i);
        }
        
        Filter::mu(avu)-=0.25*(usqr-uOldsqr);
//------------------------------------------------------------------

//sigma--------------------------------------------------------------
         for( int i=0; i<dimDomain;++i) 
         Filter::sigma(avu,i)=Filter::sigma(vu,i)-Filter::dphi(du,i);
//------------------------------------------------------------------        

         
         
         //add to result
         //wlocal+=avu*phi+diffusion*dphi
         wLocal.axpy( quadrature[ pt ], avu,diffusion );
      }
    }
    if ( ! space().continuous() )
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
         LocalFunctionType uOutLocal = u.localFunction( outside ); // local u on outisde element

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
::computeIntersection( const IntersectionType& intersection,
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

    // compute penalty factor
    const double intersectionArea = intersectionGeometry.volume();
    const double penaltyFactor =  intersectionArea / std::min( area, neighbor.geometry().volume() ); 
   
    const int quadOrderEn = uEn.order() + wLocal.order();
    const int quadOrderNb = uNb.order() + wLocal.order();
    
    FaceQuadratureType quadInside( space().gridPart(), intersection, quadOrderEn, FaceQuadratureType::INSIDE );
    FaceQuadratureType quadOutside( space().gridPart(), intersection, quadOrderNb, FaceQuadratureType::OUTSIDE );

    const size_t numQuadraturePoints = quadInside.nop();

    for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
      {
        const typename FaceQuadratureType::LocalCoordinateType &x = quadInside.localPoint( pt );
        const DomainType normal = intersection.integrationOuterNormal( x );
        const double weight = quadInside.weight( pt );
            
        RangeType valueEn, valuNb;
        JacobianRangeType dvalue,advalue;

        RangeType vuEn,vuNb,vuEnOld,vuNbOld,gLeft,gRight;
        JacobianRangeType duEn,duEnOld,duNb,duNbOld; 

        uEn.evaluate( quadInside[ pt ], vuEn );
        uEn.jacobian( quadInside[ pt ], duEn );
        uOldEn.evaluate( quadInside[ pt ], vuEnOld );
        uOldNb.jacobian( quadInside[ pt ], duEnOld );
        
        uNb.evaluate( quadOutside[ pt ], vuNb );
        uNb.jacobian( quadOutside[ pt ], duNb );
        uOldNb.evaluate( quadOutside[ pt ], vuNbOld );
        uOldNb.jacobian( quadOutside[ pt ], duNbOld );

        RangeType uMidEn{0.},uMidNb{0.};
        JacobianRangeType duMidEn{0.},duMidNb{0.};
        
        uMidEn.axpy(factorImp_,vuEn);
        uMidEn.axpy(factorExp_,vuEnOld);
        uMidNb.axpy(factorImp_,vuNb);
        uMidNb.axpy(factorExp_,vuNbOld);

        duMidEn.axpy(factorImp_,duEn);
        duMidEn.axpy(factorExp_,duEnOld);

        duMidNb.axpy(factorImp_,duNb);
        duMidNb.axpy(factorExp_,duNbOld);


        double fluxRet;

        fluxRet=flux_.numericalFlux(normal,vuEn,vuNb,vuEnOld,vuNbOld,gLeft,gRight); 
    
        RangeType value{0.},dvalue{0.};
      
        
        
        flux_.diffusionFlux(normal,penaltyFactor,uMidEn,uMidNb,duMidEn,duMidNb,value,dvalue);
        
        
        // penalty term : beta [u] [phi] = beta (u+ - u-)(phi+ - phi-)=beta (u+ - u-)phi+ 
  
        gLeft+=jump*beta*intersectionGeometry.intersectionGeometry( x );    
       // {A grad u}.[phi] = {A grad u}.phi+ n_+ = 0.5*(grad u+ + grad u-).n_+ phi+
        //  [ u ] * { grad phi_en } = -normal(u+ - u-) * 0.5 grad phi_en

        duMidEn+=duMidNb; 
        duMidEn*=-0.5;
        model_.diffusion(duMidEn,diffusionU);
        diffusionU.umv(normal,gLeft);
         
        model_.diffusion(jumpNormal,diffusionW);
        diffusionW*=switchIP_;
        wLocal.axpy(quadOutside[pt],gLeft,diffusionW);
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
    typedef typename IntersectionType::Geometry  IntersectionGeometryType;

    const IntersectionGeometryType &intersectionGeometry = intersection.geometry();
    LocalFunctionType uOldEn=uOld_.localFunction(entity); 

    // compute penalty factor
    const double intersectionArea = intersectionGeometry.volume();
    const double beta = penalty() * intersectionArea / area;
    const int quadOrder = uEn.order() + wLocal.order();

    FaceQuadratureType quadInside( space().gridPart(), intersection, quadOrder, FaceQuadratureType::INSIDE );
    const size_t numQuadraturePoints = quadInside.nop();

    for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
      {
        const typename FaceQuadratureType::LocalCoordinateType &x = quadInside.localPoint( pt );
        const DomainType normal = intersection.integrationOuterNormal( x );
        const double weight = quadInside.weight( pt );
      
        RangeType value;
        JacobianRangeType dvalue,advalue;

        RangeType vuIn,jump, vuOld, vuMid,gLeft;
        JacobianRangeType duIn, aduIn;
        uEn.evaluate( quadInside[ pt ], vuIn );
        uEn.jacobian( quadInside[ pt ], duIn );
        
        uOldEn.evaluate( quadInside[ pt ], vuOld);
        
        flux_.boundaryFlux(normal,vuIn, vuOld,gLeft);
      }
  }










#endif //DUNE_PHASEFIELD_MIXEDOPERATOR.HH

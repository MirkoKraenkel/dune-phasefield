#ifndef DUNE_PHASEFIELD_MIXEDOPERATOR_HH
#define DUNE_PHASEFIELD_MIXEDOPERATOR_HH

//globlas includes

//DUNE includes
#include <dune/common/fmatrix.hh>



//DUNE-FEM includes

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/common/operator.hh>

#include <dune/fem/function/localfunction/temporary.hh>

#include <dune/fem/operator/common/differentiableoperator.hh>
#include <dune/fem/operator/common/stencil.hh>

#include "phasefieldfilter.hh"

template<class DiscreteFunction, class Model, class Flux, class Params>
class DGPhasefieldOperator
: public virtual Dune::Fem::Operator<DiscreteFunction>
{
  typedef DiscreteFunction DiscreteFunctionType;
  typedef Model            ModelType;
  typedef Flux             NumericalFluxType;
  typedef Params           ParameterClassType;
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

  typedef  PhasefieldFilter<RangeType> U; 
 public:
   //! constructor
   DGPhasefieldOperator(const ModelType &model,
                        const DiscreteFunctionSpaceType &space,
                        const NumericalFluxType &flux,
                        const ParameterClassType &params)
   : model_(model),
     space_( space),
     flux_(flux),
     parameters_(params),
     uOld_("uOld" , space )
  {}
  // prepare the solution vector 
  template <class Function>
  void prepare( const Function &func, DiscreteFunctionType &u ) 
  { 
  }

  //! application operator 
  virtual void
  operator() ( const DiscreteFunctionType &u, DiscreteFunctionType &w ) const;

  
 protected:

  void setTime(double &time){ time_=time;}
  
  void setDeltaT( double &deltat){ deltat_=deltat;}
  
  void setPreviousStep( DiscreteFunctionType &uOld) { uOld_=uOld;} 
  
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
  
  protected:
  ModelType model_;
  const DiscreteFunctionSpaceType space_;
  const NumericalFluxType flux_;
  const ParameterClassType parameters_;
  double time_;
  double deltat_;
  DiscreteFunctionType& uOld_;


};

template<class DiscreteFunction, class Model, class Flux, class Params>
void DGPhasefieldOperator<DiscreteFunction, Model,Flux,Params>
  ::operator() ( const DiscreteFunctionType &u, DiscreteFunctionType &w ) const 
{
  // clear destination 
  w.clear();

  // get discrete function space 
 // const DiscreteFunctionSpaceType &dfSpace = w.space();

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

template<class DiscreteFunction, class Model, class Flux, class Params>
template< class ArgType, class LocalArgType,class LFDestType >
void DGPhasefieldOperator<DiscreteFunction, Model,Flux,Params>
::localOp(const EntityType& entity,
          const ArgType& u, 
          const LocalArgType& uLocal,
          LFDestType& wLocal) const
  {
//    const LocalFunctionType uLocal = u.localFunction( entity );
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

        RangeType vu,vuOld, vuMid;
        uLocal.evaluate( quadrature[ pt ], vu );
        uOldLocal.evaluate( quadrature[ pt ], vuOld); 
       
        vuMid=0.5*(vu+vuOld);

        JacobianRangeType du,duOld,duMid;
        uLocal.jacobian( quadrature[ pt ], du );
        uOldLocal.jacobian( quadrature[ pt ], duOld);
        
        duMid=0.5*(du+duOld);
        
        RangeType avu(0.);
//rho------------------------------------------------------------- 
        U::density(avu)+=U::density(vu);
        U::density(avu)-=U::density(vuOld);
        U::density(avu)/=deltat_;
        RangeFieldType div(0.),gradrhodotv(0.);
        for(int i= 0;i<dimDomain;++i)
          { 
            div+=U::dvelocity(duMid,i,i);
            gradrhodotv+=U::drho(duMid,i)*U::velocity(vuMid,i);
          }
        
        U::density(avu)+=U::density(vuMid)*div+gradrhodotv;
//---------------------------------------------------------------

//v--------------------------------------------------------------
        for(int i= 0;i<dimDomain;++i)
          {
            U::velocity(avu,i)=U::velocity(vu,i)-U::velocity(vuOld,i);
            U::velocity(avu,i)/=deltat_;
            U::velocity(avu,i)-=U::dmu(du,i);
   
            RangeFieldType sgradv(0);
            
            for(int j=0;j<dimDomain;++j)
              sgradv+=U::velocity(vuMid,j)*(U::dvelocity(duMid,i,j)-U::dvelcoity(duMid,j,i));
            
            U::velocity(avu,i)+=sgradv;
            U::velocity(avu,i)*=U::rho(vuMid);
            U::velocity(avu,i)-=U::tau(vuMid)*U::dphi(vuMid,i); 
          }
    

//------------------------------------------------------------------

//phi---------------------------------------------------------------
        U::phi(avu)=U::phi(vu)-U::phi(vuOld);
        U::phi(avu)/=deltat_;
        
        RangeFieldType transport(0.);
        
        for( int i=0; i<dimDomain;++i) 
        { 
          transport(0.)+=U::velocity(vuMid,i)*U::dphi(vuMid,i);
        }
        U::phi(avu)+=transport-U::tau(vuMid)/U::rho(vuMid);
//------------------------------------------------------------------        
       
//tau---------------------------------------------------------------
        U::tau(avu)=U::tau(vuMid)-(model_.Psi(U::rho(vuOld),U::phi(vu))
                                  -model_.Psi(U::rho(vuOld),U::phi(vuOld)))/(U::phi(vu)-U::phi(vuOld));
        RangeFieldType divsigma(0.);

        for( int i=0; i<dimDomain;++i) 
         divsigma+=U::dsigma(vuMid,i,i);
        
        U::tau(avu)+=model_.delta()*divsigma;
//-------------------------------------------------------------------

//mu-----------------------------------------------------------------
        U::mu(avu)+=U::mu(vuMid)-(model_.Psi(U::rho(vu),U::phi(vu))
                                  -model_.Psi(U::rho(vuOld),U::phi(vu)))/(U::rho(vu)-U::rho(vuOld));

        RangeFieldType   usqr(0.),uOldsqr(0.);
        for( int i=0; i<dimDomain;++i) 
        {
          usqr+=U::velocity(vu,i)*U::velocity(vu,i);
          uOldsqr+=U::velocity(vuOld,i)*U::velocity(vuOld,i);
        }
        
        U::mu(avu)-=0.25*(usqr-uOldsqr);
//------------------------------------------------------------------

//sigma--------------------------------------------------------------
         for( int i=0; i<dimDomain;++i) 
         U::sigma(avu,i)=U::sigma(vu,i)-U::dphi(du,i);
//------------------------------------------------------------------        
        wLocal.axpy( quadrature[ pt ], avu );
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
          if ( ! model().isDirichletIntersection( intersection ) )
                continue;
          computeBoundary(intersection, entity, area, uLocal, wLocal);
        }
      }
    }
  } 






#endif //DUNE_PHASEFIELD_MIXEDOPERATOR.HH

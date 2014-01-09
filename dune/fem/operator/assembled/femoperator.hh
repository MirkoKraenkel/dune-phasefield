#ifndef DUNE_PHASEFIELD_MIXEDOPERATOR_HH
#define DUNE_PHASEFIELD_MIXEDOPERATOR_HH

//globlas includes

//DUNE includes
#include <dune/common/fmatrix.hh>

#define OPCHECK 1
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
#include  "flux.hh"
template<class DiscreteFunction, class Model, class Flux>
class FemPhasefieldOperator
: public virtual Dune::Fem::Operator<DiscreteFunction,DiscreteFunction>
{
  typedef Dune::Fem::Operator<DiscreteFunction,DiscreteFunction> BaseType;
 protected:
 
  typedef DiscreteFunction DiscreteFunctionType;
  typedef Model            ModelType;
  typedef Flux             NumericalFluxType;
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType;
  #if 0
  typedef typename BaseType::RangeFunctionType RangeFunctionType;
  typedef typename BaseType::DomainFunctionType DomainFunctionType;
  typedef typename BaseType::RangeFieldType RangeFieldType;
  typedef typename BaseType::DomainFieldType DomainFieldType;
  typedef typename BaseType::DomainSpaceType DomainSpaceType;
  typedef typename BaseType::RangeSpaceType RangeSpaceType;
  typedef DomainSpaceType DiscreteFunctionSpaceType;
#endif
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
   FemPhasefieldOperator(const ModelType &model,
                        const DiscreteFunctionSpaceType &space,
                        const NumericalFluxType &flux)
   : 
     model_(model),
     space_(space),
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

  // prepare the solution ve#endifctor 
  template <class Function>
  void prepare( const Function &func, DiscreteFunctionType &u ) 
  { 
  }

  //! application operator 
  void operator() ( const DiscreteFunctionType &u, DiscreteFunctionType &w ) const;

  
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
  const double  theta_;
  double time_;
  double deltaT_;
  double factorImp_;
  double factorExp_;
  DiscreteFunctionType uOld_;
};





template<class DiscreteFunction, class Model, class Flux, class Jacobian>
class FDJacobianFemPhasefieldOperator
: public FemPhasefieldOperator<DiscreteFunction,Model,Flux>,
  public virtual Dune::Fem::AutomaticDifferenceOperator<DiscreteFunction,DiscreteFunction, Jacobian>
  {
  
  typedef FemPhasefieldOperator<DiscreteFunction,Model,Flux> MyOperatorType;
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
   FDJacobianFemPhasefieldOperator(const ModelType &model,
                                  const DiscreteFunctionSpaceType &space)
   :MyOperatorType(model,space){} 
  
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
void FemPhasefieldOperator<DiscreteFunction, Model,Flux>
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
void FemPhasefieldOperator<DiscreteFunction, Model,Flux>
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
        const double weight = quadrature.weight( pt )* geometry.integrationElement( x );

        RangeType vu,vuOld, vuMid{0};
        uLocal.evaluate( quadrature[ pt ], vu );
        uOldLocal.evaluate( quadrature[ pt ], vuOld); 
      
        //(1+theta)/2*U^n+(1-theta)/2*U^(n-1)
        // #if OPCHECK vuMid=vuOld;      
        vuMid.axpy(factorImp_,vu);
        vuMid.axpy(factorExp_,vuOld);
  
        JacobianRangeType du,duOld,duMid, diffusion;
        uLocal.jacobian( quadrature[ pt ], du );
        uOldLocal.jacobian( quadrature[ pt ], duOld);
       
        //(1+theta)/2*DU^n+(1-theta)/2*DU^(n-1)
        // #if OPCHECK vuMid=vuOld
        duMid.axpy(factorImp_,du);
        duMid.axpy(factorExp_,duOld);
       
        
        RangeType avu(0.);
//rho------------------------------------------------------------- 

#if OPCHECK 
        Filter::rho(avu)=Filter::rho(vu);
#else
        //d_t rho=(rho^n-rho^n-1)/delta t
        Filter::rho(avu)+=Filter::rho(vu);
        Filter::rho(avu)-=Filter::rho(vuOld);
        Filter::rho(avu)/=deltaT_;
#endif
        RangeFieldType div(0.),gradrhodotv(0.);
        
        //div(rho v)=rho*div v+gradrho v
        for(int i= 0;i<dimDomain;++i)
          { 
            //sum d_i v_i 
            div+=Filter::dvelocity(duMid,i,i);
            // sum d_i rho*v_i
            gradrhodotv+=Filter::drho(duMid,i)*Filter::velocity(vuMid,i);
          }
         
        Filter::rho(avu)+=Filter::rho(vuMid)*div+gradrhodotv;
         
//---------------------------------------------------------------

//v--------------------------------------------------------------

        for(int i= 0;i<dimDomain;++i)
          {
#if OPCHECK
            Filter::velocity(avu,i)=Filter::velocity(vu,i);
#else 
            Filter::velocity(avu,i)=Filter::velocity(vu,i);
            Filter::velocity(avu,i)-=Filter::velocity(vuOld,i);
            Filter::velocity(avu,i)/=deltaT_;
#endif            
            RangeFieldType sgradv(0);
            
            //sum_j v_j( d_j v_i - d_i v_j)
            for(int j=0;j<dimDomain;++j)
              sgradv+=Filter::velocity(vuMid,j)*(Filter::dvelocity(duMid,i,j)-Filter::dvelocity(duMid,j,i));
          
            Filter::velocity(avu,i)+=sgradv;
            Filter::velocity(avu,i)+=Filter::dmu(duMid,i);
            //rho*(d_t v+S(dv)v+dmu) 
            Filter::velocity(avu,i)*=Filter::rho(vuMid);
            // -tau\nabla phi 
            // check me: dphiMid,dphiOld or dphi???
            Filter::velocity(avu,i)-=Filter::tau(vuMid)*Filter::dphi(duMid,i); 
          }
          // A(dv) 
          model_.diffusion(duMid,diffusion);
//------------------------------------------------------------------

//phi---------------------------------------------------------------
#if OPCHECK
        Filter::phi(avu)=Filter::phi(vu);
#else
        Filter::phi(avu)=Filter::phi(vu)-Filter::phi(vuOld);
        Filter::phi(avu)/=deltaT_;
#endif
        RangeFieldType transport(0.);
       

        // \nabla phi\cdot v
        for( int i=0; i<dimDomain;++i) 
        { 
          transport+=Filter::velocity(vuMid,i)*Filter::dphi(duMid,i);
        }
        Filter::phi(avu)+=transport+Filter::tau(vuMid)/Filter::rho(vuMid);
     
//------------------------------------------------------------------        
       
//tau---------------------------------------------------------------
        // dF/dphi
        double dFdphi;
        model_.tauSource(Filter::phi(vuOld),Filter::phi(vu),Filter::rho(vuOld),dFdphi);
#if OPCHECK         
        Filter::tau(avu)=Filter::tau(vu);
#else
        Filter::tau(avu)=Filter::tau(vuMid);
#endif
        Filter::tau(avu)-=dFdphi;
        Filter::tau(avu)=Filter::tau(vu);
        RangeFieldType divsigma(0.);

        for( int i=0; i<dimDomain;++i) 
         divsigma+=Filter::dsigma(duOld,i,i);
         
        Filter::tau(avu)+=model_.delta()*divsigma;
//-------------------------------------------------------------------

//mu-----------------------------------------------------------------
        //dF/drho
        double dFdrho;
        model_.muSource(Filter::rho(vuOld),Filter::rho(vu),Filter::phi(vu),dFdrho);

        //mu-d_rho F
#if OPCHECK
        Filter::mu(avu)=Filter::mu(vu);
#else
        Filter::mu(avu)=Filter::mu(vuMid);
#endif     
        Filter::mu(avu)-=dFdrho;
        RangeFieldType   usqr(0.),uOldsqr(0.);
        for( int i=0; i<dimDomain;++i) 
        {
#if OPCHECK
          // |v^n|^2
          usqr+=Filter::velocity(vuOld,i)*Filter::velocity(vuOld,i);
#else
          // |v^n|^2
          usqr+=Filter::velocity(vu,i)*Filter::velocity(vu,i);
#endif
          // |v^{n-1}|^2
          uOldsqr+=Filter::velocity(vuOld,i)*Filter::velocity(vuOld,i);
        }
        // -\frac{1}{4}( |v^n|^2-|v^{n-1}|^2)
       // Filter::mu(avu)-=0.25*(usqr+uOldsqr);
//------------------------------------------------------------------

      for(int i=0;i<dimRange;i++)
          {
            assert( avu[i]==avu[i]) ;
          }
         
          avu*=weight;
          diffusion*=weight;

          //add to result
          //wlocal+=avu*phi+diffusion*dphi
          wLocal.axpy( quadrature[ pt ], avu,diffusion );
      }
    }
    
    if ( ! space().continuous() )
    {
      std::cout<<"DONT USE WITH DGSPACE!\n";
      abort();
    }
  } 












#endif //DUNE_PHASEFIELD_MIXEDOPERATOR.HH

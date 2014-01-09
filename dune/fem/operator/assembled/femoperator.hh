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
  typedef typename Dune::FieldVector<double, dimDomain> VelocityRangeType;
  


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
       
        
        RangeType avu{0.};
        JacobianRangeType advu{0.};
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
            //  rho^n*v_i^n
            Filter::drho(avu,i)=-Filter::rho(vu)*Filter::velocity(vu,i);
            // regularization Term alpha*(\nabla\rho, \nabla\psi);
            Filter::drho(avu,i)+=alpha_*Filter::drho(vu,i);
          }
//---------------------------------------------------------------

//v--------------------------------------------------------------

        for(int i= 0;i<dimDomain;++i)
          {
#if OPCHECK
            Filter::velocity(avu,i)=Filter::velocity(vu,i);
#else 
            //1/2*rho^(n-1)dt v^n+dt (rho v) = 1/delta*(v^n rho^(n-1)/2-3/2*rho^(n-1)v^(n-1)+rho^n v^n)
            Filter::velocity(avu,i)=Filter::velocity(vu,i);
            Filter::velocity(avu,i)-=3*Filter::velocity(vuOld,i);
            Filter::velocity(avu,i)*=0.5*Filter::rho(vuOld);
            Filter::velocity(avu,i)+=Filter::velocity(vu,i)*Filter::rho(vu); 
            Filter::velocity(avu,i)/=deltaT_;
#endif            
            VelocityRangeType tansportv{0.}; 
            //sum_j v_j( d_j v_i - d_i v_j)
            for( int j=0; j<dimDomain ; ++j )
              {
                //gamma(\nabla dt v^n,\ksi)
                Filter::dvelocity(avu,i,j)=Filter::dvelocity(vu,i,j);
                Filter::dvelocity(avu,i,j)-=Filter::dvelocity(vuOld,i,j);
                Filter::dvelocity(avu,i,j)*=gamma_;
                Filter::dvelocity(avu,i,j)/=deltaT_;
                //(rho^(n-1) (v^(n-1)\cdot\nabla)\ksi,v^n)=(v^n\otimes\v^(n-1),\nabla v^n)
                Filter::dvelocity(advu,i,j)-=Filter::velocity(vu,i)*Filter::velocity(vu,j)*Filter::rho(vuOld);
                //rho^(n-1)*(v\cot\nabla) v
                transportv=Filter::dvelocity(vu,i,j)*Filter::velocity(vuOld,j)*Filter::rho(vuOld);
              }
            
            //add rho^(n-1)*(v\cot\nabla) v
            Filter::velocity(avu,i)+=transportv;
           
            //rho^n\nabla\mu^n 
            Filter::velocity(avu,i)+=Filter::dmu(duMid,i)*Filter::rho(vu);
            
            // -tau^n \nabla phi^n 
            Filter::velocity(avu,i)-=Filter::tau(vuMid)*Filter::dphi(duMid,i); 
          }
          // A(dv) 
          model_.diffusion(duMid,diffusion);
          advu+=diffusion;
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
              transport+=Filter::velocity(vu,i)*Filter::dphi(du,i);
            }
        
          Filter::phi(avu)+=transport+Filter::tau(vu)/Filter::rho(vu);
     
//------------------------------------------------------------------        
       
//tau---------------------------------------------------------------
        // dF/dphi
        double dFdphi;
        model_.tauSource(Filter::phi(vuOld),Filter::phi(vu),Filter::rho(vuOld),dFdphi);
#if OPCHECK         
        Filter::tau(avu)=Filter::tau(vu);
#else
        Filter::tau(avu)=Filter::tau(vu);
#endif
        Filter::tau(avu)-=dFdphi;
        Filter::tau(avu)=Filter::tau(vu);
        RangeFieldType divsigma(0.);

        for( int i=0; i<dimDomain;++i) 
          Filter::tau(advu,i)-=model_.delta()*Filter::dphi(du,i);
//-------------------------------------------------------------------

//mu-----------------------------------------------------------------
        //dF/drho
        double dFdrho;
        model_.muSource(Filter::rho(vuOld),Filter::rho(vu),Filter::phi(vu),dFdrho);

        //mu-d_rho F
        Filter::mu(avu)=Filter::mu(vu);
        Filter::mu(avu)-=dFdrho;
//------------------------------------------------------------------

      for(int i=0;i<dimRange;i++)
          {
            assert( avu[i]==avu[i]) ;
          }
         
          avu*=weight;
          diffusion*=weight;

          //add to result
          //wlocal+=avu*phi+diffusion*dphi
          wLocal.axpy( quadrature[ pt ], avu, advu);
      }
    }
    
    if ( ! space().continuous() )
    {
      std::cout<<"DONT USE WITH DGSPACE!\n";
      abort();
    }
  } 












#endif //DUNE_PHASEFIELD_MIXEDOPERATOR.HH

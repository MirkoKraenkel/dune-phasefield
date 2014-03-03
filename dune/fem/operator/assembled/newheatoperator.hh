#ifndef DUNE_PHASEFIELD_MIXEDOPERATOR_HH
#define DUNE_PHASEFIELD_MIXEDOPERATOR_HH


#warning  "NEWHEATOP"
//globlas includes

//DUNE includes
#include <dune/common/fmatrix.hh>
#define HEATCHECK 1
#warning "DEBUGGING VERSION"

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
 static   const int order=POLORDER;
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
  typedef Dune::Fem::LagrangePointSet<GridPartType,order> LagrangePointSetType;
   
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
     uOld_("uOld" , space ),
     uOldLocal_(space),
     uOldNeighbor_(space)
  {
      assert(theta_>=0 && theta_<=1);
      factorImp_=0.5*(1+theta_);
      factorExp_=0.5*(1-theta_);
      std::cout<<"factorExp "<<factorExp_<<"\n";
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
  void setEntity( const EntityType& entity ) const
  {
    uOldLocal_.init(entity);
    uOldLocal_.assign(uOld_.localFunction( entity ));
    areaEn_=entity.geometry().volume();
  }
  
  void setNeighbor( const EntityType& entity ) const
  {
    uOldNeighbor_.init(entity);
    uOldNeighbor_.assign(uOld_.localFunction( entity ));
    areaNb_=entity.geometry().volume();
  }
  
  
  void localIntegral( size_t  pt,
                 const GeometryType& geometry,
                 const QuadratureType& quadrature,
                  RangeType& vu,
                 JacobianRangeType& du,
                 RangeType& avu, // to be added to the result local function
                 JacobianRangeType& advu) const;
  

  void intersectionIntegral( const IntersectionType& intersection,
                             const size_t pt,  
                              const FaceQuadratureType& quadInside,
                              const FaceQuadratureType& quadOutside,
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
  mutable TemporaryLocalType uOldLocal_; 
  mutable TemporaryLocalType uOldNeighbor_; 
  mutable double areaEn_;
  mutable double areaNb_;
  mutable double penalty_;
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

  using MyOperatorType::localIntegral;
  using MyOperatorType::computeBoundary;
  using MyOperatorType::setEntity;
  using MyOperatorType::setNeighbor;
  using MyOperatorType::intersectionIntegral;
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
  using MyOperatorType::uOldLocal_; 
  using MyOperatorType::uOldNeighbor_; 
  using MyOperatorType::areaEn_;
  using MyOperatorType::areaNb_;

  
  };

//#include "newheatoperator.cc"

#if VISIT
#include "fulloperator2.cc"
#else
#include "fulloperator.cc"
#endif









#endif //DUNE_PHASEFIELD_MIXEDOPERATOR.HH

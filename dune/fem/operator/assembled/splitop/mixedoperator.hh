#ifndef DUNE_PHASEFIELD_MIXEDOPERATOR_HH
#define DUNE_PHASEFIELD_MIXEDOPERATOR_HH


//globlas includes

//DUNE includes
#include <dune/common/fmatrix.hh>


//Dune::Fem includes
//#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/quadrature/intersectionquadrature.hh>

#include <dune/fem/operator/common/operator.hh>

#include <dune/fem/function/localfunction/temporary.hh>

#include <dune/fem/operator/common/differentiableoperator.hh>
//#include <dune/fem/operator/common/automaticdifferenceoperator.hh>
#include "../common/automaticdifferenceoperator.hh"
#include <dune/fem/operator/common/stencil.hh>

#include <dune/fem/gridpart/common/capabilities.hh>

//local includes
#include "phasefieldfilter.hh"


template<class DiscreteFunction, class Model, class Flux  >
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

  typedef Dune::Fem::CachingQuadrature< GridPartType, 1 > FaceQuadratureType;

  typedef Dune::Fem::CachingQuadrature< GridPartType, 0 > QuadratureType;

  typedef Dune::Fem::TemporaryLocalFunction<DiscreteFunctionSpaceType> TemporaryLocalType;

  static const int dimDomain = LocalFunctionType::dimDomain;
  static const int dimRange = LocalFunctionType::dimRange;


  typedef Dune::Fem::LagrangePointSet<GridPartType,order> LagrangePointSetType;

  typedef  PhasefieldFilter<RangeType> Filter; 

  public:
  //! constructor
  DGPhasefieldOperator( const ModelType &model,
                        const DiscreteFunctionSpaceType &space):
                        model_(model),
                        space_(space),
                        flux_(model),
                        theta_(Dune::Fem::Parameter::getValue<double>("phasefield.mixed.theta")),
                        time_(0.),
                        deltaT_(0.),
                        maxSpeed_(0.),
                        lastSpeed_(1.),
                        uOld_("uOld" , space ),
                        uOldLocal_(space),
                        uOldNeighbor_(space),
                        outflow_(Dune::Fem::Parameter::getValue<bool>("phasefield.outflow")),
                        minArea_( std::numeric_limits<double>::max() )
    {
      uOld_.clear();
      factorImp_=0.5*(1+theta_);
      factorExp_=0.5*(1-theta_);
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

  double getDeltaT() {return deltaT_;}

  double timeStepEstimate() { return std::min( minArea_/maxSpeed_,lipschitzC());}

  double maxSpeed() { return maxSpeed_; }
  
  double lipschitzC() { return model_.lipschitzC(); }

  void setPreviousTimeStep( DiscreteFunctionType& uOld)  { uOld_.assign(uOld);} 
  
  DiscreteFunctionType& getPreviousTimeStep() { return uOld_;}


  protected:
  void setEntity( const EntityType& entity ) const
  {
    uOldLocal_.init(entity);
    uOldLocal_.assign(uOld_.localFunction( entity ));
    addDataLocal_.init( entity );
    addDataLocal_.assign( addData_.localFunction( entity ));
    areaEn_=entity.geometry().volume();
    minArea_=std::min( areaEn_ , minArea_ );
  }

  void setNeighbor( const EntityType& entity ) const
  {
    uOldNeighbor_.init(entity);
    uOldNeighbor_.assign(uOld_.localFunction( entity ));
    addDataNeighbor_.init( entity );
    addDataNeighbor_.assign( addData_.localFunction( entity ));
    areaNb_=entity.geometry().volume();
    minArea_=std::min( areaNb_ , minArea_ ); 
  }


  void localIntegral( size_t  pt,
                      const GeometryType& geometry,
                      const QuadratureType& quadrature,
                      RangeType& vu,
                      JacobianRangeType& du,
                      RangeType& avu, // to be added to the result local function
                      JacobianRangeType& advu) const;

  
  template<bool conforming> 
  void computeIntersection( const IntersectionType &intersection,
                            const LocalFunctionType& uLocal,
                            const LocalFunctionType& uNeighbor,
                            LocalFunctionType& wLocal) const;



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

  template<class LocalArgType, class LFDestType>
  void computeBoundary( const IntersectionType& intersection,
                        const EntityType& entity,
                        const double area,
                        const LocalArgType& uEn, 
                        LFDestType& wlocal) const;

  const DiscreteFunctionSpaceType& space() const {return space_;}

  const ModelType& model() const{ return model_;}

  
  protected:
  ModelType model_;
  const DiscreteFunctionSpaceType &space_;
  const NumericalFluxType flux_;
  const double  theta_;
  double time_;
  double deltaT_;
  mutable double maxSpeed_;
  mutable double lastSpeed_;
  double factorImp_;
  double factorExp_;
  mutable DiscreteFunctionType uOld_;
  mutable DiscreteFunctionType addData_;
  mutable temporarylocaltype uoldlocal_; 
  mutable temporarylocaltype uoldneighbor_; 
  mutable temporarylocaltype addDatalocal_; 
  mutable temporarylocaltype addDataneighbor_; 
  const bool outflow_;
  mutable double areaEn_;
  mutable double areaNb_;
  mutable double minArea_;
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
  FDJacobianDGPhasefieldOperator( const ModelType &model,
                                  const DiscreteFunctionSpaceType &space)
    :MyOperatorType(model,space){} 

  using MyOperatorType::prepare;
  //! application operator 
  using MyOperatorType::operator();
  using MyOperatorType::setTime;
  using MyOperatorType::setDeltaT;
  using MyOperatorType::setPreviousTimeStep;
  using MyOperatorType::getPreviousTimeStep;
  using MyOperatorType::timeStepEstimate;
  using MyOperatorType::maxSpeed;
  using MyOperatorType::lipschitzC;
  protected:

  using MyOperatorType::localIntegral;
  using MyOperatorType::computeBoundary;
  using MyOperatorType::setEntity;
  using MyOperatorType::setNeighbor;
  using MyOperatorType::intersectionIntegral;
  using MyOperatorType::space;
  using MyOperatorType::model;

  protected:
  using MyOperatorType::model_;
  using MyOperatorType::space_;
  using MyOperatorType::flux_;
  using MyOperatorType::time_;
  using MyOperatorType::deltaT_;
  using MyOperatorType::maxSpeed_;
  using MyOperatorType::lastSpeed_;
  using MyOperatorType::theta_;
  using MyOperatorType::factorImp_;
  using MyOperatorType::factorExp_;
  using MyOperatorType::uOld_;
  using MyOperatorType::uOldLocal_; 
  using MyOperatorType::uOldNeighbor_; 
  using MyOperatorType::outflow_;
  using MyOperatorType::debugmu_;
  using MyOperatorType::debugtheta_;
  using MyOperatorType::areaEn_;
  using MyOperatorType::areaNb_;


};

template<class DiscreteFunction, class Model, class Flux>
void DGPhasefieldOperator<DiscreteFunction, Model,Flux>
  ::operator() ( const DiscreteFunctionType &u, DiscreteFunctionType &w ) const 
{
  double lastSpeed=maxSpeed_;
  maxSpeed_=0;

  // clear destination 
  w.clear();
  assert(deltaT_>0);
  // iterate over grid 
  const IteratorType end = space().end();
  for( IteratorType it = space().begin(); it != end; ++it )
  {
   
    //bool boundaryElement=false;
    // get entity (here element) 
    const EntityType &entity = *it;
    // get elements geometry
    const GeometryType& geometry=entity.geometry();
    // get local representation of the discrete functions 
    const LocalFunctionType uLocal = u.localFunction( entity);

    setEntity( entity );
    RangeType vu(0.),avu(0.);
    JacobianRangeType du{0.},adu{0.};
    
    LocalFunctionType wLocal = w.localFunction( entity );
    const int quadOrder =2*space().order(entity);
    
    QuadratureType quadrature( entity, quadOrder );
    const size_t numQuadraturePoints = quadrature.nop();
    for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
      {
        uLocal.evaluate( quadrature[ pt ], vu);
        uLocal.jacobian( quadrature[ pt ], du);
        
        localIntegral( pt , geometry, quadrature , vu , du , avu , adu );   
       
        //wlocal+=avu*phi+diffusion*dphi
        wLocal.axpy( quadrature[ pt ], avu, adu);
      }   
    
      const IntersectionIteratorType iitend = space().gridPart().iend( entity ); 
      for( IntersectionIteratorType iit = space().gridPart().ibegin( entity ); iit != iitend; ++iit ) // looping over intersections
      {
        const IntersectionType &intersection = *iit;

        if ( intersection.neighbor() ) 
        {
          const EntityPointerType pOutside = intersection.outside(); // pointer to outside element.
          const EntityType &neighbor = *pOutside;

          setNeighbor(neighbor);

          LocalFunctionType uNeighbor=u.localFunction(neighbor);

          if( !Dune::Fem::GridPartCapabilities::isConforming< GridPartType >::v && !intersection.conforming())
            {
              computeIntersection<false>( intersection,
                                          uLocal,
                                          uNeighbor,
                                          wLocal);
            }
          else
            {
              computeIntersection<true>(intersection,
                                        uLocal,
                                        uNeighbor,
                                        wLocal);
         
            }

          }
        else if (  intersection.boundary())
          {
            //boundaryElement=true;
            const int quadOrderEn = uLocal.order() + wLocal.order();

            FaceQuadratureType quadInside( space().gridPart(), intersection, quadOrderEn, FaceQuadratureType::INSIDE );

            const size_t numQuadraturePoints = quadInside.nop();


            for( size_t pt=0; pt < numQuadraturePoints; ++pt )
              {
                RangeType vuEn(0.),avuLeft(0.);
                JacobianRangeType duEn{0.},aduLeft{0.};
                uLocal.evaluate( quadInside[ pt ], vuEn);
                uLocal.jacobian( quadInside[ pt ], duEn);
                //const typename QuadratureType::CoordinateType &x = quadrature.point( pt );
                const double weight = quadInside.weight( pt );


                boundaryIntegral( intersection,
                                  pt,
                                  quadInside,
                                  vuEn,
                                  duEn,
                                  avuLeft,
                                  aduLeft);

                avuLeft*=weight;
                aduLeft*=weight;

                wLocal.axpy( quadInside[ pt ] , avuLeft , aduLeft );
          }
        }

      }
#if 0     
      if( boundaryElement )
      { 
        const int order=1; 
        const LagrangePointSetType lagrangePointSet( geometry.type(), order );

        const IntersectionIteratorType iitend = space().gridPart().iend( entity ); 
        for( IntersectionIteratorType iit = space().gridPart().ibegin( entity ); iit != iitend; ++iit ) // looping over intersections
        {
          const IntersectionType &intersection = *iit;

          if ( intersection.boundary())
          {
            // get face number of boundary intersection 
            const int face = intersection.indexInInside();


            typedef typename LagrangePointSetType::template Codim< 1 >:: SubEntityIteratorType
              FaceDofIteratorType;
            // get dof iterators 
            FaceDofIteratorType faceIt = lagrangePointSet.template beginSubEntity< 1 >( face );
            const FaceDofIteratorType faceEndIt = lagrangePointSet.template endSubEntity< 1 >( face );
            for( ; faceIt != faceEndIt; ++faceIt )
            {
              const int localBlock=*faceIt;
              ntersectionQuadrature< FaceQuadratureType, conforming > IntersectionQuadratureType; 
                  const int localBlockSize=DiscreteFunctionSpaceTiyplocalBlockSize;


              const int dofOffset=localBlock*dimRange;
              wLocal[dofOffset]=0;
              for( int ii=0 ; ii < dimDomain ; ++ii)
              {
                wLocal[dofOffset+1+ii]=0;
              }

            }

          }
        }
      } 
#endif
    


  }
  // communicate data (in parallel runs)
  w.communicate();

}

template<class DiscreteFunction, class Model, class Flux >
template<bool conforming>
void DGPhasefieldOperator<DiscreteFunction, Model,Flux>
::computeIntersection( const IntersectionType &intersection,
                       const LocalFunctionType& uLocal,
                       const LocalFunctionType& uNeighbor,
                       LocalFunctionType& wLocal) const
{
  const int quadOrderEn = 2*std::max(uLocal.order(), wLocal.order())+2;
  typedef Dune::Fem::IntersectionQuadrature< FaceQuadratureType, conforming > IntersectionQuadratureType; 
  typedef typename IntersectionQuadratureType::FaceQuadratureType QuadratureImp;

  IntersectionQuadratureType interQuad(space().gridPart(), intersection,quadOrderEn);
  const QuadratureImp& quadInside=interQuad.inside();
  const QuadratureImp& quadOutside=interQuad.outside();
 


  const size_t numQuadraturePoints = quadInside.nop();

  for( size_t pt=0; pt < numQuadraturePoints; ++pt )
  {
    RangeType vuEn(0.),vuNb(0.),avuLeft(0.),avuRight(0.);
    JacobianRangeType duEn{0.},duNb{0.},aduLeft{0.}, aduRight{0.};
    uLocal.evaluate( quadInside[ pt ], vuEn);
    uLocal.jacobian( quadInside[ pt ], duEn);
    uNeighbor.evaluate( quadOutside[ pt ], vuNb);
    uNeighbor.jacobian( quadOutside[ pt ], duNb);
    const double weight = quadInside.weight( pt );


    //calculate quadrature summands avu
    intersectionIntegral( intersection,
                          pt,
                          quadInside,
                          quadOutside,
                          vuEn,
                          vuNb,
                          duEn,
                          duNb,
                          avuLeft,
                          avuRight,
                          aduLeft,
                          aduRight);

    avuLeft*=weight;
    aduLeft*=weight;

    wLocal.axpy( quadInside[ pt ] , avuLeft , aduLeft );
  } 
}

#define DIFFQUOT 0 

#if NSK
  #include "fulloperatorNSK.cc"
#else
  #if IMPLICITTAU
    #include "fulloperatorimpl.cc"
  #else
    #include "fulloperator.cc"
#endif
#endif







#endif //DUNE_PHASEFIELD_MIXEDOPERATOR.HH

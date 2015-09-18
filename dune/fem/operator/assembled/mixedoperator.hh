#ifndef DUNE_PHASEFIELD_DGOPERATOR_HH
#define DUNE_PHASEFIELD_DGOPERATOR_HH


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


template<class DiscreteFunction, class Integrator>
class DGOperator
: public virtual Dune::Fem::Operator<DiscreteFunction,DiscreteFunction>
{
  typedef Dune::Fem::Operator<DiscreteFunction,DiscreteFunction> BaseType;
  static   const int order=POLORDER;
  public:

  typedef DiscreteFunction DiscreteFunctionType;
  typedef Integrator IntegratorType;
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

  public:
  //! constructor
  DGOperator( IntegratorType &integrator,
              const DiscreteFunctionSpaceType& space):
              integrator_( integrator),
              space_(space)
    {
    }

  // prepare the solution vector 
  template <class Function>
  void prepare( const Function &func, DiscreteFunctionType &u )
  {
  }

  //! application operator 
  void operator() ( const DiscreteFunctionType &u, DiscreteFunctionType &w ) const;


  void setTime (const double time)  {integrator_.setTime(time);}

  void setDeltaT( const double deltat) { integrator_.setDeltaT(deltat);}

  double timeStepEstimate() { return integrator_.timeStepEstimate();}

  double lastSpeed () { return integrator_.lastSpeed();}
  void setEntity ( const EntityType& entity ) const
  {
    integrator_.setEntity( entity );
  }

  void setNeighbor ( const EntityType& entity ) const
  {
    integrator_.setNeighbor( entity);
  }


  
  template<bool conforming> 
  void computeIntersection( const IntersectionType &intersection,
                            const LocalFunctionType& uLocal,
                            const LocalFunctionType& uNeighbor,
                            LocalFunctionType& wLocal) const;



  template<class LocalArgType, class LFDestType>
  void computeBoundary( const IntersectionType& intersection,
                        const EntityType& entity,
                        const double area,
                        const LocalArgType& uEn, 
                        LFDestType& wlocal) const;

  const DiscreteFunctionSpaceType& space() const {return space_;}

  IntegratorType& integrator() const{ return integrator_;}

   
  protected:
  IntegratorType& integrator_;
  const DiscreteFunctionSpaceType& space_;
};






template<class DiscreteFunction, class Integrator>
void DGOperator<DiscreteFunction, Integrator>
  ::operator() ( const DiscreteFunctionType &u, DiscreteFunctionType &w ) const 
{
  // clear destination 
  w.clear();
  integrator_.setSpeed(); 
  // iterate over grid 
  const IteratorType end = space().end();
  //for( IteratorType it = space().begin(); it != end; ++it )
  for( const auto& entity: space())
  {
   
    // get elements geometry
    const GeometryType& geometry=entity.geometry();
    // get local representation of the discrete functions 
    const LocalFunctionType uLocal = u.localFunction( entity);

    setEntity( entity );
    RangeType vu(0.),avu(0.);
    JacobianRangeType du(0.),adu(0.);
    
    LocalFunctionType wLocal = w.localFunction( entity );
    const int quadOrder =2*space().order(entity);
    
    QuadratureType quadrature( entity, quadOrder );
    const size_t numQuadraturePoints = quadrature.nop();
    for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
      {
        uLocal.evaluate( quadrature[ pt ], vu);
        uLocal.jacobian( quadrature[ pt ], du);
        
        integrator_.localIntegral( pt , geometry, quadrature , vu , du , avu , adu );   
       
        //wlocal+=avu*phi+diffusion*dphi
        wLocal.axpy( quadrature[ pt ], avu, adu);
      }   
    
      const IntersectionIteratorType iitend = space().gridPart().iend( entity ); 
      for( IntersectionIteratorType iit = space().gridPart().ibegin( entity ); iit != iitend; ++iit ) // looping over intersections
      {
        const IntersectionType &intersection = *iit;

        if ( intersection.neighbor() ) 
        {
          const EntityType &neighbor = intersection.outside();

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
            const int quadOrderEn = 2*uLocal.order() + 1;

            FaceQuadratureType quadInside( space().gridPart(), intersection, quadOrderEn, FaceQuadratureType::INSIDE );

            const size_t numQuadraturePoints = quadInside.nop();


            for( size_t pt=0; pt < numQuadraturePoints; ++pt )
              {
                RangeType vuEn(0.),avuLeft(0.);
                JacobianRangeType duEn(0.),aduLeft(0.);
                uLocal.evaluate( quadInside[ pt ], vuEn);
                uLocal.jacobian( quadInside[ pt ], duEn);
                //const typename QuadratureType::CoordinateType &x = quadrature.point( pt );
                const double weight = quadInside.weight( pt );


                integrator_.boundaryIntegral( intersection,
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
    


  }
  // communicate data (in parallel runs)
  w.communicate();

}

template< class DiscreteFunction, class Integrator >
template<bool conforming>
void DGOperator<DiscreteFunction, Integrator>
::computeIntersection( const IntersectionType &intersection,
                       const LocalFunctionType& uLocal,
                       const LocalFunctionType& uNeighbor,
                       LocalFunctionType& wLocal) const
{
  const int quadOrderEn = 2*std::max(uLocal.order(), wLocal.order())+1;
  typedef Dune::Fem::IntersectionQuadrature< FaceQuadratureType, conforming > IntersectionQuadratureType; 
  typedef typename IntersectionQuadratureType::FaceQuadratureType QuadratureImp;

  IntersectionQuadratureType interQuad(space().gridPart(), intersection,quadOrderEn);
  const QuadratureImp& quadInside=interQuad.inside();
  const QuadratureImp& quadOutside=interQuad.outside();
 


  const size_t numQuadraturePoints = quadInside.nop();

  for( size_t pt=0; pt < numQuadraturePoints; ++pt )
  {
    RangeType vuEn(0.),vuNb(0.),avuLeft(0.),avuRight(0.);
    JacobianRangeType duEn(0.),duNb(0.),aduLeft(0.), aduRight(0.);
    uLocal.evaluate( quadInside[ pt ], vuEn);
    uLocal.jacobian( quadInside[ pt ], duEn);
    uNeighbor.evaluate( quadOutside[ pt ], vuNb);
    uNeighbor.jacobian( quadOutside[ pt ], duNb);
    const double weight = quadInside.weight( pt );


    //calculate quadrature summands avu
    integrator_.intersectionIntegral( intersection,
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









#endif //DUNE_PHASEFIELD_MIXEDOPERATOR.HH

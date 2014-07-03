#ifndef PHASEFIELD_BOUNDARYCORRECTION_HH
#define PHASEFIELD_BOUNDARYCORRECTION_HH

//globlas includes



#include <dune/fem/space/lagrange/lagrangepoints.hh>

template<class DiscreteFunction, class Model>
class PhasefieldBoundaryCorrection
{

  static const int order=POLORDER;
protected:

  typedef DiscreteFunction DiscreteFunctionType;
  typedef Model            ModelType;
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




  static const int dimDomain = LocalFunctionType::dimDomain;
  static const int dimRange = LocalFunctionType::dimRange;

  typedef Dune::Fem::LagrangePointSet<GridPartType,order> LagrangePointSetType;
 
public:
   //! constructor
   PhasefieldBoundaryCorrection(const ModelType &model,
                        const DiscreteFunctionSpaceType &space)
   : 
     model_(model),
     space_(space)
    {
     }


  //! application operator 
  void operator() (  DiscreteFunctionType &u ) const;

  
 
  
 
 protected:

  const DiscreteFunctionSpaceType& space() const {return space_;}
  
  const ModelType& model() const{ return model_;}

  
protected:
  ModelType model_;
  const DiscreteFunctionSpaceType &space_;
};







template<class DiscreteFunction, class Model>
void PhasefieldBoundaryCorrection<DiscreteFunction, Model>
  ::operator() ( DiscreteFunctionType &u ) const 
{

  // iterate over grid 
  const IteratorType end = space().end();
  for( IteratorType it = space().begin(); it != end; ++it )
  {
    // get entity (here element) 
    const EntityType &entity = *it;
    const GeometryType& geometry=entity.geometry();
    
    const int order=POLORDER; 
    const LagrangePointSetType lagrangePointSet( geometry.type(), order );
  
    // get local representation of the discrete functions 
    LocalFunctionType uLocal = u.localFunction( entity);
    
    const IntersectionIteratorType iitend = space().gridPart().iend( entity ); 
    for( IntersectionIteratorType iit = space().gridPart().ibegin( entity ); iit != iitend; ++iit ) // looping over intersections
    {
      const  IntersectionType &intersection=*iit;   
     
      if( intersection.boundary())
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
                const int localBlockSize=DiscreteFunctionSpaceType::localBlockSize;
                
                const int dofOffset=localBlock*dimRange;

                for( int ii=0 ; ii < dimDomain ; ++ii)
                {
                  uLocal[dofOffset+1+ii]=0;
                }
            
              }


      
      }

    } 
    
  }

}











#endif //DUNE_PHASEFIELD_MIXEDOPERATOR.HH

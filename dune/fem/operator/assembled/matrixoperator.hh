#ifndef PHASEFIELD_MATRIXOPERATOR_HH
#define PHASEFIELD_MATRIXOPERATOR_HH      
// Dune::Fem includes
#include <dune/fem/function/localfunction/temporary.hh>
#include <dune/fem/operator/common/differentiableoperator.hh>
#include <dune/fem/operator/common/stencil.hh>
#include <dune/fem/operator/common/temporarylocalmatrix.hh>

// local includes
#include "../matrixhelper.hh"
template<class Operator,class Tensor, class Jacobian>
class MatrixOperator:
public Dune::Fem::DifferentiableOperator < Jacobian >
{

  using MyOperatorType=Operator;
  using TensorType=Tensor;
  using BaseType=Dune::Fem::DifferentiableOperator< Jacobian>;

  enum{ dimDomain=MyOperatorType::dimDomain};
  
  enum{ dimRange=MyOperatorType::RangeType::dimension};
  
  typedef typename BaseType::JacobianOperatorType JacobianOperatorType;
  typedef typename MyOperatorType::DiscreteFunctionType DiscreteFunctionType;
  typedef typename MyOperatorType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionSpaceType::BasisFunctionSetType BasisFunctionSetType;
  typedef typename MyOperatorType::DomainType DomainType;
  typedef typename MyOperatorType::RangeType RangeType;
  typedef typename MyOperatorType::RangeFieldType RangeFieldType;
  typedef typename MyOperatorType::JacobianRangeType JacobianRangeType;
  typedef typename MyOperatorType::IteratorType IteratorType;
  typedef typename MyOperatorType::IntersectionIteratorType IntersectionIteratorType;
  typedef typename MyOperatorType::IntersectionType IntersectionType;
  typedef typename MyOperatorType::EntityType EntityType;
  typedef typename MyOperatorType::EntityPointerType EntityPointerType;
  typedef typename MyOperatorType::GeometryType GeometryType;
  typedef typename MyOperatorType::LocalFunctionType LocalFunctionType;
  typedef typename MyOperatorType::QuadratureType QuadratureType;
  typedef typename MyOperatorType::FaceQuadratureType FaceQuadratureType;
  using IntegratorType=typename MyOperatorType::IntegratorType;
  using ModelType=typename IntegratorType::ModelType;
 
  typedef Dune::Fem::TemporaryLocalFunction<DiscreteFunctionSpaceType> TemporaryLocalFunctionType;
  typedef typename MyOperatorType::GridPartType GridPartType;
  typedef typename GridPartType::IndexSetType IndexSetType;
  typedef Dune::Fem::TemporaryLocalMatrix< DiscreteFunctionSpaceType, DiscreteFunctionSpaceType> LocalMatrixType;
  typedef typename JacobianOperatorType::LocalMatrixType RealLocalMatrixType;
  typedef Dune::Fem::MutableArray<RangeType> RangeVectorType;
  typedef Dune::Fem::MutableArray<JacobianRangeType> JacobianVectorType;
  
  public: 
  MatrixOperator(MyOperatorType& myop,
                 const DiscreteFunctionSpaceType &space):
    myOperator_(myop),
    assembler_(myOperator_.integrator()),
    indexSet_(space.gridPart().indexSet()),
    visited_(0),
    stencil_(space,space)
    {}



  typedef Dune::Fem::DiagonalAndNeighborStencil<DiscreteFunctionSpaceType,DiscreteFunctionSpaceType> StencilType;
 
 //! application operator 
  void operator() ( const DiscreteFunctionType &u, DiscreteFunctionType &w ) const
    {
      myOperator_( u , w);
    }


  void setTime(const double time)  { myOperator_.setTime(time); }
  
  double getTime(){return myOperator_.getTime(); }

  void setDeltaT( const double deltat) { myOperator_.setDeltaT(deltat); }

  double getDeltaT() {return myOperator_.getDeltaT(); }

  double timeStepEstimate()  { return myOperator_.timeStepEstimate(); }

  double maxSpeed() { return myOperator_.integrator().maxSpeed(); }
  
  double lipschitzC() { return myOperator_.integrator().lipschitzC(); }

  void setPreviousTimeStep( DiscreteFunctionType& uOld)  { myOperator_.integrator().setPreviousTimeStep( uOld ); }

  IntegratorType& integrator() const { return myOperator_.integrator();}
  //DiscreteFunctionType& getPreviousTimeStep() { return myOperator_.getPreviousTimeStep(); }
  
  const DiscreteFunctionSpaceType& space() const { return myOperator_.space();}
  
  const ModelType& model() const { return myOperator_.integrator().model(); }

 
  void jacobian(const DiscreteFunctionType &u, JacobianOperatorType &jOp) const;
   
  
 protected: 
  void computeElement ( const EntityType& entity,
                        const GeometryType& geometry,
                        const BasisFunctionSetType &baseSet,
                        const LocalFunctionType& uLocal,
                        LocalMatrixType& jLocal) const;
  
  
  template< bool conforming >
  void computeIntersection ( const IntersectionType &intersection,
                             const GeometryType& geometry,
                             const GeometryType& geometryNb,
                             const BasisFunctionSetType &baseSet,
                             const BasisFunctionSetType &baseSetNb,
                             const LocalFunctionType& uLocal,
                             const LocalFunctionType& uLocalNb,
                             LocalMatrixType&  jLocal,
                             LocalMatrixType&  jLocalNbEn,
                             LocalMatrixType&  jLocalEnNb,
                             LocalMatrixType&  jLocalNbNb)  const;

  void computeBoundary ( const IntersectionType &intersection,
                         const GeometryType& geometry,
                         const BasisFunctionSetType &baseSet,
                         const LocalFunctionType& uLocal,
                         LocalMatrixType&  jLocal)  const;



 
  private:
  MyOperatorType& myOperator_;
  TensorType& assembler_;
  const IndexSetType& indexSet_;
  mutable Dune::Fem::MutableArray<bool> visited_;
  StencilType stencil_;
};


// Implementation of LocalFDOperator
// // ------------------------------------------------


template< class Operator ,class Tensor, class Jacobian > 
template< bool comforming >
void MatrixOperator< Operator , Tensor, Jacobian >
::computeIntersection ( const IntersectionType &intersection,
                        const GeometryType& geometry,
                        const GeometryType& geometryNb,
                        const BasisFunctionSetType &baseSet,
                        const BasisFunctionSetType &baseSetNb,
                        const LocalFunctionType& uLocal,
                        const LocalFunctionType& uLocalNb,
                        LocalMatrixType&  jLocal,
                        LocalMatrixType&  jLocalNbEn,
                        LocalMatrixType&  jLocalEnNb,
                        LocalMatrixType&  jLocalNbNb)  const
{  
  //construct quadrature 
  const int quadOrder = 2*std::max(uLocal.order(),uLocalNb.order())+1;
  typedef Dune::Fem::IntersectionQuadrature< FaceQuadratureType, conforming > IntersectionQuadratureType; 
  typedef typename IntersectionQuadratureType::FaceQuadratureType QuadratureImp;

  IntersectionQuadratureType interQuad( space().gridPart(), 
                                        intersection,
                                        quadOrder);
 
  const QuadratureImp& quadInside=interQuad.inside();
  const QuadratureImp& quadOutside=interQuad.outside();
  

  const size_t numQuadraturePoints = quadInside.nop();

  //prepare assembler
  assembler_.computeCoefficientsEn( uLocal, quadInside );
  assembler_.computeCoefficientsNb( uLocalNb, quadOutside );
 
  //quad Loop
  for( size_t pt=0 ; pt < numQuadraturePoints ; ++pt )
    {
    
      assembler_.intersectionTensor( pt,
                                     intersection,
                                     geometry,
                                     geometryNb,
                                     quadInside, 
                                     quadOutside,
                                     baseSet,
                                     baseSetNb,
                                     jLocal,
                                     jLocalNbEn,
                                     jLocalEnNb,
                                     jLocalNbNb);
                                     
                                     
                                     
                                     
                                     

      }
}
template< class Operator ,class Tensor, class Jacobian > 
void MatrixOperator< Operator , Tensor, Jacobian >
::computeBoundary ( const IntersectionType &intersection,
                    const GeometryType& geometry,
                    const BasisFunctionSetType &baseSet,
                    const LocalFunctionType& uLocal,
                    LocalMatrixType&  jLocal)const
{ 
  const int quadOrder = 2*uLocal.order()+1;

  FaceQuadratureType quadInside( space().gridPart(), 
                                 intersection,
                                 quadOrder,
                                 FaceQuadratureType::INSIDE);
  

  const size_t numQuadraturePoints = quadInside.nop();

  //prepare assembler
  assembler_.computeCoefficientsEn( uLocal, quadInside );
 
  //quad Loop
  for( size_t pt=0 ; pt < numQuadraturePoints ; ++pt )
    {
    
      assembler_.boundaryTensor( pt,
                                 intersection,
                                 geometry,
                                 quadInside, 
                                 baseSet,
                                 jLocal);

    }
}
// compute the temporary local matrix of the bilinearform for the (en,en) block
template< class Operator, class Tensor, class Jacobian> void
MatrixOperator< Operator , Tensor ,Jacobian >
::computeElement( const EntityType& entity,
                  const GeometryType& geometry,
                  const BasisFunctionSetType &baseSet,
                  const LocalFunctionType& uLocal,
                  LocalMatrixType& jLocal) const
{    

  QuadratureType quadrature( entity, 2*space().order(entity)+1 );
  
  size_t numQuadraturePoints=quadrature.nop();
  assembler_.computeCoefficientsEn( uLocal,
                                    quadrature);
  RangeType vuOld(0.),vuMid(0);
  
  for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
    {
      assembler_.elementTensor( pt,
                                geometry,
                                quadrature,
                                baseSet,
                                jLocal );

     }
    
}



template<class Operator, class Tensor,class Jacobian> void
MatrixOperator< Operator , Tensor, Jacobian >
::jacobian ( const DiscreteFunctionType &u, JacobianOperatorType &jOp ) const
{
  Dune::Fem::DiagonalAndNeighborStencil<DiscreteFunctionSpaceType,DiscreteFunctionSpaceType> stencil(space(),space());
  jOp.reserve(stencil);
  jOp.clear();
  
  //intialize visited marker
  visited_.resize( indexSet_.size(0));
  const size_t indSize = visited_.size();
  for( size_t ii = 0; ii < indSize; ++ii) visited_[ii] = false;


  const DiscreteFunctionSpaceType &dfSpace = u.space();
  const GridPartType& gridPart = dfSpace.gridPart();

  const unsigned int numDofs = dfSpace.blockMapper().maxNumDofs() * 
    DiscreteFunctionSpaceType :: localBlockSize ;
  
  assembler_.resizeBaseStorageEn( numDofs );
  assembler_.resizeBaseStorageNb( numDofs ); 
  
  const IteratorType end = dfSpace.end();
  for( IteratorType it = dfSpace.begin(); it != end; ++it )
  {
    const EntityType &entity = *it;
    const GeometryType geometry = entity.geometry();
    const LocalFunctionType uLocal = u.localFunction( entity );

    //initialize uOld
    integrator().setEntity( entity );
    LocalMatrixType jLocal( dfSpace, dfSpace);//= jOp.localMatrix( entity, entity );
    jLocal.init(entity,entity);
    jLocal.clear();
    RealLocalMatrixType realLocal=jOp.localMatrix( entity, entity);

    const BasisFunctionSetType &baseSet = jLocal.domainBasisFunctionSet();
    
    computeElement( entity,
                    geometry,
                    baseSet,
                    uLocal,
                    jLocal );
 
    MatrixHelper::add(jLocal, realLocal);
    
    jLocal.clear();
    if (space().continuous() )
      continue;
    
    const IntersectionIteratorType endiit = gridPart.iend( entity );
    for( IntersectionIteratorType iit = gridPart.ibegin( entity ); iit != endiit ; ++ iit )
      {
        const IntersectionType& intersection = *iit ;

        if( intersection.neighbor() )
          {
            EntityPointerType ep = intersection.outside();
            const EntityType& neighbor = *ep ;
            const GeometryType geometryNb = neighbor.geometry();

            if( !visited_[ indexSet_.index( neighbor ) ])
              {
                //evaluate additional data
                integrator().setNeighbor( neighbor );
                //initialize all needed  tempory local matrices
                jLocal.clear();
                LocalMatrixType jLocalNbEn( dfSpace,dfSpace);
                jLocalNbEn.init( neighbor, entity);
                
                jLocalNbEn.clear();
                LocalMatrixType jLocalEnNb( dfSpace, dfSpace);
                jLocalEnNb.init( entity, neighbor);
                jLocalEnNb.clear();
                LocalMatrixType jLocalNbNb( dfSpace,dfSpace);
                jLocalNbNb.init( neighbor,neighbor);
                jLocalNbNb.clear();
                
                // get local matrix for face entries and neighbor
                RealLocalMatrixType realLocalNbEn = jOp.localMatrix( neighbor,entity );
                RealLocalMatrixType realLocalEnNb = jOp.localMatrix( entity, neighbor );
                RealLocalMatrixType realLocalNbNb = jOp.localMatrix( neighbor,neighbor); 


                const LocalFunctionType uLocalNb = u.localFunction(neighbor);
                
                // get neighbor's base function set 
                const BasisFunctionSetType &baseSetNb = jLocalNbEn.domainBasisFunctionSet();
                 
                if( !Dune::Fem::GridPartCapabilities::isConforming< GridPartType >::v
                    && !intersection.conforming())
                  {
                    
                    computeIntersection<false>( intersection,
                                                geometry,
                                                geometryNb,
                                                baseSet,
                                                baseSetNb,
                                                uLocal,
                                                uLocalNb,
                                                jLocal,
                                                jLocalNbEn,
                                                jLocalEnNb,
                                                jLocalNbNb);
                  }
                  else
                  {
                    computeIntersection<true>( intersection,
                                               geometry,
                                               geometryNb,
                                               baseSet,
                                               baseSetNb,
                                               uLocal,
                                               uLocalNb,
                                               jLocal,
                                               jLocalNbEn,
                                               jLocalEnNb,
                                               jLocalNbNb);

                  }
                  
                MatrixHelper::add( jLocal, realLocal);
                jLocal.clear();
                MatrixHelper::add( jLocalNbEn, realLocalNbEn);
                jLocalNbEn.clear();
                MatrixHelper::add( jLocalEnNb, realLocalEnNb);
                jLocalEnNb.clear();
                MatrixHelper::add( jLocalNbNb, realLocalNbNb);
                jLocalNbNb.clear();
          } 
      
        }
        else if ( intersection.boundary() )
        {
          
          computeBoundary( intersection,
                           geometry,
                           baseSet,
                           uLocal,
                           jLocal); 
          
          MatrixHelper::add( jLocal, realLocal);
          jLocal.clear();
        }
      }
    visited_[indexSet_.index( entity )]=true;
  } // end grid traversal 
  jOp.communicate();
}


#endif //PHASEFIELD_MATRIXOPERATOR_HH


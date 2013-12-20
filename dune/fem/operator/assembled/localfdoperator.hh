#ifndef LOCALFD_OPERATOR_HH
#define LOCALFD_OPERATOR_HH



#include <dune/fem/function/localfunction/temporary.hh>

#include <dune/fem/operator/common/differentiableoperator.hh>
#include <dune/fem/operator/common/stencil.hh>

#include "mixedoperator.hh"



template<class DiscreteFunction,class Model, class Flux, class Jacobian>
class LocalFDOperator
 :public Dune::Fem::DifferentiableOperator<Jacobian>,
  protected DGPhasefieldOperator<DiscreteFunction,Model,Flux>
{
 
  typedef DGPhasefieldOperator<DiscreteFunction,Model,Flux> MyOperatorType;
  
  typedef Dune::Fem::DifferentiableOperator<Jacobian> BaseType;
 

  typedef typename BaseType::JacobianOperatorType JacobianOperatorType;

  typedef typename MyOperatorType::DiscreteFunctionType DiscreteFunctionType;
  typedef typename MyOperatorType::ModelType ModelType;
  typedef typename MyOperatorType::NumericalFluxType NumericalFluxType;
  typedef typename MyOperatorType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename MyOperatorType::BasisFunctionSetType BasisFunctionSetType;
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
  
  typedef typename MyOperatorType::TemporaryLocalFunctionType TemporaryLocalFunctionType;
  typedef typename MyOperatorType::GridPartType GridPartType;
  public: 
  LocalFDOperator(const ModelType &model,
                  const DiscreteFunctionSpaceType &space,
                  const NumericalFluxType &flux)
  :MyOperatorType(model,space,flux),
    stencil_(space,space),
    epsilon_(Dune::Fem::Parameter::getValue<double>("phasefield.fdjaconian.epsilon"))
  {}

  using MyOperatorType::localOp;
  using MyOperatorType::computeIntersection;
  using MyOperatorType::operator();
  using MyOperatorType::space;
  void jacobian(const DiscreteFunctionType &u, JacobianOperatorType &jOp) const;
  
  private:
  Dune::Fem::DiagonalAndNeighborStencil<DiscreteFunctionSpaceType,DiscreteFunctionSpaceType> stencil_;
  double epsilon_;
};


// Implementation of LocalFDOperator
// // ------------------------------------------------

template<class DiscreteFunction,class Model, class Flux, class Jacobian> void
LocalFDOperator<DiscreteFunction, Model, Flux,  Jacobian>
  ::jacobian ( const DiscreteFunctionType &u, JacobianOperatorType &jOp ) const
{
  typedef typename JacobianOperatorType::LocalMatrixType LocalMatrixType;
  typedef typename DiscreteFunctionSpaceType::BasisFunctionSetType BasisFunctionSetType;

  jOp.reserve(stencil_);
  jOp.clear();
  TemporaryLocalFunctionType fu( space() );
  TemporaryLocalFunctionType fu_j( space() );
  TemporaryLocalFunctionType u_j( space() );
  

  const DiscreteFunctionSpaceType &dfSpace = u.space();
  const GridPartType& gridPart = dfSpace.gridPart();

  const unsigned int numDofs = dfSpace.blockMapper().maxNumDofs() * 
                               DiscreteFunctionSpaceType :: localBlockSize ;

  std::vector< RangeType> phi( numDofs ); 
  std::vector< JacobianRangeType > dphi( numDofs );

  std::vector< RangeType > phiNb( numDofs );
  std::vector< JacobianRangeType > dphiNb( numDofs );

  const IteratorType end = dfSpace.end();
  for( IteratorType it = dfSpace.begin(); it != end; ++it )
  {
    const EntityType &entity = *it;
    const GeometryType geometry = entity.geometry();

    const LocalFunctionType uLocal = u.localFunction( entity );
  
    fu.init(entity);
    fu_j.init(entity);
    u_j.init(entity);
    fu.clear();
    fu_j.clear();
    u_j.clear();
   
    LocalMatrixType jLocal = jOp.localMatrix( entity, entity );

    const BasisFunctionSetType &baseSet = jLocal.domainBasisFunctionSet();
    const unsigned int numBasisFunctions = baseSet.size();
    std::vector<TemporaryLocalFunctionType> ujLocal(numBasisFunctions,TemporaryLocalFunctionType( dfSpace));
    std::vector<TemporaryLocalFunctionType> fujLocal(numBasisFunctions,TemporaryLocalFunctionType( dfSpace));

    for( size_t j=0; j<numBasisFunctions; ++j)
    {
      ujLocal[j].init(entity);
      fujLocal[j].init(entity);
      ujLocal[j].clear();
      fujLocal[j].clear();
      ujLocal[j].assign(uLocal);
      ujLocal[j][j]+=epsilon_;
    }
    
    localOp(entity,u,uLocal,fu);
    RangeFieldType eps =epsilon_;//sqrt((1+sqrt(normU))*epsilon_);
     double epsinv=1./eps;
 
    for( unsigned int localCol = 0; localCol < numBasisFunctions; ++localCol )
    {
      localOp(entity,u,ujLocal[localCol],fujLocal[localCol]);
      for( unsigned int k = 0; k <numBasisFunctions; ++k)
        {
          double diff=fujLocal[localCol][k]-fu[k];
          diff*=epsinv;
          jLocal.add(k,localCol,diff);
          
        }
    }
 
    QuadratureType quadrature( entity, 2*dfSpace.order() );

   if ( dfSpace.continuous() )
      continue;


    double area = geometry.volume();
    const IntersectionIteratorType endiit = gridPart.iend( entity );
    for ( IntersectionIteratorType iit = gridPart.ibegin( entity );
          iit != endiit ; ++ iit )
    {
      const IntersectionType& intersection = *iit ;

      if( intersection.neighbor() )
      {
        EntityPointerType ep = intersection.outside();
        const EntityType& neighbor = *ep ;
        typedef typename IntersectionType::Geometry  IntersectionGeometryType;
        const IntersectionGeometryType &intersectionGeometry = intersection.geometry();
      
        // get local matrix for face entries 
        LocalMatrixType localMatNb = jOp.localMatrix( neighbor,entity );
        const LocalFunctionType uNb = u.localFunction(neighbor);
        // get neighbor's base function set 
        const BasisFunctionSetType &baseSetNb = localMatNb.domainBasisFunctionSet();
        const unsigned int numBasisFunctionsNb = baseSetNb.size();
        std::vector<TemporaryLocalFunctionType> ujNb(numBasisFunctions,TemporaryLocalFunctionType( dfSpace));
        std::vector<TemporaryLocalFunctionType> fujNb(numBasisFunctions,TemporaryLocalFunctionType( dfSpace));

        for( size_t j=0; j<numBasisFunctionsNb; ++j)
          {
            ujNb[j].init(neighbor);
            fujNb[j].init(entity);
            ujNb[j].clear();
            fujNb[j].clear();
            ujNb[j].assign( uNb );
            ujNb[j][j]+=epsilon_;
          }
        
        fu.clear();
        RangeFieldType eps=epsilon_;
        double epsinv=1./eps;
        computeIntersection(intersection,entity,neighbor,area,uLocal,uNb,fu);
        for( unsigned int localCol = 0; localCol < numBasisFunctionsNb; ++localCol )
          {
            computeIntersection(intersection,entity,neighbor,area,uLocal,ujNb[localCol],fujNb[localCol]);
        
            for( unsigned int k = 0; k <numBasisFunctions; ++k)
              {
                double diff=fujNb[localCol][k]-fu[k];
               
                diff*=epsinv;
               localMatNb.add(k,localCol,diff);
              }
          }          
      } 
   }
  } // end grid traversal 
  jOp.communicate();
}











#endif //LOCALFD_OPERATOR_HH

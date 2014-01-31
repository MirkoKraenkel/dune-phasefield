#ifndef LOCALFD_OPERATOR_HH
#define LOCALFD_OPERATOR_HH



#include <dune/fem/function/localfunction/temporary.hh>

#include <dune/fem/operator/common/differentiableoperator.hh>
#include <dune/fem/operator/common/stencil.hh>

//LocalFDOperator
//------------------------------

/** \class LocalFDOperator
 *  \brief operator providing Jacobian throug elementwise finite differences
 *  
 *  \note The user needs to provide the methods localOp and computeIntersection
 *        to calculate the local residuals
 *
 */


template<class DomainFunction,class RangeFunction=DomainFunction,
        class JacobianOperator>
class LocalFDOperatorInterface
 :public Dune::Fem::DifferentiableOperator< JacobianOperator >
{
  
  typedef Dune::Fem::DifferentiableOperator< JacobianOperatorType > 

public:
  typedef typename BaseType::RangeFunctionType RangeFunctionType;
  typedef typename BaseType::DomainFunctionType DomainFunctionType;
  typedef typename BaseType::RangeFieldType RangeFieldType;
  typedef typename BaseType::DomainFieldType DomainFieldType;

  typedef typename BaseType::JacobianOperatorType JacobianOperatorType;
  typedef typename RangeFunctionType::DiscreteFunctionSpaceType RangeSpaceType;
  typedef typename DomainFunctionType::DiscreteFunctionSpaceType DomainSpaceType;
  
  

  typedef Dune::Fem::TemporaryLocalFunction<DiscreteFunctionSpaceType> TemporaryLocalFunctionType;
 
  typedef typename MyOperatorType::GridPartType GridPartType;
  public: 
  LocalFDOperator()
  : globalEpsilon_(Dune::Fem::Parameter::getValue<double>("phasefield.fdjacobian.epsilon"))
  {}

  virtual void jacobian(const DiscreteFunctionType &u, JacobianOperatorType &jOp) const;


  // note this shoudl call jOp.reserve(stencil_) as stencil_ 
  // is only known to the implemetation class
  virtual void reserve(JacobianOperatorType &jOp) const = 0;

  private:
  //let the user choose how to get the fd epsilon(finalize in derived class)
  template<class LocalArgumentType>
  virtual double calculateEpsilon(const LocalArgumentType& u) const = 0;

  private:
  double globalEpsilon_;
};


// Implementation of LocalFDOperator
// // ------------------------------------------------

template< class Jacobian> void
LocalFDOperator<Jacobian>
  ::jacobian ( const DiscreteFunctionType &u, JacobianOperatorType &jOp ) const
{
  typedef typename JacobianOperatorType::LocalMatrixType LocalMatrixType;
  typedef typename DiscreteFunctionSpaceType::BasisFunctionSetType BasisFunctionSetType;
  


  reserve(jOp);
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
    RangeFieldType eps =calculateEpsilon(uLocal);//sqrt((1+sqrt(normU))*epsilon_);
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

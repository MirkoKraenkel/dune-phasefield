#ifndef LOCALFD_OPERATOR_HH
#define LOCALFD_OPERATOR_HH



#include <dune/fem/function/localfunction/temporary.hh>

#include <dune/fem/operator/common/differentiableoperator.hh>
#include <dune/fem/operator/common/stencil.hh>

#include "mixedoperator.hh"



template<class DiscreteFunction,class Model, class Flux, class Jacobian>
class LocalFDOperator
 :public Dune::Fem::DifferentiableOperator < Jacobian >,
  protected DGPhasefieldOperator<DiscreteFunction,Model,Flux>
{
 
  typedef DGPhasefieldOperator<DiscreteFunction,Model,Flux> MyOperatorType;
  
  typedef Dune::Fem::DifferentiableOperator< Jacobian> BaseType;
 

  typedef typename BaseType::JacobianOperatorType JacobianOperatorType;

  typedef typename MyOperatorType::DiscreteFunctionType DiscreteFunctionType;
  typedef typename MyOperatorType::ModelType ModelType;
  typedef typename MyOperatorType::NumericalFluxType NumericalFluxType;
  typedef typename MyOperatorType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
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
  typedef typename DiscreteFunctionSpaceType::BasisFunctionSetType BasisFunctionSetType;


  typedef Dune::Fem::TemporaryLocalFunction<DiscreteFunctionSpaceType> TemporaryLocalFunctionType;
  typedef typename MyOperatorType::GridPartType GridPartType;
  typedef typename GridPartType::IndexSetType IndexSetType;
  public: 
  LocalFDOperator(const ModelType &model,
      const DiscreteFunctionSpaceType &space,
      const NumericalFluxType &flux)
    :MyOperatorType(model,space,flux),
     indexSet_(space_.gridPart().indexSet()),
     visited_(0),
     stencil_(space,space),
     epsilon_(Dune::Fem::Parameter::getValue<double>("phasefield.fdjacobian.epsilon"))
  {}

  using MyOperatorType::localIntegral;
  using MyOperatorType::intersectionIntegral;
  using MyOperatorType::boundaryIntegral;
  using MyOperatorType::setEntity;
  using MyOperatorType::setNeighbor;
  using MyOperatorType::operator();
  using MyOperatorType::setTime;
  using MyOperatorType::setDeltaT;
  using MyOperatorType::setPreviousTimeStep;
  using MyOperatorType::getPreviousTimeStep; 
  using MyOperatorType::space;
  using MyOperatorType::space_;


  typedef Dune::Fem::DiagonalAndNeighborStencil<DiscreteFunctionSpaceType,DiscreteFunctionSpaceType> StencilType;


  void jacobian(const DiscreteFunctionType &u, JacobianOperatorType &jOp) const;
  template< class LocalMatrixImp >
  void computeEntity( const EntityType &entity, const LocalFunctionType &uLocal,const BasisFunctionSetType &baseSet,LocalMatrixImp &jLocal) const;

  private:
  const IndexSetType& indexSet_;
  mutable Dune::Fem::MutableArray<bool> visited_;
  StencilType stencil_;
  mutable double epsilon_;
  mutable  std::vector< RangeType> phi_; 
  mutable  std::vector< JacobianRangeType > dphi_;
  mutable std::vector< RangeType > phiNb_;
  mutable std::vector< JacobianRangeType > dphiNb_;





};
template<class DiscreteFunction,class Model, class Flux, class Jacobian>
template<class LocalMatrixImp>
void LocalFDOperator< DiscreteFunction, Model, Flux,  Jacobian>
::computeEntity( const EntityType &entity, const LocalFunctionType &uLocal, const BasisFunctionSetType& baseSet, LocalMatrixImp &jLocal ) const
{
   const unsigned int numBasisFunctions = baseSet.size();
    const GeometryType geometry = entity.geometry();



    double epsInv=1./epsilon_;

    QuadratureType quadrature( entity, 2*uLocal.order() );
    const size_t numQuadraturePoints = quadrature.nop();

    std::vector<RangeType> uValues(numQuadraturePoints);
    std::vector<JacobianRangeType> uJacobians(numQuadraturePoints);

    uLocal.evaluateQuadrature(quadrature, uValues);
    uLocal.evaluateQuadrature(quadrature,uJacobians);
    phi_.resize(baseSet.size());  
    dphi_.resize(baseSet.size());
    
    for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
    {
      baseSet.evaluateAll( quadrature[ pt ], phi_);
      baseSet.jacobianAll( quadrature[ pt ], dphi_);

      RangeType vu(0.) , fu(0.);
      JacobianRangeType dvu{0.} , fdu{0.};


      localIntegral( pt,
                     geometry,
                     quadrature,
                     uValues[pt],
                     uJacobians[pt],
                     fu,
                     fdu);

      for( size_t jj = 0; jj < numBasisFunctions ; ++jj )
      {
        RangeType ueps , fueps(0.);
        JacobianRangeType dueps{0.} , fdueps{0.};
        ueps=uValues[pt];
        ueps.axpy( epsilon_ , phi_[ jj ] );
        dueps=uJacobians[pt];
        dueps.axpy( epsilon_ , dphi_[ jj ] );

        localIntegral(pt,
            geometry,
            quadrature,
            ueps,
            dueps,
            fueps,
            fdueps);

        fueps-=fu;
        fueps*=epsInv;
        fdueps-=fdu;
        fdueps*=epsInv;
        jLocal.column( jj ).axpy( phi_ , dphi_ , fueps , fdueps );
      }
    }
}
// Implementation of LocalFDOperator
// // ------------------------------------------------

template<class DiscreteFunction,class Model, class Flux, class Jacobian> void
LocalFDOperator< DiscreteFunction, Model, Flux,  Jacobian>
::jacobian ( const DiscreteFunctionType &u, JacobianOperatorType &jOp ) const
{
  typedef typename JacobianOperatorType::LocalMatrixType LocalMatrixType;
  Dune::Fem::DiagonalAndNeighborStencil<DiscreteFunctionSpaceType,DiscreteFunctionSpaceType> stencil(space(),space());

  jOp.reserve(stencil);

  RangeFieldType normU=std::sqrt(u.scalarProductDofs(u));
  jOp.clear();
  //  visited_.resize( indexSet_.size(0));
  // const size_t indSize = visited_.size();
  //  for( size_t ii = 0; i<indSize; ++i) visited_[i] = false;
  visited_.resize( indexSet_.size(0));
  const size_t indSize = visited_.size();
  for( size_t ii = 0; ii < indSize; ++ii) visited_[ii] = false;




  if( epsilon_ <= RangeFieldType( 0 ) )
  {
    const RangeFieldType machine_eps = std::numeric_limits< RangeFieldType >::epsilon();

    epsilon_ = std::sqrt( (RangeFieldType( 1 ) + normU) * machine_eps  );
  }


  double epsInv=1./epsilon_;
  const DiscreteFunctionSpaceType &dfSpace = u.space();
  const GridPartType& gridPart = dfSpace.gridPart();

//  const unsigned int numDofs = dfSpace.blockMapper().maxNumDofs() * 
  //  DiscreteFunctionSpaceType :: localBlockSize ;

  const IteratorType end = dfSpace.end();
  for( IteratorType it = dfSpace.begin(); it != end; ++it )
  {
    const EntityType &entity = *it;
    const LocalFunctionType uLocal = u.localFunction( entity );

    setEntity( entity );
    std::cout<<"EntityIndex="<<indexSet_.index(entity)<<"\n";
       
    LocalMatrixType jLocal = jOp.localMatrix( entity, entity );
 
    const BasisFunctionSetType &baseSet = jLocal.domainBasisFunctionSet();
#if 1    
    computeEntity( entity,
                   uLocal,
                   baseSet,
                   jLocal);
#endif

    if ( space().continuous() )
      continue;

    const IntersectionIteratorType endiit = gridPart.iend( entity );
    for ( IntersectionIteratorType iit = gridPart.ibegin( entity );
        iit != endiit ; ++ iit )
    {
      const unsigned int numBasisFunctions = baseSet.size();
 
      const IntersectionType& intersection = *iit ;

      if( intersection.neighbor() )
      {
        EntityPointerType ep = intersection.outside();
        const EntityType& neighbor = *ep ;
        std::cout<<"NeighborIndex="<<indexSet_.index(neighbor)<<"\n";
        
        if(!visited_[indexSet_.index(neighbor)])
        { 
         setNeighbor( neighbor );
        typedef typename IntersectionType::Geometry  IntersectionGeometryType;
        //const IntersectionGeometryType &intersectionGeometry = intersection.geometry();

        // get local matrix for face entries 
        LocalMatrixType jLocalNb = jOp.localMatrix( neighbor,entity );


        const LocalFunctionType uLocalNb = u.localFunction(neighbor);
        LocalMatrixType jLocalEnNb = jOp.localMatrix( entity, neighbor );
        LocalMatrixType jLocalNbNb = jOp.localMatrix( neighbor,neighbor); 
 
        // get neighbor's base function set 
        const BasisFunctionSetType &baseSetNb = jLocalNb.domainBasisFunctionSet();
        //   const unsigned int numBasisFunctionsNb = baseSetNb.size();
        
        if( neighbor.partitionType()!=InteriorEntity)
        {
         std::cout<<"Neighbor is not Interior\n";
         LocalMatrixType jNbNb=jOp.localMatrix( neighbor,neighbor);
         computeEntity(neighbor,uLocalNb,baseSetNb,jNbNb);
        }

  
    

        const int quadOrderEn = 2*uLocal.order();
        const int quadOrderNb = 2*uLocalNb.order();

        FaceQuadratureType quadInside( space().gridPart(), intersection, quadOrderEn, FaceQuadratureType::INSIDE );
        FaceQuadratureType quadOutside( space().gridPart(), intersection, quadOrderNb, FaceQuadratureType::OUTSIDE );
        const size_t numQuadraturePoints = quadInside.nop();

        std::vector<RangeType> vuEn(numQuadraturePoints);
        std::vector<JacobianRangeType> duEn(numQuadraturePoints);
        std::vector<RangeType> vuNb(numQuadraturePoints);
        std::vector<JacobianRangeType> duNb(numQuadraturePoints);

        uLocal.evaluateQuadrature(quadInside,vuEn);
        uLocal.evaluateQuadrature(quadInside,duEn);

        uLocalNb.evaluateQuadrature(quadOutside,vuNb);
        uLocalNb.evaluateQuadrature(quadOutside,duNb);


        for( size_t pt=0 ; pt < numQuadraturePoints ; ++pt )
        {
          RangeType   avuLeft(0.),avuRight(0.);
          JacobianRangeType aduLeft{0.},aduRight{0.};
          const double weight=quadInside.weight( pt ); 
          const double weightOutside=quadOutside.weight(pt);
            
          phi_.resize(baseSet.size());
          dphi_.resize(baseSet.size());
     
          baseSet.evaluateAll( quadInside[ pt ] , phi_);
          baseSet.jacobianAll( quadInside[ pt ] , dphi_);
         
          phiNb_.resize(baseSetNb.size());
          dphiNb_.resize(baseSetNb.size());
          baseSetNb.evaluateAll( quadOutside[ pt ] , phiNb_ );
          baseSetNb.jacobianAll( quadOutside[ pt ] , dphiNb_ );
          intersectionIntegral( intersection,                  
              pt, 
              quadInside,   
              quadOutside, 
              vuEn[pt],
              vuNb[pt], 
              duEn[pt], 
              duNb[pt],
              avuLeft,
              avuRight,
              aduLeft,
              aduRight);
        
          for( size_t jj=0 ; jj < numBasisFunctions ; ++jj)
          {
            RangeType ueps,fueps(0.),fuepsRight(0.), uepsNb , fuepsNb(0.),fuepsNbRight(0.);
            JacobianRangeType dueps{0.} , fdueps{0.} ,fduepsRight{0.}, duepsNb{0.} , fduepsNb{0.},fduepsNbRight{0.};
            ueps=vuEn[pt];
            ueps.axpy( epsilon_ , phi_[ jj ] );
            dueps=duEn[pt];
            dueps.axpy( epsilon_, dphi_[ jj ] );
            uepsNb=vuNb[pt];
            uepsNb.axpy( epsilon_ , phiNb_[ jj ] );
            duepsNb=duNb[pt];
            duepsNb.axpy( epsilon_, dphiNb_[ jj ] );

            intersectionIntegral( intersection,
                pt,
                quadInside,
                quadOutside,
                ueps,
                vuNb[pt],
                dueps,
                duNb[pt],
                fueps,
                fuepsRight,
                fdueps,
                fduepsRight);
            intersectionIntegral( intersection,
                pt,
                quadInside,
                quadOutside,
                vuEn[pt],
                uepsNb,
                duEn[pt],
                duepsNb,
                fuepsNb,
                fuepsNbRight,
                fduepsNb,
                fduepsNbRight);

            fueps-=avuLeft;
            fueps*=epsInv;
          
            fdueps-=aduLeft;
            fdueps*=epsInv;
           
            fuepsNb-=avuLeft;
            fuepsNb*=epsInv;

            fduepsNb-=aduLeft;
            fduepsNb*=epsInv;
        
            fuepsRight-=avuRight;
            fuepsRight*=epsInv;
              
            fduepsRight-=aduLeft;
            fduepsRight*=epsInv;

            fuepsNbRight-=avuRight;
            fuepsNbRight*=epsInv;
              
            fduepsNbRight-=aduRight;
            fduepsNbRight*=epsInv;

#if 0      
            jLocal.column( jj ).axpy( phi_ , dphi_ , fueps , fdueps, weight );
            jLocalNb.column( jj ).axpy( phi_, dphi_, fuepsNb,fduepsNb,weight); 
            //g(u_Nb+eps*phi_Nb,uNb)*phi_Nb
            jLocalNbNb.column( jj).axpy( phiNb_, dphiNb_, fuepsNbRight,fduepsNbRight, weightOutside);
            //g(u_Nb,u_En+eps*phi_E)*phi_Nb
            jLocalEnNb.column( jj ).axpy( phiNb_, dphiNb_, fuepsRight, fduepsRight, weightOutside);
#endif

          }
        }
      }
      } 
      else if ( intersection.boundary() )
      {
        const int quadOrderEn = 2*uLocal.order();

        FaceQuadratureType quadInside( space().gridPart(), intersection, quadOrderEn, FaceQuadratureType::INSIDE );
        const size_t numQuadraturePoints = quadInside.nop();

        for( size_t pt=0 ; pt < numQuadraturePoints ; ++pt )
        {
          RangeType vuEn,vuNb,avuLeft(0.);
          JacobianRangeType duEn{0.},duNb{0.},aduLeft{0.};
          uLocal.evaluate( quadInside[ pt ], vuEn);
          uLocal.jacobian( quadInside[ pt ], duEn);

          baseSet.evaluateAll( quadInside[ pt ] , phi_);
          baseSet.jacobianAll( quadInside[ pt ] , dphi_);



          boundaryIntegral( intersection,                  
              pt, 
              quadInside,   
              vuEn,
              duEn, 
              avuLeft,
              aduLeft );

          for( size_t jj=0 ; jj < numBasisFunctions ; ++jj)
          {
            RangeType ueps,fueps(0.);
            JacobianRangeType dueps{0.} , fdueps{0.} , duepsNb{0.} , fduepsNb{0.};
            ueps=vuEn;
            ueps.axpy( epsilon_ , phi_[ jj ] );
            dueps=duEn;
            dueps.axpy( epsilon_, dphi_[ jj ] );

            boundaryIntegral( intersection,
                pt,
                quadInside,
                ueps,
                dueps,
                fueps,
                fdueps);


            fueps-=avuLeft;
            fueps*=epsInv;
            fdueps-=aduLeft;
            fdueps*=epsInv;


            jLocal.column( jj ).axpy( phi_ , dphi_ , fueps , fdueps );
          }
        }
      }
    }
  visited_[indexSet_.index( entity )]= true;    
  } // end grid traversal 
  jOp.communicate();
}











#endif //LOCALFD_OPERATOR_HH

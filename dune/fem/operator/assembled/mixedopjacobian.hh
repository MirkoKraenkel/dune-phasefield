#ifndef LOCALFD_OPERATOR_HH
#define LOCALFD_OPERATOR_HH



#include <dune/fem/function/localfunction/temporary.hh>

#include <dune/fem/operator/common/differentiableoperator.hh>
#include <dune/fem/operator/common/stencil.hh>

#include "newheatoperator.hh"
#include "phasefieldfilter.hh"


template<class DiscreteFunction,class Model, class Flux, class Jacobian>
class PhasefieldJacobianOperator
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
  
  typedef Dune::Fem::TemporaryLocalFunction<DiscreteFunctionSpaceType> TemporaryLocalFunctionType;
  typedef typename MyOperatorType::GridPartType GridPartType;
  typedef typename GridPartType::IndexSetType IndexSetType;
  public: 
  LocalFDOperator(const ModelType &model,
                  const DiscreteFunctionSpaceType &space,
                  const NumericalFluxType &flux)
  :MyOperatorType(model,space,flux),
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

   typedef Dune::Fem::DiagonalAndNeighborStencil<DiscreteFunctionSpaceType,DiscreteFunctionSpaceType> StencilType;
  

  void jacobian(const DiscreteFunctionType &u, JacobianOperatorType &jOp) const;
  
  private:
  StencilType stencil_;
  double epsilon_;
};


// Implementation of LocalFDOperator
// // ------------------------------------------------

template<class DiscreteFunction,class Model, class Flux, class Jacobian> void
LocalFDOperator< DiscreteFunction, Model, Flux,  Jacobian>
  ::jacobian ( const DiscreteFunctionType &u, JacobianOperatorType &jOp ) const
{
  typedef typename JacobianOperatorType::LocalMatrixType LocalMatrixType;
  typedef typename DiscreteFunctionSpaceType::BasisFunctionSetType BasisFunctionSetType;
  
  Dune::Fem::DiagonalAndNeighborStencil<DiscreteFunctionSpaceType,DiscreteFunctionSpaceType> stencil(space(),space());

  jOp.reserve(stencil);

  RangeFieldType normU=std::sqrt(u.scalarProductDofs(u));
  jOp.clear();

  double deltaInv=1./deltaT_;
  
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

    //initialize uOld
    setEntity( entity );
    
    LocalMatrixType jLocal = jOp.localMatrix( entity, entity );
    const BasisFunctionSetType &baseSet = jLocal.domainBasisFunctionSet();
    const unsigned int numBasisFunctions = baseSet.size();
 
 
    
    QuadratureType quadrature( entity, 2*dfSpace.order() );
    const size_t numQuadraturePoints = quadrature.nop();
    
    std::vector<RangeType>         uValues(numQuadraturePoints);
    std::vector<JacobianRangeType> uJacobians(numQuadraturePoints);

    std::vector<RangeType>         uOldValues(numQuadraturePoints);
    std::vector<JacobianRangeType> uOldJacobians(numQuadraturePoints);


    uLocal.evaluateQuadrature(quadrature, uValues);
    uLocal.evaluateQuadrature(quadrature,uJacobians);
    
    uOldLocal_.evaluateQuadrature( quadrature, uOldValues); 
    uOldLocal_.evaluateQuadrature( quadrature, uOldJacobians);
    
    const DomainType xgl = geometry.global(x);
    RangeType vuOld{0.},vuMid{0};
    
       
   for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
      {
        baseSet.evaluateAll( quadrature[ pt ], phi);
        baseSet.jacobianAll( quadrature[ pt ], dphi);

        RangeType vu{0.} , vuMid{0.} ,fu{0.};
        JacobianRangeType dvu{0.} , duMid{0.}, fdu{0.};
          
        //(1+theta)/2*U^n+(1-theta)/2*U^(n-1)
        vuMid.axpy(factorImp_,uValuesp[pt]);
        vuMid.axpy(factorExp_,uOldValues[pt]);
  
        JacobianRangeType duOld,duMid ;
       
        //(1+theta)/2*DU^n+(1-theta)/2*DU^(n-1)
        // #if OPCHECK vuMid=vuOld
        duMid.axpy(factorImp_,uJacobians[pt]);
        duMid.axpy(factorExp_,uOldJacobians[pt]);

         
        
        for( size_t jj = 0; jj < numBasisFunctions ; ++jj )
        {
          RangeFieldType div{0.},grad{0.};
          for( size_t ii = 0; ii < dimDomain ; ++ ii)
            {
              div+=Filter::dvelocity(duMid,ii)
            
            }

            Filter::rho(fu)=deltaInv*Filter::rho( phi[ jj ] )
           
          
          
          
          
          
          jLocal.column( jj ).axpy( phi , dphi , fu , fdu );
        }
      } 
   
      
    if ( space().continuous() )
      continue;
    
    const IntersectionIteratorType endiit = gridPart.iend( entity );
    for ( IntersectionIteratorType iit = gridPart.ibegin( entity );
          iit != endiit ; ++ iit )
    {
      const IntersectionType& intersection = *iit ;

      if( intersection.neighbor() )
      {
        EntityPointerType ep = intersection.outside();
        const EntityType& neighbor = *ep ;
        
        setNeighbor( neighbor );
        typedef typename IntersectionType::Geometry  IntersectionGeometryType;
        //const IntersectionGeometryType &intersectionGeometry = intersection.geometry();
      
        // get local matrix for face entries 
        LocalMatrixType jLocalNb = jOp.localMatrix( neighbor,entity );

  
        const LocalFunctionType uLocalNb = u.localFunction(neighbor);
        // get neighbor's base function set 
        const BasisFunctionSetType &baseSetNb = jLocalNb.domainBasisFunctionSet();
     //   const unsigned int numBasisFunctionsNb = baseSetNb.size();
          
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
       //  RangeType vuEn{0.},vuNb{0.},
        RangeType   avuLeft{0.},avuRight{0.};
        //   JacobianRangeType duEn{0.},duNb{0.},
        JacobianRangeType aduLeft{0.},aduRight{0.};
    //       uLocal.evaluate( quadInside[ pt ], vuEn);
    //      uLocal.jacobian( quadInside[ pt ], duEn);
   //        uLocalNb.evaluate( quadOutside[ pt ], vuNb);
   //        uLocalNb.jacobian( quadOutside[ pt ], duNb);
          const double weight=quadInside.weight( pt ); 
          
          baseSet.evaluateAll( quadInside[ pt ] , phi);
          baseSet.jacobianAll( quadInside[ pt ] , dphi);
              
          baseSetNb.evaluateAll( quadOutside[ pt ] , phiNb );
          baseSetNb.jacobianAll( quadOutside[ pt ] , dphiNb );
          
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
             RangeType ueps{0.},fueps{0.},fuepsRight{0.}, uepsNb{0.} , fuepsNb{0.},fuepsNbRight{0.};
             JacobianRangeType dueps{0.} , fdueps{0.} ,fduepsRight{0.}, duepsNb{0.} , fduepsNb{0.},fduepsNbRight{0.};
             ueps=vuEn[pt];
             ueps.axpy( eps , phi[ jj ] );
             dueps=duEn[pt];
             dueps.axpy( eps, dphi[ jj ] );
             uepsNb=vuNb[pt];
             uepsNb.axpy( eps , phiNb[ jj ] );
             duepsNb=duNb[pt];
             duepsNb.axpy( eps, dphiNb[ jj ] );
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
          
            jLocal.column( jj ).axpy( phi , dphi , fueps , fdueps, weight );
            jLocalNb.column( jj ).axpy( phi, dphi, fuepsNb,fduepsNb,weight); 
           
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
              RangeType vuEn{0.},vuNb{0.},avuLeft{0.};
              JacobianRangeType duEn{0.},duNb{0.},aduLeft{0.};
              uLocal.evaluate( quadInside[ pt ], vuEn);
              uLocal.jacobian( quadInside[ pt ], duEn);
           
              baseSet.evaluateAll( quadInside[ pt ] , phi);
              baseSet.jacobianAll( quadInside[ pt ] , dphi);
                
 

              boundaryIntegral( intersection,                  
                                pt, 
                                quadInside,   
                                vuEn,
                                duEn, 
                                avuLeft,
                                aduLeft );
              
              for( size_t jj=0 ; jj < numBasisFunctions ; ++jj)
                {
                  RangeType ueps{0.},fueps{0.}, uepsNb{0.} , fuepsNb{0.};
                  JacobianRangeType dueps{0.} , fdueps{0.} , duepsNb{0.} , fduepsNb{0.};
                  ueps=vuEn;
                  ueps.axpy( eps , phi[ jj ] );
                  dueps=duEn;
                  dueps.axpy( eps, dphi[ jj ] );
             
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
            

                  jLocal.column( jj ).axpy( phi , dphi , fueps , fdueps );
                }
            }
      }
    }
  } // end grid traversal 
  jOp.communicate();
}











#endif //LOCALFD_OPERATOR_HH

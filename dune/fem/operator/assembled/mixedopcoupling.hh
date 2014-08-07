#ifndef PHASEFIELDJACOBIANOPERATOR_HH
#define PHASEFIELDJACOBIANOPERATOR_HH

// Dune::Fem includes
#include <dune/fem/function/localfunction/temporary.hh>
#include <dune/fem/operator/common/differentiableoperator.hh>
#include <dune/fem/operator/common/stencil.hh>
#include <dune/fem/operator/common/temporarylocalmatrix.hh>

// local includes
#include "mixedoperator.hh"
#include "phasefieldfilter.hh"
#include "fluxes/jacobianfluxcoupling.hh"
#include "matrixhelper.hh"
template<class DiscreteFunction,class Model, class Flux, class Jacobian>
class PhasefieldJacobianOperator
:public Dune::Fem::DifferentiableOperator < Jacobian >,
  protected DGPhasefieldOperator<DiscreteFunction,Model,Flux>
{

  typedef DGPhasefieldOperator<DiscreteFunction,Model,Flux> MyOperatorType;

  typedef Dune::Fem::DifferentiableOperator< Jacobian> BaseType;

  enum{dimDomain=MyOperatorType::dimDomain};
  
  enum{ dimRange=MyOperatorType::RangeType::dimension};
  
  typedef typename BaseType::JacobianOperatorType JacobianOperatorType;
  typedef typename MyOperatorType::DiscreteFunctionType DiscreteFunctionType;
  typedef typename MyOperatorType::ModelType ModelType;
  typedef typename MyOperatorType::NumericalFluxType NumericalFluxType;
  typedef typename MyOperatorType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
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
  typedef typename MyOperatorType::Filter Filter;

  typedef  JacobianFlux<ModelType> JacobianFluxType;

  typedef Dune::Fem::TemporaryLocalFunction<DiscreteFunctionSpaceType> TemporaryLocalFunctionType;
  typedef typename MyOperatorType::GridPartType GridPartType;
  typedef typename GridPartType::IndexSetType IndexSetType;
  typedef typename Model::DiffusionTensorType DiffusionType;
#if TEMPORARY
  typedef Dune::Fem::TemporaryLocalMatrix< DiscreteFunctionSpaceType, DiscreteFunctionSpaceType> LocalMatrixType;
  typedef typename JacobianOperatorType::LocalMatrixType RealLocalMatrixType;
#else
  typedef typename JacobianOperatorType::LocalMatrixType LocalMatrixType;
#endif
  typedef Dune::Fem::MutableArray<RangeType> RangeVectorType;
  typedef Dune::Fem::MutableArray<JacobianRangeType> JacobianVectorType;
  
  typedef typename MatrixHelper::Alignment< dimRange > DofAlignmentType;


  typedef typename Dune::FieldMatrix<double,dimRange,dimRange> FluxRangeType;
  typedef typename Dune::FieldVector<double,1> ComponentRangeType;
  typedef typename Dune::FieldMatrix<double,1,dimDomain> ComponentJacobianType;
  typedef typename std::array<DomainType, dimDomain> DiffusionValueType;

  public: 
  PhasefieldJacobianOperator(const ModelType &model,
      const DiscreteFunctionSpaceType &space,
      int volQuadOrder=-1)
    :MyOperatorType(model,space),
    indexSet_(space_.gridPart().indexSet()),
    visited_(0),
    stencil_(space,space),
    jacFlux_(model),
    zerocount_(0),
    allcount_(0)
  {}

  using MyOperatorType::localIntegral;
  using MyOperatorType::intersectionIntegral;
  using MyOperatorType::boundaryIntegral;
  using MyOperatorType::setEntity;
  using MyOperatorType::setNeighbor;
  using MyOperatorType::operator();
  using MyOperatorType::setTime;
  using MyOperatorType::getTime;
  using MyOperatorType::setDeltaT;
  using MyOperatorType::setPreviousTimeStep;
  using MyOperatorType::getPreviousTimeStep; 
  using MyOperatorType::space;
  using MyOperatorType::space_;
  using MyOperatorType::deltaT_;
  using MyOperatorType::uOldLocal_;
  using MyOperatorType::uOldNeighbor_;
  using MyOperatorType::factorImp_;
  using MyOperatorType::factorExp_;
  using MyOperatorType::model_;
  using MyOperatorType::areaEn_;
  using MyOperatorType::areaNb_;
  using MyOperatorType::flux_;
  using MyOperatorType::time_;



  typedef Dune::Fem::DiagonalAndNeighborStencil<DiscreteFunctionSpaceType,DiscreteFunctionSpaceType> StencilType;
 
  template< class RangeVector, class JacobianVector>
  void myaxpy( const size_t jj,
                const RangeVector& phi,
                const JacobianVector& dphi,
                const RangeType& factor,
                const JacobianRangeType& jacobianFactor,
                const double weight,
                LocalMatrixType& jLocal) const;
   
   template< class JacobianVector,class DiffusionTensor>
   void diffusionaxpy( const size_t local_i,
                 const size_t local_j,
                 const JacobianVector& dphi,
                 const DiffusionTensor du,
                 const double weight,
                 LocalMatrixType& jLocal) const;

  
  template< class Matrix1, class Matrix2>
  void add( const Matrix1& m1, Matrix2& m2) const
  {
   for( int ii = 0 ; ii<m2.rows() ; ++ii)
    for( int jj = 0 ; jj<m2.columns() ; ++jj)
    m2.add( ii , jj, m1.get(ii , jj));
  }
  
  int globalbf( int scalarbf, int component) const
  {
    return scalarbf*dimRange+component;
  }

  void jacobian(const DiscreteFunctionType &u, JacobianOperatorType &jOp) const;
  //template< class IntersectionQuad >
  //void assembleIntersection
  
  private:
  const IndexSetType& indexSet_;
  StencilType stencil_;
  mutable Dune::Fem::MutableArray<bool> visited_;
  const JacobianFluxType jacFlux_;
  mutable RangeVectorType uEn_,uNb_,uEnOld_,uNbOld_;
  mutable JacobianVectorType duEn_,duNb_,duEnOld_,duNbOld_;
  mutable int zerocount_;
  mutable int allcount_; 
  mutable MatrixHelper::Couplings<dimDomain> couplings_;
};


// Implementation of LocalFDOperator
// // ------------------------------------------------

template<class DiscreteFunction,class Model, class Flux, class Jacobian> 
template< class RangeVector, class JacobianVector>
void PhasefieldJacobianOperator< DiscreteFunction, Model, Flux,  Jacobian>
::myaxpy( const size_t jj ,
          const RangeVector& phi,
          const JacobianVector& dphi,
          const RangeType& factor,
          const JacobianRangeType& jacobianFactor,
          const double weight,
          LocalMatrixType& jLocal) const
{
 
  const unsigned int numBasisFunctions = jLocal.rows();
  int scalarNumBasisFunctions=numBasisFunctions/RangeType::dimension;
  int dim=RangeType::dimension; 
  for( int i = 0;i < numBasisFunctions; ++i)
   {
     int range=i%dim;
     int scalarbf=i/dim;
     RangeFieldType value = factor[range]*phi[i][range];

    for( int k = 0; k < jacobianFactor.rows; ++k )
          value += jacobianFactor[ k ][range] * dphi[ i ][ k ][range];
      
     jLocal.add( i , jj , weight * value );
   }
}



template<class DiscreteFunction,class Model, class Flux, class Jacobian> 
template< class JacobianVector, class DiffusionTensor>
void PhasefieldJacobianOperator< DiscreteFunction, Model, Flux,  Jacobian>
::diffusionaxpy( const size_t local_i,
                 const size_t local_j,
                 const JacobianVector& dphi,
                 const DiffusionTensor du,
                 const double weight,
                 LocalMatrixType& jLocal) const
{   
   //for each velocomponent
  for( size_t ii = 0 ; ii < dimDomain ; ++ii )
    { 
      size_t global_i=local_i*dimRange + (1+ii);//veloDof( local_i, comp);
      assert( global_i == DofAlignmentType::vectorialIndex( 1 + ii , local_i) ); 
      for( size_t jj = 0 ; jj < dimDomain; ++jj )
        {
          size_t global_j =  local_j*dimRange+(1+jj);//veloDof( local_j, comp);
          assert( global_j == DofAlignmentType::vectorialIndex( 1 + jj , local_j) ); 
          double  value(0);
          value=du[ jj ][ ii ]*dphi[ local_i ][ 0 ];         
                
          jLocal.add( global_i , global_j ,0.5*value*weight );
        }
     }
}




template<class DiscreteFunction,class Model, class Flux, class Jacobian> void
PhasefieldJacobianOperator< DiscreteFunction, Model, Flux,  Jacobian>
::jacobian ( const DiscreteFunctionType &u, JacobianOperatorType &jOp ) const
{
 typedef typename DiscreteFunctionSpaceType::BasisFunctionSetType BasisFunctionSetType;

  Dune::Fem::DiagonalAndNeighborStencil<DiscreteFunctionSpaceType,DiscreteFunctionSpaceType> stencil(space(),space());
  jOp.reserve(stencil);

  RangeFieldType normU=std::sqrt(u.scalarProductDofs(u));
  jOp.clear();
  //intialize visited marker
  visited_.resize( indexSet_.size(0));
  const size_t indSize = visited_.size();
  for( size_t ii = 0; ii < indSize; ++ii) visited_[ii] = false;


  double deltaTInv=1./deltaT_;

  const DiscreteFunctionSpaceType &dfSpace = u.space();
  const GridPartType& gridPart = dfSpace.gridPart();

  const unsigned int numDofs = dfSpace.blockMapper().maxNumDofs() * 
    DiscreteFunctionSpaceType :: localBlockSize ;


#if VECTORIAL
  std::vector< RangeType> phi( numDofs ); 
  std::vector< JacobianRangeType > dphi( numDofs );
  std::vector< RangeType > phiNb( numDofs );
  std::vector< JacobianRangeType > dphiNb( numDofs );
#else
  std::vector< ComponentRangeType> phi( numDofs/dimRange);
  std::vector< ComponentJacobianType> dphi( numDofs/dimRange);
  std::vector< ComponentRangeType> phiNb( numDofs/dimRange);
  std::vector< ComponentJacobianType> dphiNb( numDofs/dimRange);
#endif
 
 
  const IteratorType end = dfSpace.end();
  for( IteratorType it = dfSpace.begin(); it != end; ++it )
  {
    const EntityType &entity = *it;
    const GeometryType geometry = entity.geometry();

    const LocalFunctionType uLocal = u.localFunction( entity );

    //initialize uOld
    setEntity( entity );
#if TEMPORARY     
    LocalMatrixType jLocal( dfSpace, dfSpace);//= jOp.localMatrix( entity, entity );
    jLocal.init(entity,entity);
    RealLocalMatrixType realLocal=jOp.localMatrix( entity, entity);
#else   
    LocalMatrixType jLocal = jOp.localMatrix( entity, entity );
#endif
    const BasisFunctionSetType &baseSet = jLocal.domainBasisFunctionSet();
    const unsigned int numBasisFunctions = baseSet.size();
      
    const unsigned int numScalarBf = numBasisFunctions/dimRange;

    QuadratureType quadrature( entity, 2*dfSpace.order() );
    const size_t numQuadraturePoints = quadrature.nop();
    
    uEn_.resize(numQuadraturePoints);
    uLocal.evaluateQuadrature( quadrature, uEn_);
    duEn_.resize(numQuadraturePoints); 
    uLocal.evaluateQuadrature( quadrature, duEn_);

    uEnOld_.resize(numQuadraturePoints);
    uOldLocal_.evaluateQuadrature( quadrature, uEnOld_); 
    duEnOld_.resize(numQuadraturePoints);
    uOldLocal_.evaluateQuadrature( quadrature, duEnOld_);

    RangeType vuOld(0.),vuMid(0);
     
    for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
    {
      const typename QuadratureType::CoordinateType &x = quadrature.point( pt );
      const double weight = quadrature.weight( pt )* geometry.integrationElement( x );

#if VECTORIAL
      baseSet.evaluateAll( quadrature[ pt ], phi);
      baseSet.jacobianAll( quadrature[ pt ], dphi);
#else
      MatrixHelper::evaluateScalarAll( quadrature[ pt ], baseSet.shapeFunctionSet(),phi);
      MatrixHelper::jacobianScalarAll( quadrature[ pt ], geometry, baseSet.shapeFunctionSet(),dphi);      
//      baseSet.evaluateScalarAll( quadrature[ pt ], phi);
//      baseSet.jacobianScalarAll( quadrature[ pt ], dphi);
#endif
    
      
      RangeType vu(0.) , vuMid(0.) ,fu(0.);
      JacobianRangeType dvu(0.) , duMid(0.), fdu(0.);
      FluxRangeType flux(0.);
      vuOld=uEnOld_[pt];
      vu=uEn_[ pt ];
      //(1+theta)/2*U^n+(1-theta)/2*U^(n-1)
      vuMid.axpy(factorImp_,uEn_[pt]);
      vuMid.axpy(factorExp_,uEnOld_[pt]);

      //(1+theta)/2*DU^n+(1-theta)/2*DU^(n-1)
      // #if OPCHECK vuMid=vuOld
      duMid.axpy(factorImp_,duEn_[pt]);
      duMid.axpy(factorExp_,duEnOld_[pt]);
      
      for( size_t jj=0 ;  jj < numScalarBf ; ++jj)
        {
          flux=0.;
          //evaluate Operator  for all basisFunctions
          //rhog
          flux[0][0]=deltaTInv*phi[ jj ][0];
          
          for( size_t kk = 0 ; kk < dimDomain ; ++kk )
            {
              //div u*phi_rho + v\cdot\nabla phi_rho 
              flux[0][0]+=0.5*(duMid[ 1+kk ][ kk ]*phi[ jj ][ 0 ]+vuMid[ 1+kk ]*dphi[ jj ][ 0 ][ kk ]);
              // rho*div phi_v +\nabla\rho\cdot\phi_v
              flux[ 0 ][ 1+kk ]=0.5*vuMid[0]*dphi[ jj ][0 ][kk]+duMid[0][ kk ]*phi[ jj ][0];
      
              flux[1+kk][0] = (vu[ 1 + kk ]-vuOld[ 1 + kk ])*phi[ jj ][0]*0.5;
              flux[1+kk][1+kk]=phi[ jj ][ 0 ]*vuMid[0]*deltaTInv;
   
              for( size_t ll = 0; ll<dimDomain; ++ll)
                { 
                  flux[ 1+kk ][ 0 ]+=0.5*(duMid[ 1 + kk ][ kk ] - duMid[ 1 + ll ][ kk ])*vuMid[ 1 + ll ]*phi[ jj ][ 0 ];
                  flux[ 1+kk ][ 1+ll ]+=0.5*(duMid[ 1 + kk ][ kk ] - duMid[ 1 + ll ][ kk ])* phi[ jj ][ 0 ]*vuMid[0];
                  flux[ 1+kk ][ 1+kk ]+=0.5*( dphi[ jj ][ 0 ][ ll ] )*vuMid[ 1 + ll ]*vuMid[ 0 ];
                  flux[ 1+kk ][ 1+ll ]-=0.5*( dphi[ jj ][ 0 ][ kk ] )*vuMid[ 1 + ll ]*vuMid[ 0 ];
                } 
               
               flux[ 1+kk ][ 0 ]+=0.5*phi[ jj ][ 0 ]*duMid[ 2 + dimDomain ][ kk ];
               flux[ 1+kk ][ dimDomain+2]+=0.5*vuMid[ 0 ]*dphi[ jj ][ 0 ][ kk ];
               flux[ 1+kk][ dimDomain+1]-=0.5*vuMid[ dimDomain + 3 ]*dphi[ jj ][0][ kk ];

               flux[ 1+kk][ dimDomain+3]-=0.5*phi[ jj ][ 0 ]*duMid[ dimDomain + 1 ][ kk ] ;
       

            }
          //phi
          
          //( - tau phi_rho)/ rho*rho
          flux[ dimDomain+1 ][ 0 ]-=0.5* vuMid[ dimDomain + 3 ]*phi[ jj ][ 0 ]/(vuMid[0]*vuMid[0]);
          // phi_phi/deltaT
          flux[ dimDomain+1 ][ dimDomain+1]=phi[jj][0]*deltaTInv;
          for( size_t kk = 0 ; kk < dimDomain ; ++kk )
            {
              // \nabla phi_phi\cdot v
              flux[ dimDomain+1 ][ dimDomain+1]+=0.5*dphi[ jj ][0][kk]*vuMid[ 1+ kk];
              //\nabla phi \cdot phi_v
              flux[ dimDomain+1 ][ 1+kk ] +=0.5* duMid[dimDomain+1][kk]*phi[ jj ][0];
            }
          //(phi_tau rho)/ rho*rho
          flux[ dimDomain+1 ][ dimDomain+3 ]=0.5*phi[ jj ][ 0 ]*vuMid[ 0 ]/(vuMid[0]*vuMid[0]);
          //mu
          RangeFieldType drhomu(0.),dphimu(0.);
          model_.drhomuSource( vu[ 0 ] , vu[ 0 ] , vu[ dimDomain + 1 ] , drhomu );
          model_.dphimuSource( vu[ 0 ] , vu[ 0 ] , vu[ dimDomain + 1 ] , dphimu );
          
          flux[ dimDomain+2 ][ 0 ]=-1*drhomu*phi[ jj ][ 0 ];
          
          for( size_t kk = 0 ; kk < dimDomain ; ++kk )
            flux[ dimDomain+2 ][ 1+kk]-=2*vu[ 1+kk ]*phi[ jj ][ 0 ]*0.25;
      
          flux[ dimDomain+2 ][ dimDomain+1 ]=-1*dphimu*phi[ jj ][ 0];
          flux[ dimDomain+2 ][ dimDomain+2 ]=0.5*phi[ jj ][ 0 ];
        
          //tau
          RangeFieldType dphitau;
          model_.dphitauSource( vu[dimDomain+1],vu[dimDomain+1], vu[0],dphitau);
          flux[ dimDomain+3 ][ dimDomain+1 ]-=dphitau*phi[ jj ][ 0 ];
          flux[ dimDomain+3 ][ dimDomain+3 ]+=0.5*phi[ jj ][ 0 ];
          
          for( size_t kk = 0 ; kk < dimDomain ; ++kk )
            {
              flux[ dimDomain+3 ][ dimDomain+4+kk ]+=model_.delta()*0.5*dphi[ jj ][ 0 ][ kk ];
              flux[ dimDomain+4+kk ][ dimDomain+1 ]-=dphi[ jj ][ 0 ][ kk ];
              flux[ dimDomain+4+kk ][ dimDomain+4+kk]=phi[ jj ][ 0 ];
            }  
          DiffusionType du;
        
          for(int i = 0;  i<dimDomain ; ++i)
            du[i]=0;
 
          model_.scalar2vectorialDiffusion( dphi[jj], du);  
          
          for( size_t ii = 0 ; ii < numScalarBf ; ++ii )
            {
              diffusionaxpy( ii ,jj , dphi, du , weight, jLocal);
          
              MatrixHelper::axpyElement< DofAlignmentType >( couplings_,
                                                             phi,
                                                             flux,
                                                             ii,
                                                             jj,
                                                             weight,
                                                             jLocal);
            } 
      }
    } 
#if TEMPORARY   
    add(jLocal, realLocal);
    jLocal.clear();
#endif
    if ( !space().continuous() )
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
                setNeighbor( neighbor );
                typedef typename IntersectionType::Geometry  IntersectionGeometryType;
                const IntersectionGeometryType &intersectionGeometry = intersection.geometry();
#if TEMPORARY  
                LocalMatrixType jLocalNbEn( dfSpace,dfSpace);
                jLocalNbEn.init( neighbor, entity);
                LocalMatrixType jLocalEnNb( dfSpace, dfSpace);
                jLocalEnNb.init( entity, neighbor);
                LocalMatrixType jLocalNbNb( dfSpace,dfSpace);
                jLocalNbNb.init( neighbor,neighbor);
                // get local matrix for face entries 
                RealLocalMatrixType realLocalNbEn = jOp.localMatrix( neighbor,entity );
                RealLocalMatrixType realLocalEnNb = jOp.localMatrix( entity, neighbor );
                // get local matrix on neighbor 
                RealLocalMatrixType realLocalNbNb = jOp.localMatrix( neighbor,neighbor); 


#else
                // get local matrix for face entries 
                LocalMatrixType jLocalNbEn = jOp.localMatrix( neighbor,entity );
                LocalMatrixType jLocalEnNb = jOp.localMatrix( entity, neighbor );
                // get local matrix on neighbor 
                LocalMatrixType jLocalNbNb = jOp.localMatrix( neighbor,neighbor); 

#endif
                const LocalFunctionType uLocalNb = u.localFunction(neighbor);
                // get neighbor's base function set 
                const BasisFunctionSetType &baseSetNb = jLocalNbEn.domainBasisFunctionSet();
                //const unsigned int numBasisFunctionsNb = baseSetNb.size();

                const int quadOrderEn = 2*uLocal.order();
                const int quadOrderNb = 2*uLocalNb.order();

                FaceQuadratureType quadInside( space().gridPart(), intersection, quadOrderEn, FaceQuadratureType::INSIDE );
                FaceQuadratureType quadOutside( space().gridPart(), intersection, quadOrderNb, FaceQuadratureType::OUTSIDE );
                const size_t numQuadraturePoints = quadInside.nop();
            
                uEn_.resize( numQuadraturePoints );
                uLocal.evaluateQuadrature(quadInside,uEn_);
                duEn_.resize( numQuadraturePoints );
                uLocal.evaluateQuadrature(quadInside,duEn_);
                uNb_.resize( numQuadraturePoints );
                uLocalNb.evaluateQuadrature(quadOutside,uNb_);
                duNb_.resize( numQuadraturePoints );
                uLocalNb.evaluateQuadrature(quadOutside,duNb_);

                uEnOld_.resize( numQuadraturePoints );
                uOldLocal_.evaluateQuadrature(quadInside,uEnOld_);
                duEnOld_.resize( numQuadraturePoints );
                uOldLocal_.evaluateQuadrature(quadInside,duEnOld_);
                uNbOld_.resize( numQuadraturePoints ); 
                uOldNeighbor_.evaluateQuadrature(quadOutside,uNbOld_);
                duNbOld_.resize( numQuadraturePoints );
                uOldNeighbor_.evaluateQuadrature(quadOutside,duNbOld_);


                for( size_t pt=0 ; pt < numQuadraturePoints ; ++pt )
                  {
                    RangeType    vuMidEn(0.), vuMidNb(0.);
                    JacobianRangeType aduLeft(0.),aduRight(0.),duMidNb(0.), duMidEn(0.);
                    FluxRangeType fluxLeft(0.), fluxRight(0.);
                    
                    const double weightInside=quadInside.weight( pt ); 
                    const double weightOutside=quadOutside.weight( pt ); 
#if VECTORIAL
                    baseSet.evaluateAll( quadInside[ pt ] , phi);
                    baseSet.jacobianAll( quadInside[ pt ] , dphi);

                    baseSetNb.evaluateAll( quadOutside[ pt ] , phiNb );
                    baseSetNb.jacobianAll( quadOutside[ pt ] , dphiNb );
#else
    
                    
                    //baseSet.evaluateScalarAll( quadInside[ pt ] , phi);
                    //baseSet.jacobianScalarAll( quadInside[ pt ] , dphi);
                    MatrixHelper::evaluateScalarAll( quadInside[ pt ], baseSet.shapeFunctionSet(), phi );
                    MatrixHelper::jacobianScalarAll( quadInside[ pt ], geometry, baseSet.shapeFunctionSet(), dphi);

                    //baseSetNb.evaluateScalarAll( quadOutside[ pt ] , phiNb );
                    //baseSetNb.jacobianScalarAll( quadOutside[ pt ] , dphiNb );
                    MatrixHelper::evaluateScalarAll( quadOutside[ pt ], baseSetNb.shapeFunctionSet(), phiNb);
                    MatrixHelper::jacobianScalarAll( quadOutside[ pt ], geometryNb, baseSetNb.shapeFunctionSet(), dphiNb);


#endif
                    vuMidEn.axpy( factorImp_ , uEn_[ pt ] );
                    vuMidEn.axpy( factorExp_ , uEnOld_[ pt ] );

                    vuMidNb.axpy( factorImp_ , uNb_[ pt ] );
                    vuMidNb.axpy( factorExp_ , uNbOld_[ pt ]);

                    duMidEn.axpy( factorImp_ , duEn_[ pt ] );
                    duMidEn.axpy( factorExp_ , duEnOld_[ pt ] );

                    duMidNb.axpy( factorImp_ , duNb_[ pt ] );
                    duMidNb.axpy( factorExp_ , duNbOld_[ pt ] );

                    const typename FaceQuadratureType::LocalCoordinateType &x = quadInside.localPoint( pt );

                    const DomainType normal = intersection.integrationOuterNormal( x );

                    // compute penalty factor
                    const double intersectionArea = normal.two_norm();
                    const double penaltyFactor = intersectionArea / std::min( areaEn_, areaNb_ ); 
                    const double area=std::min(areaEn_,areaNb_); 



                    jacFlux_.numericalFlux( normal,
                                            area,
                                            vuMidEn,
                                            vuMidNb,
                                            fluxLeft,
                                            fluxRight);

                    for( size_t jj=0 ; jj < numScalarBf ; ++jj)
                      {
                        RangeType avuLeft(0.), avuRight(0.), valueLeft(0.),valueRight(0.);
                        DiffusionType aLeft,aRight;
                        DiffusionValueType bLeft, bRight; 

                        jacFlux_.scalar2vectorialDiffusionFlux( normal,
                                                                penaltyFactor,
                                                                phi[ jj ],
                                                                phiNb[ jj ],
                                                                dphi[ jj ],
                                                                dphiNb[ jj ],
                                                                aLeft,
                                                                aRight,
                                                                bLeft,
                                                                bRight );


                        for( size_t ii = 0; ii < numScalarBf ; ++ii )
                          {
                            MatrixHelper::axpyIntersection< DofAlignmentType >( couplings_,
                                                                                phi,
                                                                                phi,
                                                                                fluxLeft,
                                                                                ii,
                                                                                jj,
                                                                                weightInside,
                                                                                jLocal );

                            MatrixHelper::axpyIntersection< DofAlignmentType >( couplings_,
                                                                                phi,
                                                                                phiNb,
                                                                                fluxRight,
                                                                                ii,
                                                                                jj,
                                                                                weightInside,
                                                                                jLocalNbEn );

#if 1 
                            MatrixHelper::axpyIntersection< DofAlignmentType >( couplings_,
                                                                                phiNb,
                                                                                phi,
                                                                                fluxLeft,
                                                                                ii,
                                                                                jj,
                                                                                weightOutside,
                                                                                jLocalEnNb );

                            MatrixHelper::axpyIntersection< DofAlignmentType >( couplings_,
                                                                                phiNb,
                                                                                phiNb,
                                                                                fluxRight,
                                                                                ii,
                                                                                jj,
                                                                                weightOutside,
                                                                                jLocalNbNb );
#endif  
#if 1
                        for(int i  = 0; i  < dimDomain ; ++i )
                          {
                            int global_i=ii*dimRange+1+i; 
                            for(int j  = 0 ; j  < dimDomain ; ++j )
                              {
                                int global_j= jj*dimRange+1+j;
                                double valueEn,valueNb; 
                                
                                valueEn=aLeft[ j ][ i ]*dphi[ ii ][ 0 ];
                                valueEn+=bLeft[ j ][ i ]*phi[ ii ];
                                valueNb=-1*(aRight[ j ][ i ]*dphi[ ii ][ 0 ]);
                                valueNb+=bRight[ j ][ i ]*phi[ ii ];

                                jLocal.add( global_i , global_j , weightInside*valueEn*0.5); 
                                jLocalNbEn.add( global_i , global_j , weightInside*valueNb*0.5);
                                
                                valueEn=(aLeft[ j ][ i ]*dphiNb[ ii ][ 0 ]);
                                valueEn-=bLeft[ j ][ i ]*phiNb[ ii ];
                                valueNb=-1*(aRight[ j ][ i ]*dphiNb[ ii][ 0 ]);
                                valueNb-=bRight[ j ][ i ]*phiNb[ ii ];

                                jLocalEnNb.add( global_i , global_j , weightOutside*valueEn*0.5); 
                                jLocalNbNb.add( global_i , global_j , weightOutside*valueNb*0.5);
                               
                              }
                          }
#endif
                      } 
                }
          }
#if TEMPORARY 
        add( jLocal, realLocal);
        add( jLocalEnNb, realLocalEnNb);
        add( jLocalNbEn, realLocalNbEn);
        add( jLocalNbNb, realLocalNbNb);
#endif
    
        } 
      
      }
      else if ( intersection.boundary() )
      {
        const int quadOrderEn = 2*uLocal.order();
         typedef typename IntersectionType::Geometry  IntersectionGeometryType;
        const IntersectionGeometryType &intersectionGeometry = intersection.geometry();


        FaceQuadratureType quadInside( space().gridPart(), intersection, quadOrderEn, FaceQuadratureType::INSIDE );
        const size_t numQuadraturePoints = quadInside.nop();
        std::vector<RangeType> vuEn(numQuadraturePoints);
        std::vector<JacobianRangeType> duEn(numQuadraturePoints);
        std::vector<RangeType> vuOldEn(numQuadraturePoints);
        std::vector<JacobianRangeType> duOldEn(numQuadraturePoints);
 
        uLocal.evaluateQuadrature(quadInside,vuEn);
        uLocal.evaluateQuadrature(quadInside,duEn);

        uOldLocal_.evaluateQuadrature(quadInside,vuOldEn);
        uOldLocal_.evaluateQuadrature(quadInside,duOldEn);

   
        for( size_t pt=0 ; pt < numQuadraturePoints ; ++pt )
        {
          RangeType vuMidEn(0.),avuLeft(0.),bndValue(0.);
          JacobianRangeType duMidEn(0.),aduLeft(0.);
          
          const double weight=quadInside.weight( pt ); 

          vuMidEn.axpy( factorImp_ , vuEn[ pt ] );
          vuMidEn.axpy( factorExp_ , vuOldEn[ pt ] );

          duMidEn.axpy( factorImp_ , duEn[ pt ] );
          duMidEn.axpy( factorExp_ , duOldEn[ pt ] );

         // baseSet.evaluateAll( quadInside[ pt ] , phi);
          //baseSet.jacobianAll( quadInside[ pt ] , dphi);

          const typename FaceQuadratureType::LocalCoordinateType &x = quadInside.localPoint( pt );

          const DomainType normal = intersection.integrationOuterNormal( x );
          const DomainType xgl = intersectionGeometry.global(x);
          model_.dirichletValue( time_,xgl, bndValue);


          // compute penalty factor
          const double intersectionArea = normal.two_norm();
          const double penaltyFactor = intersectionArea /  areaEn_; 
          const double area=areaEn_; 

 
          for( size_t jj=0 ; jj < numBasisFunctions ; ++jj)
          {
            RangeType avuLeft(0.), avuRight(0.),dummy(0.), valueLeft(0.),valueRight(0.);
            JacobianRangeType aduLeft(0.),aduRight(0.);
       #if  0 
            double fluxRet=jacFlux_.numericalFlux( normal,
                                                   area,
                                                   vuEn[ pt ],
                                                   bndValue, 
                                                   vuMidEn,
                                                   bndValue,
                                                   phi[ jj ],
                                                   dummy,
                                                   avuLeft,
                                                   avuRight); 

          #endif

#if 0 
            double  fluxRet=jacFlux_.boundaryFlux( normal,
                                                   area,
                                                   vuEn[ pt ],
                                                   vuMidEn,
                                                   phi[ jj ],
                                                   avuLeft ); 

#endif
#if 0 
              double fluxRet;
              fluxRet+=jacFlux_.diffusionBoundaryFlux( normal,
                                                     penaltyFactor,
                                                     phi[ jj ],
                                                     dphi[ jj ],
                                                     valueLeft,
                                                     aduLeft );                
#endif
            avuLeft+=valueLeft;
//#if MYAXPY
          myaxpy( jj, phi, dphi, avuLeft, aduLeft, weight, jLocal);
//#else
 //         jLocal.column( jj ).axpy( phi , dphi , avuLeft,aduLeft , weight );
//#endif


          }
        }
      }
    }
  visited_[indexSet_.index( entity )]=true;
  } // end grid traversal 
  jOp.communicate();
  std::cout<<"addFactor:"<<static_cast<double>(zerocount_)/allcount_<<" zeros:"<<zerocount_<<"\n";
}











#endif //PHASEFIELDJACOBIANOPERATOR_HH


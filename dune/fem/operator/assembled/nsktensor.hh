#ifndef PHASEFIELD_MIXED_TENSOR_HH
#define PHASEFIELD_MIXED_TENSOR_HH


// Dune::Fem includes
#include <dune/fem/function/localfunction/temporary.hh>
#include <dune/fem/operator/common/differentiableoperator.hh>
#include <dune/fem/operator/common/stencil.hh>
#include <dune/fem/operator/common/temporarylocalmatrix.hh>

// local includes
#include "matrixhelper.hh"

template<class DiscreteFunction,class Model, class Flux, class JacobianFlux>
class NskMixedTensor:
public NskIntegrator<DiscreteFunction,Model,Flux>
{
  
  using BaseType= NskIntegrator<DiscreteFunction,Model,Flux>;

 public:
  using DiscreteFunctionType      = typename BaseType::DiscreteFunctionType;
  using ModelType                 = typename BaseType::ModelType;
  using DiscreteFunctionSpaceType = typename BaseType::DiscreteFunctionSpaceType;
  using BasisFunctionSetType      = typename DiscreteFunctionSpaceType::BasisFunctionSetType;  
  using RangeFieldType            = typename BaseType::RangeFieldType;
  using LocalFunctionType         = typename BaseType::LocalFunctionType;
  using RangeType                 = typename BaseType::RangeType;
  using JacobianRangeType         = typename BaseType::JacobianRangeType;
  static const int dimDomain      = BaseType::dimDomain;
  static const int dimRange       = BaseType::dimRange;
 
  using IteratorType              = typename BaseType::IteratorType;
  using EntityType                = typename BaseType::EntityType;

  using GeometryType              = typename BaseType::GeometryType;

  using DomainType                = typename BaseType::DomainType; 

  using GridPartType              = typename BaseType::GridPartType;
  using IntersectionIteratorType  = typename BaseType::IntersectionIteratorType;
  using IntersectionType          = typename BaseType::IntersectionType;

  using QuadratureType            = typename BaseType::QuadratureType; 
  using FaceQuadratureType        = typename BaseType::FaceQuadratureType;

  using JacobianFluxType=JacobianFlux;

  using TemporaryLocalTypex      = typename BaseType::TemporaryLocalType;
  
  typedef Dune::Fem::TemporaryLocalMatrix<DiscreteFunctionSpaceType,
                                          DiscreteFunctionSpaceType> LocalMatrixType;
  
  typedef Dune::Fem::MutableArray<RangeType> RangeVectorType;
  typedef Dune::Fem::MutableArray<JacobianRangeType> JacobianVectorType;
  typedef typename MatrixHelper::Alignment< dimRange > DofAlignmentType;

  typedef typename Dune::FieldMatrix<double,dimRange,dimRange> TensorRangeType;
  typedef typename Dune::FieldVector<double,1> ComponentRangeType;
  typedef typename Dune::FieldMatrix<double,1,dimDomain> ComponentJacobianType;
  typedef std::vector< ComponentRangeType> BasefunctionStorage;
  typedef std::vector< ComponentJacobianType> BaseJacobianStorage;
  
  typedef typename std::array<DomainType, dimDomain> DiffusionValueType;
  using  DiffusionType=typename Model::DiffusionTensorType; 


  public: 
  NskMixedTensor(const ModelType &model,
      const DiscreteFunctionSpaceType &space)
    :BaseType(model,space),
    jacFlux_(model)
    {}

    template<class Quad>
    void computeCoefficientsEn ( const LocalFunctionType& uLocal,
                                const Quad& quadrature ) const;
    template<class Quad>
    void computeCoefficientsNb ( const LocalFunctionType& uLocal,
                               const Quad& quadrature ) const;

    void resizeBaseStorageEn ( int numDofs ) const
      {
        const unsigned int scalarDofs=numDofs/dimRange;
        phi_.resize( scalarDofs );
        dphi_.resize( scalarDofs );
      }
  
    void resizeBaseStorageNb ( int numDofs ) const
      {
        const unsigned int scalarDofs=numDofs/dimRange;
        phiNb_.resize( scalarDofs );
        dphiNb_.resize( scalarDofs );
      }
  
 
  void elementTensor ( size_t pt,
                       const GeometryType& geometry,
                       const QuadratureType& quadrature,
                       const BasisFunctionSetType& baseSet,
                       LocalMatrixType& jLocal) const;
  template<class Quad>
  void intersectionTensor ( size_t pt,
                            const IntersectionType& intersection,
                            const GeometryType& geometry,
                            const GeometryType& geometryNb,
                            const Quad& quadInside, 
                            const Quad& quadOutside,
                            const BasisFunctionSetType& baseSet,
                            const BasisFunctionSetType& baseSetNb,
                            LocalMatrixType& jLocal,
                            LocalMatrixType& jocalNbEn,
                            LocalMatrixType& jLocalEnNb,
                            LocalMatrixType& jLocalNbNb) const;
  
  void boundaryTensor( size_t pt,
                       const IntersectionType& intersection,
                       const GeometryType& geometry,
                       const FaceQuadratureType& quadInside, 
                       const BasisFunctionSetType& baseSet,
                       LocalMatrixType& jLocal) const;
  
  
  
  private:
  typedef Dune::Fem::DiagonalAndNeighborStencil<DiscreteFunctionSpaceType,DiscreteFunctionSpaceType> StencilType;
  template< class JacobianVector,class DiffusionTensor>
  void diffusionaxpy ( const size_t local_i,
                       const size_t local_j,
                       const JacobianVector& dphi,
                       const DiffusionTensor du,
                       const double weight,
                       LocalMatrixType& jLocal) const;

  
  int globalbf( int scalarbf, int component) const
  {
    return scalarbf*dimRange+component;
  }

  using BaseType::uOld_;
  using BaseType::uOldLocal_;
  using BaseType::uOldNeighbor_;
  using BaseType::model_; 
  using BaseType::time_;
  using BaseType::deltaTInv_; 
  using BaseType::lastSpeed_;
  using BaseType::factorImp_;
  using BaseType::factorExp_;
  using BaseType::areaEn_;
  using BaseType::areaNb_;
  using BaseType::outflow_;

  private:
  const JacobianFluxType jacFlux_;
  mutable RangeVectorType uEn_,uOldEn_,uNb_,uOldNb_;
  mutable JacobianVectorType duEn_,duOldEn_,duNb_,duOldNb_;
  mutable BasefunctionStorage phi_;
  mutable BaseJacobianStorage dphi_;
  mutable BasefunctionStorage phiNb_;
  mutable BaseJacobianStorage dphiNb_;
 
  mutable MatrixHelper::Couplings<dimDomain> couplings_;

};


// Implementation of NskMixedTensor
// // ------------------------------------------------

// evaluate the coefficient functions in all quadrature points
template<class Operator, class Model, class Flux,class JacFlux>
template<class Quad>
void NskMixedTensor<Operator , Model, Flux,JacFlux >
::computeCoefficientsEn( const LocalFunctionType& uLocal,
                       const Quad& quadrature ) const
{
  const size_t numQuadraturePoints = quadrature.nop();

  uEn_.resize(numQuadraturePoints);
  uLocal.evaluateQuadrature( quadrature, uEn_);

  duEn_.resize(numQuadraturePoints);
  uLocal.evaluateQuadrature( quadrature, duEn_);

  uOldEn_.resize(numQuadraturePoints);
  uOldLocal_.evaluateQuadrature( quadrature, uOldEn_);

  duOldEn_.resize(numQuadraturePoints);
  uOldLocal_.evaluateQuadrature( quadrature, duOldEn_);
}

// evaluate the coefficient functions in all quadrature points
template<class Operator, class Model, class Flux,class JacFlux>
template<class Quad>
void NskMixedTensor<Operator , Model, Flux,JacFlux >
::computeCoefficientsNb( const LocalFunctionType& uLocal,
                         const Quad& quadrature ) const
{
  const size_t numQuadraturePoints = quadrature.nop();

  uNb_.resize(numQuadraturePoints);
  uLocal.evaluateQuadrature( quadrature, uNb_);

  duNb_.resize(numQuadraturePoints); 
  uLocal.evaluateQuadrature( quadrature, duNb_);

  uOldNb_.resize(numQuadraturePoints);
  uOldNeighbor_.evaluateQuadrature( quadrature, uOldNb_);

  duOldNb_.resize(numQuadraturePoints); 
  uOldNeighbor_.evaluateQuadrature( quadrature, duOldNb_);
}


template<class DiscreteFunction,  class Model, class Flux,class JacFlux>
void NskMixedTensor<DiscreteFunction , Model, Flux , JacFlux >
::elementTensor ( size_t pt,
                  const GeometryType& geometry,
                  const QuadratureType& quadrature,
                  const BasisFunctionSetType& baseSet,
                  LocalMatrixType& jLocal ) const


{
  const typename QuadratureType::CoordinateType &x = quadrature.point( pt );
  const double weight = quadrature.weight( pt )* geometry.integrationElement( x );
  const int numScalarBf=baseSet.size()/dimRange;
 

  MatrixHelper::evaluateScalarAll( quadrature[ pt ], baseSet.shapeFunctionSet(),phi_);
  MatrixHelper::jacobianScalarAll( quadrature[ pt ], geometry, baseSet.shapeFunctionSet(),dphi_);      
    

  RangeType vu(0.) ,vuOld(0.), vuMid(0.) ,fu(0.);
  JacobianRangeType dvu(0.) ,duOld(0), duMid(0.), fdu(0.);
  TensorRangeType flux(0.);
  vuOld=uOldEn_[pt];
  vu=uEn_[ pt ];
  //(1+theta)/2*U^n+(1-theta)/2*U^(n-1)
  vuMid.axpy(factorImp_,uEn_[pt]);
  vuMid.axpy(factorExp_,uOldEn_[pt]);

  //(1+theta)/2*DU^n+(1-theta)/2*DU^(n-1)
  duMid.axpy(factorImp_,duEn_[pt]);
  duMid.axpy(factorExp_,duOldEn_[pt]);

  //mu
  RangeFieldType dphimu, drhomu, dphitau, drhotau;
  model_.drhomuSource( vu[ 0 ] , vuOld[ 0 ]  , drhomu );


  JacobianRangeType diffusionU;
  model_.diffusionprime(vuMid,duMid,diffusionU);

  for( size_t jj=0 ;  jj < numScalarBf ; ++jj)
    {
      flux*=0.;
      //(rho,rho)
      flux[0][0]=deltaTInv_*phi_[ jj ][0];
      //(v2,rho) gravity
      DomainType xgl(0.);
      RangeType baseFunc(0.);
      baseFunc[0]=phi_[jj][0];

      for( size_t kk = 0 ; kk < dimDomain ; ++kk )
        {
         //div u*phi_rho + v\cdot\nabla phi_rho
          flux[0][0]+=factorImp_*(duMid[ 1+kk ][ kk ]*phi_[ jj ][ 0 ]+vuMid[ 1+kk ]*dphi_[ jj ][ 0 ][ kk ]);

          // rho*div phi_v +\nabla\rho\cdot\phi_v
          flux[ 0 ][ 1+kk ]=factorImp_*(vuMid[0]*dphi_[ jj ][0 ][kk]+duMid[0][ kk ]*phi_[ jj ][0]);
      
          flux[1+kk][0] = (vu[ 1 + kk ]-vuOld[ 1 + kk ])*phi_[ jj ][0]*factorImp_*deltaTInv_;
          flux[1+kk][1+kk]=phi_[ jj ][ 0 ]*vuMid[0]*deltaTInv_;
   
          for( size_t ll = 0; ll<dimDomain; ++ll)
            { 
              flux[ 1+kk ][ 0 ]+=factorImp_*(duMid[ 1 + kk ][ ll ] - duMid[ 1 + ll ][ kk ])*vuMid[ 1 + ll ]*phi_[ jj ][ 0 ];
              flux[ 1+kk ][ 1+ll ]+=factorImp_*(duMid[ 1 + kk ][ ll ] - duMid[ 1 + ll ][ kk ])* phi_[ jj ][ 0 ]*vuMid[0];
              flux[ 1+kk ][ 1+kk ]+=factorImp_*( dphi_[ jj ][ 0 ][ ll ] )*vuMid[ 1 + ll ]*vuMid[ 0 ];
              flux[ 1+kk ][ 1+ll ]-=factorImp_*( dphi_[ jj ][ 0 ][ kk ] )*vuMid[ 1 + ll ]*vuMid[ 0 ];
              
            } 
              
          //(v,rho,mu*)
          flux[ 1+kk ][ 0 ]+=factorImp_*phi_[ jj ][ 0 ]*duMid[ 2 + dimDomain ][ kk ];
          //(v,rho*,mu)
          flux[ 1+kk ][ dimDomain+1]+=factorImp_*vuMid[ 0 ]*dphi_[ jj ][ 0 ][ kk ];

          //(mu,v_kk)
          flux[ dimDomain+1 ][ 1+kk]-=2*vu[ 1+kk ]*phi_[ jj ][ 0 ]*0.25;


          flux[ dimDomain+2+kk ][ 0 ]-=dphi_[ jj ][ 0 ][ kk ];
          //(sigma,sigma)
          flux[ dimDomain+2+kk ][ 0 ]=phi_[ jj ][ 0 ];

          //(muu, gradsigma)
          flux[ dimDomain+1 ][ dimDomain+2+kk ]+=model_.delta()*factorImp_*dphi_[ jj ][ 0 ][ kk ];
        }


          //<mu,rho>=(dmu/drho,rho))
          flux[ dimDomain+1 ][ 0 ]=-drhomu*phi_[ jj ][ 0 ];
          
          //(mu,mu)
          flux[ dimDomain+1 ][ dimDomain+1 ]=factorImp_*phi_[ jj ][ 0 ];
        

          DiffusionType du;

          for(int i = 0;  i<dimDomain ; ++i)
            du[i]=0;
 
          model_.scalar2vectorialDiffusion( vuMid , dphi_[jj] , du );
          
          diffusionU*=phi_[jj][0]*factorImp_;
          for( size_t ii = 0 ; ii < numScalarBf ; ++ii )
            {
              diffusionaxpy( ii ,jj , dphi_, du ,weight, jLocal);

              MatrixHelper::diffPhiAxpy< DofAlignmentType >( ii,
                                                             jj,
                                                             dimDomain,
                                                             dphi_,
                                                             diffusionU,
                                                             weight,
                                                             jLocal);
              MatrixHelper::axpyElement< DofAlignmentType >( couplings_,
                                                             phi_,
                                                             flux,
                                                             ii,
                                                             jj,
                                                             weight,
                                                             jLocal);
          } 
      }
} 






template<class DiscreteFunction,class Model, class Flux, class Jacobian> 
template< class JacobianVector, class DiffusionTensor>
void NskMixedTensor< DiscreteFunction, Model, Flux,  Jacobian>
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
      size_t global_i= DofAlignmentType::vectorialIndex( 1 + ii , local_i); 
      for( size_t jj = 0 ; jj < dimDomain; ++jj )
        {
          size_t global_j = DofAlignmentType::vectorialIndex( 1 + jj , local_j); 
          double  value(0);
          value=du[ jj ][ ii ]*dphi[ local_i ][ 0 ];         
                
          jLocal.add( global_i , global_j ,factorImp_*value*weight );
        }
     }
}

template<class DiscreteFunction,  class Model, class Flux,class JacFlux>
template<class Quad>
void NskMixedTensor<DiscreteFunction, Model, Flux,JacFlux >
::intersectionTensor ( size_t pt,
                       const IntersectionType& intersection,
                       const GeometryType& geometry,
                       const GeometryType& geometryNb,
                       const Quad& quadInside, 
                       const Quad& quadOutside,
                       const BasisFunctionSetType& baseSet,
                       const BasisFunctionSetType& baseSetNb,
                       LocalMatrixType& jLocal,
                       LocalMatrixType& jLocalNbEn,
                       LocalMatrixType& jLocalEnNb,
                       LocalMatrixType& jLocalNbNb) const

{  
  const unsigned int numScalarBf=baseSet.size()/dimRange;

    
  RangeType    vuMidEn(0.), vuMidNb(0.);
  JacobianRangeType aduLeft(0.),aduRight(0.),duMidNb(0.), duMidEn(0.);
  TensorRangeType fluxLeft(0.), fluxRight(0.);
  TensorRangeType fluxLeftNeg(0.), fluxRightNeg(0.);
                    
  const double weightInside=quadInside.weight( pt ); 
  const double weightOutside=quadOutside.weight( pt ); 

  MatrixHelper::evaluateScalarAll( quadInside[ pt ], baseSet.shapeFunctionSet(), phi_ );
  MatrixHelper::jacobianScalarAll( quadInside[ pt ], geometry, baseSet.shapeFunctionSet(), dphi_);

  MatrixHelper::evaluateScalarAll( quadOutside[ pt ], baseSetNb.shapeFunctionSet(), phiNb_);
  MatrixHelper::jacobianScalarAll( quadOutside[ pt ], geometryNb, baseSetNb.shapeFunctionSet(), dphiNb_);

  vuMidEn.axpy( factorImp_ , uEn_[ pt ] );
  vuMidEn.axpy( factorExp_ , uOldEn_[ pt ] );

  vuMidNb.axpy( factorImp_ , uNb_[ pt ] );
  vuMidNb.axpy( factorExp_ , uOldNb_[ pt ]);

  const typename FaceQuadratureType::LocalCoordinateType &x = quadInside.localPoint( pt );

  const DomainType normal = intersection.integrationOuterNormal( x );

  // compute penalty factor
  const double intersectionArea = normal.two_norm();
  DomainType unitnormal=normal;
  unitnormal/=intersectionArea;
  const double localMinArea=std::min(areaNb_,areaEn_);
  const double penaltyFactor = intersectionArea/std::min(areaEn_,areaNb_);
  const double localMaxSpeed =  std::max( model_.maxSpeed( unitnormal , uOldEn_[pt] ), model_.maxSpeed( unitnormal , uOldNb_[pt] ) );
  
  double localwidth = 0.5*localMinArea*lastSpeed_;
 
  jacFlux_.numericalFlux( normal,
                          localwidth,
                          penaltyFactor,
                          vuMidEn,
                          vuMidNb,
                          vuMidEn,
                          vuMidNb,
                          fluxLeft,
                          fluxRight);
  
  DomainType negnormal=normal;
  negnormal*=-1;
  
  jacFlux_.numericalFlux( negnormal,
                          localwidth,
                          penaltyFactor,
                          vuMidNb,
                          vuMidEn,
                          vuMidNb,
                          vuMidEn,
                          fluxLeftNeg,
                          fluxRightNeg);

  for( size_t jj=0 ; jj < numScalarBf ; ++jj)
    {
      RangeType avuLeft(0.), avuRight(0.), valueLeft(0.),valueRight(0.);
      DiffusionType aLeft,aRight;
      JacobianRangeType aphiEn,aphiNb;
      DiffusionValueType bLeft, bRight;
      RangeType bphiEn,bphiNb;

      jacFlux_.scalar2vectorialDiffusionFlux( normal,
                                              penaltyFactor,
                                              vuMidEn,
                                              vuMidNb,
                                              phi_[ jj ],
                                              phiNb_[ jj ],
                                              dphi_[ jj ],
                                              dphiNb_[ jj ],
                                              aLeft,
                                              aRight,
                                              bLeft,
                                              bRight );

      jacFlux_.diffPhiDiffusionFlux( normal,
                                     penaltyFactor,
                                     vuMidEn,
                                     vuMidNb,
                                     duMidEn,
                                     duMidNb,
                                     phi_[ jj ],
                                     phiNb_[ jj ],
                                     aphiEn,
                                     aphiNb,
                                     bphiEn,
                                     bphiNb );
 
      for( size_t ii = 0; ii < numScalarBf ; ++ii )
        {

          MatrixHelper::axpyIntersection< DofAlignmentType >( couplings_,
                                                              phi_,
                                                              phi_,
                                                              fluxLeft,
                                                              ii,
                                                              jj,
                                                              weightInside,
                                                              jLocal );

          MatrixHelper::axpyIntersection< DofAlignmentType >( couplings_,
                                                              phi_,
                                                              phiNb_,
                                                              fluxRight,
                                                              ii,
                                                              jj,
                                                              weightInside,
                                                              jLocalNbEn );

          MatrixHelper::axpyIntersection< DofAlignmentType >( couplings_,
                                                              phiNb_,
                                                              phi_,
                                                              fluxRightNeg,
                                                              ii,
                                                              jj,
                                                              weightOutside,
                                                              jLocalEnNb );

          MatrixHelper::axpyIntersection< DofAlignmentType >( couplings_,
                                                              phiNb_,
                                                              phiNb_,
                                                              fluxLeftNeg,
                                                              ii,
                                                              jj,
                                                              weightOutside,
                                                              jLocalNbNb );

            //column indexfor \phi
            int global_j_phi=jj*dimRange+(dimDomain+1);
                 // Adding DiffusionFluxTerms
            for(int i  = 0; i  < dimDomain ; ++i )
            {
                int global_i=ii*dimRange+1+i; 
                double valueEn,valueNb;
 
                valueEn=aphiEn[ i+1 ]*dphi_[ ii ][ 0 ];
                valueEn+=bphiEn[ i+1 ]*phi_[ ii ];

                valueNb=-1*(aphiNb[ i+1 ]*dphi_[ ii ][ 0 ]);
                valueNb+=bphiNb[ i+1 ]*phi_[ ii ];

                jLocal.add( global_i , global_j_phi , weightInside*valueEn*factorImp_);
                jLocalNbEn.add( global_i , global_j_phi , weightInside*valueNb*factorImp_);

                valueEn=(aphiEn[ i+1 ]*dphiNb_[ ii ][ 0 ]);
                valueEn-=bphiEn[ i+1 ]*phiNb_[ ii ];
                valueNb=-1*(aphiNb[ i+1 ]*dphiNb_[ ii][ 0 ]);
                valueNb-=bphiNb[ i+1 ]*phiNb_[ ii ];

                jLocalEnNb.add( global_i , global_j_phi , weightOutside*valueEn*factorImp_);
                jLocalNbNb.add( global_i , global_j_phi , weightOutside*valueNb*factorImp_);


              
              
              for(int j  = 0 ; j  < dimDomain ; ++j )
              {
                int global_j= jj*dimRange+1+j;
                valueEn=0;
                valueNb=0; 

                valueEn=aLeft[ j ][ i ]*dphi_[ ii ][ 0 ];
                valueEn+=bLeft[ j ][ i ]*phi_[ ii ];
                valueNb=-1*(aRight[ j ][ i ]*dphi_[ ii ][ 0 ]);
                valueNb+=bRight[ j ][ i ]*phi_[ ii ];


                jLocal.add( global_i , global_j , weightInside*valueEn*factorImp_); 
                jLocalNbEn.add( global_i , global_j , weightInside*valueNb*factorImp_);

                valueEn=(aLeft[ j ][ i ]*dphiNb_[ ii ][ 0 ]);
                valueEn-=bLeft[ j ][ i ]*phiNb_[ ii ];
                valueNb=-1*(aRight[ j ][ i ]*dphiNb_[ ii][ 0 ]);
                valueNb-=bRight[ j ][ i ]*phiNb_[ ii ];

                jLocalEnNb.add( global_i , global_j , weightOutside*valueEn*factorImp_); 
                jLocalNbNb.add( global_i , global_j , weightOutside*valueNb*factorImp_);

              
              }
            }
          } 
        }
}    
template<class DiscreteFunction, class Model, class Flux,class JacFlux>
void NskMixedTensor<DiscreteFunction ,  Model, Flux,JacFlux >
::boundaryTensor ( size_t pt,
                   const IntersectionType& intersection,
                   const GeometryType& geometry,
                   const FaceQuadratureType& quadInside, 
                   const BasisFunctionSetType& baseSet,
                   LocalMatrixType& jLocal) const
{
  const typename FaceQuadratureType::LocalCoordinateType &x = quadInside.localPoint( pt );
  const DomainType normal = intersection.integrationOuterNormal( x );
  size_t boundaryIndex=intersection.boundaryId();

  // compute penalty factor
  const double intersectionArea = normal.two_norm();
  const double penaltyFactor = intersectionArea /  areaEn_;
  const double localMaxSpeed = model_.maxSpeed( normal , uOldEn_[ pt ] );
  const double localwidth = areaEn_*lastSpeed_; 


  const double weightInside=quadInside.weight( pt );
  const int numScalarBf=baseSet.size()/dimRange;
  MatrixHelper::evaluateScalarAll( quadInside[ pt ], baseSet.shapeFunctionSet(), phi_ );
  MatrixHelper::jacobianScalarAll( quadInside[ pt ], geometry, baseSet.shapeFunctionSet(), dphi_);
  RangeType vuMidEn(0.);
  JacobianRangeType duMidEn(0.);
  vuMidEn.axpy( factorImp_ , uEn_[ pt ] );
  vuMidEn.axpy( factorExp_ , uOldEn_[ pt ] );
  duMidEn.axpy( factorImp_ , duEn_[ pt ] );
  duMidEn.axpy( factorExp_ , duOldEn_[ pt ] );



  TensorRangeType fluxLeft(0.);
  
  if( boundaryIndex==1 || !outflow_)
    {
      jacFlux_.boundaryFlux( normal,
                             localwidth,
                             vuMidEn,
                             fluxLeft);
    }
  else
    {
      jacFlux_.outFlowFlux( normal,
                            localwidth,
                            vuMidEn,
                            fluxLeft);

    }

  for( size_t jj=0 ; jj < numScalarBf ; ++jj)
    {
      RangeType bphiEn(0.), avuRight(0.),dummy(0.), valueLeft(0.),valueRight(0.);
      JacobianRangeType aphiEn(0.),aduRight(0.);
      DiffusionType aLeft;
      DiffusionValueType bLeft;
 
      if( boundaryIndex==1 || !outflow_)
        {
          jacFlux_.scalar2vectorialBoundaryFlux(normal,
                                                penaltyFactor,
                                                vuMidEn,
                                                phi_[ jj ],
                                                dphi_[ jj ],
                                                aLeft,
                                                bLeft);

          jacFlux_.diffPhiBoundaryFlux( normal,
                                        penaltyFactor,
                                        vuMidEn,
                                        duMidEn,
                                        phi_[ jj ],
                                        aphiEn,
                                        bphiEn);

          
        }


      for( size_t ii = 0; ii < numScalarBf ; ++ii )     
        {
          MatrixHelper::axpyIntersection< DofAlignmentType >( couplings_,
                                                              phi_,
                                                              phi_,
                                                              fluxLeft,
                                                              ii,
                                                              jj,
                                                              weightInside,
                                                              jLocal );
          if( boundaryIndex==1 || !outflow_)
            {
              //column indexfor \phi
              int global_j_phi=jj*dimRange+(dimDomain+1);
        
              for(int i  = 0; i  < dimDomain ; ++i )
                {
                  int global_i=ii*dimRange+1+i; 
                  double valueEn;
 
                  valueEn=aphiEn[ i+1 ]*dphi_[ ii ][ 0 ];
                  valueEn+=bphiEn[ i+1 ]*phi_[ ii ];

                  jLocal.add( global_i , global_j_phi , weightInside*valueEn*factorImp_);
                  
                  for(int j  = 0 ; j  < dimDomain ; ++j )
                    {
                      int global_j= jj*dimRange+1+j;
                      double valueEn;
                     
                      valueEn=aLeft[ j ][ i ]*dphi_[ ii ][ 0 ];
                      valueEn+=bLeft[ j ][ i ]*phi_[ ii ];

                      jLocal.add( global_i , global_j , weightInside*valueEn*factorImp_);
                    }
                }
            }  
          }
        
    }
}

 










#endif //PHASEFIELDJACOBIANOPERATOR_HH


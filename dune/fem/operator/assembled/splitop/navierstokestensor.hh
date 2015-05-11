#ifndef PHASEFIELD_NAVIERSTOKES_TENSOR_HH
#define PHASEFIELD_NAVIERSTOKES_TENSOR_HH

// Dune::Fem includes
#include <dune/fem/function/localfunction/temporary.hh>
#include <dune/fem/operator/common/differentiableoperator.hh>
#include <dune/fem/operator/common/stencil.hh>
#include <dune/fem/operator/common/temporarylocalmatrix.hh>

// local includes
#include "../matrixhelper.hh"

template<class DiscreteFunction, class AddFunction, class Model, class Flux,class JacobianFlux>
class PhasefieldNavierStokesTensor:
public PhasefieldNavierStokesIntegrator<DiscreteFunction,AddFunction,Model,Flux>
{

  typedef PhasefieldNavierStokesIntegrator<DiscreteFunction,AddFunction,Model,Flux> BaseType;
  
  public:
  using DiscreteFunctionType      = typename BaseType::DiscreteFunctionType;
  using AddFunctionType           = typename BaseType::AddFunctionType;
  using ModelType                 = typename BaseType::ModelType;
  using DiscreteFunctionSpaceType = typename BaseType::DiscreteFunctionSpaceType;
  using BasisFunctionSetType      = typename DiscreteFunctionSpaceType::BasisFunctionSetType;  
  using RangeFieldType            = typename BaseType::RangeFieldType;
  using LocalFunctionType         = typename BaseType::LocalFunctionType;
  using RangeType                 = typename BaseType::RangeType;
  using JacobianRangeType         = typename BaseType::JacobianRangeType;
  static const int dimDomain      = BaseType::dimDomain;
  static const int dimRange       = BaseType::dimRange;
 
  using AddLocalFunctionType      = typename BaseType::AddLocalFunctionType;
  using AddRangeType              = typename BaseType::AddRangeType;
  using AddJacobianRangeType      = typename BaseType::AddJacobianRangeType;
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
   PhasefieldNavierStokesTensor(const ModelType &model,
                                const DiscreteFunctionSpaceType& space ):
    BaseType(model,space),
    jacFlux_(model),
    imexFactor_(Dune::Fem::Parameter::getValue<double>("phasefield.IMEX"))
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
   
  using BaseType::addFunction_;
  using BaseType::uOldLocal_;
  using BaseType::addLocal_;
  using BaseType::addNeighbor_;
  using BaseType::model_; 
  using BaseType::deltaTInv_; 
  using BaseType::lastSpeed_;
  using BaseType::areaEn_;
  using BaseType::areaNb_;
  using BaseType::outflow_;

  private:
  const JacobianFluxType jacFlux_;
  const double imexFactor_;
  mutable RangeVectorType uEn_,uOldEn_,uNb_,addEn_,addNb_;
  mutable JacobianVectorType duEn_,duNb_,duAddEn_,duAddNb_;
  mutable BasefunctionStorage phi_;
  mutable BaseJacobianStorage dphi_;
  mutable BasefunctionStorage phiNb_;
  mutable BaseJacobianStorage dphiNb_;
 
  mutable MatrixHelper::NvStCouplings<dimDomain> couplings_;
};


// Implementation of LocalFDOperator
// // ------------------------------------------------

// evaluate the coefficient functions in all quadrature points
template<class Operator, class AddFunction, class Model, class Flux,class JacFlux>
template<class Quad>
void PhasefieldNavierStokesTensor<Operator , AddFunction, Model, Flux,JacFlux >
::computeCoefficientsEn( const LocalFunctionType& uLocal,
                       const Quad& quadrature ) const
{
  const size_t numQuadraturePoints = quadrature.nop();

  uEn_.resize(numQuadraturePoints);
  uLocal.evaluateQuadrature( quadrature, uEn_);

  uOldEn_.resize(numQuadraturePoints);
  uOldLocal_.evaluateQuadrature( quadrature, uOldEn_);

  duEn_.resize(numQuadraturePoints);
  uLocal.evaluateQuadrature( quadrature, duEn_);

  addEn_.resize(numQuadraturePoints);
  addLocal_.evaluateQuadrature( quadrature, addEn_);

  duAddEn_.resize(numQuadraturePoints);
  addLocal_.evaluateQuadrature( quadrature, duAddEn_);
}

// evaluate the coefficient functions in all quadrature points
template<class Operator, class AddFunction, class Model, class Flux,class JacFlux>
template<class Quad>
void PhasefieldNavierStokesTensor<Operator , AddFunction, Model, Flux,JacFlux >
::computeCoefficientsNb( const LocalFunctionType& uLocal,
                         const Quad& quadrature ) const
{
  const size_t numQuadraturePoints = quadrature.nop();

  uNb_.resize(numQuadraturePoints);
  uLocal.evaluateQuadrature( quadrature, uNb_);

  duNb_.resize(numQuadraturePoints);
  uLocal.evaluateQuadrature( quadrature, duNb_);

  addNb_.resize(numQuadraturePoints);
  addNeighbor_.evaluateQuadrature( quadrature, addNb_);

  duAddNb_.resize(numQuadraturePoints);
  addNeighbor_.evaluateQuadrature( quadrature, duAddNb_);
}




template<class Operator, class AddFunction, class Model, class Flux,class JacFlux>
void PhasefieldNavierStokesTensor<Operator , AddFunction, Model, Flux,JacFlux >
::elementTensor( size_t pt,
                  const GeometryType& geometry,
                  const QuadratureType& quadrature,
                  const BasisFunctionSetType& baseSet,
                  LocalMatrixType& jLocal ) const

{    
  const typename QuadratureType::CoordinateType &x = quadrature.point( pt );
  const double weight = quadrature.weight( pt )* geometry.integrationElement( x );
  //get this before application
  const int numScalarBf=baseSet.size()/dimRange;
   
  //this should be provided by the actuial basisfunctionset
  MatrixHelper::evaluateScalarAll( quadrature[ pt ], baseSet.shapeFunctionSet(),phi_);
  MatrixHelper::jacobianScalarAll( quadrature[ pt ], geometry, baseSet.shapeFunctionSet(),dphi_);      
    
  RangeType vu(0.),vuOld(0.);
  AddRangeType uAdd(0);
  JacobianRangeType dvu(0.);
  AddJacobianRangeType duAdd(0.);
      
  vu=uEn_[ pt ];
  vuOld=uOldEn_[pt];
  dvu=duEn_[ pt ]; 
  uAdd=addEn_[ pt ];
  duAdd=duAddEn_[ pt ];
  
  TensorRangeType flux;
  for( size_t jj=0 ;  jj < numScalarBf ; ++jj)
    {
      flux=0.;
      
      //(rho,rho)
      flux[0][0]=deltaTInv_*phi_[ jj ][0];
          
      for( size_t kk = 0 ; kk < dimDomain ; ++kk )
        {
          //div u*phi_rho + v\cdot\nabla phi_rho 
          //(rho,rho)
          flux[0][0]+=(dvu[ 1+kk ][ kk ]*phi_[ jj ][ 0 ]+vu[ 1+kk ]*dphi_[ jj ][ 0 ][ kk ]);
          // rho*div phi_v +\nabla\rho\cdot\phi_v
          //(rho,v)
          flux[ 0 ][ 1+kk ]=(vu[0]*dphi_[ jj ][0 ][kk]+dvu[0][ kk ]*phi_[ jj ][0]);
      
          flux[1+kk][0] = (vu[ 1 + kk ]-vuOld[ 1 + kk ])*phi_[ jj ][0]*deltaTInv_;
          flux[1+kk][1+kk]=phi_[ jj ][ 0 ]*vu[0]*deltaTInv_;
   
          for( size_t ll = 0; ll<dimDomain; ++ll)
            { 
              flux[ 1+kk ][ 0 ]+=(dvu[ 1 + kk ][ kk ] - dvu[ 1 + ll ][ kk ])*vu[ 1 + ll ]*phi_[ jj ][ 0 ];
              flux[ 1+kk ][ 1+ll ]+=(dvu[ 1 + kk ][ ll ] - dvu[ 1 + ll ][ kk ])* phi_[ jj ][ 0 ]*vu[0];
              flux[ 1+kk ][ 1+kk ]+=( dphi_[ jj ][ 0 ][ ll ] )*vu[ 1 + ll ]*vu[ 0 ];
              flux[ 1+kk ][ 1+ll ]-=( dphi_[ jj ][ 0 ][ kk ] )*vu[ 1 + ll ]*vu[ 0 ];
            } 
              
          //(v,rho,mu*)
          flux[ 1+kk ][ 0 ]+=phi_[ jj ][ 0 ]*dvu[ 1 + dimDomain ][ kk ];
          //(v,rho*,mu)
          flux[ 1+kk ][ 1+dimDomain ]+=imexFactor_*vu[ 0 ]*dphi_[ jj ][ 0 ][ kk ];
          //(mu,v_kk)
          flux[ 1+dimDomain ][ 1+kk ]-=0.5*vu[ 1+kk ]*phi_[ jj ][ 0 ];

        }


          //(mu,mu)
          flux[ 1+dimDomain ][ 1+dimDomain ]=imexFactor_*phi_[ jj ][ 0 ];
        

          DiffusionType du;

          for(int i = 0;  i<dimDomain ; ++i)
            du[i]=0;
 
          model_.scalar2vectorialDiffusion( dphi_[jj], du);
          
          for( size_t ii = 0 ; ii < numScalarBf ; ++ii )
            {
              MatrixHelper::diffusionaxpy< DofAlignmentType >(ii,
                                                              jj,
                                                              dimDomain,
                                                              dphi_,
                                                              du,
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

template<class Operator, class AddFunction, class Model, class Flux,class JacFlux>
template<class Quad>
void PhasefieldNavierStokesTensor<Operator , AddFunction, Model, Flux,JacFlux >
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
  const typename FaceQuadratureType::LocalCoordinateType &x = quadInside.localPoint( pt );
  const DomainType normal = intersection.integrationOuterNormal( x );

  // compute penalty factor
  const double intersectionArea = normal.two_norm();
  const double penaltyFactor = intersectionArea/std::min(areaEn_,areaNb_);
  const double localwidth = lastSpeed_*std::min(areaEn_,areaNb_)/intersectionArea;

  const double weightInside=quadInside.weight( pt );
  const double weightOutside=quadOutside.weight( pt );
  const int numScalarBf=baseSet.size()/dimRange;
  MatrixHelper::evaluateScalarAll( quadInside[ pt ], baseSet.shapeFunctionSet(), phi_ );
  MatrixHelper::jacobianScalarAll( quadInside[ pt ], geometry, baseSet.shapeFunctionSet(), dphi_);

  MatrixHelper::evaluateScalarAll( quadOutside[ pt ], baseSetNb.shapeFunctionSet(), phiNb_);
  MatrixHelper::jacobianScalarAll( quadOutside[ pt ], geometryNb, baseSetNb.shapeFunctionSet(), dphiNb_);


  TensorRangeType fluxLeft(0.), fluxRight(0.);
  TensorRangeType fluxLeftNeg(0.), fluxRightNeg(0.);

  jacFlux_.numericalFlux( normal,
                          localwidth,
                          penaltyFactor,
                          uEn_[pt],
                          uNb_[pt],
                          addEn_[pt],
                          addNb_[pt],
                          fluxLeft,
                          fluxRight);
  DomainType negnormal=normal;

  negnormal*=-1;

  jacFlux_.numericalFlux( negnormal,
                          localwidth,
                          penaltyFactor,
                          uNb_[pt],
                          uEn_[pt],
                          addNb_[pt],
                          addEn_[pt],
                          fluxLeftNeg,
                          fluxRightNeg);

      for( size_t jj=0 ; jj < numScalarBf ; ++jj)
        {
        
          RangeType avuLeft(0.), avuRight(0.), valueLeft(0.),valueRight(0.);
          DiffusionType aLeft,aRight;
          DiffusionValueType bLeft, bRight; 

          jacFlux_.scalar2vectorialDiffusionFlux( normal,
                                                  penaltyFactor,
                                                  phi_[ jj ],
                                                  phiNb_[ jj ],
                                                  dphi_[ jj ],
                                                  dphiNb_[ jj ],
                                                  aLeft,
                                                  aRight,
                                                  bLeft,
                                                  bRight );


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
#if 1
// Adding DiffusionFluxTerms
            for(int i  = 0; i  < dimDomain ; ++i )
            {
              int global_i=ii*dimRange+1+i; 
              for(int j  = 0 ; j  < dimDomain ; ++j )
              {
                int global_j= jj*dimRange+1+j;
                double valueEn,valueNb; 

                valueEn=aLeft[ j ][ i ]*dphi_[ ii ][ 0 ];
                valueEn+=bLeft[ j ][ i ]*phi_[ ii ];
                valueNb=-1*(aRight[ j ][ i ]*dphi_[ ii ][ 0 ]);
                valueNb+=bRight[ j ][ i ]*phi_[ ii ];

                jLocal.add( global_i , global_j , weightInside*valueEn); 
                jLocalNbEn.add( global_i , global_j , weightInside*valueNb);

                valueEn=(aLeft[ j ][ i ]*dphiNb_[ ii ][ 0 ]);
                valueEn-=bLeft[ j ][ i ]*phiNb_[ ii ];
                valueNb=-1*(aRight[ j ][ i ]*dphiNb_[ ii][ 0 ]);
                valueNb-=bRight[ j ][ i ]*phiNb_[ ii ];

                jLocalEnNb.add( global_i , global_j , weightOutside*valueEn); 
                jLocalNbNb.add( global_i , global_j , weightOutside*valueNb);

              }
            }
#endif
 

              }

          }
}
                                     
template<class Operator, class AddFunction, class Model, class Flux,class JacFlux>
void PhasefieldNavierStokesTensor<Operator , AddFunction, Model, Flux,JacFlux >
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

  const double weightInside=quadInside.weight( pt );
  const int numScalarBf=baseSet.size()/dimRange;
  MatrixHelper::evaluateScalarAll( quadInside[ pt ], baseSet.shapeFunctionSet(), phi_ );
  MatrixHelper::jacobianScalarAll( quadInside[ pt ], geometry, baseSet.shapeFunctionSet(), dphi_);


  TensorRangeType fluxLeft(0.);

  if( boundaryIndex==1 || !outflow_)
    {
      jacFlux_.boundaryFlux( normal,
                             penaltyFactor,
                             uEn_[pt],
                             fluxLeft);
    }
  else
    {
      jacFlux_.outFlowFlux( normal,
                            penaltyFactor,
                            uEn_[pt],
                            fluxLeft);
    }

  for( size_t jj=0 ; jj < numScalarBf ; ++jj)
    {
      RangeType avuLeft(0.), avuRight(0.),dummy(0.), valueLeft(0.),valueRight(0.);
      JacobianRangeType aduLeft(0.),aduRight(0.);
      DiffusionType aLeft;
      DiffusionValueType bLeft;

      if( boundaryIndex==1 || !outflow_)
        {
          jacFlux_.scalar2vectorialBoundaryFlux(normal,
                                                penaltyFactor,
                                                phi_[ jj ],
                                                dphi_[ jj ],
                                                aLeft,
                                                bLeft);
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
            for(int i  = 0; i  < dimDomain ; ++i )
              {
                int global_i=ii*dimRange+1+i;
                for(int j  = 0 ; j  < dimDomain ; ++j )
                  {
                    int global_j= jj*dimRange+1+j;
                    double valueEn;

                    valueEn=aLeft[ j ][ i ]*dphi_[ ii ][ 0 ];
                    valueEn+=bLeft[ j ][ i ]*phi_[ ii ];

                    jLocal.add( global_i , global_j , weightInside*valueEn);
                  }
                }
            }
    }

}

 



#endif //PHASEFIELD_NAVIERSTOKES_TENSOR_HH


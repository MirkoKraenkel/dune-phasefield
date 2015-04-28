#ifndef PHASEFIELD_ALLENCAHN_TENSOR_HH
#define PHASEFIELD_ALLENCAHN_TENSOR_HH

// Dune::Fem includes
#include <dune/fem/function/localfunction/temporary.hh>
#include <dune/fem/operator/common/differentiableoperator.hh>
#include <dune/fem/operator/common/stencil.hh>
#include <dune/fem/operator/common/temporarylocalmatrix.hh>

// local includes
#include "../matrixhelper.hh"

template<class DiscreteFunction, class AddFunction, class Model, class Flux,class JacobianFlux>
class PhasefieldAllenCahnTensor:
public PhasefieldAllenCahnIntegrator<DiscreteFunction,AddFunction,Model,Flux>
{
  typedef PhasefieldAllenCahnIntegrator<DiscreteFunction,AddFunction,Model,Flux> BaseType;

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
  using IntersectionGeometryType  = typename BaseType::IntersectionGeometryType;

  using QuadratureType            = typename BaseType::QuadratureType; 
  using FaceQuadratureType        = typename BaseType::FaceQuadratureType;

  using JacobianFluxType=JacobianFlux;

  using TemporaryLocalType       = typename BaseType::TemporaryLocalType;
  
  
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


  public: 
  PhasefieldAllenCahnTensor(const ModelType &model,
                            const DiscreteFunctionSpaceType& space):
    BaseType( model,space),
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
  
  void elementTensor (size_t pt,
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
                            LocalMatrixType& jocalEnNb,
                            LocalMatrixType& jLocalNbEn,
                            LocalMatrixType& jLocalNbNb) const;
                                     
       
  void boundaryTensor ( size_t pt,
                        const IntersectionType& intersection,
                        const GeometryType& geometry,
                        const FaceQuadratureType& quadInside, 
                        const BasisFunctionSetType& baseSet,
                        LocalMatrixType& jLocal) const;
                           

  using BaseType::addFunction_;
  using BaseType::addLocal_;
  using BaseType::addNeighbor_;
  using BaseType::model_; 
  using BaseType::deltaTInv_; 
  using BaseType::lastSpeed_;
  using BaseType::areaEn_;
  using BaseType::areaNb_;
  private:
  const JacobianFluxType jacFlux_;
  const double imexFactor_;
  mutable RangeVectorType uEn_,uNb_,addEn_,addNb_;
  mutable JacobianVectorType duEn_,duNb_,duAddEn_,duAddNb_;
  mutable BasefunctionStorage phi_;
  mutable BaseJacobianStorage dphi_;
  mutable BasefunctionStorage phiNb_;
  mutable BaseJacobianStorage dphiNb_;
 
  mutable MatrixHelper::Couplings<dimDomain> couplings_;
};


// Implementation of LocalFDOperator
// // ------------------------------------------------

// evaluate the coefficient functions in all quadrature points
template<class Operator, class AddFunction, class Model, class Flux,class JacFlux>
template<class Quad>
void PhasefieldAllenCahnTensor<Operator , AddFunction, Model, Flux,JacFlux >
::computeCoefficientsEn( const LocalFunctionType& uLocal,
                         const Quad& quadrature ) const
{
  const size_t numQuadraturePoints = quadrature.nop();
  uEn_.resize(numQuadraturePoints);
  uLocal.evaluateQuadrature( quadrature, uEn_);
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
void PhasefieldAllenCahnTensor<Operator , AddFunction, Model, Flux,JacFlux >
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
void PhasefieldAllenCahnTensor<Operator , AddFunction, Model, Flux,JacFlux >
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
    
  RangeType vu(0.);
  AddRangeType uAdd(0);
  JacobianRangeType dvu(0.);
  AddJacobianRangeType duAdd(0.);
      
  vu=uEn_[ pt ];
  dvu=duEn_[ pt ]; 
  uAdd=addEn_[ pt ];
  duAdd=duAddEn_[ pt ];
  
  TensorRangeType flux;
  //mu
  RangeFieldType dphimu=0.;
  model_.dphimuSource( uAdd[ 0 ] , uAdd[ 0 ] , vu[ 0 ] , dphimu );

  for( size_t jj=0 ;  jj < numScalarBf ; ++jj)
    {
      flux=0.;
      
      for( size_t kk = 0 ; kk < dimDomain ; ++kk )
        {
          //\nabla phi_phi\cdot v
          flux[ 0 ][ 0 ]+=0.5*dphi_[ jj ][0][kk]*uAdd[ 1+ kk];
 
          //(sigma,gradphi)
          flux[ 2+kk ][ 0 ]-=dphi_[ jj ][ 0 ][ kk ];
          
          //(sigma,sigma)
          flux[ 2+kk ][ 2+kk]=phi_[ jj ][ 0 ];

          //(tau, gradsigma)
          flux[ 1 ][ 2+kk ]+=model_.delta()*0.5*dphi_[ jj ][ 0 ][ kk ];
        }

        // phi_phi/deltaT
        flux[ 0 ][ 0 ]+=phi_[jj][0]*deltaTInv_;


        //(phi_tau rho)/ rho*rho
        flux[ 0 ][ 1 ]=model_.reactionFactor()*imexFactor_*phi_[ jj ][ 0 ]*uAdd[ 0 ]/(uAdd[0]*uAdd[0]);

        flux[ 1 ][ 1 ]+=imexFactor_*phi_[ jj ][ 0 ];

       for( size_t ii = 0 ; ii < numScalarBf ; ++ii )
        {      
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
void PhasefieldAllenCahnTensor<Operator , AddFunction, Model, Flux,JacFlux >
::intersectionTensor ( size_t pt,
                       const IntersectionType& intersection,
                       const GeometryType& geometry,
                       const GeometryType& geometryNb,
                       const Quad& quadInside, 
                       const Quad& quadOutside,
                       const BasisFunctionSetType& baseSet,
                       const BasisFunctionSetType& baseSetNb,
                       LocalMatrixType& jLocal,
                       LocalMatrixType& jLocalEnNb,
                       LocalMatrixType& jLocalNbEn,
                       LocalMatrixType& jLocalNbNb) const
{
  const typename FaceQuadratureType::LocalCoordinateType &x = quadInside.localPoint( pt );

  const DomainType normal = intersection.integrationOuterNormal( x );

  // compute penalty factor
  const double intersectionArea = normal.two_norm();
  const double localwidth = lastSpeed_*std::min(areaEn_,areaNb_)/intersectionArea;
  const double penaltyFactor = 1./localwidth;

  const double weightInside=quadInside.weight( pt ); 
  const double weightOutside=quadOutside.weight( pt ); 
   const int numScalarBf=baseSet.size()/dimRange;

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


              }
 
          }
}
                                     
       
template<class Operator, class AddFunction, class Model, class Flux,class JacFlux>
void PhasefieldAllenCahnTensor<Operator , AddFunction, Model, Flux,JacFlux >
::boundaryTensor ( size_t pt,
                   const IntersectionType& intersection,
                   const GeometryType& geometry,
                   const FaceQuadratureType& quadInside, 
                   const BasisFunctionSetType& baseSet,
                   LocalMatrixType& jLocal) const
{

}
 









#endif //PHASEFIELD_ALLENCAHN_TENSOR_HH


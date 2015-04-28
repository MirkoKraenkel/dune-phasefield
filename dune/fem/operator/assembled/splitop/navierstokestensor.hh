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
  using IntersectionGeometryType  = typename BaseType::IntersectionGeometryType;

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
                            LocalMatrixType& jocalEnNb,
                            LocalMatrixType& jLocalNbEn,
                            LocalMatrixType& jLocalNbNb) const;
  
  void boundaryTensor( size_t pt,
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
void PhasefieldNavierStokesTensor<Operator , AddFunction, Model, Flux,JacFlux >
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
    
  RangeType vu(0.);
  AddRangeType uAdd(0);
  JacobianRangeType dvu(0.);
  AddJacobianRangeType duAdd(0.);
      
  vu=uEn_[ pt ];
  dvu=duEn_[ pt ]; 
  uAdd=addEn_[ pt ];
  duAdd=duAddEn_[ pt ];
  
  TensorRangeType flux;
#if 0
  for( size_t jj=0 ;  jj < numScalarBf ; ++jj)
    {
      flux=0.;
      
      //(rho,rho)
      flux[0][0]=deltaTInv_*phi[ jj ][0];
          
      for( size_t kk = 0 ; kk < dimDomain ; ++kk )
        {
          //div u*phi_rho + v\cdot\nabla phi_rho 
          flux[0][0]+=0.5*(divu[ 1+kk ][ kk ]*phi[ jj ][ 0 ]+vu[ 1+kk ]*dphi[ jj ][ 0 ][ kk ]);
          // rho*div phi_v +\nabla\rho\cdot\phi_v
          flux[ 0 ][ 1+kk ]=0.5*(vuMid[0]*dphi[ jj ][0 ][kk]+duMid[0][ kk ]*phi[ jj ][0]);
      
          flux[1+kk][0] = (vu[ 1 + kk ]-vuOld[ 1 + kk ])*phi[ jj ][0]*0.5*deltaTInv;
          flux[1+kk][1+kk]=phi[ jj ][ 0 ]*vuMid[0]*deltaTInv_;
   
          for( size_t ll = 0; ll<dimDomain; ++ll)
            { 
              flux[ 1+kk ][ 0 ]+=0.5*(duMid[ 1 + kk ][ kk ] - duMid[ 1 + ll ][ kk ])*vuMid[ 1 + ll ]*phi[ jj ][ 0 ];
              flux[ 1+kk ][ 1+ll ]+=0.5*(duMid[ 1 + kk ][ ll ] - duMid[ 1 + ll ][ kk ])* phi[ jj ][ 0 ]*vuMid[0];
              flux[ 1+kk ][ 1+kk ]+=0.5*( dphi[ jj ][ 0 ][ ll ] )*vuMid[ 1 + ll ]*vuMid[ 0 ];
              flux[ 1+kk ][ 1+ll ]-=0.5*( dphi[ jj ][ 0 ][ kk ] )*vuMid[ 1 + ll ]*vuMid[ 0 ];
            } 
              
          //(v,rho,mu*)
          flux[ 1+kk ][ 0 ]+=0.5*phi[ jj ][ 0 ]*duImEx[ 2 + dimDomain ][ kk ];
          //(v,rho*,mu)
          flux[ 1+kk ][ dimDomain+2]+=imexFactor_*vuMid[ 0 ]*dphi[ jj ][ 0 ][ kk ];
          //(mu,v_kk)
          flux[ dimDomain+2 ][ 1+kk]-=2*vu[ 1+kk ]*phi[ jj ][ 0 ]*0.25;

           }


          //(mu,mu)
          flux[ dimDomain+2 ][ dimDomain+2 ]=imexFactor_*phi[ jj ][ 0 ];
        

          DiffusionType du;

          for(int i = 0;  i<dimDomain ; ++i)
            du[i]=0;
 
          model_.scalar2vectorialDiffusion( dphi[jj], du);  
          
          for( size_t ii = 0 ; ii < numScalarBf ; ++ii )
            {
              MatrixHelper::diffusionaxpy< DofAlignmentType >(ii,
                                                              jj,
                                                              dphi,
                                                              du,
                                                              weight,
                                                              jLocal);
          
              MatrixHelper::axpyElement< DofAlignmentType >( couplings_,
                                                             phi,
                                                             flux,
                                                             ii,
                                                             jj,
                                                             weight,
                                                             jLocal);
            } 
      }
#endif
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
                       LocalMatrixType& jocalEnNb,
                       LocalMatrixType& jLocalNbEn,
                       LocalMatrixType& jLocalNbNb) const
                         {}
                                     
template<class Operator, class AddFunction, class Model, class Flux,class JacFlux>
void PhasefieldNavierStokesTensor<Operator , AddFunction, Model, Flux,JacFlux >
::boundaryTensor ( size_t pt,
                   const IntersectionType& intersection,
                   const GeometryType& geometry,
                   const FaceQuadratureType& quadInside, 
                   const BasisFunctionSetType& baseSet,
                   LocalMatrixType& jLocal) const
                      {}
 



#endif //PHASEFIELD_NAVIERSTOKES_TENSOR_HH


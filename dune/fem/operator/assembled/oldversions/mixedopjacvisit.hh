#ifndef LOCALFD_OPERATOR_HH
#define LOCALFD_OPERATOR_HH



#include <dune/fem/function/localfunction/temporary.hh>

#include <dune/fem/operator/common/differentiableoperator.hh>
#include <dune/fem/operator/common/stencil.hh>
#include <dune/fem/operator/common/temporarylocalmatrix.hh>


#include "mixedoperator.hh"
#include "phasefieldfilter.hh"
#include "fluxes/jacobianflux_nofilter.hh"

template<class DiscreteFunction,class Model, class Flux, class Jacobian>
class PhasefieldJacobianOperator
:public Dune::Fem::DifferentiableOperator < Jacobian >,
  protected DGPhasefieldOperator<DiscreteFunction,Model,Flux>
{

  typedef DGPhasefieldOperator<DiscreteFunction,Model,Flux> MyOperatorType;

  typedef Dune::Fem::DifferentiableOperator< Jacobian> BaseType;

  enum{dimDomain=MyOperatorType::dimDomain};

  typedef typename BaseType::JacobianOperatorType JacobianOperatorType;
  typedef typename JacobianOperatorType::LocalMatrixType LocalMatrixType;
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
  
  typedef Dune::Fem::MutableArray<RangeType> RangeVectorType;
  typedef Dune::Fem::MutableArray<JacobianRangeType> JacobianVectorType;

  public: 
  PhasefieldJacobianOperator(const ModelType &model,
      const DiscreteFunctionSpaceType &space,
      int volQuadOrder=-1)
    :MyOperatorType(model,space),
    indexSet_(space_.gridPart().indexSet()),
    visited_(0),
    stencil_(space,space),
    jacFlux_(model)
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







template<class DiscreteFunction,class Model, class Flux, class Jacobian> void
PhasefieldJacobianOperator< DiscreteFunction, Model, Flux,  Jacobian>
::jacobian ( const DiscreteFunctionType &u, JacobianOperatorType &jOp ) const
{
 typedef typename DiscreteFunctionSpaceType::BasisFunctionSetType BasisFunctionSetType;

  Dune::Fem::DiagonalAndNeighborStencil<DiscreteFunctionSpaceType,DiscreteFunctionSpaceType> stencil(space(),space());
//  typedef typename JacobianOperatorType::LocalMatrixType LocalMatrixType;
//  typedef Dune::Fem::TemporaryLocalMatrix< DiscreteFunctionSpaceType, DiscreteFunctionSpaceType> LocalMatrixType;
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
    
    //LocalMatrixType jLocal( dfSpace, dfSpace);//= jOp.localMatrix( entity, entity );
    //jLocal.init(entity,entiy);
    const BasisFunctionSetType &baseSet = jLocal.domainBasisFunctionSet();
    const unsigned int numBasisFunctions = baseSet.size();


    QuadratureType quadrature( entity, 2*dfSpace.order() );
    const size_t numQuadraturePoints = quadrature.nop();

  //  std::vector<RangeType>         uValues(numQuadraturePoints);
  //  std::vector<JacobianRangeType> uJacobians(numQuadraturePoints);

 //   std::vector<RangeType>         uOldValues(numQuadraturePoints);
 //   std::vector<JacobianRangeType> uOldJacobians(numQuadraturePoints);
    
    uEn_.resize(numQuadraturePoints);
    uLocal.evaluateQuadrature( quadrature, uEn_);
    duEn_.resize(numQuadraturePoints); 
    uLocal.evaluateQuadrature( quadrature, duEn_);

    uEnOld_.resize(numQuadraturePoints);
    uOldLocal_.evaluateQuadrature( quadrature, uEnOld_); 
    duEnOld_.resize(numQuadraturePoints);
    uOldLocal_.evaluateQuadrature( quadrature, duEnOld_);

    //  const DomainType xgl = geometry.global(x);
    RangeType vuOld(0.),vuMid(0);

    for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
    {

      const typename QuadratureType::CoordinateType &x = quadrature.point( pt );
      const double weight = quadrature.weight( pt )* geometry.integrationElement( x );

      baseSet.evaluateAll( quadrature[ pt ], phi);
      baseSet.jacobianAll( quadrature[ pt ], dphi);

      RangeType vu(0.) , vuMid(0.) ,fu(0.);
      JacobianRangeType dvu(0.) , duMid(0.), fdu(0.);
      vuOld=uEnOld_[pt];
      vu=uEn_[ pt ];
      //(1+theta)/2*U^n+(1-theta)/2*U^(n-1)
      vuMid.axpy(factorImp_,uEn_[pt]);
      vuMid.axpy(factorExp_,uEnOld_[pt]);


      //(1+theta)/2*DU^n+(1-theta)/2*DU^(n-1)
      // #if OPCHECK vuMid=vuOld
      duMid.axpy(factorImp_,duEn_[pt]);
      duMid.axpy(factorExp_,duEnOld_[pt]);


      for( size_t jj = 0; jj < numBasisFunctions ; ++jj )
      {
        RangeFieldType div(0.),grad(0.);
        for(int ii = 0; ii < dimDomain ; ++ ii)
        {
          div+=duMid[ 1+ii ][ ii ] * phi[ jj ][ 0 ] + dphi[ jj ][ 1 + ii ][ ii ]*vuMid[ 0 ];
          grad+=duMid[ 0 ][ ii ]*phi[ jj ][ 1 + ii ] + dphi[ jj ][ 0 ][ ii ]*vuMid[ 1 + ii];
        }
        fu[ 0 ]=deltaTInv*phi[ jj ][0]+0.5*( div + grad );
        
        //////////////////////////////////////////////////////////////////////////////////////////////////////////

        for( size_t ii = 0;ii <dimDomain ; ++ii)
        {
          fu[1+ii] = (vu[ 1 + ii ]-vuOld[ 1 + ii ])*phi[ jj ][0]*0.5+phi[ jj ][ 1 + ii ]*vuMid[0]*deltaTInv;
          RangeFieldType  sgradv(0.);

          for( size_t kk = 0 ; kk < dimDomain ; ++kk )
          {
            sgradv+=(duMid[ 1 + ii ][ kk ] - duMid[ 1 + kk ][ ii ]);
            sgradv*=(vuMid[ 1 + kk ]*phi[ jj ][ 0 ] + phi[ jj ][ 1 + kk ])*vuMid[0];
            sgradv+= ( dphi[ jj ][ 1 + ii ][ kk ] - dphi[ jj ][ 1 + kk ][ ii ])*vuMid[ 1 + kk ]*vuMid[ 0 ];
          }

          sgradv+=phi[ jj ][ 0 ]*duMid[ 2 + dimDomain ][ ii ] + vuMid[ 0 ]*dphi[jj][ 2 + dimDomain][ ii];
          sgradv-=phi[ jj ][ dimDomain + 3 ]*duMid[ dimDomain + 1 ][ ii ] + vuMid[ dimDomain + 3 ]*dphi[ jj ][ dimDomain+1][ ii ];
          
          sgradv*=0.5;
          fu[ 1 + ii ]+=sgradv;
        }   
        //////////////////////////////////////////////////////////////////////////////////////////////////////////

        RangeFieldType gradphiv(0.);
        
        fu[ dimDomain + 1 ] = phi[ jj ][dimDomain + 1 ]*deltaTInv;

        for( size_t ii=0 ; ii < dimDomain ; ++ ii)
        {
          gradphiv+=dphi[ jj ][ dimDomain + 1][ii] * vuMid[ 1 + ii ] + duMid[dimDomain + 1 ][ ii ]*phi[ jj ][ 1 + ii ];
        }

        //(phi_tau rho - tau phi_rho)/ rho*rho
        RangeFieldType taurho =  phi[ jj ][ dimDomain + 3]*vuMid[ 0 ] - vuMid[ dimDomain + 3 ]*phi[ jj ][ 0 ];
        
        taurho/= vuMid[ 0 ]*vuMid[ 0 ];

        fu[ dimDomain + 1]+=0.5*( gradphiv+taurho);
        //////////////////////////////////////////////////////////////////////////////////////////////////////////
        //mu  
        fu[ dimDomain + 2] = 0.5 * phi[ jj ][ dimDomain+2];
        
        RangeFieldType drhomu(0.),dphimu(0.);
        model_.drhomuSource( vu[ 0 ] , vu[ 0 ] , vu[ dimDomain + 1 ] , drhomu ); 
        model_.dphimuSource( vu[ 0 ] , vu[ 0 ] , vu[ dimDomain + 1 ] , dphimu );

        fu[ dimDomain + 2]-=drhomu*phi[ jj ][ 0 ]; 
        fu[ dimDomain + 2]-=dphimu*phi[ jj ][ dimDomain + 1];
        
        for( size_t ii = 0 ; ii < dimDomain ; ++ ii)
          fu[ dimDomain + 2 ]-=0.5*vu[ 1 + ii ]*phi[ jj ][ 1 + ii ]*0.25;

        //////////////////////////////////////////////////////////////////////////////////////////////////////////
        //tau
        fu[ dimDomain + 3 ]=0.5*phi[ jj ][ dimDomain + 3 ];

        RangeFieldType dphitau(0.);
        
        model_.dphitauSource( vu[dimDomain+1],vu[dimDomain+1], vu[0],dphitau);
        fu[ dimDomain + 3 ]-=dphitau*phi[ jj ][ dimDomain + 1]; 
        
        
        RangeFieldType divsigma(0.);
        for( size_t ii=0 ; ii < dimDomain ; ++ii )
          divsigma+=0.5*dphi[jj][dimDomain+4+ii][ii];

        fu[ dimDomain + 3 ]+=model_.delta()*divsigma;
        //////////////////////////////////////////////////////////////////////////////////////////////////////////
        for( size_t ii=0;  ii <dimDomain ; ++ii)
        {
          fu[ dimDomain + 4 + ii ] =  phi[ jj ][ dimDomain + 4 + ii] - dphi[ jj ][ dimDomain + 1][ ii ];
        }

        model_.diffusion( dphi[jj],fdu); 
        fdu*=0.5;

        fdu*=weight; 
        fu*=weight;
#if MYAXPY        
        myaxpy( jj , phi , dphi , fu , fdu , 1 , jLocal);
#else
        jLocal.column( jj ).axpy( phi,dphi, fu,fdu  );
#endif
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
        if( !visited_[ indexSet_.index( neighbor ) ])
          {
            setNeighbor( neighbor );
            typedef typename IntersectionType::Geometry  IntersectionGeometryType;
            const IntersectionGeometryType &intersectionGeometry = intersection.geometry();
#if 1
            // get local matrix for face entries 
           LocalMatrixType jLocalNbEn = jOp.localMatrix( neighbor,entity );
           LocalMatrixType jLocalEnNb = jOp.localMatrix( entity, neighbor );
            // get local matrix on neighbor 
            LocalMatrixType jLocalNbNb = jOp.localMatrix( neighbor,neighbor); 
#else
            LocalMatrixType jLocalNbEn( dfSpace,dfSpace);
            jLocalNbEn.init( neighbor, entity);
            LocalMatrixType jLocalEnNb( dfSpace, dfSpace);
            jLocalEnNb.init( entity, neighbor);
            LocalMatrixType jLocalNbNb( dfSpace,dfSpace);
            jLocalNbNb.init( neighbor,neighbor);
#endif
       const LocalFunctionType uLocalNb = u.localFunction(neighbor);
            // get neighbor's base function set 
            const BasisFunctionSetType &baseSetNb = jLocalNbEn.domainBasisFunctionSet();
            //   const unsigned int numBasisFunctionsNb = baseSetNb.size();

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
                const double weightInside=quadInside.weight( pt ); 
                const double weightOutside=quadOutside.weight( pt ); 


                baseSet.evaluateAll( quadInside[ pt ] , phi);
                baseSet.jacobianAll( quadInside[ pt ] , dphi);

                baseSetNb.evaluateAll( quadOutside[ pt ] , phiNb );
                baseSetNb.jacobianAll( quadOutside[ pt ] , dphiNb );

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

                baseSet.evaluateAll( quadInside[ pt ] , phi);
                baseSet.jacobianAll( quadInside[ pt ] , dphi);

                baseSetNb.evaluateAll( quadOutside[ pt ] , phiNb );
                baseSetNb.jacobianAll( quadOutside[ pt ] , dphiNb );


                for( size_t jj=0 ; jj < numBasisFunctions ; ++jj)
                  {
                    RangeType avuLeft(0.), avuRight(0.), valueLeft(0.),valueRight(0.);
                    JacobianRangeType aduLeft(0.),aduRight(0.);

                    double fluxRet=jacFlux_.numericalFlux( normal,                                                        
                                                           area,
                                                           uEn_[ pt ],
                                                           uNb_[ pt ],
                                                           vuMidEn,
                                                           vuMidNb,
                                                           phi[ jj ],
                                                           phiNb[ jj ],
                                                           avuLeft,
                                                           avuRight); 

                    fluxRet+=jacFlux_.diffusionFlux( normal,
                                                     penaltyFactor,
                                                     phi[ jj ],
                                                     phiNb[ jj ],
                                                     dphi[ jj ],
                                                     dphiNb[ jj ],
                                                     valueLeft,
                                                     valueRight,
                                                     aduLeft,
                                                     aduRight);

                    avuLeft+=valueLeft;
                    avuRight+=valueRight;
#if MYAXPY
                    myaxpy( jj, phi, dphi, avuLeft, aduLeft, weightInside, jLocal);
                    myaxpy( jj, phi, dphi, avuRight, aduRight, weightInside, jLocalNbEn);
#else

                    jLocal.column( jj ).axpy( phi , dphi , avuLeft,aduLeft , weightInside );
                    jLocalNbEn.column( jj ).axpy( phi,dphi , avuRight,aduRight, weightInside); 
#endif       
                    DomainType negnormal=normal; 
                    negnormal*=-1.;
                    fluxRet=jacFlux_.numericalFlux( negnormal,
                                                    area,
                                                    uNb_[ pt ],
                                                    uEn_[ pt ],
                                                    vuMidNb,
                                                    vuMidEn,
                                                    phiNb[ jj ],
                                                    phi[ jj ],
                                                    avuLeft,
                                                    avuRight); 

                    fluxRet+=jacFlux_.diffusionFlux( negnormal,
                                                     penaltyFactor,
                                                     phiNb[ jj ],
                                                     phi[ jj ],
                                                     dphi[ jj ],
                                                     dphiNb[ jj ],
                                                     valueLeft,
                                                     valueRight,
                                                     aduLeft,
                                                     aduRight);
             
                    avuLeft+=valueLeft;
                    avuRight+=valueRight;
#if MYAXPY
                    myaxpy( jj, phiNb, dphiNb, avuLeft, aduLeft, weightOutside,jLocalNbNb);
                    myaxpy( jj, phiNb, dphiNb, avuRight,aduRight, weightOutside, jLocalEnNb);
#else
                    jLocalNbNb.column( jj ).axpy( phiNb , dphiNb , avuLeft , aduLeft , weightOutside );
                    jLocalEnNb.column( jj ).axpy( phiNb , dphiNb , avuRight, aduRight, weightOutside); 
#endif    

                }
          }
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

          baseSet.evaluateAll( quadInside[ pt ] , phi);
          baseSet.jacobianAll( quadInside[ pt ] , dphi);

          const typename FaceQuadratureType::LocalCoordinateType &x = quadInside.localPoint( pt );

          const DomainType normal = intersection.integrationOuterNormal( x );
          const DomainType xgl = intersectionGeometry.global(x);
          model_.dirichletValue( time_,xgl, bndValue);


          // compute penalty factor
          const double intersectionArea = normal.two_norm();
          const double penaltyFactor =intersectionArea /  areaEn_; 
          const double area=areaEn_; 

 
          for( size_t jj=0 ; jj < numBasisFunctions ; ++jj)
          {
            RangeType avuLeft(0.), avuRight(0.),dummy(0.), valueLeft(0.),valueRight(0.);
            JacobianRangeType aduLeft(0.),aduRight(0.);
        
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

          

#if 0 
            double  fluxRet=jacFlux_.boundaryFlux( normal,
                                                   area,
                                                   vuEn[ pt ],
                                                   vuMidEn,
                                                   phi[ jj ],
                                                   avuLeft ); 

#endif
#if 1 
            fluxRet+=jacFlux_.diffusionBoundaryFlux( normal,
                                                     penaltyFactor,
                                                     phi[ jj ],
                                                     dphi[ jj ],
                                                     valueLeft,
                                                     aduLeft );                
#endif
            avuLeft+=valueLeft;
#if MYAXPY
          myaxpy( jj, phi, dphi, avuLeft, aduLeft, weight, jLocal);
#else
          jLocal.column( jj ).axpy( phi , dphi , avuLeft,aduLeft , weight );
#endif


          }
        }
      }
    }
  visited_[indexSet_.index( entity )]=true;
  } // end grid traversal 
  jOp.communicate();
}











#endif //LOCALFD_OPERATOR_HH

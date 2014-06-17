#ifndef LOCALFD_OPERATOR_HH
#define LOCALFD_OPERATOR_HH



#include <dune/fem/function/localfunction/temporary.hh>

#include <dune/fem/operator/common/differentiableoperator.hh>
#include <dune/fem/operator/common/stencil.hh>

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
  public: 
  PhasefieldJacobianOperator(const ModelType &model,
      const DiscreteFunctionSpaceType &space,
      const NumericalFluxType &flux,
      int volQuadOrder=-1)
    :MyOperatorType(model,space,flux),
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
  using MyOperatorType::setDeltaT;
  using MyOperatorType::setPreviousTimeStep;
  using MyOperatorType::getPreviousTimeStep; 
  using MyOperatorType::space;
  using MyOperatorType::penalty;
  using MyOperatorType::deltaT_;
  using MyOperatorType::uOldLocal_;
  using MyOperatorType::uOldNeighbor_;
  using MyOperatorType::factorImp_;
  using MyOperatorType::factorExp_;
  using MyOperatorType::model_;
  using MyOperatorType::areaEn_;
  using MyOperatorType::areaNb_;
  using MyOperatorType::flux_;


  typedef Dune::Fem::DiagonalAndNeighborStencil<DiscreteFunctionSpaceType,DiscreteFunctionSpaceType> StencilType;


  void jacobian(const DiscreteFunctionType &u, JacobianOperatorType &jOp) const;

  private:
  StencilType stencil_;
  const JacobianFluxType jacFlux_;
};


// Implementation of LocalFDOperator
// // ------------------------------------------------

template<class DiscreteFunction,class Model, class Flux, class Jacobian> void
PhasefieldJacobianOperator< DiscreteFunction, Model, Flux,  Jacobian>
::jacobian ( const DiscreteFunctionType &u, JacobianOperatorType &jOp ) const
{
  typedef typename JacobianOperatorType::LocalMatrixType LocalMatrixType;
  typedef typename DiscreteFunctionSpaceType::BasisFunctionSetType BasisFunctionSetType;

  Dune::Fem::DiagonalAndNeighborStencil<DiscreteFunctionSpaceType,DiscreteFunctionSpaceType> stencil(space(),space());

  jOp.reserve(stencil);

  RangeFieldType normU=std::sqrt(u.scalarProductDofs(u));
  jOp.clear();

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
      vuOld=uOldValues[pt];
      vu=uValues[ pt ];
      //(1+theta)/2*U^n+(1-theta)/2*U^(n-1)
      vuMid.axpy(factorImp_,uValues[pt]);
      vuMid.axpy(factorExp_,uOldValues[pt]);


      //(1+theta)/2*DU^n+(1-theta)/2*DU^(n-1)
      // #if OPCHECK vuMid=vuOld
      duMid.axpy(factorImp_,uJacobians[pt]);
      duMid.axpy(factorExp_,uOldJacobians[pt]);


      for( size_t jj = 0; jj < numBasisFunctions ; ++jj )
      {
        RangeFieldType div{0.},grad{0.};
        for(int ii = 0; ii < dimDomain ; ++ ii)
        {
          //   div+=Filter::dvelocity(duMid, ii, ii)*Filter::rho( phi[ jj ])
          //   +Filter::dvelocity( dphi[jj], ii, ii)*Filter::rho( vuMid );
          div+=duMid[ 1+ii ][ ii ] * phi[ jj ][ 0 ] + dphi[ jj ][ 1 + ii ][ ii ]*vuMid[ 0 ];
          
         // grad+=Filter::drho( duMid , ii )*Filter::velocity( phi[ jj ], ii )
         //   +Filter::drho(dphi[jj], ii )*Filter::velocity( vuMid, ii);
          grad+=duMid[ 0 ][ ii ]*phi[ jj ][ 1 + ii ] + dphi[ jj ][ 0 ][ ii ]*vuMid[ 1 + ii];
        }

        //Filter::rho(fu)=deltaInv*Filter::rho( phi[ jj ] )
        //  +0.5*(div+grad);
        fu[ 0 ]=deltaTInv*phi[ jj ][0]+0.5*( div + grad );
        
        //////////////////////////////////////////////////////////////////////////////////////////////////////////

        for( size_t ii = 0;ii <dimDomain ; ++ii)
        {
          //Filter::velocity( fu , ii) = (Filter::velocity( vu, ii )-Filter::velocity( vuOld ,ii ))*Filter::rho( phi[ jj ] )*0.5
          //  +Filter::velocity( phi[ jj ] , ii)*Filter::rho( vuMid)*deltaInv;

          fu[1+ii] = (vu[ 1 + ii ]-vuOld[ 1 + ii ])*phi[ jj ][0]*0.5+phi[ jj ][ 1 + ii ]*vuMid[0]*deltaTInv;


          RangeFieldType  sgradv(0.);

          for( size_t kk = 0 ; kk < dimDomain ; ++kk )
          {
            //sgradv+=(Filter::dvelocity( duMid, ii , kk ) - Filter::dvelocity( duMid , kk , ii));
            sgradv+=(duMid[ 1 + ii ][ kk ] - duMid[ 1 + kk ][ ii ]);
            
           // sgradv*=(Filter::velocity( vuMid, kk)*Filter::rho( phi[ jj ] ) + Filter::velocity( phi[ jj ] , kk )*Filter::rho( vuMid));
            sgradv*=(vuMid[ 1 + kk ]*phi[ jj ][ 0 ] + phi[ jj ][ 1 + kk ])*vuMid[0];

            //sgradv+=(Filter::dvelocity( dphi[ jj ] , ii, kk ) - Filter::dvelocity( dphi[ jj ] , kk , ii ))*Filter::velocity( vuMid , kk )*Filter::rho( vuMid );
            sgradv+= ( dphi[ jj ][ 1 + ii ][ kk ] - dphi[ jj ][ 1 + kk ][ ii ])*vuMid[ 1 + kk ]*vuMid[ 0 ];
          }

          //sgradv+=Filter::rho( phi[ jj ]  )*Filter::dmu( duMid, ii ) + Filter::rho( vuMid )*Filter::dmu( dphi[ jj ] , ii );
          sgradv+=phi[ jj ][ 0 ]*duMid[ 2 + dimDomain ][ ii ] + vuMid[ 0 ]*dphi[jj][ 2 + dimDomain][ ii];
     //     sgradv-=Filter::tau( phi[ jj ]  )*Filter::dphi( duMid, ii) + Filter::tau( vuMid )*Filter::dphi(dphi[ jj ] , ii );
          sgradv-=phi[ jj ][ dimDomain + 3 ]*duMid[ dimDomain + 1 ][ ii ] + vuMid[ dimDomain + 3 ]*dphi[ jj ][ dimDomain+1][ ii ];
          
          sgradv*=0.5;
        //  Filter::velocity( fu , ii )+=sgradv;
          
          fu[ 1 + ii ]+=sgradv;
        }   
        //////////////////////////////////////////////////////////////////////////////////////////////////////////

        RangeFieldType gradphiv{0.};

        //    Filter::phi( fu )=Filter::phi( phi[ jj ] )*deltaInv;
        
        fu[ dimDomain + 1 ] = phi[ jj ][dimDomain + 1 ]*deltaTInv;

        for( size_t ii=0 ; ii < dimDomain ; ++ ii)
        {
          //gradphiv+=Filter::dphi( dphi[ jj ] , ii) * Filter::velocity( vuMid, ii )+ Filter::dphi( duMid, ii )*Filter::velocity( phi[ jj ], ii  );
          gradphiv+=dphi[ jj ][ dimDomain + 1][ii] * vuMid[ 1 + ii ] + duMid[dimDomain + 1 ][ ii ]*phi[ jj ][ 1 + ii ];
        }

        //(phi_tau rho - tau phi_rho)/ rho*rho
        //RangeFieldType taurho=Filter::tau( phi[ jj ])*Filter::rho( vuMid ) - Filter::tau( vuMid )*Filter::rho( phi[ jj ] );
        RangeFieldType taurho =  phi[ jj ][ dimDomain + 3]*vuMid[ 0 ] - vuMid[ dimDomain + 3 ]*phi[ jj ][ 0 ];
        
        //taurho/=Filter::rho( vuMid )*Filter::rho( vuMid );
        taurho/= vuMid[ 0 ]*vuMid[ 0 ];

       // Filter::phi( fu )+=0.5*(gradphiv+taurho);
        fu[ dimDomain + 1]+=0.5*( gradphiv+taurho);
        //////////////////////////////////////////////////////////////////////////////////////////////////////////
        //mu  
        //Filter::mu( fu )=0.5*Filter::mu( phi[ jj ] );
        fu[ dimDomain + 2] = 0.5 * phi[ jj ][ dimDomain+2];
        
        RangeFieldType drhomu(0.),dphimu(0.);
        //model_.drhomuSource( Filter::rho( vu ),  Filter::rho( vu ) , Filter::phi( vu ),drhomu);
        model_.drhomuSource( vu[ 0 ] , vu[ 0 ] , vu[ dimDomain + 1 ] , drhomu ); 
        //model_.dphimuSource( Filter::rho( vu ),  Filter::rho( vu ) , Filter::phi( vu ),dphimu);
        model_.dphimuSource( vu[ 0 ] , vu[ 0 ] , vu[ dimDomain + 1 ] , dphimu );

        //Filter::mu( fu )-=drhomu*Filter::rho( phi[ jj ] ); 
        fu[ dimDomain + 2]-=drhomu*phi[ jj ][ 0 ]; 
        //Filter::mu( fu )-=dphimu*Filter::phi( phi[ jj ] );
        fu[ dimDomain + 2]-=dphimu*phi[ jj ][ dimDomain + 1];
        
        for( size_t ii = 0 ; ii < dimDomain ; ++ ii)
          fu[ dimDomain + 2 ]-=0.5*vu[ 1 + ii ]*phi[ jj ][ 1 + ii ]*0.25;
        //Filter::mu( fu )-=0.5*Filter::velocity( vu , ii )*Filter::velocity( phi[ jj ] , ii )*0.25;  

        //////////////////////////////////////////////////////////////////////////////////////////////////////////
        //tau
        //Filter::tau( fu )=0.5*Filter::tau( phi[ jj ] );
        fu[ dimDomain + 3 ]=0.5*phi[ jj ][ dimDomain + 3 ];

        RangeFieldType dphitau(0.);
        
       // model_.dphitauSource( Filter::phi(vu),
       //     Filter::phi( vu ), 
       //     Filter::rho( vu ),
       //     dphitau);
        model_.dphitauSource( vu[dimDomain+1],vu[dimDomain+1], vu[0],dphitau);
        fu[ dimDomain + 3 ]-=dphitau*phi[ jj ][ dimDomain + 1]; 
        
     //   Filter::tau( fu )-=dphitau*Filter::phi( phi[ jj ] );
        
        RangeFieldType divsigma(0.);
        for( size_t ii=0 ; ii < dimDomain ; ++ii )
          divsigma+=0.5*dphi[jj][dimDomain+4+ii][ii];
//          divsigma+=0.5*Filter::dsigma( dphi[ jj ] ,  ii , ii);

        //Filter::tau( fu )+=model_.delta()*divsigma;
        fu[ dimDomain + 3 ]+=model_.delta()*divsigma;
        //////////////////////////////////////////////////////////////////////////////////////////////////////////
        for( size_t ii=0;  ii <dimDomain ; ++ii)
        {
          //Filter::sigma( fu , ii )=Filter::sigma( phi[ jj ], ii)-Filter::dphi( dphi[ jj ], ii);
          fu[ dimDomain + 4 + ii ] =  phi[ jj ][ dimDomain + 4 + ii] - dphi[ jj ][ dimDomain + 1][ ii ];
        }

        model_.diffusion( dphi[jj],fdu); 
        fdu*=0.5;

        fdu*=weight; 
        fu*=weight;
        jLocal.column( jj ).axpy( phi,dphi, fu,fdu  );
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
        std::vector<RangeType> vuOldEn(numQuadraturePoints);
        std::vector<JacobianRangeType> duOldEn(numQuadraturePoints);
        std::vector<RangeType> vuOldNb(numQuadraturePoints);
        std::vector<JacobianRangeType> duOldNb(numQuadraturePoints);

        uLocal.evaluateQuadrature(quadInside,vuEn);
        uLocal.evaluateQuadrature(quadInside,duEn);

        uLocalNb.evaluateQuadrature(quadOutside,vuNb);
        uLocalNb.evaluateQuadrature(quadOutside,duNb);

        uOldLocal_.evaluateQuadrature(quadInside,vuOldEn);
        uOldLocal_.evaluateQuadrature(quadInside,duOldEn);

        uOldNeighbor_.evaluateQuadrature(quadOutside,vuOldNb);
        uOldNeighbor_.evaluateQuadrature(quadOutside,duOldNb);


        for( size_t pt=0 ; pt < numQuadraturePoints ; ++pt )
        {
          RangeType    vuMidEn{0.}, vuMidNb{0.};
          JacobianRangeType aduLeft{0.},aduRight{0.},duMidNb{0.}, duMidEn{0.};
          const double weight=quadInside.weight( pt ); 

          baseSet.evaluateAll( quadInside[ pt ] , phi);
          baseSet.jacobianAll( quadInside[ pt ] , dphi);

          baseSetNb.evaluateAll( quadOutside[ pt ] , phiNb );
          baseSetNb.jacobianAll( quadOutside[ pt ] , dphiNb );

          vuMidEn.axpy( factorImp_ , vuEn[ pt ] );
          vuMidEn.axpy( factorExp_ , vuOldEn[ pt ] );

          vuMidNb.axpy( factorImp_ , vuNb[ pt ] );
          vuMidNb.axpy( factorExp_ , vuOldNb[ pt ]);

          duMidEn.axpy( factorImp_ , duEn[ pt ] );
          duMidEn.axpy( factorExp_ , duOldEn[ pt ] );

          duMidNb.axpy( factorImp_ , duNb[ pt ] );
          duMidNb.axpy( factorExp_ , duOldNb[ pt ] );

          const typename FaceQuadratureType::LocalCoordinateType &x = quadInside.localPoint( pt );

          const DomainType normal = intersection.integrationOuterNormal( x );

          // compute penalty factor
          const double intersectionArea = normal.two_norm();
          const double penaltyFactor = penalty()*intersectionArea / std::min( areaEn_, areaNb_ ); 
          const double area=std::min(areaEn_,areaNb_); 

          baseSet.evaluateAll( quadInside[ pt ] , phi);
          baseSet.jacobianAll( quadInside[ pt ] , dphi);

          baseSetNb.evaluateAll( quadOutside[ pt ] , phiNb );
          baseSetNb.jacobianAll( quadOutside[ pt ] , dphiNb );


          for( size_t jj=0 ; jj < numBasisFunctions ; ++jj)
          {
                 RangeType avuLeft{0.}, avuRight{0.}, valueLeft{0.},valueRight{0.};
            JacobianRangeType aduLeft{0.},aduRight{0.};

            double  fluxRet=jacFlux_.numericalFlux(normal,
                area,
                vuEn[ pt ],
                vuNb[ pt ],
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
#if 0
             std::cout<<"Phi="<<phi[ jj ]<<"\n";

      std::cout<<"avuLeft="<<avuLeft<<"\n";

            std::cout<<"PhiNb="<<phiNb[ jj ]<<"\n";
            std::cout<<"avuRight="<<avuRight<<"\n"; 
#endif       
            jLocal.column( jj ).axpy( phi , dphi , avuLeft,aduLeft , weight );
            jLocalNb.column( jj ).axpy( phi,dphi , avuRight,aduRight, weight); 


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

          }
        }
      }
    }
  } // end grid traversal 
  jOp.communicate();
}











#endif //LOCALFD_OPERATOR_HH

#if 0
template<class DiscreteFunction, class Model, class Flux>
void DGPhasefieldOperator<DiscreteFunction, Model,Flux>
  ::operator() ( const DiscreteFunctionType &u, DiscreteFunctionType &w ) const 
{
 
  // clear destination 
  w.clear();
  assert(deltaT_>0);
  // iterate over grid 
  const IteratorType end = space().end();
  for( IteratorType it = space().begin(); it != end; ++it )
  {
    
    bool boundaryElement=false;
    // get entity (here element) 
    const EntityType &entity = *it;
    // get elements geometry
    const GeometryType& geometry=entity.geometry();
    // get local representation of the discrete functions 
    const LocalFunctionType uLocal = u.localFunction( entity);

    setEntity( entity );
    RangeType vu(0.),avu(0.);
    JacobianRangeType du(0.),adu(0.);
    
    LocalFunctionType wLocal = w.localFunction( entity );
    const int quadOrder = uLocal.order() + wLocal.order();
    
    QuadratureType quadrature( entity, quadOrder );
    const size_t numQuadraturePoints = quadrature.nop();
    for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
      {
        uLocal.evaluate( quadrature[ pt ], vu);
        uLocal.jacobian( quadrature[ pt ], du);
        
        localIntegral( pt , geometry, quadrature , vu , du , avu , adu );   
       
        //wlocal+=avu*phi+diffusion*dphi
        wLocal.axpy( quadrature[ pt ], avu, adu);
      }   
   
    if ( !space().continuous() )
    {
      //const double area = entity.geometry().volume();
      const IntersectionIteratorType iitend = space().gridPart().iend( entity ); 
      for( IntersectionIteratorType iit = space().gridPart().ibegin( entity ); iit != iitend; ++iit ) // looping over intersections
      {
        const IntersectionType &intersection = *iit;

        if ( intersection.neighbor() ) 
        {
          const EntityPointerType pOutside = intersection.outside(); // pointer to outside element.
          const EntityType &neighbor = *pOutside;

          //evaluate additional quantities on neighbor
          //penaltyfactor
          setNeighbor(neighbor);
          // compute penalty factor


          LocalFunctionType uNeighbor=u.localFunction(neighbor);

          const int quadOrderEn = uLocal.order() + wLocal.order();
          const int quadOrderNb = uNeighbor.order() + wLocal.order();

          FaceQuadratureType quadInside( space().gridPart(), intersection, quadOrderEn, FaceQuadratureType::INSIDE );
          FaceQuadratureType quadOutside( space().gridPart(), intersection, quadOrderNb, FaceQuadratureType::OUTSIDE );

          const size_t numQuadraturePoints = quadInside.nop();

          for( size_t pt=0; pt < numQuadraturePoints; ++pt )
          {
            RangeType vuEn(0.),vuNb(0.),avuLeft(0.),avuRight(0.);
            JacobianRangeType duEn(0.),duNb(0.),aduLeft(0.), aduRight(0.);
            uLocal.evaluate( quadInside[ pt ], vuEn);
            uLocal.jacobian( quadInside[ pt ], duEn);
            uNeighbor.evaluate( quadOutside[ pt ], vuNb);
            uNeighbor.jacobian( quadOutside[ pt ], duNb);
            const typename QuadratureType::CoordinateType &x = quadrature.point( pt );
            const double weight = quadInside.weight( pt );


            //calculate quadrature summands avu
            intersectionIntegral( intersection,                  
                pt, 
                quadInside,   
                quadOutside, 
                vuEn,
                vuNb, 
                duEn, 
                duNb,
                avuLeft,
                avuRight,
                aduLeft,
                aduRight);

            avuLeft*=weight;
            aduLeft*=weight;

            wLocal.axpy( quadInside[ pt ] , avuLeft , aduLeft );
          }


        }
        else if (  intersection.boundary())
        {
          boundaryElement=true;
          const int quadOrderEn = uLocal.order() + wLocal.order();

          FaceQuadratureType quadInside( space().gridPart(), intersection, quadOrderEn, FaceQuadratureType::INSIDE );

          const size_t numQuadraturePoints = quadInside.nop();


          for( size_t pt=0; pt < numQuadraturePoints; ++pt )
          {
            RangeType vuEn(0.),avuLeft(0.);
            JacobianRangeType duEn(0.),aduLeft(0.);
            uLocal.evaluate( quadInside[ pt ], vuEn);
            uLocal.jacobian( quadInside[ pt ], duEn);
            const typename QuadratureType::CoordinateType &x = quadrature.point( pt );
            const double weight = quadInside.weight( pt );


            boundaryIntegral( intersection,
                pt,
                quadInside,
                vuEn,
                duEn,
                avuLeft,
                aduLeft);
            avuLeft*=weight;
            aduLeft*=weight;

            wLocal.axpy( quadInside[ pt ] , avuLeft , aduLeft );
          }
        }

      }
#if 0     
      if( boundaryElement )
      { 
        const int order=1; 
        const LagrangePointSetType lagrangePointSet( geometry.type(), order );

        const IntersectionIteratorType iitend = space().gridPart().iend( entity ); 
        for( IntersectionIteratorType iit = space().gridPart().ibegin( entity ); iit != iitend; ++iit ) // looping over intersections
        {
          const IntersectionType &intersection = *iit;

          if ( intersection.boundary())
          {
            // get face number of boundary intersection 
            const int face = intersection.indexInInside();


            typedef typename LagrangePointSetType::template Codim< 1 >:: SubEntityIteratorType
              FaceDofIteratorType;
            // get dof iterators 
            FaceDofIteratorType faceIt = lagrangePointSet.template beginSubEntity< 1 >( face );
            const FaceDofIteratorType faceEndIt = lagrangePointSet.template endSubEntity< 1 >( face );
            for( ; faceIt != faceEndIt; ++faceIt )
            {
              const int localBlock=*faceIt;
              const int localBlockSize=DiscreteFunctionSpaceType::localBlockSize;


              const int dofOffset=localBlock*dimRange;
              wLocal[dofOffset]=0;
              for( int ii=0 ; ii < dimDomain ; ++ii)
              {
                wLocal[dofOffset+1+ii]=0;
              }

            }

          }
        }
      } 
#endif
    }


  }
  // communicate data (in parallel runs)
  w.communicate();

}
#endif

template<class DiscreteFunction, class Model, class Flux >
void DGPhasefieldOperator<DiscreteFunction, Model,Flux>
::localIntegral( size_t  pt,
    const GeometryType& geometry,
    const QuadratureType& quadrature,
    RangeType& vu,
    JacobianRangeType& du,
    RangeType& avu, // to be added to the result local function
    JacobianRangeType& adu) const
{

  const typename QuadratureType::CoordinateType &x = quadrature.point( pt );
  const double weight = quadrature.weight( pt )* geometry.integrationElement( x );
  const DomainType xgl = geometry.global(x);
  RangeType vuOld(0.),vuMid(0);

  //this should stay instide local Integral as it is operator specific
  uOldLocal_.evaluate( quadrature[ pt ], vuOld); 

  double deltaInv=1./deltaT_;
  //(1+theta)/2*U^n+(1-theta)/2*U^(n-1)
  vuMid.axpy(factorImp_,vu);
  vuMid.axpy(factorExp_,vuOld);

  JacobianRangeType duOld,duMid ;
  uOldLocal_.jacobian( quadrature[ pt ], duOld);

  //(1+theta)/2*DU^n+(1-theta)/2*DU^(n-1)
  // #if OPCHECK vuMid=vuOld
  duMid.axpy(factorImp_,du);
  duMid.axpy(factorExp_,duOld);

  RangeType  source(0.);
  model_.systemSource(time_, xgl, source);        
  //rho------------------------------------------------------------- 
  //d_t rho=(rho^n-rho^n-1)/delta t
  Filter::rho(avu)=Filter::rho(vu);
  Filter::rho(avu)-=Filter::rho(vuOld);
  Filter::rho(avu)*=deltaInv;

  RangeFieldType div(0.),gradrhodotv(0.);
  //div(rho v)=rho*div v+gradrho v
  for(int ii = 0; ii <dimDomain ; ++ii )
  { 
    //sum d_i v_i 
    div+=Filter::dvelocity(duMid, ii , ii );
    // sum d_i rho*v_i
    gradrhodotv+=Filter::drho(duMid , ii )*Filter::velocity(vuMid, ii );
  }

  Filter::rho(avu)+=div*Filter::rho(vuMid)+gradrhodotv;
  //---------------------------------------------------------------

  //v--------------------------------------------------------------

  for( int ii = 0; ii < dimDomain ; ++ii )
  {
    Filter::velocity(avu,ii)=Filter::velocity(vu,ii);
    Filter::velocity(avu,ii)-=Filter::velocity(vuOld,ii);
    Filter::velocity(avu,ii)*=deltaInv;
    RangeFieldType sgradv(0);

    //sum_j v_j( d_j v_i - d_i v_j)
    for(int jj = 0;jj < dimDomain ; ++jj)
      sgradv+=Filter::velocity( vuMid , jj )*(Filter::dvelocity(duMid, ii , jj )
          -Filter::dvelocity( duMid, jj , ii));

    Filter::velocity( avu , ii )+=sgradv;
    Filter::velocity( avu , ii )+=Filter::dmu( duMid, ii);
    Filter::velocity( avu , ii )*=Filter::rho( vuMid);

    //-tau\nabla phi
    //Filter::velocity( avu , ii )-=Filter::tau( vu )*Filter::dphi( duMid , ii );

    Filter::velocity( avu , ii )-=Filter::tau( vuMid )*Filter::dphi( duMid , ii );
    
  }
  //Filter::velocity(avu,0)+=Filter::rho(vuMid)*0.5;
  // A(dv) 
  model_.diffusion( duMid , adu );
  //------------------------------------------------------------------

  //phi---------------------------------------------------------------

  Filter::phi( avu )=Filter::phi( vu )-Filter::phi( vuOld );
  Filter::phi( avu )*=deltaInv;

  RangeFieldType transport(0.);

  // \nabla phi\cdot v
  for( int ii = 0; ii < dimDomain ; ++ii ) 
  { 
    transport+=Filter::velocity( vuMid , ii )*Filter::dphi(duMid, ii );
  }
  Filter::phi( avu )+=transport+model_.reactionFactor()*Filter::tau( vuMid )/Filter::rho(vuMid);
//  Filter::phi( avu )+=transport+model_.reactionFactor()*Filter::tau( vu )/Filter::rho(vuMid);

  //tau---------------------------------------------------------------
  // dF/dphi
  double dFdphi;
  model_.tauSource( Filter::phi(vuOld),
      Filter::phi(vu),
      Filter::rho(vuOld),
      dFdphi);

 Filter::tau( avu )=Filter::tau( vuMid );
// Filter::tau( avu )=Filter::tau( vu );
  Filter::tau( avu )-=dFdphi;

  RangeFieldType divsigma(0.);

  for( int ii = 0 ; ii < dimDomain ; ++ii) 
    divsigma+=Filter::dsigma( duMid, ii , ii );

  Filter::tau( avu )+=model_.delta()*divsigma;
  //-------------------------------------------------------------------

  //mu-----------------------------------------------------------------
  double dFdrho;
  model_.muSource(Filter::rho(vuOld),Filter::rho(vu),Filter::phi(vu),dFdrho);


  Filter::mu(avu)=Filter::mu( vuMid );
  Filter::mu(avu)-=dFdrho;
  RangeFieldType usqr(0.),uOldsqr(0.);

  for( int ii = 0; ii < dimDomain ; ++ii) 
  {
    // |v^n|^2
    usqr+=Filter::velocity( vu , ii )*Filter::velocity( vu , ii );

    // |v^{n-1}|^2
    uOldsqr+=Filter::velocity( vuOld , ii )*Filter::velocity( vuOld , ii );
  }

  Filter::mu(avu)-=0.25*(usqr+uOldsqr);
  // Filter::mu(avu)*=deltaInv;

  //------------------------------------------------------------------

  //sigma--------------------------------------------------------------
  //\sigma-\nabla\phi
  for( int ii = 0 ; ii < dimDomain ; ++ii) 
  {
    //sigma^n
    Filter::sigma( avu , ii )=Filter::sigma( vu , ii );
    //\nabla\phi^n
    Filter::sigma( avu , ii )-=Filter::dphi( du , ii );
    //   Filter::sigma( avu, ii )*=deltaInv;
  }
  //------------------------------------------------------------------        
  for(int ii = 0; ii < dimRange ; ii++)
  {
    assert( avu[ii]==avu[ii]) ;
  }

  //avu-=source;
  avu*=weight;
  adu*=weight;

}



template<class DiscreteFunction, class Model, class Flux>
template< class IntersectionQuad>
void DGPhasefieldOperator<DiscreteFunction, Model,Flux>
::intersectionIntegral( const IntersectionType& intersection,
    const size_t pt,  
    const IntersectionQuad& quadInside,
    const IntersectionQuad& quadOutside,
    const RangeType& vuEn,
    const RangeType& vuNb, 
    const JacobianRangeType& duEn,
    const JacobianRangeType& duNb,
    RangeType& avuLeft,
    RangeType& avuRight,
    JacobianRangeType& aduLeft,
    JacobianRangeType& aduRight) const
{
///  double deltaInv=1./deltaT_;

  typedef typename IntersectionType::Geometry  IntersectionGeometryType;
  //    const IntersectionGeometryType &intersectionGeometry = intersection.geometry();




  RangeType vuOldEn(0.),vuMidEn(0.),vuOldNb(0.),vuMidNb(0.);
  JacobianRangeType duOldEn(0.),duOldNb(0.),duMidEn(0.), duMidNb(0.);



  //calc vuOldEn....
  uOldLocal_.evaluate( quadInside[ pt ] , vuOldEn );
  uOldLocal_.jacobian( quadInside[ pt ] , duOldEn );
  uOldNeighbor_.evaluate( quadOutside[ pt ] , vuOldNb );
  uOldNeighbor_.jacobian( quadOutside[ pt ] , duOldNb );

  vuMidEn.axpy( factorImp_ , vuEn );
  vuMidEn.axpy( factorExp_ , vuOldEn );

  vuMidNb.axpy( factorImp_ , vuNb );
  vuMidNb.axpy( factorExp_ , vuOldNb);

  duMidEn.axpy( factorImp_ , duEn );
  duMidEn.axpy( factorExp_ , duOldEn );

  duMidNb.axpy( factorImp_ , duNb );
  duMidNb.axpy( factorExp_ , duOldNb );


  const typename FaceQuadratureType::LocalCoordinateType &x = quadInside.localPoint( pt );

  const DomainType normal = intersection.integrationOuterNormal( x );
  // compute penalty factor
  const double intersectionArea = normal.two_norm();
  const double penaltyFactor = penalty()*intersectionArea / std::min( areaEn_, areaNb_ ); 
  const double area=std::min(areaEn_,areaNb_); 
  // const typename FaceQuadratureType::LocalCoordinateType &x = quadInside.localPoint( pt );
  //   DomainType xgl=intersectionGeometry.global(x); 
  //  const double weight = quadInside.weight( pt );

  JacobianRangeType dvalue(0.),advalue(0.);
  double fluxRet;
  //RangeType gLeft(0.);
// RangeType gRight(0.);

  fluxRet=flux_.numericalFlux(normal,
      area,
      vuEn,
      vuNb,
      vuMidEn,  
      vuMidNb, 
      avuLeft,
      avuRight); 
#if 0  
  std::cout<<"Flux===\n";
  std::cout<<"normal"<<normal<<"\n";
  std::cout<<"area="<<area<<"\n";
  std::cout<<"vuEn="<<vuEn<<"\n";
  std::cout<<"vuNb="<<vuNb<<"\n";
  std::cout<<"vuMidEn="<<vuMidEn<<"\n";
  std::cout<<"vuMidNb="<<vuMidNb<<"\n";
  std::cout<<"Fleft="<<avuLeft<<"\n";
  std::cout<<"Fright="<<avuRight<<"\n";

#endif
  

  RangeType value(0);
#if 1        
  fluxRet+=flux_.diffusionFlux(normal,
      penaltyFactor,
      vuMidEn,
      vuMidNb,
      duMidEn,
      duMidNb,
      value,
      aduLeft);
#endif
  avuLeft+=value;
  avuRight-=value; 
#if 0 
  Filter::tau(avuLeft)*=deltaInv;
  Filter::tau(avuRight)*=deltaInv;
  Filter::mu( avuLeft)*=deltaInv;
  Filter::mu( avuRight)*=deltaInv;
#endif
  //      for(int ii=0; ii<dimRange ; ++ii)
  //      {
  //      Filter::sigma( avuLeft , ii )*=deltaInv;
  //    Filter::sigma( avuRight , ii )*=deltaInv;
  // } 
  aduRight=aduLeft;
  aduRight*=-1.;  
}

//Boundary Intgral
template<class DiscreteFunction, class Model, class Flux>
void DGPhasefieldOperator<DiscreteFunction, Model,Flux>
::boundaryIntegral( const IntersectionType& intersection,
    const size_t pt,  
    const FaceQuadratureType& quadInside,
    const RangeType& vuEn,
    const JacobianRangeType& duEn,
    RangeType& avuLeft,
    JacobianRangeType& aduLeft) const
{
  typedef typename IntersectionType::Geometry  IntersectionGeometryType;
  const IntersectionGeometryType &intersectionGeometry = intersection.geometry();


  RangeType vuOldEn(0.),vuMidEn(0.);
  JacobianRangeType duOldEn(0.),duMidEn(0.);



  //calc vuOldEn....
  uOldLocal_.evaluate( quadInside[ pt ] , vuOldEn );
  uOldLocal_.jacobian( quadInside[ pt ] , duOldEn );

  vuMidEn.axpy( factorImp_ , vuEn );
  vuMidEn.axpy( factorExp_ , vuOldEn );


  duMidEn.axpy( factorImp_ , duEn );
  duMidEn.axpy( factorExp_ , duOldEn );



  // compute penalty factor
  const double intersectionArea = intersectionGeometry.volume();
  const double penaltyFactor = penalty()*intersectionArea /  areaEn_; 
  const double area=std::min(areaEn_,areaNb_); 
  const typename FaceQuadratureType::LocalCoordinateType &x = quadInside.localPoint( pt );
  const DomainType normal = intersection.integrationOuterNormal( x );

  DomainType xgl=intersectionGeometry.global(x); 
  //  const double weight = quadInside.weight( pt );

  JacobianRangeType dvalue(0.),advalue(0.);
  double fluxRet;
  RangeType gLeft(0.),gRight(0.);
  fluxRet=flux_.boundaryFlux(normal,
      area,
      vuEn,
      vuMidEn,  
      avuLeft);


  RangeType value(0.);
#if 1        
  fluxRet+=flux_.diffusionBoundaryFlux(normal,
      penaltyFactor,
      vuMidEn,
      duMidEn,    
      value,
      advalue);
#endif
  avuLeft+=value;



}

template< class DiscreteFunction,class Model, class Flux>
template< class LocalArgType, class LFDestType>
void DGPhasefieldOperator<DiscreteFunction, Model, Flux>
::computeBoundary( const IntersectionType& intersection,
    const EntityType& entity,
    const double area,
    const LocalArgType& uEn,
    LFDestType& wLocal) const
{
  abort();

}











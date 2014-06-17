
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
  RangeType vuOld(0.),vuMid{0};

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

  //RangeType  source(0.);
//  model_.systemSource(time_, xgl, source);        
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
    Filter::velocity( avu , ii )-=Filter::tau( vuMid )*Filter::dphi( duMid , ii );

  }
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
  //------------------------------------------------------------------        

  //tau---------------------------------------------------------------
  // dF/dphi
  double dFdphi;

  model_.tauSource( Filter::phi(vuOld),
      Filter::phi(vu),
      Filter::rho(vuOld),
      dFdphi);

  Filter::tau( avu )=Filter::tau( vuMid );
  Filter::tau( avu )-=dFdphi;

  RangeFieldType divsigma(0.) , gradrhosigma(0.);

  for( int ii = 0 ; ii < dimDomain ; ++ii) 
    {
      divsigma+=Filter::dalpha( duMid, ii , ii );
     // gradrhosigma+=Filter::sigma( vuMid, ii )*Filter::drho( duMid, ii);
    }
  
  //delta*\Div(\rho\sigma)=delta*\div \nabla\rho\cdot\sigma+\rho\Div\sigma
  Filter::tau( avu )+=model_.delta()*(divsigma);
  //  Filter::tau( avu)*=deltaInv;         
  //-------------------------------------------------------------------

  //mu-----------------------------------------------------------------
  double dFdrho;

  //model_.muSource(Filter::rho(vu),Filter::rho(vu),Filter::phi(vu),dFdrho);
  //old version like Paris talk
  model_.muSource(Filter::rho(vuOld),Filter::rho(vu),Filter::phi(vu),dFdrho);



  Filter::mu(avu)=Filter::mu( vuMid );
  Filter::mu(avu)-=dFdrho;
  RangeFieldType usqr(0.),uOldsqr(0.),sigmasqr(0.),sigmaOldsqr(0.);

  for( int ii = 0; ii < dimDomain ; ++ii) 
  {
    // |v^n|^2
    usqr+=Filter::velocity( vu , ii )*Filter::velocity( vu , ii );

    // |v^{n-1}|^2
    uOldsqr+=Filter::velocity( vuOld , ii )*Filter::velocity( vuOld , ii );
   
    // |sigma^n|^2
    sigmasqr+=model_.hprime(Filter::rho(vu))*Filter::sigma( vu , ii )*Filter::sigma( vu , ii );
    
    // |sigma^{n-1}|^2
    sigmaOldsqr+=model_.hprime(Filter::rho(vuOld))*Filter::sigma( vuOld , ii )*Filter::sigma( vuOld , ii );

  }

  Filter::mu(avu)-=0.25*(usqr+uOldsqr);
  Filter::mu(avu)-=0.25*model_.delta()*(sigmasqr+sigmaOldsqr);

  //------------------------------------------------------------------

  //sigma--------------------------------------------------------------
  //\sigma-\nabla\phi
  for( int ii = 0 ; ii < dimDomain ; ++ii) 
  {
    //sigma^n
    Filter::sigma( avu , ii )=Filter::sigma( vu , ii );
    //\nabla\phi^n
    Filter::sigma( avu , ii )-=Filter::dphi( du , ii );
 
    Filter::alpha( avu, ii )=Filter::alpha(vuMid,ii)-model_.h(Filter::rho(vuMid))*Filter::sigma(vuMid, ii);
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
template<class IntersectionQuad>
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
  double deltaInv=1./deltaT_;

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
  RangeType gLeft(0.),gRight(0.);
  fluxRet=flux_.numericalFlux(normal,
      area,
      vuEn,
      vuNb,
      vuMidEn,  
      vuMidNb, 
      avuLeft,
      avuRight); 

  RangeType value(0.);
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











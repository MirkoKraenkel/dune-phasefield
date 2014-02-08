
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
    // get entity (here element) 
    const EntityType &entity = *it;
    // get elements geometry
    const GeometryType& geometry=entity.geometry();
    // get local representation of the discrete functions 
    const LocalFunctionType uLocal = u.localFunction( entity);

    setEntity( entity );
    RangeType vu{0.},avu{0.};
    JacobianRangeType du{0.},adu{0.};
    
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
           RangeType vuEn{0.},vuNb{0.},avuLeft{0.};
           JacobianRangeType duEn{0.},duNb{0.},aduLeft{0.};
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
                                 aduLeft );
 
            avuLeft*=weight;
            aduLeft*=weight;
        
            wLocal.axpy( quadInside[ pt ] , avuLeft , aduLeft );
          }
        
            
        }
        else if ( intersection.boundary() )
        {
#if 0
          boundaryIntegral( intersection,
                            pt,
                            quadInside,
                            vuEn,
                            duEn,
                            avuLeft,
                            aduLeft);

#endif   
        }
      }
    }
      
    
  }
  // communicate data (in parallel runs)
  w.communicate();

}

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
    RangeType vuOld{0.},vuMid{0};
    
    //this should stay instide local Integral as it is perator specific
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
    RangeFieldType div{0.},gradrhodotv{0.};
#if 1 
    //div(rho v)=rho*div v+gradrho v
    for(int ii = 0; ii <dimDomain ; ++ii )
      { 
        //sum d_i v_i 
        div+=Filter::dvelocity(duMid, ii , ii );
        // sum d_i rho*v_i
        gradrhodotv+=Filter::drho(duMid , ii )*Filter::velocity(vuMid, ii );
      }

    Filter::rho(avu)+=Filter::rho(vuMid)*div+gradrhodotv;
#endif
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
       // Filter::phi( avu )+=transport+Filter::tau( vuMid )/Filter::rho( vuMid );
        Filter::phi( avu )+=transport;
        Filter::phi( avu )+=Filter::tau( vuMid );
//------------------------------------------------------------------        
       
//tau---------------------------------------------------------------
        // dF/dphi
        double dFdphi;
        model_.tauSource( Filter::phi(vuOld),
                          Filter::phi(vu),
                          Filter::rho(vuOld),
                          dFdphi);

        Filter::tau(avu)=Filter::tau( vuMid );
        Filter::tau( avu )-=dFdphi;

        RangeFieldType divsigma(0.);

        for( int ii = 0 ; ii < dimDomain ; ++ii) 
          divsigma+=Filter::dsigma( duMid, ii , ii );
         
        Filter::tau( avu )+=model_.delta()*divsigma;
     //-------------------------------------------------------------------

//mu-----------------------------------------------------------------

     //   Filter::mu(avu)=Filter::mu(vu)-Filter::mu(vuOld);
        Filter::mu(avu)=Filter::mu( vuMid );
        
        RangeFieldType usqr{0.},uOldsqr{0.};
       
        for( int ii = 0; ii < dimDomain ; ++ii) 
        {
         // |v^n|^2
          usqr+=Filter::velocity( vu , ii )*Filter::velocity( vu , ii );

          // |v^{n-1}|^2
          uOldsqr+=Filter::velocity( vuOld , ii )*Filter::velocity( vuOld , ii );
        }
     
        Filter::mu(avu)-=0.25*(usqr+uOldsqr);
      
        
        //------------------------------------------------------------------

//sigma--------------------------------------------------------------
        //\sigma-\nabla\phi
        for( int ii = 0 ; ii < dimDomain ; ++ii) 
          {
            //sigma^n
           Filter::sigma( avu , ii )=Filter::sigma( vu , ii );
           //\nabla\phi^n
           Filter::sigma( avu , ii )-=Filter::dphi( du , ii );
          }
          //------------------------------------------------------------------        
          for(int ii = 0; ii < dimRange ; ii++)
          {
            assert( avu[ii]==avu[ii]) ;
          }
          
          avu-=source;
          avu*=weight;
          adu*=weight;
       
      }
    


template<class DiscreteFunction, class Model, class Flux>
void DGPhasefieldOperator<DiscreteFunction, Model,Flux>
::intersectionIntegral( const IntersectionType& intersection,
                        const size_t pt,  
                        const FaceQuadratureType& quadInside,
                        const FaceQuadratureType& quadOutside,
                        const RangeType& vuEn,
                        const RangeType& vuNb, 
                        const JacobianRangeType& duEn,
                        const JacobianRangeType& duNb,
                        RangeType& avuLeft,
                        JacobianRangeType& aduLeft) const
  {
    typedef typename IntersectionType::Geometry  IntersectionGeometryType;
    const IntersectionGeometryType &intersectionGeometry = intersection.geometry();
  
    
    
    
    
    RangeType vuOldEn{0.},vuMidEn{0.},vuOldNb{0.},vuMidNb{0.};
    JacobianRangeType duOldEn{0.},duOldNb{0.},duMidEn{0.}, duMidNb{0.};


    
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


    // compute penalty factor
    const double intersectionArea = intersectionGeometry.volume();
    const double penaltyFactor = penalty()*intersectionArea / std::min( areaEn_, areaNb_ ); 
    const double area=std::min(areaEn_,areaNb_); 
    const typename FaceQuadratureType::LocalCoordinateType &x = quadInside.localPoint( pt );
    const DomainType normal = intersection.integrationOuterNormal( x );

    DomainType xgl=intersectionGeometry.global(x); 
  //  const double weight = quadInside.weight( pt );
             
     JacobianRangeType dvalue{0.},advalue{0.};
     double fluxRet;
     RangeType gLeft{0.},gRight{0.};
     fluxRet=flux_.numericalFlux(normal,
                                 area,
                                 vuEn,
                                 vuNb,
                                 vuMidEn,  
                                 vuMidNb, 
                                 avuLeft,
                                 gRight); 
   
        RangeType value{0.};
#if 1        
        fluxRet+=flux_.diffusionFlux(normal,
                                    penaltyFactor,
                                    vuMidEn,
                                    vuMidNb,
                                    duMidEn,
                                    duMidNb,
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
    typedef typename IntersectionType::Geometry  IntersectionGeometryType;

    const IntersectionGeometryType &intersectionGeometry = intersection.geometry();
    LocalFunctionType uOldEn=uOld_.localFunction(entity); 
    
    TemporaryLocalType ufMidEn(space());

    ufMidEn.init( entity );
    ufMidEn.clear();

    ufMidEn.axpy( factorImp_, uEn );
    ufMidEn.axpy( factorExp_, uOldEn );



    DomainType xglobal{0};

    // compute penalty factor
    const double intersectionArea = intersectionGeometry.volume();
    const double penaltyFactor=penalty()*intersectionArea / area;
    const int quadOrder = uEn.order() + wLocal.order();

    FaceQuadratureType quadInside( space().gridPart(), intersection, quadOrder, FaceQuadratureType::INSIDE );
    const size_t numQuadraturePoints = quadInside.nop();
    double fluxRet=0.;
    std::vector<RangeType> valuesEn( numQuadraturePoints );
    std::vector<RangeType> midValuesEn( numQuadraturePoints );
    std::vector<JacobianRangeType> midJacobiansEn( numQuadraturePoints );
    
    uEn.evaluateQuadrature( quadInside, valuesEn );
    ufMidEn.evaluateQuadrature( quadInside,midValuesEn );
    ufMidEn.evaluateQuadrature( quadInside,midJacobiansEn );
 
    for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
      {
        const typename FaceQuadratureType::LocalCoordinateType &x = quadInside.localPoint( pt );
        xglobal = intersectionGeometry.global(x);
      

        const DomainType normal = intersection.integrationOuterNormal( x );
        const double weight = quadInside.weight( pt );
      
        RangeType value, valueBnd,valueBndOld;
        JacobianRangeType dvalue{0.},advalue{0.};

        RangeType vuIn,jump, vuOld, vuMid,gLeft;
        JacobianRangeType duIn, aduIn, duMidNb{0.};

        model_.dirichletValue(time_,xglobal,valueBnd); 
        model_.dirichletValue(time_-deltaT_, xglobal,valueBndOld);
      
        valueBnd*=factorImp_;
        valueBnd.axpy( factorExp_,valueBndOld );


        
        flux_.boundaryFlux(normal,midValuesEn[ pt ],gLeft);

#warning "DIFFUSIONFLUX MISSING AT BOUNDARY"

#if 1     
        fluxRet+=flux_.diffusionFlux( normal,
                                      penaltyFactor,
                                      midValuesEn[ pt ],
                                      valueBnd,
                                      midJacobiansEn[ pt ],
                                      duMidNb,//=0.
                                      value,
                                      advalue);
        gLeft+=value;
        gLeft*=weight;

        advalue*=weight;
#endif 
        wLocal.axpy(quadInside[pt],gLeft,advalue);
     }
  }











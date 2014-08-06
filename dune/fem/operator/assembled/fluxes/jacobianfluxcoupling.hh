#ifndef PHASEFIELD_JACOBIAN_FLUX_HH
#define PHASEFIELD_JACBIAN_FLUX_HH

#include <dune/common/fvector.hh>

#include "../phasefieldfilter.hh"
#warning "JACOBIANFLUX"
template<class Model>
class JacobianFlux
{
  //typedefs

 typedef Model ModelType;
 
  enum{ dimDomain=Model::dimDomain};
  enum{ dimRange=Model::dimRange};
 
  typedef typename Dune::FieldVector<double,dimRange> RangeType;
  typedef typename Dune::FieldVector<double,dimDomain> DomainType;
  typedef typename Dune::FieldMatrix<double,dimRange,dimDomain> JacobianRangeType;
  typedef typename Dune::FieldMatrix<double,dimRange,dimRange> FluxRangeType;
  typedef  PhasefieldFilter<RangeType> Filter; 

  
public:
  JacobianFlux(const ModelType& model):
    model_(model),
    penalty_( Dune::Fem::Parameter::getValue< double >( "phasefield.penalty" )),
    acpenalty_( Dune::Fem::Parameter::getValue< double >( "phasefield.acpenalty" )),
    switchIP_(Dune::Fem::Parameter::getValue<int>("phasefield.ipswitch",1)),
    numVisc_(Dune::Fem::Parameter::getValue<double>("phasefield.addvisc",0))  
    {
    }


  void  numericalFlux ( const DomainType& normal,
                        const double penaltyFactor,
                         const RangeType& vuEn,
                         const RangeType& vuNb,
                         FluxRangeType& gLeft,
                         FluxRangeType& gRight) const;

  template< class ScalarType, class JacobianType, class DiffusionTensorType, class DiffusionVectorType>
  void scalar2vectorialDiffusionFlux( const DomainType& normal,
                                  const double penaltyFactor,
                                  const ScalarType& phiEn,
                                  const ScalarType& phiNb,
                                  const JacobianType& dphiEn,
                                  const JacobianType& dphiNb,
                                  DiffusionTensorType& aLeft,
                                  DiffusionTensorType& aRight,
                                  DiffusionVectorType& bLeft,
                                  DiffusionVectorType& bRight) const;

  double diffusionFlux ( const DomainType& normal,
                         const double penaltyFactor,
                         const RangeType& uEn,
                         const RangeType& uNb,
                         const JacobianRangeType& duEn,
                         const JacobianRangeType& duNb,
                         RangeType& valueLeft,
                         RangeType& valueRight,
                         JacobianRangeType& dvalueLeft,
                         JacobianRangeType& dvalueRight) const;
  
  void boundaryFlux( const DomainType& normal, 
                       const double penaltyFactor,                
                       const RangeType& vuEn,  
                       const RangeType& vuMidEn,
                       RangeType& phiEn,
                       RangeType& gLeft) const;

  double diffusionBoundaryFlux( const DomainType& normal,
                                const double penaltyFactor,
                                const RangeType& uEn,
                                const JacobianRangeType& duEn,
                                RangeType& value,
                                JacobianRangeType& dvalue) const;

private:
  const ModelType& model_;
  double penalty_;
  double acpenalty_;
  const int switchIP_; 
  const double numVisc_;
};


template< class Model >
void JacobianFlux<Model>
::boundaryFlux(const DomainType& normal,
               const double penaltyFactor,
               const RangeType& vuEn,
               const RangeType& vuMidEn,
               RangeType& phiEn,
               RangeType& gLeft) const
  {
    RangeType midEn ;

    midEn=vuMidEn;

    double vNormalEn(0.);

    //rho-------------------------------------------------------------
  
    for(int i = 0; i<dimDomain;++i)
     {
         vNormalEn+=Filter::velocity(midEn,i)*normal[i];
     }
  
    Filter::rho(gLeft)=-1*vNormalEn*Filter::rho(midEn);
    //Filter::rho(gLeft)*=-0.5;
  
    //----------------------------------------------------------------
    
    //v---------------------------------------------------------------
  
    for(int i = 0; i<dimDomain;++i)
      {
        Filter::velocity(gLeft,i)=0;
      } 
    //----------------------------------------------------------------

    //phi-------------------------------------------------------------
  
    double laplaceFlux(0.);

    for(int i = 0; i<dimDomain;++i)
      {
        laplaceFlux-=Filter::sigma(midEn,i)*normal[i];
      }  
  
    //----------------------------------------------------------------

    //tau-------------------------------------------------------------

    Filter::tau(gLeft)=model_.delta()*laplaceFlux;

    //----------------------------------------------------------------

    //sigma-----------------------------------------------------------

    for(int i = 0; i<dimDomain;++i)
      {   
        Filter::sigma(gLeft,i)=0;
      } 
    //----------------------------------------------------------------
  
    return 0.;
  }




template< class Model >
void  JacobianFlux<Model>
::numericalFlux( const DomainType& normal,
                const double penaltyFactor,              
                const RangeType& vuEnMid,
                const RangeType& vuNbMid,
                FluxRangeType& gLeft,
                FluxRangeType& gRight) const
  {
      RangeType valEn,valNb,midEn(0.), midNb(0.),jump,mean, jumpPhi(0.);
      
      gLeft=0.;
      gRight=0.;

      midEn=vuEnMid;
      midNb=vuNbMid;
      jump = midEn;
      jump-= midNb;
      mean = midEn ;
      mean+= midNb;
      mean*=0.5;
      
      
      //rho-------------------------------------------------------------
 
      double vNormalEn(0),testNormalEn(0.),vNormalNb(0.),testNormalNb(0.);
    
      for(int ii = 0; ii < dimDomain ; ++ii )
        {
          //v^+\cdot n^+
          vNormalEn+=midEn[ 1 + ii ]*normal[ ii ];
          //v^-\cdot n^+
          vNormalNb+=midNb[ 1 + ii ]*normal[ ii ];
          gLeft[ 0 ][ 1 + ii ]  = -0.25*normal[ ii ]*midEn[ 0 ];
          gRight[ 0 ][ 1 + ii ] =  0.25*normal[ ii ]*midNb[ 0 ];
        }
  

      //F_1=-0.5*( \rho^+*v^+\cdot n^+ -\rho^-*v-\cdot n^+)  
      gLeft[ 0 ][ 0 ]=-0.25*vNormalEn;
      //from linerarization  vuMid=0.5(vu+vuOld), d vuMid/d vu=0.5 
      gRight[ 0 ][ 0 ]=0.25*vNormalNb;

      //----------------------------------------------------------------
    
      //v---------------------------------------------------------------
    
      for(int ii = 0; ii < dimDomain ; ++ii )
        { 
          //(\rho_j,\v_i)
          gLeft[ 1 + ii ][ 0 ]=-0.5*jump[ dimDomain + 2 ]*normal[ ii ];
          //(\mu_j,\v_i)
          gLeft[ 1 + ii ][ dimDomain + 2 ]=-0.5*normal[ ii ]*midEn[ 0 ];
          gRight[ 1 + ii ][ dimDomain + 2 ]=0.5*normal[ ii ]*midEn[ 0 ]; 
         //(\tau_j,v_i)
          gLeft[ 1 + ii ][ dimDomain + 3 ]=jump[ dimDomain + 1 ]*normal[ ii ]*0.5;
          //(\phi_j,v_i)
          gLeft[ 1 + ii ][ dimDomain + 1] =0.5*normal[ ii ]*midEn[ dimDomain + 3 ];
          gRight[ 1 + ii ][ dimDomain + 1 ]=-0.5*normal[ ii ]*midEn[ dimDomain + 3 ];
          
          gLeft[ 1 + ii ]*=0.5;
          gRight[ 1 + ii ]*=0.5;

        } 
    
      //----------------------------------------------------------------
      double laplaceFluxLeft(0.),laplaceFluxRight(0.);

      //phi-------------------------------------------------------------
      for( int ii = 0 ; ii<dimDomain ; ++ii )
        {
          //F_{3.1}
       
          //-(\phi^+-\phi^-)*n[i]*v[i]^+*0.5 
          gLeft[ dimDomain + 1 ][ dimDomain + 1 ]+=normal[ ii ]*midEn[ 1 + ii ]*-0.25;
          gRight[ dimDomain + 1 ][ dimDomain +1 ]+=normal[ ii ]*midEn[ 1 + ii ]*0.25;

          gLeft[ dimDomain + 1 ][ 1 + ii ]=jump[ dimDomain + 1 ]*normal[ ii ]*-0.25;
              //tau
          //F_{3.2}
          //(\sigma^+-\sigma^-)\cdot n * 0.5
          laplaceFluxLeft=normal[ ii ]*0.25;
          laplaceFluxRight=normal[ ii ]*-0.25;
       //tau-------------------------------------------------------------
      gLeft[ dimDomain + 3 ][ dimDomain + 4+ii ]=-model_.delta()*laplaceFluxLeft;
      gRight[ dimDomain + 3 ][dimDomain + 4+ii ]=-model_.delta()*laplaceFluxRight;
        } 
    
      //----------------------------------------------------------------
     //----------------------------------------------------------------

      //sigma-----------------------------------------------------------
      for(int ii = 0; ii < dimDomain ; ++ii )
        {  
          //F_4
          //(\phi^+-\phi^-)\normal*0.5
          gLeft[ dimDomain + 4 + ii ][ dimDomain + 1 ]=normal[ ii ]*0.5;
          gRight[ dimDomain + 4 + ii ][ dimDomain + 1 ]=normal[ ii ]*-0.5;
  
        } 
      //----------------------------------------------------------------
  }
template< class Model >
template< class ScalarType , class JacobianType ,class DiffusionTensorType, class DiffusionVectorType>
void JacobianFlux<Model>
:: scalar2vectorialDiffusionFlux( const DomainType& normal,
                                  const double penaltyFactor,
                                  const ScalarType& phiEn,
                                  const ScalarType& phiNb,
                                  const JacobianType& dphiEn,
                                  const JacobianType& dphiNb,
                                  DiffusionTensorType& aLeft,
                                  DiffusionTensorType& aRight,
                                  DiffusionVectorType& bLeft,
                                  DiffusionVectorType& bRight) const
{
  double integrationElement=normal.two_norm();
  //[\phiEn]\otimes n ; [\phiNb]\otimes n
  JacobianType jumpNormalLeft(0.),jumpNormalRight(0.);

  for(int j=0; j<dimDomain; ++j)
    { 
      jumpNormalLeft[0][j]=-0.5*phiEn[0]*normal[j];
      jumpNormalRight[0][j]=-0.5*phiNb[0]*normal[j];
    }

  DiffusionTensorType bbLeft,bbRight;
  jumpNormalLeft*=switchIP_;
  jumpNormalRight*=switchIP_;
  
  model_.scalar2vectorialDiffusion(jumpNormalLeft ,aLeft);
  model_.scalar2vectorialDiffusion(jumpNormalRight,aRight);
  
  model_.scalar2vectorialDiffusion( dphiEn, bbLeft);
  model_.scalar2vectorialDiffusion( dphiNb, bbRight);
  
  //loop over components
  for( int ii = 0 ; ii < dimDomain ; ++ ii)
  {
    bbLeft[ ii ]*=-0.5;
    bbLeft[ ii ].mv( normal , bLeft[ ii ]);
    bbLeft[ii]=0;
    bbRight[ ii ]*=-0.5;
    bbRight[ ii ].mv( normal , bRight[ ii ]);
   // for( int jj = 0 ; jj < dimDomain ;++jj)
    {
      bRight[ii][ii]-=penalty_*penaltyFactor*integrationElement*phiNb[ 0 ];
      bLeft[ii][ii]+=penalty_*penaltyFactor*integrationElement*phiEn[ 0 ];
    }
  }
  
}






template< class Model >
double JacobianFlux<Model>
:: diffusionFlux( const DomainType& normal,
                  const double penaltyFactor,
                  const RangeType& uEn,
                  const RangeType& uNb,
                  const JacobianRangeType& duEn,
                  const JacobianRangeType& duNb,
                  RangeType& valueLeft,
                  RangeType& valueRight,
                  JacobianRangeType& dvalueLeft,
                  JacobianRangeType& dvalueRight) const
{
  
  RangeType jump{0};
  jump=uEn;
  jump-=uNb;  
  RangeType phiEn=uEn;
  RangeType phiNb=uNb;
  JacobianRangeType aduEn(0.), aduNb(0.); 
  double integrationElement=normal.two_norm();
  
  
  for( int i=0; i<dimDomain;++i)
    {
      Filter::velocity(valueLeft,i)=penalty_*penaltyFactor*Filter::velocity(phiEn,i)*integrationElement*0.5;
      Filter::velocity(valueRight,i)=-penalty_*penaltyFactor*Filter::velocity(phiNb,i)*integrationElement*0.5;
    }
  JacobianRangeType jumpNormalLeft(0.),jumpNormalRight(0.);
 
  // [u]\otimes n
  for(int i=0; i<dimDomain; ++i)
    for(int j=0; j<dimDomain; ++j)
    { 
      jumpNormalLeft[i+1][j]=-0.5*uEn[i+1]*normal[j];
      jumpNormalRight[i+1][j]=0.5*uNb[i+1]*normal[j];
    }

  jumpNormalLeft*=switchIP_;
  jumpNormalRight*=switchIP_;
  model_.diffusion(jumpNormalLeft,dvalueLeft);
  model_.diffusion(jumpNormalRight,dvalueRight);
   
  JacobianRangeType mean(0.), Amean(0.);
  mean=duEn;
  mean*=-0.25;
  model_.diffusion(mean,Amean);

  Amean.umv(normal,valueLeft);
  Amean=0;
  mean=duNb;
  mean*=-0.25;
  model_.diffusion(mean,Amean);

  Amean.umv(normal,valueRight);
 
  
  
  return penalty_;
}
template< class Model >
double JacobianFlux<Model>
:: diffusionBoundaryFlux( const DomainType& normal,
                  const double penaltyFactor,
                  const RangeType& uEn,
                  const JacobianRangeType& duEn,
                  RangeType& value,
                  JacobianRangeType& dvalue) const
{
  
  RangeType jump;
  jump=uEn;
  JacobianRangeType aduEn(0.), aduNb(0.); 
  double integrationElement=normal.two_norm();
  
  for( int i=0; i<dimDomain;++i)
    {
      Filter::velocity(value,i)=penalty_*penaltyFactor*Filter::velocity(jump,i)*integrationElement;
    }
  JacobianRangeType jumpNormal(0.);
 
  // [u]\otimes n
  for(int i=0; i<dimDomain; ++i)
    for(int j=0; j<dimDomain; ++j)
      jumpNormal[i+1][j]=-0.5*jump[i+1]*normal[j];
 

  jumpNormal*=switchIP_;
  model_.diffusion(jumpNormal,dvalue);
   
  JacobianRangeType mean(0.), Amean(0.);
  mean=duEn;
  //mean+=duNb;
  //mean*=-0.5;
  model_.diffusion(mean,Amean);

  Amean.umv(normal,value);
  return penalty_;
}

#endif

#ifndef PHASEFIELD_JACOBIAN_SPLITFLUX_HH
#define PHASEFIELD_JACOBIAN_SPLITFLUX_HH

#include <dune/common/fvector.hh>

template<class Model>
class AllenCahnJacobianFlux
{
  //typedefs

 typedef Model ModelType;
 
  enum{ dimDomain=Model::dimDomain};
  enum{ dimRange=Model::dimAcRange};
  
  static const int phi=0;
  static const int tau=1;
  static const int sigma=2;

  typedef typename Dune::FieldVector<double,dimRange> RangeType;
  typedef typename Dune::FieldVector<double,dimDomain> DomainType;
  typedef typename Dune::FieldMatrix<double,dimRange,dimDomain> JacobianRangeType;
  typedef typename Dune::FieldMatrix<double,dimRange,dimRange> FluxRangeType;
  
public:
  AllenCahnJacobianFlux(const ModelType& model):
    model_(model),
    acpenalty_( Dune::Fem::Parameter::getValue< double >( "phasefield.acpenalty" ))
    {
    }


  void  numericalFlux ( const DomainType& normal,
                        const double area,
                        const double penaltyFactor,
                        const RangeType& vuEn,
                        const RangeType& vuNb,
                        const RangeType& vuEnAdd,
                        const RangeType& vuNbAdd,
                        FluxRangeType& gLeft,
                        FluxRangeType& gRight) const;

 
  void boundaryFlux( const DomainType& normal, 
                     const double area,
                     const RangeType& midEn,
                     FluxRangeType& gLeft) const;


private:
  const ModelType& model_;
  double acpenalty_;
};


template< class Model >
void AllenCahnJacobianFlux<Model>
::boundaryFlux(const DomainType& normal,
               const double area,
               const RangeType& midEn,
               FluxRangeType& gLeft) const
  {
    double laplaceFlux(0.);
    gLeft=0;
    for(int ii = 0; ii < dimDomain ; ++ii )
      {
        laplaceFlux=normal[ ii ];//*0.5;
        //std::cout<<"normal["<<ii<<"]="<<laplaceFlux<<"\n";
        // F_{3.3}(\sigma^+)\cdot n
        gLeft[ tau ][ sigma + ii ]=-model_.delta()*laplaceFlux;
      } 
  }



template< class Model >
void AllenCahnJacobianFlux<Model>
::numericalFlux( const DomainType& normal,
                const double area,
                const double penaltyFactor,
                const RangeType& vuEn,
                const RangeType& vuNb,
                const RangeType& addEn,
                const RangeType& addNb,
                FluxRangeType& gLeft,
                FluxRangeType& gRight) const
  {
      RangeType jump ;
      double integrationElement=normal.two_norm();
  
      gLeft=0.;
      gRight=0.;

      jump = vuEn;
      jump-= vuNb;
 
      //rho-------------------------------------------------------------
 
      double laplaceFluxLeft(0.),laplaceFluxRight(0.);

      for(int ii = 0; ii < dimDomain ; ++ii )
        {
          //F_{3.1}
          //-(\phi^+-\phi^-)*n[i]*v[i]^+*0.5
          gLeft[ phi ][ phi ]-=normal[ ii ]*addEn[ 1 + ii ]*0.5;
          gRight[ phi ][ phi ]+=normal[ ii ]*addEn[ 1 + ii ]*0.5;
          
          //tau
          //F_{3.2}
          //(\sigma^+-\sigma^-)\cdot n * 0.5
          laplaceFluxLeft=normal[ ii ]*0.5;
          laplaceFluxRight=normal[ ii ]*-0.5;
          
          double penaltyTerm=acpenalty_;
          penaltyTerm*=penaltyFactor;
          penaltyTerm*=integrationElement;
   
          //tau-------------------------------------------------------------
          gLeft[ tau ][ sigma+ii ]=-model_.delta()*laplaceFluxLeft;
          gRight[ tau ][ sigma+ii ]=-model_.delta()*laplaceFluxRight;
          //factor 0.5 comes from linerarization of 1/2(phi^(n+1)+phi^n)
          gLeft[ tau ][ phi ]=-model_.delta()*penaltyTerm;
          gRight[ tau ][ phi ]=model_.delta()*penaltyTerm; 
          
          //F_4
          //(\phi^+-\phi^-)\normal*0.5
          gLeft[ sigma + ii ][ phi ]=normal[ ii ]*0.5;
          gRight[ sigma + ii ][ phi ]=normal[ ii ]*-0.5;
        }
  }


template<class Model>
class NavierStokesJacobianFlux
{
  //typedefs

 typedef Model ModelType;
 
  enum{ dimDomain=Model::dimDomain};
  enum{ dimRange=Model::dimNvStRange};
 
  typedef typename Dune::FieldVector<double,dimRange> RangeType;
  typedef typename Dune::FieldVector<double,dimDomain> DomainType;
  typedef typename Dune::FieldMatrix<double,dimRange,dimDomain> JacobianRangeType;
  typedef typename Dune::FieldMatrix<double,dimRange,dimRange> FluxRangeType;
  typedef NvStFilter<RangeType> Filter; 
  typedef AcFilter<RangeType> AddFilter;

  
public:
  NavierStokesJacobianFlux(const ModelType& model):
    model_(model),
    penalty_( Dune::Fem::Parameter::getValue< double >( "phasefield.penalty" )),
    acpenalty_( Dune::Fem::Parameter::getValue< double >( "phasefield.acpenalty" )),
    switchIP_(Dune::Fem::Parameter::getValue<int>("phasefield.ipswitch",1)),
    numVisc_(Dune::Fem::Parameter::getValue<double>("phasefield.addvisc",0)),
    numViscMu_(Dune::Fem::Parameter::getValue<double>("phasefield.muvisc",0)),
    imexFactor_(Dune::Fem::Parameter::getValue<double>("phasefield.IMEX"))
    {
    }


  void  numericalFlux ( const DomainType& normal,
                        const double area,
                        const double penaltyFactor,
                        const RangeType& vuEn,
                        const RangeType& vuNb,
                        const RangeType& vuEnAdd,
                        const RangeType& vuNbAdd,
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

  template< class ScalarType , class JacobianType ,class DiffusionTensorType, class DiffusionVectorType>
  void  scalar2vectorialBoundaryFlux( const DomainType& normal,
                                      const double penaltyFactor,
                                      const ScalarType& phiEn,
                                      const JacobianType& dphiEn,
                                      DiffusionTensorType& aLeft,
                                      DiffusionVectorType& bLeft) const;

  void diffusionFlux ( const DomainType& normal,
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
                     const double area,
                     const RangeType& midEn,
                     FluxRangeType& gLeft) const;

  void outFlowFlux( const DomainType& normal,
                    const double area,
                    const RangeType& midEn,
                    FluxRangeType& gLeft) const;

  void diffusionBoundaryFlux( const DomainType& normal,
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
  const double numViscMu_;
  double imexFactor_;
};


template< class Model >
void NavierStokesJacobianFlux<Model>
::boundaryFlux(const DomainType& normal,
               const double area,
               const RangeType& midEn,
               FluxRangeType& gLeft) const
  {
    double vNormalEn(0.);
    for(int ii = 0; ii < dimDomain ; ++ii )
      {
        vNormalEn+=midEn[1 + ii]*normal[ ii ];
        gLeft[ 0 ][ 1 + ii ] = -1*normal[ii]*midEn[0];
      } 
    gLeft[ 0 ][ 0 ]=-1*vNormalEn;
  }

template< class Model >
void NavierStokesJacobianFlux<Model>
::outFlowFlux(const DomainType& normal,
              const double area,
              const RangeType& midEn,
              FluxRangeType& gLeft) const
  {
      
  }


template< class Model >
void  NavierStokesJacobianFlux<Model>
::numericalFlux( const DomainType& normal,
                const double area,
                const double penaltyFactor,
                const RangeType& vuEn,
                const RangeType& vuNb,
                const RangeType& vuEnAdd,
                const RangeType& vuNbAdd,
                FluxRangeType& gLeft,
                FluxRangeType& gRight) const
  {
      RangeType valEn,valNb,jump,mean;
  
      gLeft=0.;
      gRight=0.;

      jump = vuEn;
      jump-= vuNb;
 
      
      
      //rho-------------------------------------------------------------
 
      double vNormalEn(0),vNormalNb(0.);

      for(int ii = 0; ii < dimDomain ; ++ii )
        {
          //v^+\cdot n^+
          vNormalEn+=vuEn[ 1 + ii ]*normal[ ii ];
          //v^-\cdot n^+
          vNormalNb+=vuNb[ 1 + ii ]*normal[ ii ];
          gLeft[ 0 ][ 1 + ii ]  = -0.5*normal[ ii ]*vuEn[ 0 ];
          gRight[ 0 ][ 1 + ii ] =  0.5*normal[ ii ]*vuNb[ 0 ];

 
           //(\rho_j,\mu*,\v_i)
          gLeft[ 1 + ii ][ 0 ]=-1*jump[ dimDomain + 1 ]*normal[ ii ];

          //(\rho*,\mu_j,\v_i)
          gLeft[ 1 + ii ][ dimDomain + 1 ]=-1*normal[ ii ]*vuEn[ 0 ];
          gRight[ 1 + ii ][ dimDomain + 1 ]=normal[ ii ]*vuEn[ 0 ]; 

          gLeft[ 1 + ii ]*=0.5;
          gRight[ 1 + ii ]*=0.5;

 
 
       }
  
 
      //F_1=-0.5*( \rho^+*v^+\cdot n^+ -\rho^-*v-\cdot n^+)
      gLeft[ 0 ][ 0 ]=-0.5*vNormalEn;
      //from linerarization  vuMid=0.5(vu+vuOld), d vuMid/d vu=0.5
      gRight[ 0 ][ 0 ]=0.5*vNormalNb;
      //[[mu]]
      gLeft[ 0 ][ dimDomain+1 ]=0.5*numViscMu_*area;
      gRight[ 0 ][ dimDomain+1 ]=-0.5*numViscMu_*area;
      //[[rho]]
      gLeft[ 0 ][ 0 ]+=0.5*numVisc_*area;
      gRight[ 0 ][0]+=-0.5*numVisc_*area;
      //----------------------------------------------------------------

  }
template< class Model >
template< class ScalarType , class JacobianType ,class DiffusionTensorType, class DiffusionVectorType>
void NavierStokesJacobianFlux<Model>
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

      bLeft[ii][ii]+=penalty_*penaltyFactor*integrationElement*phiEn[ 0 ];
      bRight[ii][ii]-=penalty_*penaltyFactor*integrationElement*phiNb[ 0 ];
    }
  
}



template< class Model >
template< class ScalarType , class JacobianType ,class DiffusionTensorType, class DiffusionVectorType>
void NavierStokesJacobianFlux<Model>
:: scalar2vectorialBoundaryFlux( const DomainType& normal,
                                 const double penaltyFactor,
                                  const ScalarType& phiEn,
                                  const JacobianType& dphiEn,
                                  DiffusionTensorType& aLeft,
                                  DiffusionVectorType& bLeft) const
{
  double integrationElement=normal.two_norm();
  //[\phiEn]\otimes n ; [\phiNb]\otimes n
  JacobianType jumpNormalLeft(0.);

  for(int j=0; j<dimDomain; ++j)
    {
      jumpNormalLeft[0][j]=-1*phiEn[0]*normal[j];
    }

  DiffusionTensorType bbLeft;
  jumpNormalLeft*=switchIP_;

  model_.scalar2vectorialDiffusion(jumpNormalLeft ,aLeft);
  model_.scalar2vectorialDiffusion( dphiEn, bbLeft);

  //loop over components
  for( int ii = 0 ; ii < dimDomain ; ++ ii)
    {
      bbLeft[ ii ]*=-1;
      bbLeft[ ii ].mv( normal , bLeft[ ii ]);
      bLeft[ii][ii]+=penalty_*penaltyFactor*integrationElement*phiEn[ 0 ];
    }
}




template< class Model >
void NavierStokesJacobianFlux<Model>
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
  abort();
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
void NavierStokesJacobianFlux<Model>
:: diffusionBoundaryFlux( const DomainType& normal,
                  const double penaltyFactor,
                  const RangeType& uEn,
                  const JacobianRangeType& duEn,
                  RangeType& value,
                  JacobianRangeType& dvalue) const
{
  abort(); 
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


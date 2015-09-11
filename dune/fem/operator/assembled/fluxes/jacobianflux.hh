#ifndef PHASEFIELD_JACOBIAN_FLUX_HH
#define PHASEFIELD_JACBIAN_FLUX_HH

#include <dune/common/fvector.hh>

#include "../phasefieldfilter.hh"
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
    numVisc_(Dune::Fem::Parameter::getValue<double>("phasefield.addvisc",0)),
    numViscMu_(Dune::Fem::Parameter::getValue<double>("phasefield.muvisc",0)),
    theta_(Dune::Fem::Parameter::getValue<double>("phasefield.mixed.theta")),
    factorImp_(0.5*(1+theta_))
    {

    }


  void  numericalFlux ( const DomainType& normal,
                        const double area,
                        const double penaltyFactor,
                        const RangeType& vuEn,
                        const RangeType& vuNb,
                        const RangeType& vuEnImEx,
                        const RangeType& vuNbImEx,
                        FluxRangeType& gLeft,
                        FluxRangeType& gRight) const;

  template< class ScalarType, class JacobianType, class DiffusionTensorType, class DiffusionVectorType>
  void scalar2vectorialDiffusionFlux( const DomainType& normal,
                                      const double penaltyFactor,
                                      const RangeType& vuEn,
                                      const RangeType& vuNb,
                                      const ScalarType& phiEn,
                                      const ScalarType& phiNb,
                                      const JacobianType& dphiEn,
                                      const JacobianType& dphiNb,
                                      DiffusionTensorType& aLeft,
                                      DiffusionTensorType& aRight,
                                      DiffusionVectorType& bLeft,
                                      DiffusionVectorType& bRight) const;
  
  template< class ScalarType >
  void diffPhiDiffusionFlux( const DomainType& normal,
                             const double penaltyFactor,
                             const RangeType& vuEn,
                             const RangeType& vuNb,
                             const JacobianRangeType& duEn,
                             const JacobianRangeType& duNb,
                             const ScalarType& phiEn,
                             const ScalarType& phiNb,
                             JacobianRangeType& aLeft,
                             JacobianRangeType& aRight,
                             RangeType& bLeft,
                             RangeType& bRight) const;


  template< class ScalarType , class JacobianType ,class DiffusionTensorType, class DiffusionVectorType>
  void  scalar2vectorialBoundaryFlux( const DomainType& normal,
                                      const double penaltyFactor,
                                      const RangeType& vuEn,
                                      const ScalarType& phiEn,
                                      const JacobianType& dphiEn,
                                      DiffusionTensorType& aLeft,
                                      DiffusionVectorType& bLeft) const;
  template< class ScalarType >
  void diffPhiBoundaryFlux( const DomainType& normal,
                            const double penaltyFactor,
                            const RangeType& vuEn,
                            const JacobianRangeType& duEn,
                            const ScalarType& phiEn,
                            JacobianRangeType& aLeft,
                            RangeType& bLeft) const;


    void boundaryFlux( const DomainType& normal, 
                     const double area,
                     const RangeType& midEn,
                     FluxRangeType& gLeft) const;

  void outFlowFlux( const DomainType& normal,
                    const double area,
                    const RangeType& midEn,
                    FluxRangeType& gLeft) const;



private:
  const ModelType& model_;
  double penalty_;
  double acpenalty_;
  const int switchIP_; 
  const double numVisc_;
  const double numViscMu_;
  const double theta_;
  const double factorImp_;
};


template< class Model >
void JacobianFlux<Model>
::boundaryFlux(const DomainType& normal,
               const double area,
               const RangeType& midEn,
               FluxRangeType& gLeft) const
  {
    double vNormalEn(0.);
    double laplaceFlux(0.);

    for(int ii = 0; ii < dimDomain ; ++ii )
      {
        vNormalEn+=midEn[1 + ii]*normal[ ii ];
        gLeft[ 0 ][ 1 + ii ] = -factorImp_*normal[ii]*midEn[0];
        laplaceFlux=normal[ ii ]*factorImp_;
#if LAMBDASCHEME
        // F_{3.3}(\lambda^+)\cdot n
        gLeft[ dimDomain + 3 ][ 2*dimDomain + 4 +ii ]=-model_.delta()*laplaceFlux;
#else
        // F_{3.3}(\sigma^+)\cdot n
        gLeft[ dimDomain + 3 ][ dimDomain + 4 +ii ]=-model_.delta()*laplaceFlux;
#endif
      } 
    gLeft[ 0 ][ 0 ]=-factorImp_*vNormalEn;
  }

template< class Model >
void JacobianFlux<Model>
::outFlowFlux(const DomainType& normal,
              const double area,
              const RangeType& midEn,
              FluxRangeType& gLeft) const
  {
    double laplaceFlux(0.);
  
    for(int ii = 0; ii < dimDomain ; ++ii )
      {
        //gLeft[ 0 ][ 1 + ii ] = -0.5*normal[ii]*midEn[0];
        laplaceFlux=normal[ ii ]*0.5;
#if LAMBDASCHEME
        // F_{3.3}(\lambda^+)\cdot n
        gLeft[ dimDomain + 3 ][ 2*dimDomain + 4 +ii ]=-model_.delta()*laplaceFlux;
#else
        // F_{3.3}(\sigma^+)\cdot n
        gLeft[ dimDomain + 3 ][ dimDomain + 4 +ii ]=-model_.delta()*laplaceFlux;
#endif
      }
  }


template< class Model >
void  JacobianFlux<Model>
::numericalFlux( const DomainType& normal,
                const double area,
                const double penaltyFactor,
                const RangeType& vuEnMid,
                const RangeType& vuNbMid,
                const RangeType& vuEnImEx,
                const RangeType& vuNbImEx,
                FluxRangeType& gLeft,
                FluxRangeType& gRight) const
  {
      RangeType valEn,valNb,midEn(0.), midNb(0.),jump,jumpImEx,mean, jumpPhi(0.);
      double integrationElement=normal.two_norm();
  
      gLeft=0.;
      gRight=0.;

      midEn=vuEnMid;
      midNb=vuNbMid;
      jump = midEn;
      jump-= midNb;
      jumpImEx = vuEnImEx;
      jumpImEx-= vuNbImEx;
 
      mean = midEn ;
      mean+= midNb;
      mean*=0.5;
      
      
      //rho-------------------------------------------------------------
 
      double vNormalEn(0),vNormalNb(0.);
      double laplaceFluxLeft(0.),laplaceFluxRight(0.);

      for(int ii = 0; ii < dimDomain ; ++ii )
        {
          //v^+\cdot n^+
          vNormalEn+=midEn[ 1 + ii ]*normal[ ii ];
          //v^-\cdot n^+
          vNormalNb+=midNb[ 1 + ii ]*normal[ ii ];
          gLeft[ 0 ][ 1 + ii ]  = -0.5*factorImp_*normal[ ii ]*midEn[ 0 ];
          gRight[ 0 ][ 1 + ii ] =  0.5*factorImp_*normal[ ii ]*midNb[ 0 ];


           //(\rho_j,\mu*,\v_i)
          gLeft[ 1 + ii ][ 0 ]=-factorImp_*jump[ dimDomain + 2 ]*normal[ ii ];
          //(\rho*,\mu_j,\v_i)
          gLeft[ 1 + ii ][ dimDomain + 2 ]=-factorImp_*normal[ ii ]*midEn[ 0 ];
          gRight[ 1 + ii ][ dimDomain + 2 ]=factorImp_*normal[ ii ]*midEn[ 0 ];
          //(\phi*,\tau_j,v_i)
          gLeft[ 1 + ii ][ dimDomain + 3 ]=jump[ dimDomain + 1 ]*normal[ ii ]*factorImp_;
          //(\phi_j,\tau*,v_i)
          gLeft[ 1 + ii ][ dimDomain + 1] =factorImp_*normal[ ii ]*vuEnImEx[ dimDomain + 3 ];
          gRight[ 1 + ii ][ dimDomain + 1 ]=-factorImp_*normal[ ii ]*vuEnImEx[ dimDomain + 3 ];

          gLeft[ 1 + ii ]*=0.5;
          gRight[ 1 + ii ]*=0.5;

          //F_{3.1}

          //-(\phi^+-\phi^-)*n[i]*v[i]^+*0.5
          gLeft[ dimDomain + 1 ][ dimDomain + 1 ]+=normal[ ii ]*midEn[ 1 + ii ]*-0.5*factorImp_;
          gRight[ dimDomain + 1 ][ dimDomain +1 ]+=normal[ ii ]*midEn[ 1 + ii ]*0.5*factorImp_;

          gLeft[ dimDomain + 1 ][ 1 + ii ]=jump[ dimDomain + 1 ]*normal[ ii ]*-0.5*factorImp_;
          //tau
          //F_{3.2}
          //(\sigma^+-\sigma^-)\cdot n * 0.5
          laplaceFluxLeft=normal[ ii ]*0.5*factorImp_;
          laplaceFluxRight=normal[ ii ]*-0.5*factorImp_;
          
          double penaltyTerm=acpenalty_;
          penaltyTerm*=penaltyFactor;
          penaltyTerm*=integrationElement;
   
          //tau-------------------------------------------------------------
#if LAMBDASCHEME
          gLeft[ dimDomain + 3 ][2*dimDomain + 4 + ii ]=-model_.delta()*laplaceFluxLeft;
          gRight[ dimDomain + 3 ][2*dimDomain + 4 + ii ]=-model_.delta()*laplaceFluxRight;
#else
          gLeft[ dimDomain + 3 ][ dimDomain + 4+ii ]=-model_.delta()*laplaceFluxLeft;
          gRight[ dimDomain + 3 ][dimDomain + 4+ii ]=-model_.delta()*laplaceFluxRight;
#endif
          gLeft[ dimDomain + 3 ][ dimDomain + 1]=-model_.delta()*factorImp_*penaltyTerm;
          gRight[ dimDomain + 3 ][ dimDomain + 1]=model_.delta()*factorImp_*penaltyTerm; 

          //F_4
          //(\phi^+-\phi^-)\normal*0.5
          gLeft[ dimDomain + 4 + ii ][ dimDomain + 1 ]=normal[ ii ]*0.5;
          gRight[ dimDomain + 4 + ii ][ dimDomain + 1 ]=normal[ ii ]*-0.5;
       }
  

      //F_1=-0.5*( \rho^+*v^+\cdot n^+ -\rho^-*v-\cdot n^+)
      gLeft[ 0 ][ 0 ]=-0.5*factorImp_*vNormalEn;
      gRight[ 0 ][ 0 ]=0.5*factorImp_*vNormalNb;
       //[[mu]]
      gLeft[ 0 ][ dimDomain+2 ]=factorImp_*numViscMu_*area;
      gRight[ 0 ][ dimDomain+2 ]=-factorImp_*numViscMu_*area;
      //[[rho]]
      gLeft[ 0 ][ 0 ]+=factorImp_*numVisc_*area;
      gRight[ 0 ][ 0 ]-=factorImp_*numVisc_*area;

      //----------------------------------------------------------------

  }
template< class Model >
template< class ScalarType , class JacobianType ,class DiffusionTensorType, class DiffusionVectorType>
void JacobianFlux<Model>
:: scalar2vectorialDiffusionFlux( const DomainType& normal,
                                  const double penaltyFactor,
                                  const RangeType& vuEn,
                                  const RangeType& vuNb,
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
  JacobianType jumpNormalEn(0.),jumpNormalNb(0.);


  for(int j=0; j<dimDomain; ++j)
    { 
      jumpNormalEn[0][j]=-0.5*phiEn[0]*normal[j];
      jumpNormalNb[0][j]=-0.5*phiNb[0]*normal[j];
    }

  DiffusionTensorType bbLeft,bbRight;
  jumpNormalEn*=switchIP_;
  jumpNormalNb*=switchIP_;
  
  model_.scalar2vectorialDiffusion( vuEn , jumpNormalEn , aLeft  );
  model_.scalar2vectorialDiffusion( vuEn , jumpNormalNb , aRight );
  
  model_.scalar2vectorialDiffusion( vuEn , dphiEn, bbLeft);
  model_.scalar2vectorialDiffusion( vuNb , dphiNb, bbRight);
  
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
template< class ScalarType >
void JacobianFlux<Model>::diffPhiDiffusionFlux ( const DomainType& normal,
                                                const double penaltyFactor,
                                                const RangeType& vuEn,
                                                const RangeType& vuNb,
                                                const JacobianRangeType& duEn,
                                                const JacobianRangeType& duNb,
                                                const ScalarType& phiEn,
                                                const ScalarType& phiNb,
                                                JacobianRangeType& aLeft,
                                                JacobianRangeType& aRight,
                                                RangeType& bLeft,
                                                RangeType& bRight) const
{
  double integrationElement=normal.two_norm();
  //[u]\otimes n 
  JacobianRangeType jumpNormal(0);
  auto jump=vuEn;
  jump-=vuNb;
  for( int i = 0 ; i < dimDomain ; ++i)
    for(  int j = 0 ; j<dimDomain ; ++j )
    {
      jumpNormal[i+1][j]=-0.5*jump[i+1]*normal[j];
    }

  //DiffusionTensorType bbLeft,bbRight;
  jumpNormal*=switchIP_;

  JacobianRangeType aaLeft(0.),bbLeft(0.),bbRight(0.);
  //dD/dphi(\bar{\phi^+},[u]\otimes\n})*\psi_phi^+
  model_.diffusionprime( vuEn , jumpNormal , aLeft  );
  aLeft*=phiEn[0];
  aRight=0.;

  model_.diffusionprime( vuEn , duEn, bbLeft);
  model_.diffusionprime( vuNb , duNb, bbRight);
  bbLeft*=-0.5*phiEn[0];
  bbLeft.mv( normal , bLeft);
  bbRight*=-0.5*phiNb[0];
  bbRight.mv( normal , bRight);

}

template< class Model >
template< class ScalarType , class JacobianType ,class DiffusionTensorType, class DiffusionVectorType>
void JacobianFlux<Model>
:: scalar2vectorialBoundaryFlux( const DomainType& normal,
                                 const double penaltyFactor,
                                 const RangeType& vuEn,
                                 const ScalarType& phiEn,
                                 const JacobianType& dphiEn,
                                 DiffusionTensorType& aLeft,
                                 DiffusionVectorType& bLeft ) const
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

  model_.scalar2vectorialDiffusion( vuEn , jumpNormalLeft , aLeft );
  model_.scalar2vectorialDiffusion( vuEn , dphiEn, bbLeft);

  //loop over components
  for( int ii = 0 ; ii < dimDomain ; ++ ii)
    {
      bbLeft[ ii ]*=-1;
      bbLeft[ ii ].mv( normal , bLeft[ ii ]);
      bLeft[ii][ii]+=penalty_*penaltyFactor*integrationElement*phiEn[ 0 ];
    }
}



template< class Model >
template< class ScalarType >
void JacobianFlux<Model>::diffPhiBoundaryFlux ( const DomainType& normal,
                                                const double penaltyFactor,
                                                const RangeType& vuEn,
                                                const JacobianRangeType& duEn,
                                                const ScalarType& phiEn,
                                                JacobianRangeType& aLeft,
                                                RangeType& bLeft) const
{
  double integrationElement=normal.two_norm();
  //[u]\otimes n 
  JacobianRangeType jumpNormal(0);
  auto jump=vuEn;
  for( int i = 0 ; i < dimDomain ; ++i)
    for(  int j = 0 ; j<dimDomain ; ++j )
    {
      jumpNormal[i+1][j]=1*jump[i+1]*normal[j];
    }

  //DiffusionTensorType bbLeft,bbRight;
  jumpNormal*=switchIP_;

  JacobianRangeType aaLeft(0.),bbLeft(0.),bbRight(0.);
  //dD/dphi(\bar{\phi^+},[u]\otimes\n})*\psi_phi^+
  model_.diffusionprime( vuEn , jumpNormal , aLeft  );
  aLeft*=phiEn[0];

  model_.diffusionprime( vuEn , duEn, bbLeft);
  bbLeft*=-1*phiEn[0];
  bbLeft.mv( normal , bLeft);

}



#endif


#ifndef PHASEFIELD_SPLITFLUX_HH
#define PHASEFIELD_SPLITFLUX_HH

#include <dune/common/fvector.hh>

#include "../splitfilter.hh"

template<class Model>
class AcMixedFlux
{
  //typedefs

 typedef Model ModelType;
 
  enum{ dimDomain=Model::dimDomain};
  enum{ dimRange=Model::dimAcRange};
 
  typedef typename Dune::FieldVector<double,dimRange> RangeType;
  typedef typename Dune::FieldVector<double,dimDomain> DomainType;
  typedef typename Dune::FieldMatrix<double,dimRange,dimDomain> JacobianRangeType;
  typedef  AcFilter<RangeType> Filter; 
  typedef  NvStFilter<RangeType> AddFilter; 
  
public:
  AcMixedFlux(const ModelType& model):
    model_(model),
    acpenalty_( Dune::Fem::Parameter::getValue< double >( "phasefield.acpenalty" ))
    {
    }


  void numericalFlux( const DomainType& normal, 
                        const double area,
                        const double penaltyFactor,
                        const RangeType& vuEn,
                        const RangeType& vuN,
                        const RangeType& addEn,
                        const RangeType& addNb,
                        RangeType& gLeft,
                        RangeType& gRight) const;



  void boundaryFlux( const DomainType& normal, 
                     const double penaltyFactor,                
                     const RangeType& vuEn,  
                     const RangeType& vuMidEn,
                     const RangeType& addEn,
                      RangeType& gLeft) const;

  void outFlowFlux( const DomainType& normal,
                    const double penaltyFactor,
                    const RangeType& vuEn,
                    const RangeType& vuMidEn,
                    const RangeType& addEn,
                    RangeType& gLeft) const;

private:
  const ModelType& model_;
  double acpenalty_;
};

////////////Implementation of AcMixedFlux///////////////
template< class Model >
void AcMixedFlux<Model>
::boundaryFlux(const DomainType& normal,
               const double penaltyFactor,
               const RangeType& vuEn,
               const RangeType& vuMidEn,
               const RangeType& addEn,
               RangeType& gLeft) const
  {
    RangeType valEn,addValEn;
    valEn=vuEn;
    addValEn=addEn; 

    gLeft=0.;

   //phi-------------------------------------------------------------

    double laplaceFlux(0.);

    for(int i = 0; i<dimDomain;++i)
      {
#if RHOMODEL
#if LAMBDASCHEME
        //(\lambda^+-\lambda^-)\cdot n * 0.5
        laplaceFlux-=(Filter::alpha(valEn,i))*normal[i];

#else
        //(h2(rho^+)\sigma^+-h2(rho^-)\sigma^-)\cdot n * 0.5
        laplaceFlux-=Filter::sigma(valEn,i)*model_.h2(AddFilter::rho(addValEn))*normal[i];
#endif

#else
        //F_{3.2}
        //(\sigma^+-\sigma^-)\cdot n * 0.5
        laplaceFlux-=Filter::sigma(valEn,i)*normal[i];
#endif
      }

    //----------------------------------------------------------------

    //tau-------------------------------------------------------------

    Filter::tau(gLeft)=model_.delta()*laplaceFlux;

    //----------------------------------------------------------------

    //sigma-----------------------------------------------------------

    for(int i = 0; i<dimDomain;++i)
      {
          //(\phi^+-\phi^-)\normal*0.5
          Filter::sigma(gLeft,i)=0;//(Filter::phi(valEn))*normal[i]*0.5;
      }
    //----------------------------------------------------------------

  }

template< class Model >
void AcMixedFlux<Model>
::outFlowFlux(const DomainType& normal,
               const double penaltyFactor,
               const RangeType& vuEn,
               const RangeType& vuMidEn,
               const RangeType& addEn,
               RangeType& gLeft) const
  {
    RangeType valEn,addValEn;
    valEn=vuEn;
    addValEn=addEn;

    gLeft=0.;

    //phi-------------------------------------------------------------
  
    double laplaceFlux(0.);

    for(int i = 0; i<dimDomain;++i)
      {
#if RHOMODEL
#if LAMBDASCHEME
        //(\lambda^+-\lambda^-)\cdot n * 0.5
        laplaceFlux-=(Filter::alpha(valEn,i))*normal[i];
        
#else
        //(h2(rho^+)\sigma^+-h2(rho^-)\sigma^-)\cdot n * 0.5
        laplaceFlux-=Filter::sigma(valEn,i)*model_.h2(AddFilter::rho(addValEn))*normal[i];
#endif

#else
        //F_{3.2}
        //(\sigma^+-\sigma^-)\cdot n * 0.5
        laplaceFlux-=Filter::sigma(valEn,i)*normal[i];
#endif
      }  
  
    //----------------------------------------------------------------

    //tau-------------------------------------------------------------

    Filter::tau(gLeft)=model_.delta()*laplaceFlux;

    //----------------------------------------------------------------

    //sigma-----------------------------------------------------------

    for(int i = 0; i<dimDomain;++i)
      {   
          //(\phi^+-\phi^-)\normal*0.5
          Filter::sigma(gLeft,i)=0;//(Filter::phi(valEn))*normal[i]*0.5;
      } 
    //----------------------------------------------------------------
  
  }







template< class Model >
void AcMixedFlux<Model>
::numericalFlux( const DomainType& normal,
                const double area,              
                const double penaltyFactor,
                const RangeType& vuEn, // needed for calculation of sigma which is fully implicit
                const RangeType& vuNb, // needed for calculation of sigma which is fully implicit
                const RangeType& addEn,
                const RangeType& addNb,
                RangeType& gLeft,
                RangeType& gRight) const
  {
    RangeType valEn,valNb,jump,addValEn,addValNb;
    valEn=vuEn;
    valNb=vuNb;
    addValEn=addEn;
    addValNb=addNb;
    double integrationElement=normal.two_norm();
 
    gLeft=0.;
    gRight=0.;

    jump = valEn;
    jump-= valNb;
   
    double laplaceFlux(0.);
    //phi-------------------------------------------------------------
    for(int i = 0; i<dimDomain;++i)
      {
        //F_{3.1}
     
        //-(\phi^+-\phi^-)*n[i]*v[i]*0.5 
        Filter::phi(gLeft)-=Filter::phi(jump)*normal[i]*AddFilter::velocity(addValEn,i)*0.5;
        Filter::phi(gRight)-=Filter::phi(jump)*normal[i]*AddFilter::velocity(addValNb,i)*0.5;
          
        //tau 
        //F_{3.2}
#if RHOMODEL
#if LAMBDASCHEME
        //(\lambda^+-\lambda^-)\cdot n * 0.5
        laplaceFlux+=(Filter::alpha(valEn,i)-Filter::alpha(valNb,i))*normal[i]*0.5;
#else
        //(h2(rho^+)\sigma^+-h2(rho^-)\sigma^-)\cdot n * 0.5
        laplaceFlux+=(Filter::sigma(valEn,i)*model_.h2(AddFilter::rho(addValEn))-Filter::sigma(valNb,i)*model_.h2(AddFilter::rho(addValNb)))*normal[i]*0.5;
#endif
#else
        //F_{3.2}
        //(\sigma^+-\sigma^-)\cdot n * 0.5
        laplaceFlux+=Filter::sigma(jump,i)*normal[i]*0.5;
#endif        
      }  
    //----------------------------------------------------------------
    //tau-------------------------------------------------------------
 
    double penaltyTerm=acpenalty_;
    penaltyTerm*=penaltyFactor;
    penaltyTerm*=integrationElement;
    penaltyTerm*=Filter::phi(jump);

    Filter::tau(gLeft)-=model_.delta()*(laplaceFlux+penaltyTerm);
    Filter::tau(gRight)=Filter::tau( gLeft );
   
    //----------------------------------------------------------------

    //sigma-----------------------------------------------------------
    for(int i = 0; i<dimDomain;++i)
      {  
        //F_4
        //(\phi^+-\phi^-)\normal*0.5
        Filter::sigma(gLeft,i)=(Filter::phi(valEn)-Filter::phi(valNb))*normal[i]*0.5;
        Filter::sigma(gRight,i)=Filter::sigma( gLeft,i);

      } 
    //----------------------------------------------------------------
}
////////////Implementation of AcMixedFlux///////////////



template<class Model>
class NvStMixedFlux
{
  //typedefs

 typedef Model ModelType;
 
  enum{ dimDomain=Model::dimDomain};
  enum{ dimRange=Model::dimNvStRange};
 
  typedef typename Dune::FieldVector<double,dimRange> RangeType;
  typedef typename Dune::FieldVector<double,dimDomain> DomainType;
  typedef typename Dune::FieldMatrix<double,dimRange,dimDomain> JacobianRangeType;
  typedef NvStFilter<RangeType> Filter; 
  typedef AcFilter<RangeType> AddFilter; 
  
public:
  NvStMixedFlux(const ModelType& model):
    model_(model),
    penalty_( Dune::Fem::Parameter::getValue< double >( "phasefield.penalty" )),
    switchIP_(Dune::Fem::Parameter::getValue<int>("phasefield.ipswitch",1)),
    numVisc_(Dune::Fem::Parameter::getValue<double>("phasefield.addvisc",0)),
    numViscMu_(Dune::Fem::Parameter::getValue<double>("phasefield.muvisc",0))
    {
    }


   void numericalFlux( const DomainType& normal, 
                        const double area,
                        const double penaltyFactor,
                        const RangeType& vuEn,
                        const RangeType& vuN,
                        const RangeType& addEn,
                        const RangeType& addNb,
                        RangeType& gLeft,
                        RangeType& gRight) const;


   void diffusionFlux( const DomainType& normal,
                        const double penaltyFactor,
                        const RangeType& uEn,
                        const RangeType& uNb,
                        const JacobianRangeType& duEn,
                        const JacobianRangeType& duNb,
                        RangeType& value,
                        JacobianRangeType& dvalue) const;

  void boundaryFlux( const DomainType& normal, 
                       const double penaltyFactor,                
                       const RangeType& vuEn,  
                       const RangeType& vuMidEn,
                       RangeType& gLeft) const;

  void outFlowFlux( const DomainType& normal,
                      const double penaltyFactor,
                      const RangeType& vuEn,
                      const RangeType& vuMidEn,
                      RangeType& gLeft) const;

  void diffusionBoundaryFlux( const DomainType& normal,
                                const double penaltyFactor,
                                const RangeType& uEn,
                                const JacobianRangeType& duEn,
                                RangeType& value,
                                JacobianRangeType& dvalue) const;

  void diffusionOutFlowFlux( RangeType& value ,JacobianRangeType& dvalue) const;

private:
  const ModelType& model_;
  double penalty_;
  int switchIP_;
  const double numVisc_;
  const double numViscMu_;
};


template< class Model >
void NvStMixedFlux<Model>
::boundaryFlux(const DomainType& normal,
               const double penaltyFactor,
               const RangeType& vuEn,
               const RangeType& vuMidEn,
               RangeType& gLeft) const
  {
    RangeType valEn;
    valEn=vuEn;


    gLeft=0.;

    double vNormalEn(0.);

    //rho-------------------------------------------------------------
  
    for(int i = 0; i<dimDomain;++i)
      {
        vNormalEn+=Filter::velocity(valEn,i)*normal[i];
      }

    Filter::rho(gLeft)=-1*vNormalEn*Filter::rho(valEn);
  
    //----------------------------------------------------------------
    
    //v---------------------------------------------------------------
  
    for(int i = 0; i<dimDomain;++i)
      {
        Filter::velocity(gLeft,i)=0;
      } 
    //----------------------------------------------------------------


  }

template< class Model >
void NvStMixedFlux<Model>
::outFlowFlux(const DomainType& normal,
               const double penaltyFactor,
               const RangeType& vuEn,
               const RangeType& vuMidEn,
               RangeType& gLeft) const
  {
    gLeft=0.;
  }







template< class Model >
void NvStMixedFlux<Model>
::numericalFlux( const DomainType& normal,
                const double area,              
                const double penaltyFactor,
                const RangeType& vuEn, // needed for calculation of sigma which is fully implicit
                const RangeType& vuNb, // needed for calculation of sigma which is fully implicit
                const RangeType& addEn,
                const RangeType& addNb,
                RangeType& gLeft,
                RangeType& gRight) const
  {
    RangeType valEn,valNb,jump,mean,jumpNew,jumpAdd,addValEn,addValNb;
    valEn=vuEn;
    valNb=vuNb;
    addValEn=addEn;
    addValNb=addNb;
 
    gLeft=0.;
    gRight=0.;


    jump = vuEn;
    jump-= vuNb;
    jumpAdd=addEn;
    jumpAdd-=addNb;
   
    //rho-------------------------------------------------------------
 
    double vNormalEn(0.),vNormalNb(0.);
    
    for(int i = 0; i<dimDomain;++i)
      {
        //v^+\cdot n^+
        vNormalEn+=Filter::velocity(valEn,i)*normal[i];
        //v^-\cdot n^+
        vNormalNb+=Filter::velocity(valNb,i)*normal[i];
      }
    
    //F_1=-0.5*( \rho^+*v^+\cdot n^+ -\rho^-*v-\cdot n^+)  
    //Filter::rho(gLeft)=vNormalEn-vNormalNb;
    Filter::rho(gLeft)=vNormalEn*Filter::rho(valEn)-vNormalNb*Filter::rho(valNb);
    
    Filter::rho(gLeft)*=-0.5;
     
    double viscrho=Filter::rho(jump);
    double viscmu=Filter::mu( jump );
    
    Filter::rho( gLeft )+=numVisc_*area*viscrho;
    Filter::rho( gRight )=Filter::rho( gLeft );
    
    Filter::rho( gLeft )+=numViscMu_*area*viscmu;
    Filter::rho( gRight )=Filter::mu( gLeft );
    //----------------------------------------------------------------
    
    //v---------------------------------------------------------------
    
    for(int i = 0; i<dimDomain;++i)
      { 
        //F_2=F_{2.1}+F_{2.2}
        //F_{2.1}=-(\mu^+-\mu^-)*n[i]*\rho^+*0.5;
        Filter::velocity(gLeft,i)-=Filter::mu(jump)*normal[i]*Filter::rho(valEn)*0.5;
        Filter::velocity(gRight,i)-=Filter::mu(jump)*normal[i]*Filter::rho(valNb)*0.5;

        //F_{2.2}=+(\phi^+-\phi^-)*n[i]*\tau
        Filter::velocity(gLeft,i)+=AddFilter::phi(jumpAdd)*normal[i]*AddFilter::tau(addValEn)*0.5;
        Filter::velocity(gRight,i)+=AddFilter::phi(jumpAdd)*normal[i]*AddFilter::tau(addValNb)*0.5;
      } 
    
    //----------------------------------------------------------------
}

template< class Model>
void NvStMixedFlux<Model>
::diffusionOutFlowFlux(RangeType& value,JacobianRangeType& dvalue) const
{
 value=0;
 dvalue=0;
}

template< class Model >
void NvStMixedFlux<Model>
:: diffusionFlux( const DomainType& normal,
                  const double penaltyFactor,
                  const RangeType& uEn,
                  const RangeType& uNb,
                  const JacobianRangeType& duEn,
                  const JacobianRangeType& duNb,
                  RangeType& value,
                  JacobianRangeType& dvalue) const
{
  
  RangeType jump;
  jump=uEn;
  jump-=uNb;  
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
  mean+=duNb;
  mean*=-0.5;
  model_.diffusion(mean,Amean);

  Amean.umv(normal,value);
}
template< class Model >
void NvStMixedFlux<Model>
:: diffusionBoundaryFlux( const DomainType& normal,
                  const double penaltyFactor,
                  const RangeType& uEn,
                  const JacobianRangeType& duEn,
                  RangeType& value,
                  JacobianRangeType& dvalue) const
{
  
  RangeType jump;
  jump=uEn;
  double integrationElement=normal.two_norm();
  
  for( int i=0; i<dimDomain;++i)
    {
      Filter::velocity(value,i)=penalty_*penaltyFactor*Filter::velocity(jump,i)*integrationElement;
    }
  JacobianRangeType jumpNormal(0.);
 
  // [u]\otimes n
  for(int i=0; i<dimDomain; ++i)
    for(int j=0; j<dimDomain; ++j)
      jumpNormal[ i+1 ][ j ]=-1*jump[i+1]*normal[j];
 

  jumpNormal*=switchIP_;
  model_.diffusion(jumpNormal,dvalue);
   
  JacobianRangeType mean(0.), Amean(0.);
  mean=duEn;
  mean*=-1;
  model_.diffusion(mean,Amean);

  Amean.umv(normal,value);
}

#endif


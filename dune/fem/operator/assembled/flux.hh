#ifndef PHASEFIELD_ASSEMBLED_FLUX_HH
#define PHASEFIELD_ASSEMBLED_FLUX_HH

#include <dune/common/fvector.hh>

#include "phasefieldfilter.hh"

template<class Model>
class MixedFlux
{
  //typedefs

 typedef Model ModelType;
 
  enum{ dimDomain=Model::dimDomain};
  enum{ dimRange=Model::dimRange};
 
  typedef typename Dune::FieldVector<double,dimRange> RangeType;
  typedef typename Dune::FieldVector<double,dimDomain> DomainType;
  typedef typename Dune::FieldMatrix<double,dimRange,dimDomain> JacobianRangeType;
  typedef  PhasefieldFilter<RangeType> Filter; 

  
public:
  MixedFlux(const ModelType& model,double penalty):
    model_(model),
    beta_(penalty),
    switchIP_(Dune::Fem::Parameter::getValue<int>("phasefield.ipswitch",1))
    {
    }


  double numericalFlux( const DomainType& normal, 
                        const double penaltyFactor,
                        const RangeType& vuEn,
                        const RangeType& vuN,
                        const RangeType& vuEnOld,
                        const RangeType& vuNbOld,
                        RangeType& gLeft,
                        RangeType& gRight) const;


  double diffusionFlux( const DomainType& normal,
                        const double penaltyFactor,
                        const RangeType& uEn,
                        const RangeType& uNb,
                        const JacobianRangeType& duEn,
                        const JacobianRangeType& duNb,
                        RangeType& value,
                        JacobianRangeType& dvalue) const;

  double boundaryFlux( const DomainType& normal, 
                       //const RangeType& vuEn,
                       const RangeType& vuMidEn,
                       RangeType& gLeft) const;


private:
  const ModelType& model_;
  double beta_;
  const int switchIP_; 
};


template< class Model >
double MixedFlux<Model>
::boundaryFlux(const DomainType& normal,
               const RangeType& vuMidEn,
             //  const RangeType& vuOld,
               RangeType& gLeft) const
  {
    RangeType midEn ;
//    valEn=vuEn;

 //   midEn = vuEn;
 //   midEn+= vuOld ;
  //  midEn*= 0.5;
    midEn=vuMidEn;

    double vNormalEn(0.);

    //rho-------------------------------------------------------------
  
    for(int i = 0; i<dimDomain;++i)
      {
        vNormalEn+=Filter::velocity(midEn,i)*normal[i];
      }
   
    Filter::rho(gLeft)=vNormalEn*Filter::rho(midEn);
    Filter::rho(gLeft)*=-0.5;
  
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
        laplaceFlux-=Filter::sigma(midEn,i)*normal[i]*0.5;
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
double MixedFlux<Model>
::numericalFlux( const DomainType& normal,
                const double penaltyFactor,              
                const RangeType& vuEn, // needed for calculation of sigma which is fully implicit
                const RangeType& vuNb, // needed for calculation of sigma which is fully implicit
                const RangeType& vuEnMid,
                const RangeType& vuNbMid,
                RangeType& gLeft,
                RangeType& gRight) const
  {
      RangeType valEn,valNb,midEn{0.}, midNb{0.},jump,mean;
      valEn=vuEn;
      valNb=vuNb;

      midEn=vuEnMid;
      midNb=vuNbMid;
      jump = midEn;
      jump-= midNb;
      mean = midEn ;
      mean+= midNb;
      mean*=0.5;

      //rho-------------------------------------------------------------
 
      double vNormalEn(0.),vNormalNb(0.);
    
      for(int i = 0; i<dimDomain;++i)
        {
          //v^+\cdot n^+
          vNormalEn+=Filter::velocity(midEn,i)*normal[i];
          //v^-\cdot n^+
          vNormalNb+=Filter::velocity(midNb,i)*normal[i];
        }
    
      //F_1=( \rho^+*v^+\cdot n^+ -\rho^-*v-\cdot n^+)*-0.5  
      Filter::rho(gLeft)=vNormalEn*Filter::rho(midEn)-vNormalNb*Filter::rho(midNb);
      Filter::rho(gLeft)*=-0.5;
      Filter::rho(gLeft)=0; 
      //----------------------------------------------------------------
    
      //v---------------------------------------------------------------
    
      for(int i = 0; i<dimDomain;++i)
        { 
          //F_2=F_{2.1}+F_{2.2}
          //F_{2.1}=-(\mu^+-\mu^-)*n[i]*\rho^+*0.5;
//          Filter::velocity(gLeft,i)=-1*Filter::mu(jump)*normal[i]*Filter::rho(midEn)*0.5;
          //F_{2.2}=+(\phi^+-\phi^-)*n[i]*\tau
          Filter::velocity(gLeft,i)+= Filter::phi(jump)*normal[i]*Filter::tau(midEn)*0.5;
          Filter::velocity(gLeft,i)=0.;
       } 
    
      //----------------------------------------------------------------
      double laplaceFlux(0.);
      //phi-------------------------------------------------------------
      for(int i = 0; i<dimDomain;++i)
        {
          //F_{3.1}
       
          //-(\phi^+-\phi^-)*n[i]*v[i]*0.5 
          Filter::phi(gLeft)-=Filter::phi(jump)*normal[i]*Filter::velocity(midEn,i)*0.5;
          //tau
          //F_{3.2}
          //(\sigma^+-\sigma^-)\cdot n * 0.5
          laplaceFlux+=Filter::sigma(jump,i)*normal[i]*0.5;
        
        } 
     // Filter::phi( gLeft )+=laplaceFlux;
    
      //----------------------------------------------------------------
      //tau-------------------------------------------------------------
      Filter::tau(gLeft)-=model_.delta()*laplaceFlux;
      
      //----------------------------------------------------------------

      //sigma-----------------------------------------------------------
      for(int i = 0; i<dimDomain;++i)
        {  
          //F_4
          //(\phi^+-\phi^-)\normal*0.5
          Filter::sigma(gLeft,i)=(Filter::phi(valEn)-Filter::phi(valNb))*normal[i]*0.5;
        } 
      //----------------------------------------------------------------
      return 0.;
  }

template< class Model >
double MixedFlux<Model>
:: diffusionFlux( const DomainType& normal,
                  const double penaltyFactor,
                  const RangeType& uEn,
                  const RangeType& uNb,
                  const JacobianRangeType& duEn,
                  const JacobianRangeType& duNb,
                  RangeType& value,
                  JacobianRangeType& dvalue) const
{
  
  RangeType jump{0};
  jump=uEn;
  jump-=uNb;  
  JacobianRangeType aduEn{0.}, aduNb{0.}; 
  double integrationElement=normal.two_norm();
  
  
  for( int i=0; i<dimDomain;++i)
    {
      Filter::velocity(value,i)=beta_*penaltyFactor*Filter::velocity(jump,i)*integrationElement;
    }
  JacobianRangeType jumpNormal{0.};
 
  // [u]\otimes n
  for(int i=0; i<dimDomain; ++i)
    for(int j=0; j<dimDomain; ++j)
      jumpNormal[i+1][j]=-0.5*jump[i+1]*normal[j];
 

  jumpNormal*=switchIP_;
  model_.diffusion(jumpNormal,dvalue);
   
  JacobianRangeType mean{0.}, Amean{0.};
  mean=duEn;
  mean+=duNb;
  mean*=-0.5;
  model_.diffusion(mean,Amean);

  Amean.umv(normal,value);
  return beta_;
}

#endif


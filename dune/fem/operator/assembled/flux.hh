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
  typedef  PhasefieldFilter<RangeType> Filter; 

  
public:
  MixedFlux(const ModelType& model,double penalty):
    model_(model),
    penalty_(penalty)
    {}


  double numericalFlux( const DomainType& normal, 
                        const RangeType& vuEn,
                        const RangeType& vuNb,
                        const RangeType& vuEnOld,
                        const RangeType& vuNbOld,
                        RangeType gLeft,
                        RangeType gRight) const;



  double boundaryFlux( const DomainType& normal, 
                        const RangeType& vuEn,
                        const RangeType& vuEnOld,
                        RangeType gLeft) const;


private:
  const ModelType& model_;
  double penalty_;

};


template< class Model >
double MixedFlux<Model>
::boundaryFlux(const DomainType& normal,
               const RangeType& vuEn,
               const RangeType& vuOld,
               RangeType gLeft) const
  {
    RangeType valEn,midEn ;
    valEn=vuEn;

    midEn = vuEn;
    midEn+= vuOld ;
    midEn*= 0.5;
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
        laplaceFlux+=Filter::sigma(midEn,i)*normal[i]*0.5;
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
                      const RangeType& vuEn,
                      const RangeType& vuNb,
                      const RangeType& vuEnOld,
                      const RangeType& vuNbOld,
                      RangeType gLeft,
                      RangeType gRight) const
  {
      RangeType valEn,valNb,midEn, midNb,jump,mean;
      valEn=vuEn;
      valNb=vuNb;

      midEn =  vuEn;
      midEn+= vuEnOld ;
      midEn*=0.5;

      midNb = vuNb ;
      midNb+= vuNbOld ;
      midNb*=0.5;
    
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
    
      //----------------------------------------------------------------
    
      //v---------------------------------------------------------------
    
      for(int i = 0; i<dimDomain;++i)
        { 
          //F_2=F_{2.1}+F_{2.2}
          //F_{2.1}=-(\mu^+-\mu^-)*n[i]*\rho^+*0.5;
          Filter::velocity(gLeft,i)=-1*Filter::mu(jump)*normal[i]*Filter::rho(midEn)*0.5;
          //F_{2.2}=+(\phi^+-\phi^-)*n[i]*\tau
          Filter::velocity(gLeft,i)+= Filter::phi(jump)*normal[i]*Filter::tau(midEn)*0.5;
        } 
    
      //----------------------------------------------------------------
      double laplaceFlux(0.);
      //phi-------------------------------------------------------------
      for(int i = 0; i<dimDomain;++i)
        {
          //F_{3.1}
          //-(\phi^+-\phi^-)*n[i]*v[i]*0.5 
          Filter::phi(gLeft)+=Filter::phi(jump)*normal[i]*Filter::velocity(midEn,i)*0.5;
          //tau
          //F_{3.2}
          //(\sigma^+-\sigma^-)\cdot n * 0.5
          laplaceFlux+=Filter::sigma(jump,i)*normal[i]*0.5;
        }  
      //----------------------------------------------------------------

      //tau-------------------------------------------------------------
      Filter::tau(gLeft)+=model_.delta()*laplaceFlux;
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

#endif


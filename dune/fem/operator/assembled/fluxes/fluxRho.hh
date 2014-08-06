#ifndef PHASEFIELD_ASSEMBLED_FLUX_HH
#define PHASEFIELD_ASSEMBLED_FLUX_HH

#include <dune/common/fvector.hh>

#include "../phasefieldfilter.hh"
#warning "OLDFLUX"
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
  MixedFlux(const ModelType& model):
    model_(model),
    penalty_( Dune::Fem::Parameter::getValue< double >( "phasefield.penalty" )),
    acpenalty_( Dune::Fem::Parameter::getValue< double >( "phasefield.acpenalty" )),
    switchIP_(Dune::Fem::Parameter::getValue<int>("phasefield.ipswitch",1)),
    numVisc_(Dune::Fem::Parameter::getValue<double>("phasefield.addvisc",0)), 
    numViscOld_( Dune::Fem::Parameter::getValue< double >( "phasefield.oldvisc" ,0))
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
                       const double penaltyFactor,                
                       const RangeType& vuEn,  
                       const RangeType& vuMidEn,
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
  const double numViscOld_;
};


template< class Model >
double MixedFlux<Model>
::boundaryFlux(const DomainType& normal,
               const double penaltyFactor,
               const RangeType& vuEn,
               const RangeType& vuMidEn,
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
double MixedFlux<Model>
::numericalFlux( const DomainType& normal,
                const double area,              
                const RangeType& vuEn, // needed for calculation of sigma which is fully implicit
                const RangeType& vuNb, // needed for calculation of sigma which is fully implicit
                const RangeType& vuEnMid,
                const RangeType& vuNbMid,
                RangeType& gLeft,
                RangeType& gRight) const
  {
      RangeType valEn,valNb,midEn(0.), midNb(0.),jump,mean,jumpOld;
      valEn=vuEn;
      valNb=vuNb;

      midEn=vuEnMid;
      midNb=vuNbMid;
      jump = midEn;
      jump-= midNb;
      mean = midEn ;
      mean+= midNb;
      mean*=0.5;
      jumpOld=valEn;
      jumpOld-=valNb;
      //rho-------------------------------------------------------------
 
      double vNormalEn(0.),vNormalNb(0.);
    
      for(int i = 0; i<dimDomain;++i)
        {
          //v^+\cdot n^+
          vNormalEn+=Filter::velocity(midEn,i)*normal[i];
          //v^-\cdot n^+
          vNormalNb+=Filter::velocity(midNb,i)*normal[i];
        }
    
      //F_1=-0.5*( \rho^+*v^+\cdot n^+ -\rho^-*v-\cdot n^+)  
      //Filter::rho(gLeft)=vNormalEn-vNormalNb;
      Filter::rho(gLeft)=vNormalEn*Filter::rho(midEn)-vNormalNb*Filter::rho(midNb);
      Filter::rho(gLeft)*=-0.5;
     
      //double visc=Filter::rho(jump);
      //Filter::rho(gLeft)+=numVisc_*area*visc;
      Filter::rho( gRight )=Filter::rho( gLeft );
      //----------------------------------------------------------------
    
      //v---------------------------------------------------------------
    
      for(int i = 0; i<dimDomain;++i)
        { 
          //F_2=F_{2.1}+F_{2.2}
          //F_{2.1}=-(\mu^+-\mu^-)*n[i]*\rho^+*0.5;
          Filter::velocity(gLeft,i)-=Filter::mu(jump)*normal[i]*Filter::rho(midEn)*0.5;
          Filter::velocity(gRight,i)-=Filter::mu(jump)*normal[i]*Filter::rho(midNb)*0.5;
          
           //F_{2.2}=+(\phi^+-\phi^-)*n[i]*\tau
          Filter::velocity(gLeft,i)+= Filter::phi(jump)*normal[i]*Filter::tau(midEn)*0.5;
          Filter::velocity(gRight,i)+= Filter::phi(jump)*normal[i]*Filter::tau(midNb)*0.5;

        } 
    
      //----------------------------------------------------------------
      double laplaceFlux(0.);
      //phi-------------------------------------------------------------
      for(int i = 0; i<dimDomain;++i)
        {
          //F_{3.1}
       
          //-(\phi^+-\phi^-)*n[i]*v[i]*0.5 
          Filter::phi(gLeft)-=Filter::phi(jump)*normal[i]*Filter::velocity(midEn,i)*0.5;
          Filter::phi(gRight)-=Filter::phi(jump)*normal[i]*Filter::velocity(midNb,i)*0.5;
       
          //tau
          //F_{3.2}
          //(h2(rho^+)\sigma^+-h2(rho^-)\sigma^-)\cdot n * 0.5
          laplaceFlux+=(Filter::sigma(midEn,i)*model_.h2(Filter::rho(midEn))-Filter::sigma(midNb,i)*model_.h2(Filter::rho(midNb)))*normal[i]*0.5;
       //   laplaceFlux+=(Filter::alpha(midEn,i)-Filter::alpha(midNb,i))*normal[i]*0.5;
        

        } 
    
      //----------------------------------------------------------------
      //tau-------------------------------------------------------------
      Filter::tau(gLeft)-=model_.delta()*laplaceFlux;
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
  mean+=duNb;
  mean*=-0.5;
  model_.diffusion(mean,Amean);

  Amean.umv(normal,value);
  return penalty_;
}
template< class Model >
double MixedFlux<Model>
:: diffusionBoundaryFlux( const DomainType& normal,
                  const double penaltyFactor,
                  const RangeType& uEn,
                  const JacobianRangeType& duEn,
                  RangeType& value,
                  JacobianRangeType& dvalue) const
{
  
  RangeType jump{0};
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


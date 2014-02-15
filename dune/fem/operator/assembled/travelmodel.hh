#ifndef PHASEFIELD_MIXED_MODEL_HH
#define PHASEFIELD_MIXED_MODEL_HH

#include<dune/common/fvector.hh>
#include<dune/common/fmatrix.hh>

#include<dune/fem/io/parameter.hh>

#include "phasefieldfilter.hh"

template<class Grid, class Problem>
class HeatModel
{
  public:
  typedef Problem ProblemType;
  typedef Grid GridType;
  enum{ dimDomain = GridType::dimensionworld };
 
  enum{ dimRange = ProblemType::dimRange };
  
  typedef typename ProblemType::ThermodynamicsType ThermodynamicsType;
  typedef double RangeFieldType; 
  typedef typename Dune::FieldVector<RangeFieldType,dimRange> RangeType;
  typedef typename Dune::FieldVector<RangeFieldType,dimDomain> DomainType; 
  typedef typename Dune::FieldMatrix<RangeFieldType,dimRange,dimDomain> JacobianRangeType;

  typedef PhasefieldFilter<RangeType> Filter;


//contructor
 public:
 HeatModel( const ProblemType& problem):
   problem_(problem),
   penalty_(Dune::Fem::Parameter::getValue<double>("phasefield.penalty"))
  {}


  inline void totalEnergy( const DomainType& xgl,
                           RangeType& vu,
                           double& kin,
                           double& therm,
                           double& total) const;
                       
  // additional Source for the whole syten eg. for 
  // generatring exact solutions
  inline void systemSource( const double time,
                            const DomainType& xgl,
                            RangeType& s) const;



  inline void  muSource( const RangeFieldType rho1,
                         const RangeFieldType rho2,
                         const RangeFieldType phi,
                         RangeFieldType& mu) const;
  
  inline void tauSource( const RangeFieldType phi1,
                         const RangeFieldType phi2,
                         const RangeFieldType rho,
                         RangeFieldType& tau) const;
 
  inline void dirichletValue( const double time,
                              const DomainType& xglobal,
                              RangeType& g) const;
  
  inline  double penalty() const
    {
      return penalty_;
    }
  
  inline void diffusion(JacobianRangeType& vu,
                        JacobianRangeType& diffusion) const;


  inline double delta() const
    {
      return problem_.thermodynamics().delta();
    }

  private:
    const ProblemType& problem_;
    double penalty_; 

};

template< class Grid, class Problem>
inline void HeatModel< Grid,Problem>
::totalEnergy(const DomainType& xgl,
              RangeType& vu,
              double& kin,
              double& therm,
              double& total) const
  {
    double rho=Filter::rho(vu);
    double phi=Filter::phi(vu);
    double kineticEnergy{0.},surfaceEnergy{0.};
    for(int i=0; i<dimDomain; i++)
    {
      kineticEnergy+=Filter::velocity(vu,i)*Filter::velocity(vu,i);
#if DGSCHEME
      surfaceEnergy+=Filter::sigma(vu,i)*Filter::sigma(vu,i);
#elif FEMSCHEME
#warning "TOTAL ENERGY NEEDS JACOBIANRANGETYPE - NOT IMPLEMENTED"
#endif
    }
  
    kin=rho*0.5*kineticEnergy;
    surfaceEnergy*=0.5;
    surfaceEnergy*=problem_.thermodynamics().delta();
  
    therm=problem_.thermodynamics().helmholtz(rho,phi);
    therm+=surfaceEnergy;
   
    total=therm+kin;
  }
  
template< class Grid, class Problem>
inline void HeatModel< Grid,Problem>
::systemSource( const double time,
                const DomainType& xgl,
                RangeType& s) const
  {
    s=0.;
#if 1 
    double cosx=std::cos(2*M_PI*xgl[0]);
    double cost=std::cos(M_PI*time);
    double sinx=std::sin(2*M_PI*xgl[0]);
    double sint=std::sin(M_PI*time);
    double phi=0.5*cosx*cost+0.5;
    double delta=problem_.thermodynamics().delta();
    double deltaInv=problem_.thermodynamics().deltaInv() ;
    double v=problem_.thermodynamics().velo();

    //f = 2*(d_t phi-\Delta phi)
    double f=M_PI*( 4*M_PI*std::cos( M_PI*time ) - std::sin( M_PI*time ) );
    f*=cosx;
  

    Filter::rho( s )=0;
   
    //tension= -\nabla\phi tau
    double tension=-sinx*M_PI*cost*(-1.*cosx*cost+cosx*cosx*cosx*cost*cost*cost+2*cosx*M_PI*M_PI*cost);   


    
     Filter::velocity( s , 0 )=-tension;

     if(dimDomain==2)
       Filter::velocity(s, 1)=0.;

    //f=2(d_t\phi-\Delta\phi)
    Filter::phi( s )=f*0.5;
  
    // transportphi=v\cdot\nabla\phi
    double transportphi=-M_PI*sinx*cost*v;

    // dFdphi=dF(rho,phi)/dphi
    
    double dFdphi{0.};
    tauSource(1,phi,phi,dFdphi);
    Filter::phi(s)+=dFdphi+transportphi;
      
#endif


     
  }
  
template<class Grid, class Problem > 
inline void HeatModel< Grid, Problem >
::muSource( RangeFieldType rho1,
            RangeFieldType rho2,
            RangeFieldType phi,
            RangeFieldType& mu) const
  {
    mu=problem_.thermodynamics().chemicalPotential(rho1,phi);
    mu=0;
  }

template< class Grid, class Problem > 
inline void HeatModel< Grid, Problem>
::tauSource(RangeFieldType phi1,
            RangeFieldType phi2,
            RangeFieldType rho,
            RangeFieldType& tau) const
  {
    tau=problem_.thermodynamics().reactionSource(rho,phi2);
   // tau=0;
  }

template< class Grid, class Problem>
inline void HeatModel< Grid,Problem>
::dirichletValue(const double time, const DomainType& xglobal, RangeType& g) const
{
  problem_.evaluate(time,xglobal,g);
}

template< class Grid, class Problem > 
inline void HeatModel< Grid, Problem>
::diffusion( JacobianRangeType& dvu,
            JacobianRangeType& diffusion) const
  {
    diffusion=0;
    double mu1=problem_.thermodynamics().mu1();
 //   double mu2=problem_.thermodynamics().mu2();
  
  //  diffusion*=mu1;
    

    for(int ii = 0 ; ii < dimDomain ; ++ii )
     {
       for(int jj = 0; jj < dimDomain ; ++ jj)
        Filter::dvelocity(diffusion,ii,jj )=mu1*Filter::dvelocity(dvu,ii,jj);
   
  //     for(int j=0; j<dimDomain ; ++j )
    //      {
      //      Filter::dvelocity(diffusion,i,j)+=mu2*0.5*(Filter::dvelocity(dvu,i,j)+Filter::dvelocity(dvu,j,i));
        // }
      }

}

#endif


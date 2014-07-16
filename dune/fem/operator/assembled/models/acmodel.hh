#ifndef PHASEFIELD_MIXED_MODEL_HH
#define PHASEFIELD_MIXED_MODEL_HH

#include<dune/common/fvector.hh>
#include<dune/common/fmatrix.hh>

#include<dune/fem/io/parameter.hh>

#include "../phasefieldfilter.hh"

template<class Grid, class Problem>
class PhasefieldModel
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
    PhasefieldModel( const ProblemType& problem):
      problem_(problem)
  {}


    inline void totalEnergy ( const DomainType& xgl,
                              RangeType& vu,
                              double& kin,
                              double& therm,
                              double& total) const;

    // additional Source for the whole syten eg. for 
    // generatring exact solutions
    inline void systemSource ( const double time,
                               const DomainType& xgl,
                               RangeType& s) const;


    inline void  muSource ( const RangeFieldType rho,
                            const RangeFieldType rhoOld,
                            const RangeFieldType phi,
                            RangeFieldType& mu) const;

    inline void  drhomuSource ( const RangeFieldType rho,
                                const RangeFieldType rhoOld,
                                const RangeFieldType phi,
                                RangeFieldType& mu) const;

    inline void  dphimuSource ( const RangeFieldType rho,
                                const RangeFieldType rhoOld,
                                const RangeFieldType phi,
                                RangeFieldType& mu) const;

    inline void tauSource ( const RangeFieldType phi,
                            const RangeFieldType phiOld,
                            const RangeFieldType rho,
                            RangeFieldType& tau) const;

    inline void dphitauSource ( const RangeFieldType phi,
                                const RangeFieldType phiOld,
                                const RangeFieldType rho,
                                RangeFieldType& tau) const;

    inline void dirichletValue ( const double time,
                                 const DomainType& xglobal,
                                 RangeType& g) const;


    inline void diffusion ( JacobianRangeType& vu,
                            JacobianRangeType& diffusion) const;
    
    inline RangeFieldType h2 ( double rho ) const
    {
      return 1./rho;
    }

    inline RangeFieldType h2prime( double rho ) const
    {
      return -1./(rho*rho);
    }

    inline double reactionFactor () const
    { 
      return problem_.thermodynamics().reactionFactor();
    }

    inline double delta () const
    {
      return problem_.thermodynamics().delta();
    }
    inline double deltaInv () const
    {
      return problem_.thermodynamics().deltaInv();
    }


  private:
    const ProblemType& problem_;

};

template< class Grid, class Problem>
inline void PhasefieldModel< Grid,Problem>
::totalEnergy ( const DomainType& xgl,
                RangeType& vu,
                double& kin,
                double& therm,
                double& total ) const
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
#if RHOMODEL
  surfaceEnergy*=h2(rho);
#else
#endif
  kin=rho*0.5*kineticEnergy;
  surfaceEnergy*=0.5;
  surfaceEnergy*=problem_.thermodynamics().delta();

  therm=problem_.thermodynamics().helmholtz(rho,phi);
  therm+=surfaceEnergy;

  total=therm+kin;
}

template< class Grid, class Problem>
inline void PhasefieldModel< Grid,Problem>
::systemSource ( const double time,
                 const DomainType& xgl,
                 RangeType& s ) const
{
  s=0.;
}

template<class Grid, class Problem > 
inline void PhasefieldModel< Grid, Problem >
::muSource ( RangeFieldType rho,
             RangeFieldType rhoOld,
             RangeFieldType phi,
             RangeFieldType& mu) const
{

  mu=problem_.thermodynamics().chemicalPotential(rhoOld,phi);
}
template<class Grid, class Problem > 
inline void PhasefieldModel< Grid, Problem >
::drhomuSource ( RangeFieldType rho,
                 RangeFieldType rhoOld,
                 RangeFieldType phi,
                 RangeFieldType& mu ) const

{
  mu=problem_.thermodynamics().drhochemicalPotential(rhoOld,phi);
}
template<class Grid, class Problem > 
inline void PhasefieldModel< Grid, Problem >
::dphimuSource ( RangeFieldType rho,
                 RangeFieldType rhoOld,
                 RangeFieldType phi,
                 RangeFieldType& mu ) const
{
  mu=problem_.thermodynamics().dphichemicalPotential(rhoOld,phi);
}


template< class Grid, class Problem > 
inline void PhasefieldModel< Grid, Problem>
::tauSource ( RangeFieldType phi,
              RangeFieldType phiOld,
              RangeFieldType rho,
              RangeFieldType& tau) const
{
  tau=problem_.thermodynamics().reactionSource(rho,phiOld);
}
template< class Grid, class Problem > 
inline void PhasefieldModel< Grid, Problem>
::dphitauSource ( RangeFieldType phi,
                  RangeFieldType phiOld,
                  RangeFieldType rho,
                  RangeFieldType& tau) const
{
  tau=problem_.thermodynamics().dphireactionSource(rho,phiOld);
}


template< class Grid, class Problem>
inline void PhasefieldModel< Grid,Problem>
::dirichletValue(const double time, const DomainType& xglobal, RangeType& g) const
{
  problem_.evaluate(time,xglobal,g);
}

template< class Grid, class Problem > 
inline void PhasefieldModel< Grid, Problem>
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


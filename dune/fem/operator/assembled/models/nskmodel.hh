#ifndef NSK_MIXED_MODEL_HH
#define NSK_MIXED_MODEL_HH

#include<dune/common/fvector.hh>
#include<dune/common/fmatrix.hh>

#include<dune/fem/io/parameter.hh>

#include "../nskfilter.hh"
#define LAPLACE  1
template<class Grid, class Problem>
class NSKModel
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
    typedef typename Dune::FieldMatrix<RangeFieldType, dimDomain , dimDomain> ComponentDiffusionType;
    typedef typename std::array< ComponentDiffusionType,dimDomain> DiffusionTensorType;
    typedef PhasefieldFilter<RangeType> Filter;


    //contructor
  public:
    NSKModel( const ProblemType& problem):
      problem_(problem)
  {}


    inline void totalEnergy ( const DomainType& xgl,
                              RangeType& vu,
                              double& kin,
                              double& therm,
                              double& total,
                              double& surf ) const;

    // additional Source for the whole syten eg. for 
    // generatring exact solutions
    inline void systemSource ( const double time,
                               const DomainType& xgl,
                               RangeType& s) const;


    inline void  muSource ( const RangeFieldType rho,
                            const RangeFieldType rhoOld,
                            RangeFieldType& mu) const;

    inline void  drhomuSource ( const RangeFieldType rho,
                                const RangeFieldType rhoOld,
                                RangeFieldType& mu) const;


    inline void dirichletValue ( const double time,
                                 const DomainType& xglobal,
                                 RangeType& g) const;
    inline double diffusion () const
    { 
      return problem_.thermodynamics().mu1();
    }

    template< class JacobianVector>
    inline void scalar2vectorialDiffusion ( const JacobianVector& dphi , DiffusionTensorType& diffusion ) const; 
    inline void diffusion ( JacobianRangeType& vu,
                            JacobianRangeType& diffusion) const;
    

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
inline void NSKModel< Grid,Problem>
::totalEnergy ( const DomainType& xgl,
                RangeType& vu,
                double& kin,
                double& therm,
                double& total,
                double& surf ) const
{
  double rho=Filter::rho(vu);
  double kineticEnergy{0.},surfaceEnergy{0.};
  for(int i=0; i<dimDomain; i++)
  {
    kineticEnergy+=Filter::velocity(vu,i)*Filter::velocity(vu,i);
    surfaceEnergy+=Filter::sigma(vu,i)*Filter::sigma(vu,i);
  }
  
  
  kin=rho*0.5*kineticEnergy;
  surfaceEnergy*=0.5;
  surfaceEnergy*=problem_.thermodynamics().delta();

  therm=problem_.thermodynamics().helmholtz(rho);
  therm+=surfaceEnergy;

  total=therm+kin;
  surf=surfaceEnergy;
}

template< class Grid, class Problem>
inline void NSKModel< Grid,Problem>
::systemSource ( const double time,
                 const DomainType& xgl,
                 RangeType& s ) const
{
  s=0.;
}

template<class Grid, class Problem > 
inline void NSKModel< Grid, Problem >
::muSource ( RangeFieldType rho,
             RangeFieldType rhoOld,
             RangeFieldType& mu) const
{

  mu=problem_.thermodynamics().chemicalPotential(rhoOld);
}
template<class Grid, class Problem > 
inline void NSKModel< Grid, Problem >
::drhomuSource ( RangeFieldType rho,
                 RangeFieldType rhoOld,
                 RangeFieldType& mu ) const

{
  mu=problem_.thermodynamics().drhochemicalPotential(rhoOld);
}


template< class Grid, class Problem>
inline void NSKModel< Grid,Problem>
::dirichletValue(const double time, const DomainType& xglobal, RangeType& g) const
{
  problem_.evaluate(time,xglobal,g);
}


template< class Grid, class Problem > 
template< class JacobianVector>
inline void NSKModel< Grid, Problem>
::scalar2vectorialDiffusion( const JacobianVector& dphi,DiffusionTensorType& du) const
{

  double mu1=0.5*problem_.thermodynamics().mu1();
  double mu2=problem_.thermodynamics().mu2();
  for( int ii=0 ; ii < dimDomain  ; ++ii)
    {
#if  LAPLACE
      for(int jj=0 ; jj < dimDomain ; ++jj )
        du[ ii ][ ii ][ jj ] = mu2*dphi[ 0 ][ jj ];
#else
// full diffusion
      for( int jj = 0; jj < dimDomain; ++ jj )
        {
          du[ ii ][ jj ][ ii ]=mu1*dphi[ 0 ][ jj ];
          du[ ii ][ ii ][ jj ]=mu1*dphi[ 0 ][ jj ];
        }
      for( int jj=0 ; jj < dimDomain ; ++jj ) 
        {
          du[ ii ][ jj ][ jj ]+=mu2*dphi[ 0 ][ ii ];
        }
     du[ ii ][ ii ][ ii ]+=mu2*dphi[ 0 ] [ ii ];
#endif
}
} 

template< class Grid, class Problem > 
inline void NSKModel< Grid, Problem>
::diffusion( JacobianRangeType& dvu,
    JacobianRangeType& diffusion) const
{
  diffusion=0;
 double mu2=problem_.thermodynamics().mu2();

#if LAPLACE 
   for(int ii = 0 ; ii < dimDomain ; ++ii )
    for(int jj = 0; jj < dimDomain ; ++ jj)
     Filter::dvelocity(diffusion,ii,jj )=mu2*Filter::dvelocity(dvu,ii,jj);
#else
  double mu1=problem_.thermodynamics().mu1();
 
  for(int ii = 0 ; ii < dimDomain ; ++ii )
    {
      for(int jj = 0; jj < dimDomain ; ++ jj)
        {
          Filter::dvelocity(diffusion,ii,ii)+=mu2*Filter::dvelocity(dvu,jj,jj);
        }

      for(int jj=0; jj<dimDomain ; ++jj )
        {
          Filter::dvelocity(diffusion,ii,jj)+=mu1*0.5*(Filter::dvelocity(dvu,ii,jj)+Filter::dvelocity(dvu,jj,ii));
        } 
      }
#endif
}

#endif


#ifndef PHASEFIELD_MIXED_MODEL_HH
#define PHASEFIELD_MIXED_MODEL_HH

#include<dune/common/fvector.hh>
#include<dune/common/fmatrix.hh>

#include<dune/fem/io/parameter.hh>

#include "../phasefieldfilter.hh"
#define LAPLACE  1
#define DIFFQUOT 0 
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
    typedef typename Dune::FieldMatrix<RangeFieldType, dimDomain , dimDomain> ComponentDiffusionType;
    typedef typename std::array< ComponentDiffusionType,dimDomain> DiffusionTensorType;
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
                              double& total,
                              double& surf ) const;

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

    inline double maxSpeed( const DomainType& normal,
                             const RangeType& u) const ;
    
    // this gives an estimation of the lipschitz constant of the rhs of the phasefield equation
    inline double allenCahnConstant() const
    {
      return problem_.thermodynamics().lipschitzC();
    }

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

    inline RangeFieldType pressure (double rho , double phi) const
    {
      return problem_.thermodynamics().pressure( rho , phi );
    }
    inline RangeFieldType h2 ( double rho ) const
    {
      return  problem_.thermodynamics().h2( rho );
    }

    inline RangeFieldType h2prime( double rho ) const
    {
      return problem_.thermodynamics().h2prime( rho );
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
                double& total,
                double& surf ) const
{
  double rho=Filter::rho(vu);
  double phi=Filter::phi(vu);
  double kineticEnergy{0.},surfaceEnergy{0.};
  for(int i=0; i<dimDomain; i++)
  {
    kineticEnergy+=Filter::velocity(vu,i)*Filter::velocity(vu,i);
    surfaceEnergy+=Filter::sigma(vu,i)*Filter::sigma(vu,i);
  }
  
  surfaceEnergy*=h2(rho);
  
  kin=rho*0.5*kineticEnergy;
  surfaceEnergy*=0.5;
  surfaceEnergy*=problem_.thermodynamics().delta();

  therm=problem_.thermodynamics().helmholtz(rho,phi);
  //therm+=surfaceEnergy;

  total=therm+kin+surfaceEnergy;
  surf=surfaceEnergy;
}

template< class Grid, class Problem>
inline void PhasefieldModel< Grid,Problem>
::systemSource ( const double time,
                 const DomainType& xgl,
                 RangeType& s ) const
{
 s=0.;
#if PROBLEM==6 
  double x=xgl[0];
  s[1]=problem_.veloSource(x);
  s[2]=problem_.phiSource(x);
#endif
}

template<class Grid, class Problem > 
inline void PhasefieldModel< Grid, Problem >
::muSource ( RangeFieldType rho,
             RangeFieldType rhoOld,
             RangeFieldType phi,
             RangeFieldType& mu) const
{
#if DIFFQUOT
  double diffrho=rho-rhoOld;
  if( std::abs(diffrho)<1e-9)
    mu=problem_.thermodynamics().chemicalPotential(rhoOld,phi);
  else
    {
      double fnew,fold;
      fnew=problem_.thermodynamics().helmholtz(rho,phi);
      fold=problem_.thermodynamics().helmholtz(rhoOld,phi);
      mu=(fnew-fold)/(rho-rhoOld);
    }
#else
    mu=problem_.thermodynamics().chemicalPotential(rhoOld,phi);
#endif

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
              RangeFieldType rhoOld,
              RangeFieldType& tau) const
{
#if DIFFQUOT
  double diffphi=phi-phiOld;
  if( std::abs(diffphi)<1e-9)
    tau=problem_.thermodynamics().reactionSource(rhoOld,phiOld);
  else
  {
    double fnew,fold;
    fnew=problem_.thermodynamics().helmholtz(rhoOld,phi);
    fold=problem_.thermodynamics().helmholtz(rhoOld,phiOld);
    tau=(fnew-fold)/(phi-phiOld);
  }
#else
   tau=problem_.thermodynamics().reactionSource(rhoOld,phiOld);
#endif 
}
template< class Grid, class Problem > 
inline void PhasefieldModel< Grid, Problem>
::dphitauSource ( RangeFieldType phi,
                  RangeFieldType phiOld,
                  RangeFieldType rho,
                  RangeFieldType& tau) const
{
 //abort();
 tau=problem_.thermodynamics().dphireactionSource(rho,phiOld);
}

template< class Grid, class Problem >
inline double PhasefieldModel< Grid, Problem>
::maxSpeed( const DomainType& normal,
            const RangeType& u) const
{
    double unormal(0.);
    for( int ii = 0 ; ii<dimDomain ; ++ii)
      unormal+=u[1+ii]*normal[ii];
    double c=problem_.thermodynamics().a(u[0],u[dimDomain+1]);

    return std::abs(unormal)+std::sqrt(c);
}

template< class Grid, class Problem>
inline void PhasefieldModel< Grid,Problem>
::dirichletValue(const double time, const DomainType& xglobal, RangeType& g) const
{
  problem_.evaluate(time,xglobal,g);
}


template< class Grid, class Problem > 
template< class JacobianVector>
inline void PhasefieldModel< Grid, Problem>
::scalar2vectorialDiffusion( const JacobianVector& dphi,DiffusionTensorType& du) const
{
#if  LAPLACE
#else
  double mu1=0.5*problem_.thermodynamics().mu1();
#endif
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
inline void PhasefieldModel< Grid, Problem>
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


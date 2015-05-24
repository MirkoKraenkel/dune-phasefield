#ifndef PHASEFIELD_MIXED_MODEL_HH
#define PHASEFIELD_MIXED_MODEL_HH

#include<dune/common/fvector.hh>
#include<dune/common/fmatrix.hh>

#include<dune/fem/io/parameter.hh>

#include "../phasefieldfilter.hh"
#define LAPLACE  1
template<class Grid, class Problem>
class PhasefieldModel
{
  public:
    typedef Problem ProblemType;
    typedef Grid GridType;
    enum{ dimDomain = GridType::dimensionworld };
    enum{ dimAcRange = ProblemType::dimRange/2};
    enum{ dimNvStRange=ProblemType::dimRange/2};
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
      problem_(problem),
      diffquotthresh_(Dune::Fem::Parameter::getValue<double>("phasefield.diffquotthresh"))
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
                               const RangeType& vu,
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
    inline double lipschitzC() const
    {
      return problem_.thermodynamics().lipschitzC();
    }

    inline void dirichletValue ( const double time,
                                 const DomainType& xglobal,
                                 RangeType& g) const;


    template< class JacobianVector>
    inline void scalar2vectorialDiffusion ( const RangeType& vu,
                                            const JacobianVector& dphi,
                                            DiffusionTensorType& diffusion ) const;

    template< class JacobianRange >
    inline void diffusion ( const RangeType& vu,
                            const JacobianRange& dvu,
                            JacobianRange& diffusion) const;

    template< class JacobianRange >
    inline void diffusionprime ( const RangeType& vu,
                                 const JacobianRange& dvu,
                                 JacobianRange& diffusion) const;




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
    inline void  diffFactors (const RangeFieldType phi, double mu1, double mu2 ) const
    {
#warning "CORRECT THERMODYNAMICS"
      mu2= phi*problem_.thermodynamics().mu1()+(1-phi)*problem_.thermodynamics().mu2();
      mu1= 0;
    }

    inline void diffFactorsPrime (const RangeFieldType phi, double mu1 , double mu2 ) const
    {
      mu2=problem_.thermodynamics().mu1()-problem_.thermodynamics().mu2();
      mu1=0;
    }

    template< class JacobianRange >
    inline void diffusion ( const RangeFieldType& mu1,
                            const RangeFieldType& mu2,
                            const RangeType& vu,
                            const JacobianRange& dvu,
                            JacobianRange& diffusion) const;

    template< class JacobianVector>
    inline void scalar2vectorialDiffusion ( const RangeFieldType& mu1,
                                            const RangeFieldType& mu2,
                                            const JacobianVector& dphi,
                                            DiffusionTensorType& diffusion ) const;
 
    const ProblemType& problem_;
    const double diffquotthresh_;
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
                 const RangeType& vu,
                 const DomainType& xgl,
                 RangeType& s ) const
{
 s=0.;
 s[2]=vu[0]*-9.81;
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
#if TAYLOR
  mu=problem_.thermodynamics().chemicalPotential(rho,phi,rhoOld);
#elif DIFFQUOT
  double diffrho=rho-rhoOld;
  if( std::abs(diffrho)<diffquotthresh_)
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
#if TAYLOR
  mu=problem_.thermodynamics().drhochemicalPotential(rho,phi,rhoOld);
#else
  mu=problem_.thermodynamics().drhochemicalPotential(rhoOld,phi);
#endif
}
template<class Grid, class Problem > 
inline void PhasefieldModel< Grid, Problem >
::dphimuSource ( RangeFieldType rho,
                 RangeFieldType rhoOld,
                 RangeFieldType phi,
                 RangeFieldType& mu ) const
{
#if TAYLOR
  mu=problem_.thermodynamics().dphichemicalPotential(rho,phi,rhoOld);
#else
  mu=problem_.thermodynamics().dphichemicalPotential(rhoOld,phi);
#endif
} 


template< class Grid, class Problem > 
inline void PhasefieldModel< Grid, Problem>
::tauSource ( RangeFieldType phi,
              RangeFieldType phiOld,
              RangeFieldType rhoOld,
              RangeFieldType& tau) const
{
#if TAYLOR
  tau=problem_.thermodynamics().reactionSource(rhoOld,phi,phiOld);
#elif DIFFQUOT
  double diffphi=phi-phiOld;
  if( std::abs(diffphi)<diffquotthresh_)
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
#if TAYLOR
  tau=problem_.thermodynamics().dphireactionSource(rho,phi,phiOld);
#else
  tau=problem_.thermodynamics().dphireactionSource(rho,phiOld);
#endif
}

template< class Grid, class Problem >
inline double PhasefieldModel< Grid, Problem>
::maxSpeed( const DomainType& normal,
            const RangeType& u) const
{
    double unormal(0.);
    for( int ii = 0 ; ii<dimDomain ; ++ii)
      unormal+=u[1+ii]*normal[ii];
    double c=problem_.thermodynamics().a(u[0],u[dimDomain+1])*normal.two_norm2();

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
::scalar2vectorialDiffusion( const RangeType& vu ,const JacobianVector& dphi,DiffusionTensorType& du) const
{
  double mu1,mu2;
  diffFactors( vu[dimDomain+1], mu1, mu2);
  scalar2vectorialDiffusion( mu1, mu2, dphi, du);
}

template< class Grid, class Problem >
template< class JacobianVector>
inline void PhasefieldModel< Grid, Problem>
::scalar2vectorialDiffusion ( const RangeFieldType& mu1,
                              const RangeFieldType& mu2,
                              const JacobianVector& dphi,
                              DiffusionTensorType& du) const
{
  //du[i][j][k]
  
#if  LAPLACE
  for( int ii=0 ; ii < dimDomain  ; ++ii)
    {
      for(int jj=0 ; jj < dimDomain ; ++jj )
        du[ ii ][ ii ][ jj ] = mu2*dphi[ 0 ][ jj ];
    }
#else
  for( int ii=0 ; ii < dimDomain  ; ++ii)
    {
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
    }
#endif
} 

template< class Grid, class Problem > 
template< class JacobianRange >
inline void PhasefieldModel< Grid, Problem>
::diffusion ( const RangeType& vu,
              const JacobianRange& dvu,
              JacobianRange& diff) const
{
  double mu1,mu2;
  diffFactors( vu[dimDomain+1], mu1, mu2);
  diffusion( mu1 , mu2 , vu , dvu , diff );
}

template< class Grid, class Problem >
template< class JacobianRange >
inline void PhasefieldModel< Grid, Problem>
::diffusionprime ( const RangeType& vu,
                   const JacobianRange& dvu,
                   JacobianRange& diff) const
{
  double mu1,mu2;
  diffFactorsPrime( vu[dimDomain+1] , mu1 , mu2 );
  diffusion( mu1 , mu2 , vu , dvu , diff );
}


template< class Grid, class Problem >
template< class JacobianRange >
inline void PhasefieldModel< Grid, Problem>
::diffusion ( const RangeFieldType& mu1,
              const RangeFieldType& mu2,
              const RangeType& vu,
              const JacobianRange& dvu,
              JacobianRange& diffusion) const
{
  diffusion=0;
#if LAPLACE 
   for(int ii = 0 ; ii < dimDomain ; ++ii )
    for(int jj = 0; jj < dimDomain ; ++ jj)
      diffusion[1+ii][jj]=mu2*dvu[1+ii][jj];
      
#else
  abort();
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


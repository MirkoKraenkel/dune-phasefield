#ifndef PHASEFIELD_MIXED_MODEL_HH
#define PHASEFIELD_MIXED_MODEL_HH

#include<dune/common/fvector.hh>
#include<dune/common/fmatrix.hh>


#include<dune/fem/io/parameter.hh>

#include "phasefieldfilter.hh"

template<class Grid, class Problem>
class MixedModel
{
  public:
  typedef Problem ProblemType;
  typedef Grid GridType;
  enum{ dimDomain = GridType::dimensionworld };
  enum{ dimRange = 2*dimDomain+4 };
  
  typedef typename ProblemType::ThermodynamicsType ThermodynamicsType;
  typedef double RangeFieldType; 
  typedef typename Dune::FieldVector<RangeFieldType,dimRange> RangeType;
  typedef typename Dune::FieldMatrix<RangeFieldType,dimRange,dimDomain> JacobianRangeType;

  typedef PhasefieldFilter<RangeType> Filter;


//contructor
  public:
 MixedModel( const ProblemType& problem):
   problem_(problem),
   penalty_(Dune::Fem::Parameter::getValue<double>("phasefield.penalty"))
  {}

  
  inline void  muSource(const RangeFieldType rho1,
                        const RangeFieldType rho2,
                        const RangeFieldType phi,
                        RangeFieldType& mu) const;
  
  inline void tauSource(const RangeFieldType phi1,
                        const RangeFieldType phi2,
                        const RangeFieldType rho,
                        RangeFieldType& tau) const;
 
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

template<class Grid, class Problem > 
inline void MixedModel< Grid, Problem >
::muSource( RangeFieldType rho1,
            RangeFieldType rho2,
            RangeFieldType phi,
            RangeFieldType& mu) const
  {
    mu=problem_.thermodynamics().chemicalPotential(rho1,phi);
  }

template< class Grid, class Problem > 
inline void MixedModel< Grid, Problem>
::tauSource(RangeFieldType phi1,
            RangeFieldType phi2,
            RangeFieldType rho,
            RangeFieldType& tau) const
  {
    tau=problem_.thermodynamics().reactionSource(rho,phi1);
  }

template< class Grid, class Problem > 
inline void MixedModel< Grid, Problem>
::diffusion( JacobianRangeType& dvu,
            JacobianRangeType& diffusion) const
  {
    diffusion=0;
    double mu1=problem_.thermodynamics().mu1();
    double mu2=problem_.thermodynamics().mu2();
  
    diffusion=dvu;
    diffusion*=mu2;


#if 0 
    for(int i=0; i<dimDomain ; ++i )
      {
        Filter::dvelocity(diffusion,i,i)=mu2*Filter::dvelocity(dvu,i,i);
        for(int j=0; j<dimDomain ; ++j )
          {
            Filter::dvelocity(diffusion,i,j)+=mu1*0.5*(Filter::dvelocity(dvu,i,j)+Filter::dvelocity(dvu,j,i));
          }
      }
#endif

}

#endif

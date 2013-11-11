#ifndef PHASEFIELD_MIXED_MODEL_HH
#define PHASEFIELD_MIXED_MODEL_HH

#include<dune/fem/io/parameter.hh>

#include "phasefieldfilter.hh"

template<class Problem>
class MixedModel
{
 typedef Problem ProblemType;
 typedef typename ProblemType::ThermodynamicsType ThermodynamicsType;
 typedef double RangeFieldType; 
 //contructor
  public:
 MixedModel( const ProblemType& problem):
   problem_(problem),
   penalty_(Dune::Fem::Parameter::getValue<double>("phasefield.penalty"))
  {}

  
  inline void  muSource( const RangeFieldType rho1,
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
  
  private:
    const ProblemType& problem_;
    double penalty_; 

};

template<class Problem> 
inline void MixedModel<Problem>
::muSource( const RangeFieldType rho1,
            const RangeFieldType rho2,
            const RangeFieldType phi,
            RangeFieldType& mu) const
{
  return problem_.thermoDynamics().chemaicalPotential(rho1,phi);
}



template<class Problem> 
inline void MixedModel<Problem>
::tauSource( const RangeFieldType phi1,
            const RangeFieldType phi2,
            const RangeFieldType rho,
            RangeFieldType& tau) const
{
  return problem_.thermoDynamics().reactionSource(rho,phi1);
}




#endif


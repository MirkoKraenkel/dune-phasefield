#ifndef DUNE_PHASEFIELD_SOURCETERMS_HH
#define DUNE_PHASEFIELD_SOURCETERMS_HH

#include <cmath>

struct SourceTerms
{


  template<class DomainType,class RangeType>
  static void nstkSource(const DomainType& x,
                        const double t, 
                        const double delta,
                        const double c,
                        RangeType& res) 
  {
    res=0;
  }   

  template<class DomainType, class RangeType>
  static void acSource(const DomainType x,
                        const double t,
                        const double delta,
                        const double c,
                        RangeType& res) 
  {
    
    res=0;
  }

};
#endif

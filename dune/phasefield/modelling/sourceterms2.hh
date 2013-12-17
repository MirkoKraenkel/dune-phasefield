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
                                                 
    double t22;
    double t3;
    double t4;
    double t5;
    double t6;
    double t7;
    double t8;
    
    t3 = 1/delta;
    t4 = (x+c*t)*t3;
    t5 = sinh(t4);
    t6 = cosh(t4);
    t7 = t6*t6;
    t8 = t7*t7;
    t22 = t8*t8;
    res[1]=(-0.625E-3*(192.0*t5*t8*t6+540.0*t5*t6-704.0*t8*t7-8.0*t8-7112.0*t7-21.0+120.0*t5*t7*t6)/t22*t3);

    res[0]=0;
    for(int i=2;i<RangeType::dimension;++i)
      res[i]=0;
  }   

  template<class DomainType, class RangeType>
  static void acSource(const DomainType x,
                        const double t,
                        const double delta,
                        const double c,
                        RangeType& res) 
  {
    
    double t12;
    double t4;
    double t5;
    double t6;
    double t8;
              
    t4 = (x+c*t)/delta;
    t5 = cosh(t4);
    t6 = t5*t5;
    t8 = sinh(t4);
    t12 = t6*t6;
    res[RangeType::dimension-1]=(0.1875E-1*(-26.0*t6+1.0+30.0*t8*t5)/t12/t6);
    
    for(int i=0;i<RangeType::dimension-1;++i)
      res[i]=0.;
  }

};
#endif

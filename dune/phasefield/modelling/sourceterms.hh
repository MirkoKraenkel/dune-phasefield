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
    double t1;
    double t10;
    double t11;
    double t21;
    double t4;
    double t5;
    double t6;
    t1 = 1/delta;
    t4 = (x+c*t)*t1;
    t5 = sinh(t4);
    t6 = cosh(t4);
    t10 = t6*t6;
    t11 = t10*t10;
    t21 = t11*t11;
    res[0]=0;
    res[1]=(-0.125E-2*t1*(t5+15.0*t6)*(8.0*t5*t11+4.0*t5*t10+18.0*t5-24.0*t11*t6-225.0*t6)/t21);
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
    for(int i=0;i<RangeType::dimension;++i)
      res[i]=0.;
 
    res[RangeType::dimension-1]=(0.1875E-1*(-26.0*t6+1.0+30.0*t8*t5)/t12/t6);
                          
//return  0.;
  }









};
#endif

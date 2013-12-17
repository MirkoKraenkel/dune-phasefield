#include<cmath>
#include<iostream>

#include<dune/common/fvector.hh>
#include<dune/common/fmatrix.hh>

#include "phasefieldfilter.hh"


int main( int arc, const char** argv)
{
  enum{ dimension=2};
  enum{ dimRange=2*dimension+4};
  typedef Dune::FieldVector<double,dimRange> RangeType;
  typedef PhasefieldFilter<RangeType> FilterType;
    
  RangeType u;
  for(int i=0;i<dimRange;i++)
    u[i]=i*i;
  
  for(int i=0;i<dimRange;i++)
    std::cout<<"u "<<u[i]<<"\n";

  for(int i=0; i<dimension;i++)
    std::cout<<"sigma["<< i<<"]="<<FilterType::sigma(u,i)<<"\n";
  for(int i=0; i<dimension;i++)
    std::cout<<"v["<< i<<"]="<<FilterType::velocity(u,i)<<"\n";
  FilterType::tau(u)+=FilterType::phi(u);

  std::cout<<"phi="<<FilterType::phi(u)<<"\n";
  std::cout<<"rho="<<FilterType::rho(u)<<"\n";
  std::cout<<"tau="<<FilterType::tau(u)<<"\n";
  std::cout<<"mu="<<FilterType::mu(u)<<"\n";
  


}

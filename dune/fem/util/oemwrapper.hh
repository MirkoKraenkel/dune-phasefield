#ifndef OEMWRAPPER_HH
#define OEMWRAPPER_HH

#include<dune/fem/operator/common/operator.hh>
#include<dune/fem/operator/common/automaticdifferenceoperator.hh>

#if DGSCHEME
//#include<dune/fem/operator/assembled/mixedoperator.hh>
#elif FEMSCHEME
#include<dune/fem/operator/assembled/femoperator.hh>
#endif

template<class DomainFunction,class RangeFunction=DomainFunction>
class OEMWrapper
:public Dune::Fem::AutomaticDifferenceLinearOperator<DomainFunction ,RangeFunction >
{

   typedef OEMWrapper<DomainFunction,RangeFunction> ThisType;
  public:
   typedef Dune::Fem::AutomaticDifferenceLinearOperator<DomainFunction ,RangeFunction> BaseType;
 
  
   typedef typename BaseType::DomainFunctionType DomainFunctionType;
   typedef typename BaseType::RangeFunctionType RangeFunctionType;
   typedef typename BaseType::RangeFieldType RangeFieldType;
   typedef typename BaseType::DomainFieldType DomainFieldType;
  
   typedef typename DomainFunctionType::DiscreteFunctionSpaceType DomainSpaceType;
   typedef typename RangeFunctionType::DiscreteFunctionSpaceType RangeSpaceType;

   friend  class Dune::Fem::AutomaticDifferenceOperator<DomainFunctionType,RangeFunctionType,ThisType>;
   

    typedef  ThisType MatrixType;
    typedef Dune::Fem::Operator<DomainFunction,RangeFunction> OperatorType; 

  public:
    OEMWrapper(const std::string &name, 
                const DomainSpaceType& dSpace,
                const RangeSpaceType& rSpace):
    BaseType(name,dSpace,rSpace),
   // name_(name),
    //op_(nullptr),
    domainspc_(dSpace),
    rangespc_(rSpace)
  {}

  void multOEM(const double* a, double* b) const;

  double ddotOEM(const double* a, const double* b) const;


  using BaseType::operator();
//  inline void operator()(const DomainFunctionType& u, RangeFunctionType& w) const 
  //{
    //(*op_)(u,w);
  //}

 MatrixType& matrix(){ return *this;}
 MatrixType& systemMatrix(){ return *this;}


  private:
  using BaseType::op_;
  using BaseType::name_;
  //const std::string &name_;
  //const OperatorType* op_;
   const DomainSpaceType& domainspc_; 
  const RangeSpaceType& rangespc_;
};

template<class DomainFunction,class RangeFunction>
void OEMWrapper<DomainFunction,RangeFunction>::multOEM(const double* a, double* b) const
{
  const DomainFunctionType arg("OEMWrapper::Arg",domainspc_,a);
  RangeFunctionType dest("OEMWrapper::Dest",rangespc_,b);
  (*this)(arg,dest);
}
template<class DomainFunction, class RangeFunction>
double OEMWrapper<DomainFunction,RangeFunction>::ddotOEM(const double* a ,const double* b) const
{
  const RangeFunctionType aFunc("OEMWrapper::a",rangespc_,a);
  const RangeFunctionType bFunc("OEMWrapper::b",rangespc_,b);
  return aFunc.scalarProductDofs(bFunc);
}
#endif// OEMWRAPPER_HH

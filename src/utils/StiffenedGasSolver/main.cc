#include<cmath>
#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include<valarray>
//#include "grid.hh"
class Function
{

  public:
    Function( double a1, double b1,double p1,double a2,double b2,double p2)
    :a1_(a1),
    b1_(b1),
    p1_(p1),
    a2_(a2),
    b2_(b2),
    p2_(p2){}

    inline void evaluate( double* x, double* fx)
    {
      fx[0]=a1_*log(x[0])-a2_*log(x[1])+b1_-b2_;
      fx[1]=a1_*x[0]-a2_*x[1]-p1_+p2_;
    }
 
    inline void jacobianInverse( double* x, double (&dfx)[4])
    {
      dfx[0]=-x[0]*x[1]/(a1_*(x[0]-x[1]));
      dfx[1]=x[0]/(a1_*(x[0]-x[1]));
      dfx[2]=-x[0]*x[1]/(a2_*(x[0]-x[1]));
      dfx[3]=x[1]/(a2_*(x[0]-x[1]));
      

    }


private:
  double a1_,b1_,p1_,a2_,b2_,p2_ ;
  
};


static double norm(double* u)
{
  return std::sqrt(u[0]*u[0]+u[1]*u[1]);
}
#if 1 

template<class FunctionImp>
void  newtonSolve( FunctionImp f,double epsilon,double* start, double* result)
{
  
  double update[2]; 
  double dfInv[4];
  
  f.evaluate(start,result);
  int iterations=0; 
   std::cout<<"norm="<<norm(result)<<" iter "<<iterations<<"\n";
 
  while( norm(result)>epsilon && iterations<100000)
  {
     f.jacobianInverse(start,dfInv);
 


    update[0]=dfInv[0]*result[0]+dfInv[1]*result[1];
    update[1]=dfInv[2]*result[0]+dfInv[3]*result[1];
    
    start[0]-=update[0];
    start[1]-=update[1];
    ++iterations;
    f.evaluate(start,result);
 
    std::cout<<"norm="<<norm(result)<<" iter "<<iterations<<"\n";

    
  }
  result[0]=start[0];
  result[1]=start[1];
}
#endif
int main(int arc, const char** argv)
{

  if( arc!=10)
  {
    std::cerr<<arc<<" USAGE: a1 b1 p1 a2 b2 p2 rho1 rho2 epsilon\n";
    abort();
  }
   
  
  double a1=atof(argv[1]);
  double b1=atof(argv[2]);
  double p1=atof(argv[3]);
  double a2=atof(argv[4]);
  double b2=atof(argv[5]);
  double p2=atof(argv[6]);
 double rho1=atof(argv[7]);
  double rho2=atof(argv[8]);


  double epsilon=atof(argv[9]);
  
  Function fstiff(a1,b1,p1,a2,b2,p2);
  
  double x[2]={rho1,rho2};
  double fx[2]={0.,0.};

  fstiff.evaluate(x,fx); 
   std::cout<<"F(x,y)=(" <<fx[0]<<" , "<<fx[1]<<")\n";
 
  newtonSolve(fstiff,epsilon,x,fx);
  std::cout<<"Result=(" <<fx[0]<<" , "<<fx[1]<<")\n";
  fstiff.evaluate( fx, x);
  std::cout<<"Residuum=(" <<x[0]<<" , "<<x[1]<<")\n";
 
  return 0;
}  



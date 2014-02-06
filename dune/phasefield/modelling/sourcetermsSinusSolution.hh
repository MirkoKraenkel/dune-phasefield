#ifndef DUNE_PHASEFIELD_SOURCETERMS_HH
#define DUNE_PHASEFIELD_SOURCETERMS_HH

#include <cmath>

struct SourceTerms
{


  template< class DomainType , class RangeType >
  static systemSource( const DomainType x,
                       const double t,
                       RangeType& res )
  {
    enum{dimension=RangeType::dimension};
    double cosx=std::cos(2*M_PI*x);
    double cost=std::cos(M_PI*t);
    double sinx=std::sin(2*M_PI*x);
    double sint=std::sin(M_PI*t);
    double phi=0.5*cosx*cost+0.5;
    
    //double rho=1.;
    //double dtrho=0.;
   // double drho=0.;
    double rho=0.5*cosx*cost+1;
    double drho=-M_PI*sinx*cost;
    double dtrho=-0.5*M_PI*cosx*sint;
    double v=sinx*cost;
    double dtv=-M_PI*sinx*sint;
    double dv=2*M_PI*cosx*cost;
    double laplacev=4*M_PI*M_PI*cost*sinx;

    //f = 2*(d_t phi-\Delta phi)
    double f=M_PI*( 4*M_PI*std::cos( M_PI*t ) - std::sin( M_PI*t ) );
    f*=cosx;

    // rhof=d_trho+div(rho v)
     double rhof=dtrho+drho*v+rho*dv;
    res[0]=rhof;

    //lapv= d_t v-\Delta v
    double lapv=M_PI*( 4*M_PI*std::cos( M_PI*t ) - std::sin( M_PI*t ) );
    lapv*=sinx;
   
    //tension= -\nabla\phi tau
    double tension=-sinx*M_PI*cost*(-1.*cosx*cost+cosx*cosx*cosx*cost*cost*cost+2*cosx*M_PI*M_PI*cost);   

    // rhodmu=\nabla mu
     double dmu=2*M_PI*cosx*sinx*cost*cost; 

     res[1]=(dtv+dmu)*rho+laplacev-tension;
     double dFdphi{0.};
     // transportphi=v\cdot\nabla\phi
     double transportphi=-M_PI*sinx*sinx*cost*cost;
     tauSource(1,phi,phi,dFdphi);
     
     double phiSource=f*0.5+dFdphi+transportphi;

     if(dimension==2)
     { 
       res[2]=0.;
       res[3]=phiSource; 
     }     
    
     res[2]=phiSource;
  
  }

  template<class DomainType,class RangeType>
  static void nstkSource(const DomainType& x,
                        const double t, 
                        const double delta,
                        const double c,
                        RangeType& res) 
  {
    
  }   

  template<class DomainType, class RangeType>
  static void acSource(const DomainType x,
                        const double t,
                        const double delta,
                        const double c,
                        RangeType& res) 
  {
  }

	void  tauSource(double rho,double phi,double phi, double& dFdphi ) const
	{ 
    //    return phi*phi;
   dFdphi=4*(2*phi*phi*phi-3*phi*phi+phi);
   // return 0;
  }
  








};
#endif

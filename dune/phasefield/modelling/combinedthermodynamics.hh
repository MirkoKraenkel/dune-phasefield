#ifndef COMBINEDTHERMODYNAMICS_HH
#define COMBINEDTHERMODYNAMICS_HH
// system include
#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>

#include "thermodynamicsinterface.hh"


template< class Phase1,class Phase2>
class CombinedThermodynamics
	:public Thermodynamics< CombinedThermodynamics< Phase1, Phase2> >
  {
    typedef Thermodynamics< CombinedThermodynamics< Phase1, Phase2 > > BaseType;
    typedef Phase1 Phase1Type;
    typedef Phase2 Phase2Type;
public:
  CombinedThermodynamics(Phase1Type phase1 ,Phase2Type phase2 ):
    phase1_(phase1),
    phase2_(phase2)
    {
     delta_=1.;
	   deltainv_=1.;
     deltainv_/=delta_;
    }

  inline void init() const
  {
  
    delta_=(Fem::Parameter::getValue<double>( "phasefield.delta" ));
    epsilon_=(Fem::Parameter::getValue<double>( "phasefield.epsilon",0.));
    deltainv_=1;
    deltainv_/=delta_; 
    rate_=(Fem::Parameter::getValue<double>( "phasefield.reactionrate",0.));
    std::cout<<"1/delta="<<deltainv_<<"\n";
    std::cout<<"epsilon="<<epsilon_<<"\n";
    std::cout<<"rate_"<<rate_<<"\n";
   

       
  }
  //free EnergyPart without gradients used for monitoring the free energy
  inline double helmholtz(double rho,double phi) const
  { 
    return deltainv_*rho*W(phi)+nu(phi)*phase1_.helmholtz(rho,rho)+nu(1-phi)*phase2_.helmholtz(rho,phi);
  }
  //derivative of f wrt \rho 
  inline double chemicalPotential(double rho,double phi) const
  {

    return deltainv_*W(phi)+nu(phi)*phase1_.chemicalPotential(rho,phi)+nu(1-phi)*phase2_.chemicalPotential(rho,phi);
  }



	 // derivative of f wrt to \phi
  inline double reactionSource(double rho,double phi) const
  { 
    std::cout<<"Reaction="<<rate_<<"\n";
    return rate_*rho*Wprime(phi)+nuprime(phi)*phase1_.helmholtz(rho,phi)-nuprime(1-rho)*phase2_.helmholtz(rho,phi);
  }
  // thermodynamic pressure -f+\rho*f|_\rho 
  inline double  pressure( double rho, double phi) const
  {
    double ret;
   ret= nu(phi)*phase1_.pressure(rho,phi)+nu(1-phi)*phase2_.pressure(rho,phi);
    std::cout<<"pressure"<<ret<<" rho "<<rho<<" phi "<<phi<<"\n";
    return ret; }

  inline double a(double rho, double phi) const
  {  
    return nu(phi)*phase1_.a(rho,0.)+nu(1-phi)*phase2_.a(rho,0);
  }
//not in the interface
  inline double W(double phi) const
  { 
    return (phi-1)*(phi-1)*phi*phi;
  }
  inline double Wprime(double phi) const
  {
    return 4*std::pow(phi,3)-6*std::pow(phi,2)+2*phi;
  }

  inline double nu(double phi) const
  { 
    double t1;
    double t2;
    t1 = phi*phi;
    t2 = t1*t1;
    return(6.0*t2*phi-15.0*t2+10.0*t1*phi);
  }
  inline double nuprime(double phi) const
  {
    double t2;
    t2 = phi*phi;
    return(2.0*phi-6.0*t2+4.0*t2*phi);
  }
 
  inline double delta()const
  {
    return delta_;
  }

  inline double mu1() const
  {
    std::cout<<"combinedtherm mu1="<<epsilon_<<"\n";
    return epsilon_;
  }

private:
  mutable double delta_;  
  mutable double deltainv_;
  mutable double rate_;
  mutable double epsilon_;
  Phase1Type phase1_;
  Phase2Type phase2_;

};

class LiquidDynamics
: public Thermodynamics<LiquidDynamics>
{
    typedef Thermodynamics<LiquidDynamics> BaseType;
  public:
  
  //free EnergyPart without gradients used for monitoring the free energy
  inline double helmholtz(double rho,double phi) const
  { 
    double t2;
   t2 = log(rho);
       return(1.0-0.4*rho+0.14E1*rho*t2);

  }
  //derivative of f wrt \rho 
  inline double chemicalPotential(double rho,double phi) const
  {
    double t1;
    t1 = log(rho);
    return(0.1E1+0.14E1*t1);
  }
  // derivative of f wrt to \phi
  inline double reactionSource(double rho,double phi) const
  { 
    std::cout<<"reactionSource makes no Sense for pure phase\n";
    abort();
  }
  // thermodynamic pressure -f+\rho*f|_\rho 
  inline double  pressure( double rho, double phi) const
  {
  return(-1.0+0.14E1*rho);
    }
  //dp/drho
  inline double a(double rho,double phi) const
  {
    return 0.4;
  } 
};
class VapourDynamics
: public Thermodynamics<VapourDynamics>
{
    typedef Thermodynamics<VapourDynamics> BaseType;
public:
  
  //free EnergyPart without gradients used for monitoring the free energy
  inline double helmholtz(double rho,double phi) const
  { 
    double t1;
    t1 = log(rho);
    return(0.16E1*rho*t1+0.15E1-0.15E1*rho);
  }
 //derivative of f wrt \rho 
  inline double chemicalPotential(double rho,double phi) const
  {
     double t1;
     t1 = log(rho);
     return(0.16E1*t1+0.1);
  }
  // derivative of f wrt to \phi
  inline double reactionSource(double rho,double phi) const
  { 
    std::cout<<"reactionSource makes no Sense for pure phase\n";
    abort();
  }
  // thermodynamic pressure -f+\rho*f|_\rho 
  inline double  pressure( double rho, double phi) const
  {
  return(-2.0+0.16E1*rho);

  }
 //dp/drho
  inline double a( double rho, double phi) const
  {
    return 1.6;
  }

};


#endif //  THERMODYNAMICS_HH


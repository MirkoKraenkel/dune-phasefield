#ifndef DUNEPHASEFIELD_KORTEWEGPHYSICS_HH
#define DUNEPHASEFIELD_KORTEWEGPHYSICS_HH
#define MIN_RHO 1e-10
//#warning PHYSICS INCLUDED
// system include
#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>

// dune-fem include
#include <dune/fem/io/parameter.hh>

#include "thermodynamicsinterface.hh"

namespace Dune{


template<int dimDomain, class Thermodynamics>
class PhasefieldPhysics
{
   typedef Thermodynamics ThermodynamicsType;
 
  public:
    enum { dimRange = dimDomain + 1 };
    enum { dimThetaRange =  1 };
    enum { dimThetaGradRange = dimThetaRange*dimDomain };
    enum { dimGradRange = dimRange * dimDomain };
    
    typedef double RangeFieldType;

    typedef FieldVector< double, dimDomain >                  DomainType;
    typedef FieldVector< double, dimDomain - 1 >              FaceDomainType;
    typedef FieldVector< double, dimRange >                   RangeType;
    typedef FieldVector< double, dimThetaRange >              ThetaRangeType;
    typedef FieldVector< double, dimGradRange >               GradientType;
    typedef FieldVector< double, dimThetaGradRange >          ThetaGradientRangeType;
    typedef FieldMatrix< double, dimRange, dimDomain >        JacobianRangeType;                          
    typedef FieldMatrix< double, dimRange, dimDomain >        FluxRangeType;
    typedef FieldVector< double, dimGradRange >               GradientRangeType;

   typedef FieldMatrix< double, dimThetaRange, dimDomain >    ThetaJacobianRangeType;
   typedef FieldMatrix< double, dimGradRange, dimDomain >    JacobianFluxRangeType;

  protected:
    const ThermodynamicsType& thermoDynamics_;
  public:
  PhasefieldPhysics(const ThermodynamicsType& thermodyn):
    thermoDynamics_(thermodyn),
    delta_(Dune::Fem::Parameter::getValue<double>("phasefield.delta")),
    deltaInv_(1./delta_)
  {
  }
 
  inline void conservativeToPrimitive( const RangeType& cons, RangeType& prim ) const;
 
  template< class JacobianRangeImp >
  inline void totalEnergy( const RangeType& cons, 
                           const JacobianRangeImp& grad,
                           double& kin,
                           double& therm,
                           double& total,
                           double& surf) const;
  
  inline void chemicalPotential( const RangeType& cons, 
										   double& mu ) const;

	inline void pressure( const RangeType& cons, 
																	 double& p ) const;
  
  inline void analyticalFlux( const RangeType& u, JacobianRangeType& f ) const;
  
  inline void jacobian( const RangeType& u, JacobianFluxRangeType& a) const;

  inline double maxSpeed( const DomainType& n, const RangeType& u ) const;
  
  inline double stiffSource(const DomainType& x,
                            const double time,
                            const RangeType& u,
								            const GradientRangeType& du,
								            const ThetaRangeType& theta,
								            const ThetaJacobianRangeType& dtheta,
								            const JacobianRangeType& jacU,
                            RangeType& f) const;
 
  inline double stiffSource(const DomainType& x,
                            const double time,
                            const RangeType& u,
								            const GradientRangeType& du,
                            RangeType& f) const;
 
  template< class JacobianRangeImp >
	inline void diffusion( const RangeType& u,
												 const JacobianRangeImp& du,
												 JacobianRangeType& f ) const;
  template< class JacobianRangeImp >
	inline void boundarydiffusion( const RangeType& u,
												 const JacobianRangeImp& du,
												 JacobianRangeType& f ) const;
  
  //f|phi-div(f|nabla phi)
  template< class JacobianRangeImp >
	inline void korteweg( const RangeType& u,
												 const JacobianRangeImp& du,
												 ThetaJacobianRangeType& f ) const;
   //f|phi-div(f|nabla phi)
  template< class JacobianRangeImp >
	inline void boundarykortewegn( const RangeType& u,
												 const JacobianRangeImp& du,
												 ThetaJacobianRangeType& f ) const;
 

  public:

	inline double delta()    const { return delta_;}
	inline double deltaInv() const { return deltaInv_;}
	inline double mu1()      const { abort(); return 1.;}
	inline double mu2()      const { abort(); return 1.;}
  
  protected:
	const double delta_; 
	double deltaInv_;
 };
}

#if NSK 
#include "physics/physicsnsk_inline1d.hh"
//#include "physicsnsk_inline2d.hh"
#else
#error
#endif
	//end namspace DUNE

#endif // DUNEPHASEFIELD_PHYSICS_HH

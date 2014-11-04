#ifndef DUNEPHASEFIELD_PHYSICS_HH
#define DUNEPHASEFIELD_PHYSICS_HH
#define MIN_RHO 1e-10
//#warning PHYSICS INCLUDED
// system include
#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>

// dune-fem include
#include <dune/fem/io/parameter.hh>


namespace Dune{


template<int dimDomain, class Thermodynamics>
class PhasefieldPhysics
{
   typedef Thermodynamics ThermodynamicsType;
 
  public:
    enum{ phaseId = dimDomain+1} ;
    enum { dimRange = dimDomain + 2 };
    enum { dimThetaRange =  2 };
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
                           double& surf,
                           double& total) const;
  
  inline void chemPotAndReaction( const RangeType& cons, 
																	const JacobianRangeType& du,
                                  double& mu,
																	double& reaction ) const;

	inline void pressureAndReaction( const RangeType& cons, 
																	 double& p,
																	 double& reaction ) const;
  
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
	inline void allenCahn( const RangeType& u,
												 const JacobianRangeImp& du,
												 ThetaJacobianRangeType& f ) const;
   //f|phi-div(f|nabla phi)
  template< class JacobianRangeImp >
	inline void boundaryallenCahn( const RangeType& u,
												 const JacobianRangeImp& du,
												 ThetaJacobianRangeType& f ) const;
 
  template< class JacobianRangeImp>
  inline void tension( const RangeType& u,
                       const JacobianRangeImp& du,
                       GradientRangeType& tens) const;

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

#if WELLBALANCED
#if NONCONTRANS
#include "physicswb2_inline1d.hh"
#else
#include "physicswb_inline1d.hh"
#endif
#include "physicswb_inline2d.hh"
#else
#include "physics_inline1d.hh"
#include "physics_inline2d.hh"
#endif
	//end namspace DUNE

#endif // DUNEPHASEFIELD_PHYSICS_HH

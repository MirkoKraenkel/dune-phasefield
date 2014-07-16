#ifndef DUNEPHASEFIELD_THERMODYNAMICSINTERFACE_HH
#define DUNEPHASEFIELD_THERMODYNAMICSINTERFACE_HH
// system include
#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>

#include <dune/fem/io/parameter.hh>

using namespace Dune;

template<class  Impl >
class Thermodynamics
{

public:
  Thermodynamics()
  {
  }

  //factor in front of the double well
  inline double h1(  double rho ) const
  {
    return asImp.h( rho );
  }
  
  inline double h1prime(  double rho ) const
  {
    return asImp.h( rho );
  }
  
  //factor in front of the gradient term
  inline double h2( double rho ) const
  {
    return asImp.h2( rho );
  }
  
  inline double h2prime( double rho ) const
  {
    return asImp.h2prime( rho );
  }
  
  inline double reactionFactor( double rho) const
  {
    return asImp.reationFacor();
  }

  
  //free EnergyPart without gradients used for monitoring the free energy
  inline double helmholtz( double rho , double phi ) const
  { 
		return asImp.helmholtz( rho , phi );
	}

  //derivative of f wrt \rho 
  inline double chemicalPotential( double rho , double phi ) const
  {
    return asImp.chemicalPotential( rho , phi );
  }
  // derivative of f wrt to \phi
  inline double reactionSource( double rho , double phi ) const
  { 
    return asImp().reactionSource( rho , phi );
  }
  // thermodynamic pressurw -f+\rho*f|_\rho 
  inline double  pressure( double rho , double phi ) const
  {
	  return asImp().pressure( rho , phi ); 
	}

  inline double a( double rho , double phi ) const
  {
    return asImp().a( rho , phi );
  }

  Impl& asImp()
  {
    return static_cast<Impl>(*this);
  }
	protected:


 
  
};





#endif //  DUNEPHASEFIELD_THERMODYNAMICS_HH


#ifndef DUNEPHASEFIELD_THERMODYNAMICS_HH
#define DUNEPHASEFIELD_THERMODYNAMICS_HH
// system include
#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>

#include <dune/fem/io/parameter.hh>

using namespace Dune;

class Thermodynamics
{

public:
  Thermodynamics()
  {
	}

  //factor in front of the double well
  inline double h( double& rho ) const
  {
    return 1.;
  }
  //factor in front of the gradient term
  inline double h2( double& rho ) const
  {
    return 1.;
  }
  
  inline double reactionFactor( ) const
  {
    return 1.;
  }

  //free EnergyPart without gradients used for monitoring the free energy
  inline double helmholtz(double& rho,double& phi) const
  { 
		return 0.;
	}
  //derivative of f wrt \rho 
  inline double chemicalPotential(double& rho,double& phi) const
  {
    return 0.;
  }
  // derivative of f wrt to \phi
  inline double reactionSource(double& rho,double& phi) const
  { 
    return 0.;
  }
  // thermodynamic pressurw -f+\rho*f|_\rho 
  inline double  pressure( double& rho, double& phi) const
  {
	  return 0.; 
	}
};





#endif //  DUNEPHASEFIELD_THERMODYNAMICS_HH


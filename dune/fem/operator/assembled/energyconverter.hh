#ifndef MIXEDSCHEME_ENEREGYCONVERTER_HH
#define MIXEDSCHEME_ENEREGYCONVERTER_HH

// system include
#include <sstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <time.h>

#include <dune/fem/quadrature/cachingquadrature.hh>
#include "phasefieldfilter.hh"

/** 
 *  \param[in] dF The discrete function of conservative variables
 *  \param[in] model The analytical model
 */
template< class DiscreteFunctionType, class ModelType, class EnergyFunctionType >
double energyconverter( const DiscreteFunctionType& dF, 
          	            const ModelType& model,
                        EnergyFunctionType& energyDF,
                        double& kineticEnergy,
                        double& thermodynamicEnergy, 
                        double& surfaceEnergy)
{
  typedef typename DiscreteFunctionType::Traits::DiscreteFunctionSpaceType 
    DiscreteFunctionSpaceType;
  typedef typename EnergyFunctionType::Traits::DiscreteFunctionSpaceType 
    EnergyFunctionSpaceType;
  typedef typename DiscreteFunctionSpaceType::Traits::GridPartType GridPartType;
  typedef typename GridPartType::GridType GridType;
  typedef typename GridType :: template Codim<0> :: Entity   Entity;
  typedef typename GridType :: template Codim<0> :: Geometry Geometry;
  typedef typename GridPartType::template Codim<0>::IteratorType Iterator;
  
  typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
  typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
  typedef typename EnergyFunctionSpaceType::RangeType EnergyRangeType;
  
  const DiscreteFunctionSpaceType& space =  dF.space();
	energyDF.clear();
  
  typedef typename DiscreteFunctionType::LocalFunctionType LocalFuncType;
  typedef typename EnergyFunctionType::LocalFunctionType EnergyLocalFuncType;


  RangeType vu(0.0);
  EnergyRangeType total(0.0);  
  EnergyRangeType kin(0.0);
  EnergyRangeType therm(0.0); 
  EnergyRangeType surf(0.0);
  kineticEnergy=0.;
  thermodynamicEnergy=0;
  double integral=0.;
  
  Iterator it    = space.begin();
  Iterator endit = space.end();
 
	// if empty grid, do nothing 
  if( it == endit ) return -42. ;
  
  for( ; it != endit ; ++it) 
  {
    // get entity 
    const Entity& entity = *it ;
    const Geometry& geo = entity.geometry(); 
   
		const double volume = geo.volume();
    // Get quadrature rule for L2 projection
    Dune::Fem::CachingQuadrature< GridPartType, 0 > quad( entity, 2*space.order()+1 );
  
    LocalFuncType   localU   = dF.localFunction( entity );
		EnergyLocalFuncType energyLF = energyDF.localFunction( entity );

    const int quadNop = quad.nop();
    for(int qP = 0; qP < quadNop; ++qP) 
    {
      const DomainType& xgl = geo.global( quad.point(qP) );
      // evaluate  variables 
      localU.evaluate( quad[qP], vu );

			model.totalEnergy(xgl,vu,kin[0],therm[0],total[0],surf[0]);
      total*= quad.weight(qP);
			kin*=quad.weight(qP);
      therm*=quad.weight(qP); 
      energyLF.axpy(quad[qP],total);
      kineticEnergy+=kin*volume;
      thermodynamicEnergy+=therm*volume;
      integral+=total*volume;
       
    }
  
  }
	return integral;
    
}


#endif

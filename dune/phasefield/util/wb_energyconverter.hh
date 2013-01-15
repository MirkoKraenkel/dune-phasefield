#ifndef WBENERGY_HH
#define WBENERGY_HH

// system include
#include <sstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <time.h>


/** 
 *  \param[in] consDF The discrete function of conservative variables
 *  \param[in] model The analytical model
 *  \param[out] primDF The discrete function of primitive variables
 */
template< class ConsDiscreteFunctionType,class GradientFunctionType, class ModelType, class EnergyFunctionType >
void energyconverter( const ConsDiscreteFunctionType& consDF, 
                      const GradientFunctionType& gradDF,
          	          const ModelType& model,
                      EnergyFunctionType& energyDF) 
{
  typedef typename ConsDiscreteFunctionType::Traits::DiscreteFunctionSpaceType 
    ConsDiscreteFunctionSpaceType;
  typedef typename GradientFunctionType::Traits::DiscreteFunctionSpaceType 
    GradDiscreteFunctionSpaceType;

  typedef typename EnergyFunctionType::Traits::DiscreteFunctionSpaceType 
    EnergyFunctionSpaceType;
  typedef typename ConsDiscreteFunctionSpaceType::Traits::GridPartType GridPartType;
  typedef typename GridPartType::GridType GridType;
  typedef typename GridType :: template Codim<0> :: Entity   Entity;
  typedef typename GridType :: template Codim<0> :: Geometry Geometry;
  typedef typename GridPartType::template Codim<0>::IteratorType Iterator;
  
  typedef typename ConsDiscreteFunctionSpaceType::DomainType DomainType;
  typedef typename ConsDiscreteFunctionSpaceType::RangeType ConsRangeType;
  typedef typename GradDiscreteFunctionSpaceType::RangeType GradRangeType;
  typedef typename EnergyFunctionSpaceType::RangeType EnergyRangeType;
  
  const ConsDiscreteFunctionSpaceType& space =  consDF.space();
	energyDF.clear();
  
  typedef typename ConsDiscreteFunctionType::LocalFunctionType ConsLocalFuncType;
  typedef typename GradientFunctionType::LocalFunctionType GradLocalFuncType;
  typedef typename EnergyFunctionType::LocalFunctionType EnergyLocalFuncType;

  ConsRangeType cons(0.0); 
  GradRangeType grad(0.0);
  EnergyRangeType energy(0.0);
 
  Iterator it    = space.begin();
  Iterator endit = space.end();
 
	// if empty grid, do nothing 
  if( it == endit ) return ;
  
  for( ; it != endit ; ++it) 
  {
    // get entity 
    const Entity& entity = *it ;
    const Geometry& geo = entity.geometry(); 
   
		const double volume = geo.volume();
    // Get quadrature rule for L2 projection
    Dune::Fem::CachingQuadrature< GridPartType, 0 > quad( entity, 2*space.order()+1 );
  
    ConsLocalFuncType   consLF   = consDF.localFunction( entity );
    GradLocalFuncType   gradLF   = gradDF.localFunction( entity );
		EnergyLocalFuncType energyLF = energyDF.localFunction( entity );

    const int quadNop = quad.nop();
    for(int qP = 0; qP < quadNop; ++qP) 
    {
      const DomainType& xgl = geo.global( quad.point(qP) );
      // evaluate conservative variables and gradients
      consLF.evaluate( quad[qP], cons );
      gradLF.evaluate( quad[qP], grad );
			
			model.totalEnergy(xgl, cons, grad, energy);
      energy*=  quad.weight(qP)*volume;
			energyLF.axpy(quad[qP],energy);
    }
  
}
	
    
}

#endif

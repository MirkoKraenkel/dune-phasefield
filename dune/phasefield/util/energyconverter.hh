#ifndef ENERGY_HH
#define ENERGY_HH

// system include
#include <sstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <time.h>


/** \brief 
 *  \param[in] consDF The discrete function of conservative variables
 *  \param[in] gradDF DG-Gradient of variables.
 *  \param[in] model The analytical model
 *  \param[out] energDF discrete total energy
 *
 **/
template< class ConsDiscreteFunctionType,class GradientFunctionType, class ModelType, class EnergyFunctionType >
double energyconverter( const ConsDiscreteFunctionType& consDF, 
                        const GradientFunctionType& gradDF,
          	            const ModelType& model,
                        EnergyFunctionType& energyDF,
                        double& kineticEnergy,
                        double& surfaceEnergy)
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
  typedef typename GridPartType :: template Codim<0>::IteratorType Iterator;
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
  EnergyRangeType kin(0.0);
  EnergyRangeType total(0.0);
  EnergyRangeType surf(0.0);
  kineticEnergy=0.;
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
  
    ConsLocalFuncType consLF = consDF.localFunction( entity );
    GradLocalFuncType  gradLF = gradDF.localFunction(entity);
    EnergyLocalFuncType energyLF = energyDF.localFunction( entity );

    const int quadNop = quad.nop();
    for(int qP = 0; qP < quadNop; ++qP) 
    {
      const DomainType& xgl = geo.global( quad.point(qP) );
      // evaluate conservative variables
      consLF.evaluate( quad[qP], cons );
      gradLF.evaluate( quad[qP], grad );
       
      model.totalEnergy(xgl,cons,grad,kin,total,surf);
      
      total*= quad.weight(qP);
      kin*=quad.weight(qP);
      energyLF.axpy(quad[qP],total);
      integral+=total*volume;
      kineticEnergy+=kin*volume;
      surfaceEnergy+=surf*quad.weight(qP)*volume;
    }
  
  }
  return integral;  
}

#endif

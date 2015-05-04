#ifndef DUNE_PHASEFIELD_SPLITENERGYCONVERTER_HH
#define DUNE_PHASEFIELD_SPLITENERGYCONVERTER_HH
// system include
#include <sstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <time.h>

#include <dune/fem/quadrature/cachingquadrature.hh>

/** 
 *  \param[in] dF The discrete function of conservative variables
 *  \param[in] model The analytical model
 */
template< class AcDiscreteFunctionType,class NvStDiscreteFunctionType, class ModelType, class ScalarFunctionType >
double splitenergyconverter(   const NvStDiscreteFunctionType& nvstF,
                        const AcDiscreteFunctionType& acF, 
          	            const ModelType& model,
                        ScalarFunctionType& energyDF,
                        ScalarFunctionType& pressureDF,
                        double& kineticEnergy,
                        double& thermodynamicEnergy, 
                        double& surfaceEnergy)
{
 
  typedef typename AcDiscreteFunctionType::Traits::DiscreteFunctionSpaceType 
    DiscreteFunctionSpaceType;
  typedef typename ScalarFunctionType::Traits::DiscreteFunctionSpaceType
    ScalarFunctionSpaceType;
  typedef typename DiscreteFunctionSpaceType::Traits::GridPartType GridPartType;
  typedef typename GridPartType::GridType GridType;
  typedef typename GridType :: template Codim<0> :: Entity   Entity;
  typedef typename GridType :: template Codim<0> :: Geometry Geometry;
  typedef typename GridPartType::template Codim<0>::IteratorType Iterator;
  
  typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
  typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
  typedef typename ScalarFunctionSpaceType::RangeType ScalarRangeType;
  static const int dim=GridType::dimensionworld;
  static const int  dimRange=DiscreteFunctionSpaceType::dimRange;

  const DiscreteFunctionSpaceType& space =  acF.space();
	energyDF.clear();
  pressureDF.clear();

  typedef typename AcDiscreteFunctionType::LocalFunctionType AcLocalFuncType;
  typedef typename NvStDiscreteFunctionType::LocalFunctionType NvStLocalFuncType;
  
  typedef typename ScalarFunctionType::LocalFunctionType ScalarLocalFuncType;

  RangeType vuAc(0.);
  RangeType vuNvSt(0.);
  Dune::FieldVector<double,2*dimRange> vuComb(0.);
  ScalarRangeType total(0.0);
  ScalarRangeType kin(0.0);
  ScalarRangeType therm(0.0);
  ScalarRangeType surf(0.0);
  ScalarRangeType press(0.);
  kineticEnergy=0.;
  thermodynamicEnergy=0;
  surfaceEnergy=0;
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
  
    AcLocalFuncType   localAc   = acF.localFunction( entity );
		NvStLocalFuncType localNvSt = nvstF.localFunction( entity );
    ScalarLocalFuncType energyLF = energyDF.localFunction( entity );
    ScalarLocalFuncType pressureLF= pressureDF.localFunction( entity );
    const int quadNop = quad.nop();
    for(int qP = 0; qP < quadNop; ++qP) 
    {
      const DomainType& xgl = geo.global( quad.point(qP) );
      // evaluate  variables 
      localAc.evaluate( quad[qP], vuAc );
      localNvSt.evaluate( quad[qP], vuNvSt);
      vuComb[0]=vuNvSt[0];
      vuComb[dim+1]=vuAc[0];
      vuComb[dim+2]=vuNvSt[dim+1];
      vuComb[dim+3]=vuAc[1];
      for( int ii=0;ii<dim;++ii)
        {
          vuComb[ 1+ii ]=vuNvSt[1+ii];
          vuComb[ dim+4+ii]=vuAc[2+ii];
        }
      model.totalEnergy(xgl,vuComb,kin[0],therm[0],total[0],surf[0]);
      press[0] =model.pressure(vuComb[0],vuComb[dim+1]);

      total*= quad.weight(qP);
			press*= quad.weight(qP);
      kin*=quad.weight(qP);
      therm*=quad.weight(qP); 
      energyLF.axpy(quad[qP],total);
      pressureLF.axpy(quad[qP],press);
      kineticEnergy+=kin*volume;
      thermodynamicEnergy+=therm*volume;
      integral+=total*volume;
      surfaceEnergy+=surf*quad.weight(qP)*volume;       
    }
  
  }
	return integral;
    
}


#endif

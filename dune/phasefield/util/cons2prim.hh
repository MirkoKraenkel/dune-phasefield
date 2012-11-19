#ifndef CONS2PRIM_HH
#define CONS2PRIM_HH

// system include
#include <sstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <time.h>


/** \brief converts a discrete function of conservative variables to
 *    a discrete function of primitive variables
 *  
 *  \param[in] consDF The discrete function of conservative variables
 *  \param[in] model The analytical model
 *  \param[out] primDF The discrete function of primitive variables
 */
template< class ConsDiscreteFunctionType, class ModelType, class PrimDiscreteFunctionType >
void setupAdditionalVariables( const ConsDiscreteFunctionType& consDF, 
                               const ModelType& model,
                               PrimDiscreteFunctionType& primDF ) 
{
  typedef typename ConsDiscreteFunctionType::Traits::DiscreteFunctionSpaceType 
    ConsDiscreteFunctionSpaceType;
  typedef typename PrimDiscreteFunctionType::Traits::DiscreteFunctionSpaceType 
    PrimDiscreteFunctionSpaceType;
  typedef typename ConsDiscreteFunctionSpaceType::Traits::GridPartType GridPartType;
  typedef typename ConsDiscreteFunctionSpaceType::Traits::GridType GridType;
  typedef typename GridType :: template Codim<0> :: Entity   Entity;
  typedef typename GridType :: template Codim<0> :: Geometry Geometry;
  typedef typename ConsDiscreteFunctionSpaceType::Traits::IteratorType Iterator;
  typedef typename ConsDiscreteFunctionSpaceType::DomainType DomainType;
  typedef typename ConsDiscreteFunctionSpaceType::RangeType ConsRangeType;
  typedef typename PrimDiscreteFunctionSpaceType::RangeType PrimRangeType;
  
  const ConsDiscreteFunctionSpaceType& space =  consDF.space();
  
  primDF.clear();
  
  typedef typename ConsDiscreteFunctionType::LocalFunctionType ConsLocalFuncType;
  typedef typename PrimDiscreteFunctionType::LocalFunctionType PrimLocalFuncType;
  
  ConsRangeType cons(0.0);
  PrimRangeType prim(0.0);
  
  Iterator it    = space.begin();
  Iterator endit = space.end();

  // if empty grid, do nothing 
  if( it == endit ) return ;
  
  for( ; it != endit ; ++it) 
  {
    // get entity 
    const Entity& entity = *it ;
    const Geometry& geo = entity.geometry(); 

    // Get quadrature rule for L2 projection
    Dune::Fem::CachingQuadrature< GridPartType, 0 > quad( entity, 2*space.order()+1 );
  
    ConsLocalFuncType consLF = consDF.localFunction( entity );
    PrimLocalFuncType primLF = primDF.localFunction( entity );

    const int quadNop = quad.nop();
    for(int qP = 0; qP < quadNop; ++qP) 
    {
      const DomainType& xgl = geo.global( quad.point(qP) );
      // evaluate conservative variables
      consLF.evaluate( quad[qP], cons );
			model.conservativeToPrimitive( xgl, cons, prim );
   
      prim *=  quad.weight(qP);
      primLF.axpy( quad[qP] , prim );
    }
  }
}

#endif

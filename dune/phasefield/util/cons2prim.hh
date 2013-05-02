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
template< class ConsDiscreteFunctionType,class GradientDiscreteFunctionType, class ModelType, class PrimDiscreteFunctionType >
void setupAdditionalVariables( const ConsDiscreteFunctionType& consDF, 
                               const GradientDiscreteFunctionType& sigmaDF,   
                               const ModelType& model,
                               PrimDiscreteFunctionType& primDF ) 
{
  typedef typename ConsDiscreteFunctionType::Traits::DiscreteFunctionSpaceType 
    ConsDiscreteFunctionSpaceType;
  typedef typename PrimDiscreteFunctionType::Traits::DiscreteFunctionSpaceType 
    PrimDiscreteFunctionSpaceType;
   typedef typename GradientDiscreteFunctionType::Traits::DiscreteFunctionSpaceType 
    GradientDiscreteFunctionSpaceType;

   enum{ dimDomain=ModelType::dimDomain};

  typedef typename ConsDiscreteFunctionSpaceType::Traits::GridPartType GridPartType;
  typedef typename GridPartType::GridType GridType;
  typedef typename GridType :: template Codim<0> :: Entity   Entity;
  typedef typename GridType :: template Codim<0> :: Geometry Geometry;
  typedef typename GridPartType::template Codim<0>::IteratorType Iterator;
  
  typedef typename ConsDiscreteFunctionSpaceType::DomainType DomainType;
  typedef typename ConsDiscreteFunctionSpaceType::RangeType ConsRangeType;
  typedef typename PrimDiscreteFunctionSpaceType::RangeType PrimRangeType;
  typedef typename GradientDiscreteFunctionSpaceType::RangeType GradRangeType;  
typedef typename ModelType::JacobianRangeType JacobianRangeType;
  const ConsDiscreteFunctionSpaceType& space =  consDF.space();
  
  primDF.clear();
  
  typedef typename ConsDiscreteFunctionType::LocalFunctionType ConsLocalFuncType;
  typedef typename PrimDiscreteFunctionType::LocalFunctionType PrimLocalFuncType;
  typedef typename GradientDiscreteFunctionType::LocalFunctionType GradLocalFuncType; 
  
  ConsRangeType cons(0.0);
  PrimRangeType prim(0.0);
  GradRangeType grad(0.0);  
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
    GradLocalFuncType gradLF = sigmaDF.localFunction(entity);    
    const int quadNop = quad.nop();
    for(int qP = 0; qP < quadNop; ++qP) 
    {
      const DomainType& xgl = geo.global( quad.point(qP) );
      // evaluate conservative variables
      consLF.evaluate( quad[qP], cons );
	    gradLF.evaluate( quad[qP], grad );
      Fem::FieldMatrixConverter< GradRangeType, JacobianRangeType> jac( grad); 
      
      
      double rho = cons[0];
      double rho_inv = 1. /rho;
      double phi = cons[2];
      phi*=rho_inv;
      double gradphi=jac[2][0];
#if WELLBALANCED

#else      
      gradphi-=phi*jac[0][0];
      gradphi*=rho_inv;
#endif

      model.conservativeToPrimitive( xgl, cons, prim );
      prim[dimDomain]+=model.delta()*0.5*gradphi*gradphi; 
      prim *=  quad.weight(qP);
      primLF.axpy( quad[qP] , prim );
    }
  }
}

#endif

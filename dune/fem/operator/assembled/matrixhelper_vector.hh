#ifndef DUNE_PHASEFIELD_MATRIXHELPER_HH
#define DUNE_PHASEFIELD_MATRIXHELPER_HH
// C++ includes
#include <cassert>
#include <cstddef>

// dune-common includes
#include <dune/common/nullptr.hh>

// dune-geometry includes
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

// dune-fem includes
#include <dune/fem/space/basisfunctionset/functor.hh>
#include <dune/fem/space/basisfunctionset/transformation.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/quadrature/quadrature.hh>
#include <dune/fem/version.hh>

namespace MatrixHelper
{

template< class Point, class ShapeFunctionSet, class RangeArray>
void evaluateScalarAll( const Point &x ,
                        const ShapeFunctionSet& shapeFunctionSet,
                        RangeArray &values )
{
   Dune::Fem::AssignFunctor< RangeArray > f( values );
   shapeFunctionSet.scalarShapeFunctionSet().evaluateEach( x, f );
}

template< class Point, class Geometry , class ShapeFunctionSet , class JacobianRangeArray>
void jacobianScalarAll ( const Point &x, 
                         const Geometry& geo,
                        const ShapeFunctionSet& shapeFunctionSet, 
                        JacobianRangeArray &jacobians ) 
{
  typedef Dune::Fem::JacobianTransformation< Geometry> Transformation;
  Transformation transformation( geo, coordinate( x ) );
  Dune::Fem::AssignFunctor< JacobianRangeArray, Transformation > f( jacobians, transformation );
  shapeFunctionSet.scalarShapeFunctionSet().jacobianEach( x, f );
}


template< int Dimension>
class Alignment
{
  enum{ dimRange=Dimension};

public:
  Alignment( )
  {}


  static size_t vectorialIndex(size_t component, size_t scalarIndex)
  {
    return scalarIndex*dimRange+component;
  }

private:
};


template< class Alignment ,class Couplings, class LocalMatrix,class BfRangeVector, class FluxRangeType>
void axpyIntersection ( Couplings& couplings,
                        const BfRangeVector& phiRow,
                        const BfRangeVector& phiCol,
                        const FluxRangeType& fluxValue,
                        const size_t local_i,
                        const size_t local_j,
                        const double factor,
                        LocalMatrix& localMatrix)                    
{
  for( auto coupling : couplings.intersectionCouplings())
    { 
      auto rowcomp=coupling.first;
      auto colcomp=coupling.second;

      auto i = Alignment::vectorialIndex( rowcomp , local_i );
      auto j = Alignment::vectorialIndex( colcomp , local_j );
      localMatrix.add( i,j,fluxValue[ rowcomp ][ colcomp ]*phiRow[ local_i ][0]*phiCol[ local_j ][0]*factor);
    }
}

template<class Alignment ,class Couplings, class LocalMatrix,class BfRangeVector, class FluxRangeType>
void axpyElement ( Couplings& couplings,
                   const BfRangeVector& phiRow,
                   const FluxRangeType& fluxValue,
                   const size_t local_i,
                   const size_t local_j,
                   const double factor,
                   LocalMatrix& localMatrix)                    
{
  for( auto coupling : couplings.elementCouplings())
    { 
      auto rowcomp=coupling.first;
      auto colcomp=coupling.second;
      
      auto i = Alignment::vectorialIndex( rowcomp , local_i );
      auto j = Alignment::vectorialIndex( colcomp , local_j );
      localMatrix.add( i,j, fluxValue[ rowcomp ][ colcomp ]*phiRow[ local_i ][ 0 ]*factor);
    }
}


template< class Alignment,class JacobianVector,class DiffusionTensor, class LocalMatrix>
void diffusionaxpy( const size_t local_i,
                    const size_t local_j,
                    const size_t dimDomain,
                    const JacobianVector& dphi,
                    const DiffusionTensor du,
                    const double weight,
                    LocalMatrix& jLocal)
{
   //for each velocomponent
  for( size_t ii = 0 ; ii < dimDomain ; ++ii )
    {
      size_t global_i= Alignment::vectorialIndex( 1 + ii , local_i);
      for( size_t jj = 0 ; jj < dimDomain; ++jj )
        {
          size_t global_j = Alignment::vectorialIndex( 1 + jj , local_j);
          double  value(0);
          value=du[ jj ][ ii ]*dphi[ local_i ][ 0 ];

          jLocal.add( global_i , global_j ,0.5*value*weight );
        }
     }
}







template< class Matrix1, class Matrix2>
void add( const Matrix1& m1, Matrix2& m2)
{
  for( int ii = 0 ; ii<m2.rows() ; ++ii)
    for( int jj = 0 ; jj<m2.columns() ; ++jj)
      m2.add( ii , jj, m1.get(ii , jj));
}

template< int dimRange> 
class Couplings
{
#if LAMBDASCHEME
  enum{ couplingSize = 14*dimRange+dimRange*dimRange+10};
#else
  enum{ couplingSize = 10*dimRange+dimRange*dimRange+10};
#endif
  enum{ intersectionCouplingSize = 8*dimRange+dimRange*dimRange+4 };
  
  typedef std::array< std::pair< int , int > ,couplingSize> ElementCouplingType;
  typedef std::array< std::pair< int , int >, intersectionCouplingSize >IntersectionCouplingType;
 


  public:
    Couplings()
    {
      makeElementCouplings();
      makeIntersectionCouplings();
    }
  
    void makeElementCouplings (); 
    
    void makeIntersectionCouplings (); 
  
    ElementCouplingType& elementCouplings() 
      {
        return elementCouplings_;
      }
  
    IntersectionCouplingType& intersectionCouplings() 
      {
        return intersectionCouplings_;
      }
 
  private: 
    ElementCouplingType elementCouplings_;
    IntersectionCouplingType intersectionCouplings_;
};

template<>
void Couplings< 1 >::makeElementCouplings () 
  {   
#if LAMBDASCHEME
#include "./Codegen/elementCouplings1d_lambda_CODEGEN.c"
#else
#include "./Codegen/elementCouplings1d_CODEGEN.c"
#endif
}

template<>
void Couplings< 1 >::makeIntersectionCouplings () 
 { 
#if LAMBDASCHEME
#include "./Codegen/intersectionCouplings1d_lambda_CODEGEN.c"
#else
#include "./Codegen/intersectionCouplings1d_CODEGEN.c"
#endif
}

template<>
void Couplings< 2 >::makeElementCouplings () 
  {   
#if LAMBDASCHEME
#include "./Codegen/elementCouplings2d_lambda_CODEGEN.c"
#else
#include "./Codegen/elementCouplings2d_CODEGEN.c"
#endif
}

template<>
void Couplings< 2 >::makeIntersectionCouplings () 
  {

#if LAMBDASCHEME
#include "./Codegen/intersectionCouplings2d_lambda_CODEGEN.c"
#else
#include "./Codegen/intersectionCouplings2d_CODEGEN.c"
#endif
}
//////////////////////////////////////////////////////////////////////////////////////////////
template< int dimRange>
class AcCouplings
{
  static const int couplingSize = 2*dimRange+dimRange*dimRange+4;
  static const int intersectionCouplingSize = dimRange+dimRange*dimRange+2 ;

  typedef std::vector< std::pair< int , int > > ElementCouplingType;
  typedef std::vector< std::pair< int , int > > IntersectionCouplingType;



  public:
    AcCouplings()
    {
      makeElementCouplings();
      makeIntersectionCouplings();
    }

    void makeElementCouplings ();

    void makeIntersectionCouplings ();

    ElementCouplingType& elementCouplings()
      {
        return elementCouplings_;
      }

    IntersectionCouplingType& intersectionCouplings()
      {
        return intersectionCouplings_;
      }
 
  private:
    ElementCouplingType elementCouplings_;
    IntersectionCouplingType intersectionCouplings_;
};
template<>
void AcCouplings< 1 >::makeElementCouplings ()
  {
    #include "./Codegen/acelCouplings1d_CODEGEN.c"
  }

template<>
void AcCouplings< 1 >::makeIntersectionCouplings ()
  {
    #include "./Codegen/acisCouplings1d_CODEGEN.c"
  }

template<>
void AcCouplings< 2 >::makeElementCouplings ()
  {
    #include "./Codegen/acelCouplings2d_CODEGEN.c"
  }

template<>
void AcCouplings< 2 >::makeIntersectionCouplings ()
  {
    #include "./Codegen/acisCouplings2d_CODEGEN.c"
  }
////////////////////////////////////////////////////////////////////////////////////

template< int dimRange>
class NvStCouplings
{
  static const int couplingSize = 4*dimRange+dimRange*dimRange+3;
  static const int intersectionCouplingSize = 3*dimRange+dimRange*dimRange+2;

  typedef std::vector< std::pair< int , int > > ElementCouplingType;
  typedef std::vector< std::pair< int , int > >IntersectionCouplingType;



  public:
    NvStCouplings()
    {
      makeElementCouplings();
      makeIntersectionCouplings();
    }

    void makeElementCouplings ();

    void makeIntersectionCouplings ();

    ElementCouplingType& elementCouplings()
      {
        return elementCouplings_;
      }

    IntersectionCouplingType& intersectionCouplings()
      {
        return intersectionCouplings_;
      }

  private:
    ElementCouplingType elementCouplings_;
    IntersectionCouplingType intersectionCouplings_;
};
template<>
void NvStCouplings< 1 >::makeElementCouplings ()
  {
    #include "./Codegen/nvstelCouplings1d_CODEGEN.c"
  }

template<>
void NvStCouplings< 1 >::makeIntersectionCouplings ()
  {
    #include "./Codegen/nvstisCouplings1d_CODEGEN.c"
  }

template<>
void NvStCouplings< 2 >::makeElementCouplings ()
  {
    #include "./Codegen/nvstelCouplings2d_CODEGEN.c"
  }

template<>
void NvStCouplings< 2 >::makeIntersectionCouplings ()
  {
    #include "./Codegen/nvstisCouplings2d_CODEGEN.c"
  }
}




#endif

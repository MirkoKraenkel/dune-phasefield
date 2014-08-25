#ifndef MATRIXHELPER_HH
#define MATRIXHELPER_HH
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

  static int vectorialIndex(size_t component, size_t scalarIndex)
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

#if VECTORIAL
     localMatrix.add( i,j,fluxValue[ rowcomp ][ colcomp ]*phiRow[ i ][i rowcomp ]*phiCol[ j ][ colcomp]*factor);
#else
     localMatrix.add( i,j,fluxValue[ rowcomp ][ colcomp ]*phiRow[ local_i ][0]*phiCol[ local_j ][0]*factor);
#endif
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
#if VECTORIAL
      localMatrix.add( i,j, fluxValue[ rowcomp ][ colcomp ]*phiRow[ i ][ rowcomp ]*factor);
#else
      localMatrix.add( i,j, fluxValue[ rowcomp ][ colcomp ]*phiRow[ local_i ][ 0 ]*factor);
#endif
    }
}



template< int dimRange> 
class Couplings
{
  enum{ couplingSize = 10*dimRange+dimRange*dimRange+10};
  enum{ intersectionCouplingSize = 9*dimRange+dimRange*dimRange+3 };
  
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
    elementCouplings_[  0 ] = std::pair< int , int >( 0 , 0 );// rho rho
    elementCouplings_[  1 ] = std::pair< int , int >( 0 , 1 );// rho v
    elementCouplings_[  2 ] = std::pair< int , int >( 1 , 0 );// v rho
    elementCouplings_[  3 ] = std::pair< int , int >( 1 , 1 );// v v
    elementCouplings_[  4 ] = std::pair< int , int >( 1 , 2 );// v phi
    elementCouplings_[  5 ] = std::pair< int , int >( 1 , 3 );// v mu
    elementCouplings_[  6 ] = std::pair< int , int >( 1 , 4 );// v tau
    elementCouplings_[  7 ] = std::pair< int , int >( 2 , 0 );// phi rho
    elementCouplings_[  8 ] = std::pair< int , int >( 2 , 1 );// phi v
    elementCouplings_[  9 ] = std::pair< int , int >( 2 , 2 );// phi phi
    elementCouplings_[ 10 ] = std::pair< int , int >( 2 , 4 );// phi tau
    elementCouplings_[ 11 ] = std::pair< int , int >( 3 , 0 );// mu rho
    elementCouplings_[ 12 ] = std::pair< int , int >( 3 , 1 );// mu v
    elementCouplings_[ 13 ] = std::pair< int , int >( 3 , 2 );// mu phi
    elementCouplings_[ 14 ] = std::pair< int , int >( 3 , 3 );// mu mu
    elementCouplings_[ 15 ] = std::pair< int , int >( 4 , 0 );// tau rho 
    elementCouplings_[ 16 ] = std::pair< int , int >( 4 , 2 );// tau phi
    elementCouplings_[ 17 ] = std::pair< int , int >( 4 , 4 );// tau tau
    elementCouplings_[ 18 ] = std::pair< int , int >( 4 , 5 );// tau sigma 
    elementCouplings_[ 19 ] = std::pair< int , int >( 5 , 2 );// sigma phiv
    elementCouplings_[ 20 ] = std::pair< int , int >( 5 , 5 );// sigma sigma 
  }

template<>
void Couplings< 1 >::makeIntersectionCouplings () 
 { 
   intersectionCouplings_[  0 ] = std::pair< int , int >( 0 , 0 );// rho rho
   intersectionCouplings_[  1 ] = std::pair< int , int >( 0 , 1 );// rho v
   intersectionCouplings_[  2 ] = std::pair< int , int >( 1 , 0 );// v rho
   intersectionCouplings_[  3 ] = std::pair< int , int >( 1 , 1 );// v v
   intersectionCouplings_[  4 ] = std::pair< int , int >( 1 , 2 );// v phi
   intersectionCouplings_[  5 ] = std::pair< int , int >( 1 , 3 );// v mu
   intersectionCouplings_[  6 ] = std::pair< int , int >( 1 , 4 );// v tau
   intersectionCouplings_[  7 ] = std::pair< int , int >( 2 , 1 );// phi v
   intersectionCouplings_[  8 ] = std::pair< int , int >( 2 , 2 );// phi phi
   intersectionCouplings_[  9 ] = std::pair< int , int >( 4 , 4 );// tau tau
   intersectionCouplings_[ 10 ] = std::pair< int , int >( 4 , 5 );// tau sigma 
   intersectionCouplings_[ 11 ] = std::pair< int , int >( 5 , 2 );// sigma phiv
   intersectionCouplings_[ 12 ] = std::pair< int , int >( 5 , 5 );// sigma sigma 
 }

template<>
void Couplings< 2 >::makeElementCouplings () 
  {   
    elementCouplings_[  0 ] = std::pair< int , int >( 0 , 0 );// rho rho
    elementCouplings_[  1 ] = std::pair< int , int >( 0 , 1 );// rho v
    elementCouplings_[  2 ] = std::pair< int , int >( 0 , 2 );// rho u
    elementCouplings_[  3 ] = std::pair< int , int >( 1 , 0 );// v rho
    elementCouplings_[  4 ] = std::pair< int , int >( 1 , 1 );// v v
    elementCouplings_[  5 ] = std::pair< int , int >( 1 , 2 );// v u
    elementCouplings_[  6 ] = std::pair< int , int >( 1 , 3 );// v phi
    elementCouplings_[  7 ] = std::pair< int , int >( 1 , 4 );// v mu
    elementCouplings_[  8 ] = std::pair< int , int >( 1 , 5 );// v tau
    elementCouplings_[  9 ] = std::pair< int , int >( 2 , 0 );// u rho
    elementCouplings_[ 10 ] = std::pair< int , int >( 2 , 1 );// u v
    elementCouplings_[ 11 ] = std::pair< int , int >( 2 , 2 );// u u
    elementCouplings_[ 12 ] = std::pair< int , int >( 2 , 3 );// u phi
    elementCouplings_[ 13 ] = std::pair< int , int >( 2 , 4 );// u mu
    elementCouplings_[ 14 ] = std::pair< int , int >( 2 , 5 );// u tau
    elementCouplings_[ 15 ] = std::pair< int , int >( 3 , 0 );// phi rho 
    elementCouplings_[ 16 ] = std::pair< int , int >( 3 , 1 );// phi v
    elementCouplings_[ 17 ] = std::pair< int , int >( 3 , 2 );// phi u
    elementCouplings_[ 18 ] = std::pair< int , int >( 3 , 3 );// phi phi 
    elementCouplings_[ 19 ] = std::pair< int , int >( 3 , 5 );// phi tau
    elementCouplings_[ 20 ] = std::pair< int , int >( 4 , 0 );// mu rho 
    elementCouplings_[ 21 ] = std::pair< int , int >( 4 , 1 );// mu u
    elementCouplings_[ 22 ] = std::pair< int , int >( 4 , 2 );// mu v
    elementCouplings_[ 23 ] = std::pair< int , int >( 4 , 3 );// mu phi
    elementCouplings_[ 24 ] = std::pair< int , int >( 4 , 4 );// mu mu
    elementCouplings_[ 25 ] = std::pair< int , int >( 5 , 0 );// tau rho 
    elementCouplings_[ 26 ] = std::pair< int , int >( 5 , 3 );// tau phi
    elementCouplings_[ 27 ] = std::pair< int , int >( 5 , 5 );// tau tau
    elementCouplings_[ 28 ] = std::pair< int , int >( 5 , 6 );// tau sigma_x
    elementCouplings_[ 29 ] = std::pair< int , int >( 5 , 7 );// tau sigma_y
    elementCouplings_[ 30 ] = std::pair< int , int >( 6 , 3 );// sigma_x phi 
    elementCouplings_[ 31 ] = std::pair< int , int >( 6 , 6 );// sigma_x sigma_x
    elementCouplings_[ 32 ] = std::pair< int , int >( 7 , 3 );// sigma_y phi 
    elementCouplings_[ 33 ] = std::pair< int , int >( 7 , 7 );// sigma_y sigma_y 
  }

template<>
void Couplings< 2 >::makeIntersectionCouplings () 
  {
    intersectionCouplings_[  0 ] = std::pair< int , int >( 0 , 0 );// rho rho
    intersectionCouplings_[  1 ] = std::pair< int , int >( 0 , 1 );// rho v
    intersectionCouplings_[  2 ] = std::pair< int , int >( 0 , 2 );// rho u
    intersectionCouplings_[  3 ] = std::pair< int , int >( 1 , 0 );// v rho
    intersectionCouplings_[  4 ] = std::pair< int , int >( 1 , 1 );// v v
    intersectionCouplings_[  5 ] = std::pair< int , int >( 1 , 2 );// v u
    intersectionCouplings_[  6 ] = std::pair< int , int >( 1 , 3 );// v phi
    intersectionCouplings_[  7 ] = std::pair< int , int >( 1 , 4 );// v mu
    intersectionCouplings_[  8 ] = std::pair< int , int >( 1 , 5 );// v tau
    intersectionCouplings_[  9 ] = std::pair< int , int >( 2 , 0 );// u rho
    intersectionCouplings_[ 10 ] = std::pair< int , int >( 2 , 1 );// u v
    intersectionCouplings_[ 11 ] = std::pair< int , int >( 2 , 2 );// u u
    intersectionCouplings_[ 12 ] = std::pair< int , int >( 2 , 3 );// u phi
    intersectionCouplings_[ 13 ] = std::pair< int , int >( 2 , 4 );// u mu
    intersectionCouplings_[ 14 ] = std::pair< int , int >( 2 , 5 );// u tau
    intersectionCouplings_[ 15 ] = std::pair< int , int >( 3 , 1 );// phi v
    intersectionCouplings_[ 16 ] = std::pair< int , int >( 3 , 2 );// phi u 
    intersectionCouplings_[ 17 ] = std::pair< int , int >( 3 , 3 );// phi phi
    intersectionCouplings_[ 18 ] = std::pair< int , int >( 5 , 5 );// tau tau
    intersectionCouplings_[ 19 ] = std::pair< int , int >( 5 , 6 );// tau sigma_x
    intersectionCouplings_[ 20 ] = std::pair< int , int >( 5 , 7 );// tau sigma_y
    intersectionCouplings_[ 21 ] = std::pair< int , int >( 6 , 3 );// sigma_x phi
    intersectionCouplings_[ 22 ] = std::pair< int , int >( 7 , 3 );// sigma_y phi 
    intersectionCouplings_[ 23 ] = std::pair< int , int >( 6 , 6 );// sigma_x sigma_x
    intersectionCouplings_[ 24 ] = std::pair< int , int >( 7 , 7 );// sigma_y sigma_y
  }
 



}




#endif

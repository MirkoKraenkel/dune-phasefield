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
  enum{ couplingSize = 7*dimRange+dimRange*dimRange+3};
  enum{ intersectionCouplingSize = 6*dimRange+dimRange*dimRange+2 };
  
  typedef std::array< std::makepair ,couplingSize> ElementCouplingType;
  typedef std::array< std::makepair, intersectionCouplingSize >IntersectionCouplingType;
 


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
    elementCouplings_[  0 ] = std::makepair( 0 , 0 );// rho rho
    elementCouplings_[  1 ] = std::makepair( 0 , 1 );// rho v
    elementCouplings_[  2 ] = std::makepair( 1 , 0 );// v rho
    elementCouplings_[  3 ] = std::makepair( 1 , 1 );// v v
    elementCouplings_[  4 ] = std::makepair( 1 , 3 );// v mu
    elementCouplings_[  5 ] = std::makepair( 2 , 0 );// mu rho
    elementCouplings_[  6 ] = std::makepair( 2 , 1 );// mu v
    elementCouplings_[  7 ] = std::makepair( 2 , 2 );// mu mu
    elementCouplings_[  8 ] = std::makepair( 2 , 3 );// mu sigma 
    elementCouplings_[  9 ] = std::makepair( 3 , 0 );// sigma rho 
h   elementCouplings_[ 10 ] = std::makepair( 3 , 3 );// sigma sigma 
  }

template<>
void Couplings< 1 >::makeIntersectionCouplings () 
 { 
   intersectionCouplings_[  0 ] = std::makepair( 0 , 0 );// rho rho
   intersectionCouplings_[  1 ] = std::makepair( 0 , 1 );// rho v
   intersectionCouplings_[  2 ] = std::makepair( 1 , 0 );// v rho
   intersectionCouplings_[  3 ] = std::makepair( 1 , 1 );// v v
   intersectionCouplings_[  4 ] = std::makepair( 1 , 2 );// v mu
   intersectionCouplings_[  5 ] = std::makepair( 2 , 3 );// mu sigma 
   intersectionCouplings_[  6 ] = std::makepair( 3 , 0 );// sigma rhp
   intersectionCouplings_[  7 ] = std::makepair( 3 , 3 );// sigma sigma 
 }

template<>
void Couplings< 2 >::makeElementCouplings () 
  {   
    elementCouplings_[  0 ] = std::makepair( 0 , 0 );// rho rho
    elementCouplings_[  1 ] = std::makepair( 0 , 1 );// rho v
    elementCouplings_[  2 ] = std::makepair( 0 , 2 );// rho u
    elementCouplings_[  3 ] = std::makepair( 1 , 0 );// v rho
    elementCouplings_[  4 ] = std::makepair( 1 , 1 );// v v
    elementCouplings_[  5 ] = std::makepair( 1 , 2 );// v u
    elementCouplings_[  6 ] = std::makepair( 1 , 3 );// v mu
    elementCouplings_[  7 ] = std::makepair( 2 , 0 );// u rho
    elementCouplings_[  8 ] = std::makepair( 2 , 1 );// u v
    elementCouplings_[  9 ] = std::makepair( 2 , 2 );// u u
    elementCouplings_[ 10 ] = std::makepair( 2 , 3 );// u mu
    elementCouplings_[ 11 ] = std::makepair( 3 , 0 );// mu rho 
    elementCouplings_[ 12 ] = std::makepair( 3 , 1 );// mu u
    elementCouplings_[ 13 ] = std::makepair( 3 , 2 );// mu v
    elementCouplings_[ 14 ] = std::makepair( 3 , 3 );// mu mu
    elementCouplings_[ 15 ] = std::makepair( 3 , 4 );// mu sigma_x
    elementCouplings_[ 16 ] = std::makepair( 3 , 5 );// mu sigma_y
    elementCouplings_[ 17 ] = std::makepair( 4 , 0 );// sigma_x rho 
    elementCouplings_[ 18 ] = std::makepair( 4 , 4 );// sigma_x sigma_x
    elementCouplings_[ 19 ] = std::makepair( 5 , 0 );// sigma_y rho 
    elementCouplings_[ 20 ] = std::makepair( 5 , 5 );// sigma_y sigma_y 
  }

template<>
void Couplings< 2 >::makeIntersectionCouplings () 
  {
    intersectionCouplings_[  0 ] = std::makepair( 0 , 0 );// rho rho
    intersectionCouplings_[  1 ] = std::makepair( 0 , 1 );// rho v
    intersectionCouplings_[  2 ] = std::makepair( 0 , 2 );// rho u
    intersectionCouplings_[  3 ] = std::makepair( 1 , 0 );// v rho
    intersectionCouplings_[  4 ] = std::makepair( 1 , 1 );// v v
    intersectionCouplings_[  5 ] = std::makepair( 1 , 2 );// v u
    intersectionCouplings_[  6 ] = std::makepair( 1 , 3 );// v mu
    intersectionCouplings_[  7 ] = std::makepair( 2 , 0 );// u rho
    intersectionCouplings_[  8 ] = std::makepair( 2 , 1 );// u v
    intersectionCouplings_[  9 ] = std::makepair( 2 , 2 );// u u
    intersectionCouplings_[ 10 ] = std::makepair( 2 , 3 );// u mu
    intersectionCouplings_[ 11 ] = std::makepair( 3 , 5 );// mu mu
    intersectionCouplings_[ 12 ] = std::makepair( 3 , 4 );// mu sigma_x
    intersectionCouplings_[ 13 ] = std::makepair( 3 , 5 );// mu sigma_y
    intersectionCouplings_[ 14 ] = std::makepair( 4 , 0 );// sigma_x rho
    intersectionCouplings_[ 15 ] = std::makepair( 5 , 0 );// sigma_y rho 
    intersectionCouplings_[ 16 ] = std::makepair( 4 , 4 );// sigma_x sigma_x
    intersectionCouplings_[ 17 ] = std::makepair( 5 , 5 );// sigma_y sigma_y
  }
 



}




#endif

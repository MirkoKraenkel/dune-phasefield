#ifndef COMPONENT_L2_NORM_HH
#define COMPONENT_L2_NORM_HH

#include <dune/fem/misc/l2norm.hh>

namespace Dune
{

  namespace Fem 
  {
      
      template< class GridPart >
      class ComponentL2Norm
      : public LPNormBase< GridPart,ComponentL2Norm< GridPart > >
      {
        typedef LPNormBase< GridPart,ComponentL2Norm< GridPart > > BaseType;
        typedef ComponentL2Norm< GridPart > ThisType;
      public:
        typedef GridPart GridPartType;

        using BaseType::gridPart;
        using BaseType::comm;

        template< class Function >
        struct FunctionComponentSquare;
        template< class UFunction, class VFunction >
        struct FunctionDistance;

      protected:
        typedef typename GridPartType::template Codim< 0 >::IteratorType GridIteratorType;
        typedef typename GridIteratorType::Entity EntityType;
        typedef CachingQuadrature< GridPartType, 0 > QuadratureType;
        typedef Integrator< QuadratureType > IntegratorType;

        const unsigned int order_;

      public:
        explicit ComponentL2Norm( const GridPartType &gridPart, 
                                  std::vector<unsigned int> components,
                                  unsigned int order=0);

        template< class DiscreteFunctionType >
        typename DiscreteFunctionType::RangeFieldType
        norm ( const DiscreteFunctionType &u ) const;
      
        template< class UDiscreteFunctionType, class VDiscreteFunctionType >
        typename UDiscreteFunctionType::RangeFieldType
        distance ( const UDiscreteFunctionType &u, const VDiscreteFunctionType &v ) const;


        template< class DiscreteFunctionType, class ReturnType >
        void normLocal ( const EntityType &entity, const unsigned int order,
              const DiscreteFunctionType &u, 
              ReturnType& sum ) const;
      
        template< class UDiscreteFunctionType, class VDiscreteFunctionType, class ReturnType >
        void distanceLocal ( const EntityType &entity, const unsigned int order,  
                         const UDiscreteFunctionType &u,
                         const VDiscreteFunctionType &v,
                         ReturnType& sum ) const;
        protected:
        const std::vector< unsigned int  > components_;
 
      };

   
    template< class GridPart >
    inline ComponentL2Norm< GridPart >
    ::ComponentL2Norm( const GridPartType &gridPart,
                       const std::vector<unsigned int> components,
                       const unsigned int order) 
    : BaseType( gridPart),
      order_( order ),
      components_( components )
      {
      }

    template< class GridPart >
    template< class DiscreteFunctionType >
    inline typename DiscreteFunctionType::RangeFieldType
    ComponentL2Norm< GridPart >::norm ( const DiscreteFunctionType &u ) const
    {
      typedef typename DiscreteFunctionType::RangeFieldType RangeFieldType;
      typedef FieldVector< RangeFieldType, 1 > ReturnType ;

      // calculate integral over each element
      ReturnType sum = BaseType :: forEach( u, ReturnType(0), order_ );

      // return result, e.g. sqrt of calculated sum
      return sqrt( comm().sum( sum[ 0 ] ) );
    }


    template< class GridPart >
    template< class UDiscreteFunctionType, class VDiscreteFunctionType >
    inline typename UDiscreteFunctionType::RangeFieldType
    ComponentL2Norm< GridPart >
      ::distance ( const UDiscreteFunctionType &u, const VDiscreteFunctionType &v ) const
    {
      typedef typename UDiscreteFunctionType::RangeFieldType RangeFieldType;
      typedef FieldVector< RangeFieldType, 1 > ReturnType ;

      // calculate integral over each element
      ReturnType sum = BaseType :: forEach( u, v, ReturnType(0), order_ );

      // return result, e.g. sqrt of calculated sum
      return sqrt( comm().sum( sum[ 0 ] ) );
    }


    template< class GridPart >
    template< class DiscreteFunctionType, class ReturnType >
    inline void 
    ComponentL2Norm< GridPart >
    ::normLocal ( const EntityType& entity, const unsigned int order, 
                                    const DiscreteFunctionType &u,
                                    ReturnType& sum ) const
    {
      typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;

      // evaluate norm locally 
      IntegratorType integrator( order );

      LocalFunctionType ulocal = u.localFunction( entity );
      FunctionComponentSquare< LocalFunctionType > ulocal2( ulocal,components_ );

      integrator.integrateAdd( entity, ulocal2, sum );
    }

    template< class GridPart >
    template< class UDiscreteFunctionType, 
              class VDiscreteFunctionType,
              class ReturnType >
    inline void 
    ComponentL2Norm< GridPart >::distanceLocal ( const EntityType& entity, const unsigned int order, 
                                        const UDiscreteFunctionType &u,
                                        const VDiscreteFunctionType &v,
                                        ReturnType& sum ) const
    {
      typedef typename UDiscreteFunctionType::LocalFunctionType ULocalFunctionType;
      typedef typename VDiscreteFunctionType::LocalFunctionType VLocalFunctionType;

      // evaluate norm locally 
      IntegratorType integrator( order );

      ULocalFunctionType ulocal = u.localFunction( entity );
      VLocalFunctionType vlocal = v.localFunction( entity );

      typedef FunctionDistance< ULocalFunctionType, VLocalFunctionType > LocalDistanceType;

      LocalDistanceType dist( ulocal, vlocal );
      FunctionComponentSquare< LocalDistanceType > dist2( dist,components_ );
       
      integrator.integrateAdd( entity, dist2, sum );
    }

    


    template< class GridPart >
    template< class Function >
    struct ComponentL2Norm< GridPart >::FunctionComponentSquare
    {
      typedef Function FunctionType;

      typedef typename FunctionType::RangeFieldType RangeFieldType;
      typedef FieldVector< RangeFieldType, 1 > RangeType;

      explicit FunctionComponentSquare ( const FunctionType &function ,const std::vector<unsigned int >& components)
      : function_( function ),
        components_( components) 
      {}
      
      template< class Point >
      void evaluate ( const Point &x, RangeType &ret ) const
      {
        typename FunctionType::RangeType phi;
        function_.evaluate( x, phi );
        //for( auto i : components_)
          ret =phi[0]*phi[0];
      }
      
    private:
      const FunctionType &function_;
      const std::vector<unsigned int> components_;
    };

    template< class GridPart >
    template< class UFunction, class VFunction >
    struct ComponentL2Norm< GridPart >::FunctionDistance
    {
      typedef UFunction UFunctionType;
      typedef VFunction VFunctionType;

      typedef typename UFunctionType::RangeFieldType RangeFieldType;
      typedef typename UFunctionType::RangeType RangeType;
      typedef typename UFunctionType::JacobianRangeType JacobianRangeType;

      FunctionDistance ( const UFunctionType &u, const VFunctionType &v )
      : u_( u ), v_( v )
      {}
 
      template< class Point >
      void evaluate ( const Point &x, RangeType &ret ) const
      {
        RangeType phi;
        u_.evaluate( x, ret );
        v_.evaluate( x, phi );
        ret -= phi;
      }

      template< class Point >
      void jacobian ( const Point &x, JacobianRangeType &ret ) const
      {
        JacobianRangeType phi;
        u_.jacobian( x, ret );
        v_.jacobian( x, phi );
        ret -= phi;
      }

    private:
      const UFunctionType &u_;
      const VFunctionType &v_;
    };






 }//end namespace Fem

}//end namespace

#endif//OMPONENT_L2_NORM_HH

#ifndef COMPONENT_L2_NORM_HH
#define COMPONENT_L2_NORM_HH

#include <dune/fem/misc/l2norm.hh>

namespace Dune
{

  namespace Fem 
  {
      
      template< class GridPart >
      class ComponentL2Norm
      : public L2Norm< GridPart>
      {
        typedef L2Norm< GridPart > BaseType;
        typedef ComponentL2Norm< GridPart > ThisType;
        typedef GridPart GridPartType;
      
      protected:
        typedef typename BaseType::GridIteratorType GridIteratorType;
        typedef typename BaseType::IntegratorType IntegratorType;
        typedef typename GridIteratorType::Entity EntityType;
        
        using BaseType::gridPart;
        using BaseType::comm;
        //using BaseType::FunctionDistance;
        template< class Function >
        struct FunctionComponentSquare;

      public:
        using BaseType::norm;
        using BaseType::distance;
        explicit ComponentL2Norm( const GridPartType &gridPart, 
                                  std::vector<unsigned int> components,
                                  unsigned int order=0);
      
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
   : BaseType( gridPart, order),
      components_(components)
      {
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

      typedef typename BaseType::template FunctionDistance< ULocalFunctionType, VLocalFunctionType > LocalDistanceType;

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
        for( auto i : components_)
          ret += phi[i]*phi[i];
      }
      
    private:
      const FunctionType &function_;
      const std::vector<unsigned int> components_;
    };




 }//end namespace Fem

}//end namespace

#endif//OMPONENT_L2_NORM_HH

#ifndef GRADERROR_ESTIMATOR_HH 
#define GRADERROR_ESTIMATOR_HH 

//- Dune includes
// #include <dune/grid/common/referenceelements.hh>

#include <cmath>

//- Dune-fem includes 
#include <dune/fem/quadrature/caching/twistutility.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/common/spaceoperatorif.hh> 
#include <dune/fem/operator/matrix/blockmatrix.hh>
#include <dune/fem/space/discontinuousgalerkin.hh>
#include <dune/fem/space/combinedspace.hh>
#include <dune/fem/quadrature/intersectionquadrature.hh>
#include <dune/fem/misc/h1norm.hh>

namespace Dune
{
               
  template<class UFunction>
  class Estimator1
  {
   typedef Estimator1< UFunction> ThisType;

  public:
    typedef UFunction DiscreteFunctionType;

    typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionType :: LocalFunctionType LocalFunctionType;
  
    typedef typename DiscreteFunctionSpaceType :: DomainFieldType DomainFieldType; 
    typedef typename DiscreteFunctionSpaceType :: RangeFieldType RangeFieldType; 
    typedef typename DiscreteFunctionSpaceType :: DomainType DomainType; 
    typedef typename DiscreteFunctionSpaceType :: RangeType RangeType; 
    typedef typename DiscreteFunctionSpaceType :: JacobianRangeType JacobianRangeType; 
    typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType; 
    typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;

    typedef typename GridPartType :: GridType GridType;
    typedef typename GridPartType :: IndexSetType IndexSetType;
    typedef typename GridPartType :: IntersectionIteratorType IntersectionIteratorType;

    typedef typename IntersectionIteratorType :: Intersection IntersectionType;

    typedef typename GridType :: template Codim< 0 > :: Entity ElementType;
    typedef typename GridType :: template Codim< 0 > :: EntityPointer ElementPointerType;
    typedef typename ElementType::Geometry GeometryType;

    static const int dimension = GridType :: dimension;

    typedef Fem::CachingQuadrature< GridPartType, 0 > ElementQuadratureType;
    typedef Fem::CachingQuadrature< GridPartType, 1 > FaceQuadratureType;
    typedef std :: vector< double > ErrorIndicatorType;

  protected:
    const DiscreteFunctionType &uh_;
    const double &beta_;
    const DiscreteFunctionSpaceType &dfSpace_;
    const GridPartType &gridPart_;
    const IndexSetType &indexSet_;
    GridType &grid_;
    ErrorIndicatorType indicator_;
    double totalIndicator2_,maxIndicator_;
    const double theta_;
    int maxLevel_;
    const double coarsen_;
      // const ProblemType &problem_;

  public:
    explicit Estimator1 (const DiscreteFunctionType &uh,
												 GridType &grid)
      :  uh_( uh ),
         beta_(1.),
				 dfSpace_( uh.space() ),
				 gridPart_( dfSpace_.gridPart() ),
				 indexSet_( gridPart_.indexSet() ),
				 grid_( grid ),
				 indicator_( indexSet_.size( 0 )),
				 totalIndicator2_(0),
				 maxIndicator_(0),
				 theta_( Dune::Fem::Parameter::getValue("phasefield.adaptive.theta",0.) ),
				 maxLevel_(Dune::Fem::Parameter::getValue<int>("fem.adaptation.finestLevel")),
				 coarsen_(Dune::Fem::Parameter::getValue<double>("fem.adaptation.coarsenPercent",1e-50))
    {
			clear();
    }
    
    // make this class behave as required for a LocalFunctionAdapter
    typedef Fem::FunctionSpace< double, double, dimension, 1 > FunctionSpaceType;
    void init(const ElementType &en)
    {
      enIndex_ = indexSet_.index(en);
      entity_ = &en;
    }
    template< class PointType >
    void evaluate(const PointType& x,FieldVector<double,1> &ret)
    {
      ret[0] = indicator_[enIndex_];
    }
  private:
    const ElementType *entity_;
    int enIndex_;

  public:
    void clear ()
    {
      indicator_.resize( indexSet_.size( 0 ));
      std::fill( indicator_.begin(), indicator_.end(), 0);
    }
    
    double estimate ( )
    {
      clear();
      const IteratorType end = dfSpace_.end();
      for( IteratorType it = dfSpace_.begin(); it != end; ++it )
      {
        const ElementType &entity = *it;
        const LocalFunctionType uLocal = uh_.localFunction( entity );
  
        estimateLocal(entity, uLocal );
#if 0 
        IntersectionIteratorType end = gridPart_.iend( entity );
         for( IntersectionIteratorType inter = gridPart_.ibegin( entity ); inter != end; ++inter )
         {
           const IntersectionType &intersection = *inter;
           if( intersection.neighbor() )
             estimateIntersection( intersection, entity, uLocal );
           else
             estimateBoundary( intersection,entity,uLocal);
         }
#endif
 
      }
      return computeIndicator();
    }
    
    double computeIndicator()
    {
      totalIndicator2_ = 0.0;
  
      for (unsigned int i=0;i<indicator_.size();++i)
      {
        totalIndicator2_ += indicator_[ i ];
        maxIndicator_ = std::max(maxIndicator_,indicator_[i]);
      }
	
    	return sqrt( totalIndicator2_ );
    }


    //! mark all elements due to given tolerance 
    bool mark ( const double tolerance ) const
    {
      int marked = 0;
      if (tolerance < 0)
      {
        const IteratorType end = dfSpace_.end();
        for( IteratorType it = dfSpace_.begin(); it != end; ++it )
        {
          const ElementType &entity = *it;
           grid_.mark( 1, entity );
          ++marked;
  	    }
      }
      else
      {
      	const double localTol2=tolerance;
   
        // loop over all elements 
        const IteratorType end = dfSpace_.end();
        for( IteratorType it = dfSpace_.begin(); it != end; ++it )
        {
          const ElementType &entity = *it;
          // check local error indicator 
          if( (indicator_[ indexSet_.index( entity ) ] > localTol2 && theta_ == 0) ||
              (indicator_[ indexSet_.index( entity ) ] == -1 && theta_ > 0) )
	          { 
	            
              if(entity.level()<maxLevel_)
              {
		            grid_.mark( 1, entity );
		            
                IntersectionIteratorType end = gridPart_.iend( entity );
		  
		            for( IntersectionIteratorType inter = gridPart_.ibegin( entity ); inter != end; ++inter )
                {
                  const IntersectionType &intersection = *inter;
           
                  if( intersection.neighbor() )
                  {
			              const ElementPointerType pOutside = intersection.outside();
			              
                    const ElementType &outside = *pOutside;  
			              
                    if(outside.level()<maxLevel_)
			              {
			              
                      grid_.mark( 1, outside );
			             
                      ++marked;
			              
                    }
			            }
		            }
                
                ++marked;
		          
              }
	          }
  	        else if(indicator_[ indexSet_.index( entity ) ] < coarsen_ )
            {  
	            
              if(entity.level()>0)
		           grid_.mark(-1,entity);
	          
            }
	          else
            {
	            grid_.mark(0,entity);
	          }
        }
      }
      
      return (marked > 0);
    
    }
    
    bool estimateAndMark(double tolerance)
    {
      double esti=estimate();
      //  std::cout<<"EstimateValue="<<esti<<std::endl;
      return mark(tolerance);
    }
    

protected:
    //! caclulate error on element 
    void estimateLocal ( const ElementType &entity,
                         const LocalFunctionType &uLocal )
    {
      const typename ElementType :: Geometry &geometry = entity.geometry();

      const double volume = geometry.volume();
      double h2 = (dimension == 2 ? volume : std :: pow( volume, 2.0 / (double)dimension ));  
      const int index = indexSet_.index( entity );
     
      ElementQuadratureType quad( entity, 2*(dfSpace_.order() )+1 );
      const int numQuadraturePoints = quad.nop();
      for( int qp = 0; qp < numQuadraturePoints; ++qp )
      {
	      
        JacobianRangeType gradient;

        uLocal.jacobian(quad[qp],gradient);
	
        RangeType y(0.),tmp(0.);

	      const DomainType global = geometry.global(quad.point(qp));
	
	      const double weight = quad.weight( qp ) * geometry.integrationElement( quad.point( qp ) );
	
  //    	indicator_[ index ] += h2*weight * (gradient[dimension+1] * gradient[dimension+1]);
      }

    }

    
    
    void estimateBoundary ( const IntersectionType &intersection,
			                      const ElementType &inside,
			                      const LocalFunctionType &uInside)
    {

      // WARNING: correct boundary treatment still missing!
      return;

      const int quadOrder = 2*(dfSpace_.order()+1);
      const int insideIndex = indexSet_.index( inside );    
      
      double error;
      error=0;

      const  FaceQuadratureType quadInside( gridPart_, intersection, quadOrder, FaceQuadratureType::INSIDE );
      const typename ElementType :: Geometry &geometry = inside.geometry();

      const int numQuadraturePoints = quadInside.nop();
      for( int qp = 0; qp < numQuadraturePoints; ++qp )
      {
        RangeType dirichlet(0.0); 
        DomainType global(0.);
        DomainType integrationNormal
           = intersection.integrationOuterNormal( quadInside.localPoint( qp ) );
        const double integrationElement = integrationNormal.two_norm();
         
        global=geometry.global(quadInside.point(qp));


        RangeType uvalInside;
      
        uInside.evaluate(quadInside[qp],uvalInside);
         
        uvalInside-=dirichlet;
         
        JacobianRangeType jumpNormal;

        error+=quadInside.weight(qp)*(uvalInside*uvalInside)*integrationElement;

        if( error > 0.0 )
        {
          const double volume = (inside.geometry().volume() );
          const double h = (dimension == 1 ? volume : std::pow( volume, 1.0 / (double)dimension ));;
             
          error/=h;
             
          indicator_[ insideIndex ] +=  error;
           
        }
      }
    }
    //! caclulate error on boundary intersections 
    void estimateIntersection ( const IntersectionType &intersection,
                                const ElementType &inside,
                                const LocalFunctionType &uInside )
    {
      const ElementPointerType pOutside = intersection.outside();
      const ElementType &outside = *pOutside;

      const int quadOrder = 2*(dfSpace_.order()+1);

      const int insideIndex = indexSet_.index( inside );
      const int outsideIndex = indexSet_.index( outside );

      const bool isOutsideInterior = (outside.partitionType() == InteriorEntity);
      if( !isOutsideInterior || (insideIndex < outsideIndex) )
      {
        const LocalFunctionType uOutside = uh_.localFunction( outside );
      
      	// const double volume = std::max( inside.geometry().volume() , outside.geometry().volume() );
        // const double h = 2.*volume / intersection.geometry().volume();
      	const double volume = ( inside.geometry().volume() + outside.geometry().volume() );
 	      const double h = (dimension == 1 ? volume : std::pow(0.5* volume, 1.0 / (double)dimension ));

        double errorInside, errorOutside;
        if( ! intersection.conforming() )
        {
          estimateIntersection< false > ( intersection, inside, outside, quadOrder, insideIndex, outsideIndex, 
                                          uInside, uOutside,
                                          h,volume, 
                                          errorInside, errorOutside);
        }
        else 
        {
          estimateIntersection< true > ( intersection, inside, outside, quadOrder, insideIndex, outsideIndex, 
                                         uInside, uOutside,                                       h,volume,
                                         errorInside, errorOutside);
        }

        if( errorInside > 0.0 )
	  {
	    indicator_[ insideIndex ] +=  errorInside;
	   
	    if( isOutsideInterior )
	      indicator_[ outsideIndex ] +=  errorOutside;
	  }
      }
    }

    //! caclulate error on element intersections 
    template< bool conforming >
    void estimateIntersection ( const IntersectionType &intersection,
                                const ElementType &inside, const ElementType &outside,
                                const int quadOrder, 
                                const int insideIndex,
                                const int outsideIndex,
                                const LocalFunctionType &uInside,
                                const LocalFunctionType &uOutside, 
				const double h, double volume,
                                double &errorInside, double &errorOutside)
    {
      // use IntersectionQuadrature to create appropriate face quadratures 
      typedef Fem::IntersectionQuadrature< FaceQuadratureType, conforming > IntersectionQuadratureType;
      typedef typename IntersectionQuadratureType :: FaceQuadratureType QuadratureImp;

      // create intersection quadrature 
      IntersectionQuadratureType interQuad( gridPart_, intersection, quadOrder );

      // get appropriate references 
      const QuadratureImp &quadInside  = interQuad.inside();
      const QuadratureImp &quadOutside = interQuad.outside();
      const int numQuadraturePoints = quadInside.nop();

      // obtain all required function values on intersection
      std::vector< RangeType > uValuesEn( numQuadraturePoints );
      std::vector< RangeType > uValuesNb( numQuadraturePoints );
      uInside.evaluateQuadrature( quadInside, uValuesEn );
      uOutside.evaluateQuadrature( quadOutside, uValuesNb );
     
      
     
      errorInside = 0.0;
      errorOutside = 0.0;
      double faceVol= 0.; 
      faceVol = intersection.geometry().volume(); 
       
      volume/=3*faceVol;

      for( int qp = 0; qp < numQuadraturePoints; ++qp )
      {
	      DomainType unitNormal = intersection.integrationOuterNormal( quadInside.localPoint( qp ) );
	      const double integrationElement = unitNormal.two_norm();
	  
	      unitNormal/=integrationElement;
	  
	  
    	  RangeType jump;
	      jump=uValuesEn[qp];
	      jump-=uValuesNb[qp];
	  
    	  errorInside  += quadInside.weight( qp ) *1./h* (jump * jump) *integrationElement;
	      errorOutside += quadOutside.weight( qp ) *1./h* (jump * jump) *integrationElement;
	      
	    }
    }
     
  
  };
}
#endif

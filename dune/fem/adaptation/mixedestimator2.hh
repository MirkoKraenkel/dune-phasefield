#ifndef MIXEDESTIMATOR_HH 
#define MIXEDESTIMATOR_HH

//- Dune includes
// #include <dune/grid/common/referenceelements.hh>

#include <cmath>

//- Dune-fem includes 
#include <dune/fem/quadrature/caching/twistutility.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/common/spaceoperatorif.hh> 
#include <dune/fem/operator/matrix/blockmatrix.hh>
#include <dune/fem/space/discontinuousgalerkin.hh>
#include <dune/fem/quadrature/intersectionquadrature.hh>
#include <dune/fem/misc/h1norm.hh>

#include <dune/fem/operator/assembled/phasefieldfilter.hh>

namespace Dune
{
                 
  template<class UFunction, class Model>
  class MixedEstimator
  {
    typedef MixedEstimator< UFunction, Model> ThisType;

  public:
    typedef UFunction DiscreteFunctionType;
    typedef Model ModelType;
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
    mutable ErrorIndicatorType refined_;
   
    const  ModelType& model_;
    double totalIndicator2_,maxIndicator_;
    const double theta_;
    int maxLevel_;
    int minLevel_;
    const double coarsen_;
    const double tolfactor_;
      // const ProblemType &problem_;

  public:
    explicit MixedEstimator (const DiscreteFunctionType &uh , GridType &grid,const  ModelType& model):
      uh_( uh ),
      beta_(1.),
      dfSpace_( uh.space() ),
      gridPart_( dfSpace_.gridPart() ),
	    indexSet_( gridPart_.indexSet() ),
	    grid_( grid ),
	    indicator_( indexSet_.size( 0 )),
	    refined_( indexSet_.size( 0 )), 
      model_(model),
      totalIndicator2_(0),
	    maxIndicator_(0),
	    theta_( Dune::Fem::Parameter::getValue("phasefield.adaptive.theta",0.) ),
	    maxLevel_(Dune::Fem::Parameter::getValue<int>("fem.adaptation.finestLevel")), 
	    minLevel_(Dune::Fem::Parameter::getValue<int>("fem.adaptation.coarsestLevel")),
	    coarsen_(Dune::Fem::Parameter::getValue<double>("fem.adaptation.coarsenPercent",0.1)),
      tolfactor_( Dune::Fem::Parameter::getValue<double>("phasefield.adaptive.tolfactor",1) ) 
      {
        clear();
      }  
    
    // make this class behave as required for a LocalFunctionAdapter
    typedef Fem::FunctionSpace< double, double, dimension, 2 > FunctionSpaceType;
    void init(const ElementType &en)
    {
      enIndex_ = indexSet_.index(en);
      entity_ = &en;
    }
  
    template< class PointType >
    void evaluate(const PointType& x,FieldVector<double,2> &ret)
    {
    
      ret[0] = indicator_[enIndex_];
      ret[1] = refined_[enIndex_];
    }
  private:
    const ElementType *entity_;
    int enIndex_;

  public:
    void clear ()
    {
      indicator_.resize( indexSet_.size( 0 ));
      refined_.resize( indexSet_.size( 0 ));
      std::fill( indicator_.begin(), indicator_.end(), 0.);
      std::fill( refined_.begin(), refined_.end(), std::numeric_limits<int>::min());
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

        IntersectionIteratorType end = gridPart_.iend( entity );
#if 1        
        for( IntersectionIteratorType inter = gridPart_.ibegin( entity ); inter != end; ++inter )
         {
        
           const IntersectionType &intersection = *inter;
            
           if( intersection.neighbor() )
             estimateIntersection( intersection, entity, uLocal );
     //      else
  //  ..       estimateBoundary( intersection,entity,uLocal);

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
      	// loop over all elements
        const IteratorType end = dfSpace_.end();
        for( IteratorType it = dfSpace_.begin(); it != end; ++it )
        {
      
          const ElementType &entity = *it;
         int index=indexSet_.index(entity);
         if( std::abs(indicator_[indexSet_.index(entity)]) > tolerance)
	       {
	          if(entity.level()<maxLevel_)
		        {
 
              grid_.mark( 1, entity );
              refined_[ index ] = 1;	
  #if 1	         
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
			              refined_[indexSet_.index(outside)]=1;
                    ++marked;
			            
                  }
              }
		          }
#endif
            ++marked;
		       }
	        }
          else if( indicator_[indexSet_.index(entity)] < coarsen_*tolerance )
	        {
            if(entity.level()>minLevel_ /*&& !(refined_[indexSet_.index(entity)]>std::numeric_limits<int>::min())*/ )
            { 
              
		          grid_.mark(-1,entity);
            }
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
      //std::cout<<"MaxIndicator="<<(*std::max_element(indicator_.begin(),indicator_.end()));
      //std::cout<<"\nMin h="<<(*std::min_element(refined_.begin(),refined_.end()));
      return mark(tolerance);
    }
    

  protected:
    //! caclulate error on element 
    void estimateLocal ( const ElementType &entity,
                         const LocalFunctionType &uLocal )
    {
      const typename ElementType :: Geometry &geometry = entity.geometry();
       const Dune::ReferenceElement< double, dimension > &refElement =
       Dune::ReferenceElements< double, dimension >::general( entity.type() );
 
      const double volume = geometry.volume();
     double h2 = (dimension == 2 ? volume : std :: pow( volume, 2.0 / (double)dimension ));  
      //double h2=std::sqrt(volume);
      const int index = indexSet_.index( entity );
     
      ElementQuadratureType quad( entity, 2*(dfSpace_.order() )+1 );
      const int numQuadraturePoints = quad.nop();
      double sigmasquared=0.;
      double maxsigma=0; 
      RangeType range;
      JacobianRangeType gradient;
  //    uLocal.evaluate( refElement.position(0 , 0),range );
  //    for( int i=0 ; i < dimension ; ++ i)
      
#if 1 
  for( int qp = 0; qp < numQuadraturePoints; ++qp )
        {
          JacobianRangeType gradient;
	        RangeType range;
          double weight = quad.weight(qp) * geometry.integrationElement( quad.point( qp )) ;
	  
	        uLocal.evaluate(quad[qp],range);
          
          for(int i=0 ; i<dimension; ++i)
            sigmasquared+=PhasefieldFilter<RangeType>::sigma(range,i)*PhasefieldFilter<RangeType>::sigma(range,i);

         sigmasquared*=weight;
        // maxsigma=std::max(sigmasquared,maxsigma);

        }
#endif
      //L2-Norm Sigma
      //double normsigma=std::sqrt(sigmasquared);
      //double normsigma=std::sqrt(maxsigma);
      //normsigma*=h2;

      indicator_[ index ]=sigmasquared;
      refined_[ index ] = h2;
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
                                         uInside, uOutside,
                                         h,volume,
                                         errorInside, errorOutside);
        }

        if( errorInside < 0.0 )
	      {
	       // refined_[ insideIndex ] +=  errorInside;
            
	       // if( isOutsideInterior )
	      //   refined_[ outsideIndex ] +=  errorOutside;
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
	  
        double sigmaNormalEn=0;
        double sigmaNormalNb=0;
        for(int i=0; i<dimension ; ++i)
         { 
           sigmaNormalEn+=PhasefieldFilter<RangeType>::sigma(uValuesEn[qp],i)*unitNormal[i];
           sigmaNormalNb-=PhasefieldFilter<RangeType>::sigma(uValuesNb[qp],i)*unitNormal[i];
   
         }

          
        errorInside += quadInside.weight( qp )*sigmaNormalEn*PhasefieldFilter<RangeType>::phi( jump );
	      errorOutside += quadOutside.weight( qp )*sigmaNormalNb*PhasefieldFilter<RangeType>::phi( jump );
 
	      
	    }
    } 
     
  
  };
}
#endif

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
#include <dune/fem/space/combinedspace.hh>
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
    mutable ErrorIndicatorType refinedandcoarsened_; 
    const  ModelType& model_;
    double totalIndicator2_,maxIndicator_;
    const double theta_;
    int maxLevel_;
    int minLevel_;
    const double coarsen_;
    const double tolfactor_;
    const double maxH_;
    const bool verbose_;
      // const ProblemType &problem_;

  public:
    explicit MixedEstimator ( const DiscreteFunctionType &uh ,
                              GridType &grid,
                              const  ModelType& model):
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
      tolfactor_( Dune::Fem::Parameter::getValue<double>("phasefield.adaptive.tolfactor",1) ),
      maxH_( Dune::Fem::Parameter::getValue<double>("phasefield.adaptive.maxh",0.02 )),
      verbose_( Dune::Fem::Parameter::getValue<bool>("phasefield.adaptive.verbose",false) )
      {
        clear();
      }  
    
    // make this class behave as required for a LocalFunctionAdapter
    typedef Fem::FunctionSpace< double, double, dimension, 3 > FunctionSpaceType;
    void init(const ElementType &en)
    {
      enIndex_ = indexSet_.index(en);
      entity_ = &en;
    }
  
    template< class PointType >
    void evaluate(const PointType& x,FieldVector<double,3> &ret)
    {
    
      ret[0] = indicator_[enIndex_];
      ret[1] = refined_[enIndex_];
      ret[2] = grid_.getMark(*entity_);
   }
  private:
    const ElementType *entity_;
    int enIndex_;

  public:
    void clear ()
    {
      std::cout<<"clear\n";
      indicator_.resize( indexSet_.size( 0 ));
      refined_.resize( indexSet_.size( 0 ));
      refinedandcoarsened_.resize( indexSet_.size( 0 ));
      std::fill( indicator_.begin(), indicator_.end(), 0.);
      std::fill( refined_.begin(), refined_.end(), 0.);// std::numeric_limits<int>::min());
      std::fill( refined_.begin(), refined_.end(), 0.);
     }
    
    double estimate ( )
    {
      clear();
 
      const IteratorType end = dfSpace_.end();
      
      for( IteratorType it = dfSpace_.begin(); it != end; ++it )
      {
      
        const ElementType &entity = *it;
        int index=indexSet_.index(entity);
       
        std::vector<double> localIndicator; 
        double indicator;
        const LocalFunctionType uLocal = uh_.localFunction( entity );
  
        indicator=estimateLocal(entity, uLocal );
        localIndicator.push_back( indicator);
        IntersectionIteratorType end = gridPart_.iend( entity );
#if 0        
        for( IntersectionIteratorType inter = gridPart_.ibegin( entity ); inter != end; ++inter )
         {
        
           const IntersectionType &intersection = *inter;
            
           if( intersection.neighbor() )
            {
              const ElementPointerType pOutside = intersection.outside();
              const ElementType &outside = *pOutside;
              const LocalFunctionType uLocalNb = uh_.localFunction( outside);
              indicator=estimateLocal( entity, uLocal );
              localIndicator.push_back(indicator);
            }
         auto maxIndicator=std::max_element(localIndicator.begin(), localIndicator.end()); 
           indicator_[ index ]=*maxIndicator;
 
   //   estimateIntersection( intersection, entity, uLocal );
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
        std::cout<<"MaxIndicator="<<maxIndicator_<<"\n";	
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
          double gridFactor=indicator_[index];
          if(verbose_ && (gridFactor > 1.001))
            std::cout<<"gridFactor="<<refined_[index]<<"/ "<<indicator_[index]<<"="<<gridFactor<<"\n";
         
         if( std::abs(gridFactor)> tolerance)
	       {
	        
          if(entity.level()<maxLevel_)
		        {
              //if( refined_[index]!=1)
                {      
                  grid_.mark( 1, entity );
                  refined_[ index ] += 1;	
                }
              IntersectionIteratorType end = gridPart_.iend( entity );
		          for( IntersectionIteratorType inter = gridPart_.ibegin( entity ); inter != end; ++inter )
		          {
                const IntersectionType &intersection = *inter;
                if( intersection.neighbor() )
			          {
			            const ElementPointerType pOutside = intersection.outside();
			            const ElementType &outside = *pOutside;  
                  int outsideIndex=indexSet_.index(outside);
                  if(outside.level()<maxLevel_ );
			            {
               //     if(refined_[outsideIndex]<1)
                      { 
                        grid_.mark( 1, outside );
                        refined_[outsideIndex ] += 1;	
                      }
                    ++marked;
                  }
                }
		          }
            ++marked;
		       }
	        }
          else if( gridFactor < coarsen_ )
	        {
            //std::cout<<"GF="<<gridFactor<<"->Coarsen?\n";
            if( entity.level()>minLevel_ /* &&  (refined_[indexSet_.index(entity)]!=1)*/ )
            { 
		          if(refined_[index]==1)
               refinedandcoarsened_[index]=2.;
              else
               refinedandcoarsened_[index]=1.;
              
              grid_.mark(-1,entity);
            }
           }
	        else
          {
	        
          }
	  

	      }
      }
      
      return (marked > 0);
    
    }
    
    bool estimateAndMark(double tolerance)
    {
      double esti=estimate();
      return mark(tolerance);
    }
    

  protected:
    //! caclulate error on element 
    double estimateLocal ( const ElementType &entity,
                           const LocalFunctionType &uLocal )
    {
      const typename ElementType :: Geometry &geometry = entity.geometry();
      const Dune::ReferenceElement< double, dimension > &refElement =
      Dune::ReferenceElements< double, dimension >::general( entity.type() );
 
      const double volume = geometry.volume();
      double h2=std::sqrt(volume);
      const int index = indexSet_.index( entity );
     
      ElementQuadratureType quad( entity, 2*(dfSpace_.order() )+1 );
      const int numQuadraturePoints = quad.nop();
      double sigmasquared=0.;
      double maxsigma=0; 
      RangeType range;
      JacobianRangeType gradient;
      
      uLocal.evaluate( refElement.position(0 , 0),range );
      for( int i=0 ; i < dimension ; ++ i)
       sigmasquared+=PhasefieldFilter<RangeType>::sigma(range,i)*PhasefieldFilter<RangeType>::sigma(range,i);


      
#if 0 
  for( int qp = 0; qp < numQuadraturePoints; ++qp )
        {
          JacobianRangeType gradient;
	        RangeType range;
          double weight = quad.weight(qp);// * geometry.integrationElement( quad.point( qp )) ;
	  
	        uLocal.evaluate(quad[qp],range);
          
          for(int i=0 ; i<dimension; ++i)
            sigmasquared+=PhasefieldFilter<RangeType>::sigma(range,i)*PhasefieldFilter<RangeType>::sigma(range,i);

         sigmasquared*=weight;
        // maxsigma=std::max(sigmasquared,maxsigma);

        }
#endif
      //L2-Norm Sigma
      double normsigma=std::sqrt(sigmasquared);
      normsigma*=tolfactor_;
      indicator_[ index ]=h2*(normsigma+(1./maxH_));
      return normsigma;
    }

    
    
   
     
  
  };
}
#endif

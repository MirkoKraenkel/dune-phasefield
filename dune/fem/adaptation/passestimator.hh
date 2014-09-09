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
//#include <dune/fem/space/combinedspace.hh>
#include <dune/fem/quadrature/intersectionquadrature.hh>
#include <dune/fem/misc/h1norm.hh>
#include <dune/fem/misc/fmatrixconverter.hh>
namespace Dune
{
               
  template<class UFunction , class UGrad >
  class PassEstimator 
  {
   typedef PassEstimator< UFunction , UGrad > ThisType;
    
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
  
    typedef UGrad    DiscreteGradientType;
    typedef typename DiscreteGradientType :: DiscreteFunctionSpaceType DiscreteGradientSpaceType;
    typedef typename DiscreteGradientType :: LocalFunctionType LocalGradientType;
   
    typedef typename GridPartType :: GridType GridType;
    typedef typename GridPartType :: IndexSetType IndexSetType;
    typedef typename GridPartType :: IntersectionIteratorType IntersectionIteratorType;

    typedef typename IntersectionIteratorType :: Intersection IntersectionType;

    typedef typename GridType :: template Codim< 0 > :: Entity ElementType;
    typedef typename GridType :: template Codim< 0 > :: EntityPointer ElementPointerType;
    typedef typename ElementType::Geometry GeometryType;

    static const int dimension = GridType :: dimension;
    enum{ phaseId = dimension+1}
    typedef Fem::CachingQuadrature< GridPartType, 0 > ElementQuadratureType;
    typedef Fem::CachingQuadrature< GridPartType, 1 > FaceQuadratureType;
    typedef std :: vector< double > ErrorIndicatorType;

  protected:
    const DiscreteFunctionType &uh_;
    const DiscreteGradientType &sigmah_;
    const DiscreteFunctionSpaceType &dfSpace_;
    const GridPartType &gridPart_;
    const IndexSetType &indexSet_;
    GridType &grid_;
    ErrorIndicatorType indicator_;
    double totalIndicator2_,maxIndicator_;
    const double theta_;
    int maxLevel_;
    const double coarsen_;
    const double tolfactor_;
    const double maxH_;

  public:
    explicit PassEstimator (const DiscreteFunctionType &uh,
												    const DicsreteGradientType &sigmah, 
                            GridType &grid ):
      uh_( uh ),
      sigmah_(sigmah),
			dfSpace_( uh.space() ),
			gridPart_( dfSpace_.gridPart() ),
			indexSet_( gridPart_.indexSet() ),
			grid_( grid ),
			indicator_( indexSet_.size( 0 )),
			totalIndicator2_(0),
			maxIndicator_(0),
			theta_( Dune::Fem::Parameter::getValue("phasefield.adaptive.theta",0.) ),
			maxLevel_(Dune::Fem::Parameter::getValue<int>("fem.adaptation.finestLevel")),
			coarsen_(Dune::Fem::Parameter::getValue<double>("fem.adaptation.coarsenPercent",0))
      tolfactor_( Dune::Fem::Parameter::getValue<double>("phasefield.adaptive.tolfactor",1) ),
      maxH_( Dune::Fem::Parameter::getValue<double>("phasefield.adaptive.maxh",0.03125 )),
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
        const LocalGradientType sigmaLocal = sigmah_.localFunction( entity ); 
        estimateLocal(entity, uLocal , sigmaLocal );
       }
 
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
          double gridFactor=indicator_[index];
         
          if( std::abs(gridFactor)> tolerance)
	        {
            if(entity.level()<maxLevel_)
		          {
                grid_.mark( 1, entity );
                refined_[ index ] += 1;	
              
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
                            grid_.mark( 1, outside );
                            refined_[outsideIndex ] += 1;	
                            ++marked;
                          }
                      }     
		              }
            ++marked;
		       }
	        }
          else if( gridFactor < coarsen_*tolerance )
	        {
            //std::cout<<"GF="<<gridFactor<<"->Coarsen?\n";
            if( entity.level()>minLevel_  /*&&  (refined_[indexSet_.index(entity)]!=1)*/ )
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
    void estimateLocal ( const ElementType &entity,
                         const LocalFunctionType &uLocal,
                         const LocalGradientType &simgaLocal)
    {
      const typename ElementType :: Geometry &geometry = entity.geometry();
      const Dune::ReferenceElement< double, dimension > &refElement =
      Dune::ReferenceElements< double, dimension >::general( entity.type() );
 
      const double volume = geometry.volume();
      double h2=std::sqrt(volume);
      const int index = indexSet_.index( entity );
      GradientRangeType sigmaval; 
      JacobianRangeType graduval;
      double sigmasquared=0.
      sigmaLocal.evaluate( refElement.position( 0 , 0 ), sigmavla );
      Dune::Fem::FieldMatrixConverter< GradientRangeType , JacobianRangeType > graduval(sigma(val)); 
      for( int i = 0 ; i < dimension ; ++ i)
         sigmasquared+=graduval[phaseId][i]*graduval[phaseId][i]; 
      
      double normsigma=std::sqrt( sigmasquared );
      normsigma*=tolfactor_;
      indicator_[ index ] = h2*(normsigma+(1./maxH_));

    }

    
    
     

  
  };
}
#endif

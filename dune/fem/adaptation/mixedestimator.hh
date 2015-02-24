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
    mutable ErrorIndicatorType localsize_;
    mutable ErrorIndicatorType mark_;
    double totalIndicator_,maxIndicator_;
    const double theta_;
    int maxLevel_;
    int minLevel_;
    const double coarsen_;
    const double ifelements_;
    const double maxH_;
    const double localsizeFactor_;
    const bool verbose_;
    const bool secondNb_;
    const bool maxSigma_;
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
      localsize_( indexSet_.size( 0 )),
      mark_( indexSet_.size( 0 )),
      totalIndicator_(0),
      maxIndicator_(0),
      theta_( Dune::Fem::Parameter::getValue("phasefield.adaptive.theta",0.) ),
      maxLevel_(Dune::Fem::Parameter::getValue<int>("fem.adaptation.finestLevel")),
      minLevel_(Dune::Fem::Parameter::getValue<int>("fem.adaptation.coarsestLevel")),
      coarsen_(Dune::Fem::Parameter::getValue<double>("fem.adaptation.coarsenPercent",0.1)),
      ifelements_( Dune::Fem::Parameter::getValue<double>("phasefield.adaptive.tolfactor",1) ),
      maxH_( Dune::Fem::Parameter::getValue<double>("phasefield.adaptive.maxh",0.02 )),
      localsizeFactor_( Dune::Fem::Parameter::getValue<double>( "phasefield.adaptive.localsizefactor",0.4)),
      verbose_( Dune::Fem::Parameter::getValue<bool>("phasefield.adaptive.verbose",false) ),
      secondNb_( Dune::Fem::Parameter::getValue<bool>("phasefield.adaptive.secondneighbors",false)),
      maxSigma_( Dune::Fem::Parameter::getValue<bool>("phasefield.adaptive.maxsigma",false))
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
    
      ret[0] = indicator_[ enIndex_ ];
      ret[1] = localsize_[ enIndex_ ];
      ret[2] = mark_[ enIndex_ ];
   }
  private:
    const ElementType *entity_;
    int enIndex_;

  public:
    void clear ()
    {
      indicator_.resize( indexSet_.size( 0 ));
      localsize_.resize( indexSet_.size( 0 ));
      mark_.resize( indexSet_.size( 0 ));
      std::fill( indicator_.begin(), indicator_.end(), 0.);
      std::fill( localsize_.begin(), localsize_.end(), 0.);
      std::fill( mark_.begin(), mark_.end(), 0.);
     }
    
    void estimate ( )
    {
      clear();
 
      const IteratorType end = dfSpace_.end();
      
      for( IteratorType it = dfSpace_.begin(); it != end; ++it )
      {
      
        const ElementType &entity = *it;
       
        std::vector<double> localIndicator; 
        double indicator;
        const LocalFunctionType uLocal = uh_.localFunction( entity );
  
        indicator=estimateLocal(entity, uLocal );
        localIndicator.push_back( indicator);
        IntersectionIteratorType end = gridPart_.iend( entity );
      }
    }
    
    double computeIndicator()
    {
      totalIndicator_ = 0.0;
  
        for (unsigned int i=0;i<indicator_.size();++i)
        {
          totalIndicator_ += indicator_[ i ];
          maxIndicator_ = std::max(maxIndicator_,indicator_[i]);
        }
	      return sqrt( totalIndicator_ );
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
          double gridFactor = indicator_[index];

          if( gridFactor < localsizeFactor_*localsize_[index])
          {
	        
            if(entity.level()<maxLevel_)
		        {
              grid_.mark( 1, entity );
              mark_[index]=1;
            }

            IntersectionIteratorType end = gridPart_.iend( entity );
		          
            for( IntersectionIteratorType inter = gridPart_.ibegin( entity ); inter != end; ++inter )
		        {
              const IntersectionType &intersection = *inter;
              if( intersection.neighbor() )
			        {
			          const ElementPointerType poutside = intersection.outside();
			          const ElementType &outside = *poutside;
                int indexnb=indexSet_.index(outside);
                if(outside.level()<maxLevel_ );
			          {
                  grid_.mark( 1, outside );
                  mark_[ indexnb ]=1;
                  ++marked;
                }

                IntersectionIteratorType end = gridPart_.iend( outside );
                if(secondNb_)
                  for(IntersectionIteratorType inter2 = gridPart_.ibegin( outside ); inter2!=end; ++inter2)
                  {
                    const IntersectionType &intersection2=*inter2;

                    if(intersection2.neighbor())
                    {
                      const ElementPointerType poutside2 = intersection2.outside();
			                const ElementType &outside2 = *poutside2;
                      int indexnb2 = indexSet_.index(outside2);
                      if(outside.level()<maxLevel_ );
			                {
                        grid_.mark( 1, outside2 );
                        mark_[ indexnb2 ]=1;
                        ++marked;
                      }
                    }
                  }
		            }
              }
              ++marked;
		        }
            else if ( gridFactor > localsize_[index])
	          {
              if( entity.level()>minLevel_ )
              {
                grid_.mark(-1,entity);
                mark_[ index ]=-1;
              }
            }
	          else
            {
	            mark_[ index ]=0;
            }
	  

	        }
      }
      
      return (marked > 0);
    
    }
    
    bool estimateAndMark(double tolerance)
    {
      estimate();
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
     
      ElementQuadratureType quad( entity, 2*(dfSpace_.order() ) );
      const int numQuadraturePoints = quad.nop();
      double sigmasquared=0.;
      double maxsigma=0;
      RangeType range;
      JacobianRangeType gradient;
      
      uLocal.evaluate( refElement.position(0 , 0),range );

      
    for( int qp = 0; qp < numQuadraturePoints; ++qp )
        {
          JacobianRangeType gradient;
	        RangeType range;
          double weight = quad.weight(qp);
	  
	        uLocal.evaluate(quad[qp],range);
          if( maxSigma_)
          {
            sigmasquared=0;
            for(int i=0 ; i<dimension; ++i)
              sigmasquared+=PhasefieldFilter<RangeType>::sigma(range,i)
                            *PhasefieldFilter<RangeType>::sigma(range,i);
            maxsigma=std::max(sigmasquared,maxsigma);
          }
          else
          {
            for(int i=0 ; i<dimension; ++i)
              maxsigma+=PhasefieldFilter<RangeType>::sigma(range,i)
                        *PhasefieldFilter<RangeType>::sigma(range,i)*weight;
          }


        }
      //L2-Norm Sigma
      double normsigma=std::sqrt(maxsigma);
      normsigma*=ifelements_;
      double indicatedsize=1/(normsigma+(1./maxH_)); 
      localsize_[ index ]=h2;
      indicator_[ index ]=indicatedsize;
      return normsigma;
    }

    
    
   
     
  
  };
}
#endif

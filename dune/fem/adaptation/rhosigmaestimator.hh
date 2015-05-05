#ifndef RHOSIGMAESTIMATOR_HH 
#define RHOSIGMAESTIMATOR_HH

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
    
  private:
    const DiscreteFunctionType &uh_;
    const DiscreteFunctionSpaceType &dfSpace_;
    const GridPartType &gridPart_;
    const IndexSetType &indexSet_;
    GridType &grid_;
    ErrorIndicatorType sigmaIndicator_;
    ErrorIndicatorType rhoIndicator_;
    mutable ErrorIndicatorType localsize_;
    mutable ErrorIndicatorType mark_;
    double totalIndicator_,maxIndicator_;
    int maxLevel_;
    int minLevel_;
    const double coarsen_;
    const double ifelements_;
    const double maxH_;
    const double localsizeFactor_;
    const double rhoTol_;
    const bool verbose_;
    const bool secondNb_;
    const bool maxSigma_;
    const ElementType *entity_;
    int enIndex_;

 public:
    explicit MixedEstimator ( const DiscreteFunctionType &uh ,
                              GridType &grid,
                              const  ModelType& model):
      uh_( uh ),
      dfSpace_( uh.space() ),
      gridPart_( dfSpace_.gridPart() ),
      indexSet_( gridPart_.indexSet() ),
      grid_( grid ),
      sigmaIndicator_( indexSet_.size( 0 )),
      rhoIndicator_( indexSet_.size( 0 )),
      localsize_( indexSet_.size( 0 )),
      mark_( indexSet_.size( 0 )),
      totalIndicator_(0),
      maxIndicator_(0),
      maxLevel_(Dune::Fem::Parameter::getValue<int>("fem.adaptation.finestLevel")),
      minLevel_(Dune::Fem::Parameter::getValue<int>("fem.adaptation.coarsestLevel")),
      coarsen_(Dune::Fem::Parameter::getValue<double>("fem.adaptation.coarsenPercent",0.1)),
      ifelements_( Dune::Fem::Parameter::getValue<double>("phasefield.adaptive.tolfactor",1) ),
      maxH_( Dune::Fem::Parameter::getValue<double>("phasefield.adaptive.maxh",0.02 )),
      localsizeFactor_( Dune::Fem::Parameter::getValue<double>( "phasefield.adaptive.localsizefactor",0.4)),
      rhoTol_( Dune::Fem::Parameter::getValue<double>("phasefield.adaptive.rhoTol",0.1)),
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
    
      ret[0] = sigmaIndicator_[ enIndex_ ];
      ret[1] = rhoIndicator_[ enIndex_ ];
      ret[2] = mark_[ enIndex_ ];
    }
  public:
    void clear ()
    {
      sigmaIndicator_.resize( indexSet_.size( 0 ));
      rhoIndicator_.resize( indexSet_.size(0));
      localsize_.resize( indexSet_.size( 0 ));
      mark_.resize( indexSet_.size( 0 ));
      std::fill( sigmaIndicator_.begin(), sigmaIndicator_.end(), 0.);
      std::fill( localsize_.begin(), localsize_.end(), 0.);
      std::fill( rhoIndicator_.begin(), rhoIndicator_.end(), 0. );
      std::fill( mark_.begin(), mark_.end(), 0.);
     }
    
    void estimate ( )
    {
      clear();
      EstimateFunctor f(*this);
      dfSpace_.forEach( f );
    }
   
    
    bool refine( int index ) const
    {
      if( rhoTol_>0)      
      {
        return  ( sigmaIndicator_[ index ] < localsizeFactor_*localsize_[ index ]
                || rhoIndicator_[ index ] > rhoTol_);
      }
      else
      {
        return sigmaIndicator_[ index ] < localsizeFactor_*localsize_[ index ];
      }
    }
    
    bool coarsen( int index ) const
    {
      return ( sigmaIndicator_[ index ] > localsize_[ index ]
              && rhoIndicator_[ index ] < coarsen_*rhoTol_);
    }
    double computeIndicator()
    {
      totalIndicator_ = 0.0;
  
        for (unsigned int i=0 ; i<sigmaIndicator_.size() ; ++i)
        {
          totalIndicator_ += sigmaIndicator_[ i ];
          maxIndicator_ = std::max(maxIndicator_,sigmaIndicator_[i]);
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
          //double gridFactor = Sigmaindicator_[index];
          
          if( refine( index ) )
          {
	        
            if(entity.level()<maxLevel_)
		        {
              grid_.mark( 1 , entity );
              mark_[ index ] = 1;
              ++marked;
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
                  grid_.mark( 1 , outside );
                  mark_[ indexnb ] = 1;
                  ++marked;
                }
// loop over second neigbors
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
                        grid_.mark( 1 , outside2 );
                        mark_[ indexnb2 ] = 1;
                        ++marked;
                      }
                    }
                  }
// loop over second neighbors
                }
              }
		        }
            else if ( coarsen( index ) )
	          {
              if( entity.level()>minLevel_ )
              {
                grid_.mark( -1 , entity );
                mark_[ index ] = -1;
              }
            }
	          else
            {
	            grid_.mark( 0 , entity );
              mark_[ index ] = 0;
            }
	        }
        }
     return (marked > 0);
    }
    
    bool estimateAndMark ( double tolerance )
    {
      estimate();
      return mark(tolerance);
    }
 
    

  protected:
    
    
    //! caclulate error on element 
    double estimateLocal ( const ElementType &entity, const LocalFunctionType &uLocal )
    {
      const int index = indexSet_.index( entity );
   
      const typename ElementType :: Geometry &geometry = entity.geometry();
      const Dune::ReferenceElement< double, dimension > &refElement =
      Dune::ReferenceElements< double, dimension >::general( entity.type() );
      
      //const LocalFunctionType uLocal = uh_.localFunction( entity );
      
      const double volume = geometry.volume();
      //double h2 = (dimension == 2 ? volume : std :: pow( volume, 2.0 / (double)dimension ));  
      double localh = std::sqrt(volume);
    
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
        if( maxSigma_ )
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
      localsize_[ index ]=localh;
      sigmaIndicator_[ index ]=indicatedsize;
      return normsigma;
    }
  
    //jumEstimator for density
    void estimaterhojumpLocal ( const ElementType &entity , const LocalFunctionType &uLocal )
    {
      const int insideIndex = indexSet_.index( entity );    
 
      ElementQuadratureType quad( entity, 2*(dfSpace_.order() )+1 );
      const int numQuadraturePoints = quad.nop();
      RangeType range( 0. );
      RangeType rangeNb( 0. );

      double indLocMax=-1E100;
      double indLocMin= 1E100;
      
    for( int qp = 0; qp < numQuadraturePoints; ++qp )
    {
      DomainType xEn=quad.point( qp );
      uLocal.evaluate(quad[qp],range);
      const double rhoqp=range[0];
      indLocMax=std::max( indLocMax , rhoqp );
      indLocMin=std::min( indLocMin , rhoqp );
	    
    }

    const IntersectionIteratorType iend = gridPart_.iend( entity );
    
    for ( IntersectionIteratorType nb = gridPart_.ibegin( entity);nb!=iend; ++nb)
    {
      const IntersectionType& intersection = *nb;
      if( intersection.neighbor() )
      {
        ElementPointerType outp = intersection.outside();
        const ElementType& outside = *outp;
        const int outsideIndex = indexSet_.index( outside );

        if( outside.partitionType() == Dune::GhostEntity ||
            entity.level() > outside.level() ||
            (entity.level() == outside.level() && insideIndex < outsideIndex ) )
        {  
          const LocalFunctionType uNbLocal = uh_.localFunction( outside );
          ElementQuadratureType quadNb( entity, 2*(dfSpace_.order() )+1);
          const int numQuadraturePointsNb = quadNb.nop();

          double indLocNbMax=-1E100;
          double indLocNbMin= 1E100;
            
          for( int qpNb=0; qpNb<numQuadraturePointsNb; ++qpNb)
            {
              DomainType xNb=quadNb.point( qpNb );
              uNbLocal.evaluate( xNb, rangeNb );
              const double rhoNb=rangeNb[0];
              indLocNbMax=std::max( indLocNbMax , rhoNb );
              indLocNbMin=std::min( indLocNbMin , rhoNb );
            }

            double indLocalDiff = std::abs( indLocMax - indLocNbMin );
            indLocalDiff = std::max( indLocalDiff, std::abs( indLocMin -indLocNbMax) );
            rhoIndicator_[ insideIndex ]  = std::max( rhoIndicator_[ insideIndex ]  , indLocalDiff );
            rhoIndicator_[ outsideIndex ] = std::max( rhoIndicator_[ outsideIndex ] , indLocalDiff );
          } 
        }
      }
    }

  private:
    class EstimateFunctor
    {
      ThisType& parent_;
       
    public:
      EstimateFunctor(ThisType& parent):parent_(parent)
      {}
    
      void operator()( const ElementType& e ) const
      {
        const LocalFunctionType& uLocal=parent_.uh_.localFunction( e );
        parent_.estimateLocal( e , uLocal );
        parent_.estimaterhojumpLocal( e , uLocal );
      }
     
    
    };
    
    
    
   
     
  
  };
}
#endif

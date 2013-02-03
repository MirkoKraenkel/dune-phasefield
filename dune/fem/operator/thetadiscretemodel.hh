#ifndef DUNE_PHASEFIELD_THETADISCREMODEL_HH
#define DUNE_PHASEFIELD_THETADISCREMODEL_HH


// Dune includes
#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>

// Dune-Fem includes
#include <dune/fem/space/dgspace.hh>
#include <dune/fem/pass/dgdiscretemodel.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/misc/boundaryidentifier.hh>
#include <dune/fem/misc/fmatrixconverter.hh>

// local includes
#include "projectiondiscretemodelcommon.hh"
#include <dune/fem/fluxes/ldgfluxwellbalanced.hh>
namespace Dune {
  /*
   *  DiscreteModel for calculation the chemical potential and Allen-Cahn Resisduum
   *  res[0]=\mu(\rho,\phi), res[1]=f_\phi-\nabla\
   */

	    
  template< class Model,int polOrd, int passUId, int passGradId> class ThetaModel;
	
	template <class Model,int polOrd, int passUId, int passGradId>
	struct ThetaTraits
	{
		typedef typename Model :: Traits                                 ModelTraits;
		typedef typename ModelTraits :: GridType                         GridType;

		enum { dimRange = 2 };
		enum { dimDomain = ModelTraits::dimDomain };

		typedef MyPassTraits< Model, dimRange, polOrd >                    Traits;
		//2 dimensional Space
		typedef typename Traits :: FunctionSpaceType                     FunctionSpaceType;
		
		typedef typename Traits :: VolumeQuadratureType                  VolumeQuadratureType;
		typedef typename Traits :: FaceQuadratureType                    FaceQuadratureType;
		typedef typename Traits :: DiscreteFunctionSpaceType             DiscreteFunctionSpaceType;
		
		typedef typename Traits :: GridPartType                          GridPartType;
		
// 		typedef typename GridPartType :: IntersectionIteratorType          IntersectionIterator;
//     typedef typename IntersectionIterator :: Intersection              Intersection;
//     typedef typename DiscreteFunctionSpaceType :: EntityType           EntityType;
 		typedef typename ModelTraits::FaceDomainType              FaceDomainType;

		typedef typename Traits :: DestinationType                       DestinationType;
		typedef DestinationType                                          DiscreteFunctionType;
		typedef typename Traits :: IndicatorType                         IndicatorType;

		typedef typename DestinationType :: DomainType                   DomainType;
		//fieldvector<2,ctype>
		typedef typename DestinationType :: RangeType                    RangeType;
		typedef typename DestinationType :: RangeFieldType               RangeFieldType;
		typedef typename DestinationType :: DomainFieldType              DomainFieldType;

		//fieldmatrix<2,dimdomain,ctype>
		typedef typename DestinationType :: JacobianRangeType            JacobianRangeType;

		typedef typename Traits :: AdaptationHandlerType  AdaptationHandlerType ;

		typedef ThetaModel< Model,polOrd, passUId, passGradId >       DGDiscreteModelType;
	};



	template<class Model,int polOrd, int passUId,int passGradId>
	class ThetaModel:
		public Fem::DGDiscreteModelDefaultWithInsideOutside
		<ThetaTraits<Model,polOrd,passUId,passGradId>,passUId, passGradId>
	{
 			public:

		typedef Fem::DGDiscreteModelDefaultWithInsideOutside
		<ThetaTraits<Model,polOrd,passUId,passGradId>,passUId, passGradId> BaseType;

		using BaseType :: inside;
		using BaseType :: outside;

		integral_constant< int, passUId > uVar;
		integral_constant< int, passGradId > sigmaVar;
	public:
		typedef ThetaTraits< Model, polOrd, passUId,passGradId >         Traits;
		enum{dimDomain=Traits::dimDomain};

		typedef FieldVector< double, Traits :: dimDomain >               DomainType;
		typedef typename Traits :: RangeType                             RangeType;
		typedef typename Traits :: JacobianRangeType                     JacobianRangeType;
		typedef typename Traits :: GridType                              GridType; 
		
		
		typedef typename GridType :: template Codim< 0 > :: Entity       EntityType;
		typedef typename Traits::FaceDomainType FaceDomainType;

		typedef typename Traits :: GridPartType
		:: IntersectionIteratorType                            IntersectionIterator;
		typedef typename IntersectionIterator :: Intersection            Intersection;
		enum { dimRange = Traits::dimRange };
		enum { evaluateJacobian = false };
		static const bool ApplyInverseMassOperator = true;

	public:
    
		ThetaModel(const Model& mod):
			acpenalty_(Fem::Parameter::getValue<double>("phasefield.penalty")),
			model_(mod){}
      
		bool hasSource() const { return true; } 
		bool hasFlux() const { return true; } 
  
		template <class ArgumentTuple, class JacobianTuple>
		double source( const EntityType& en, 
									 const double time,
									 const DomainType& x,
									 const ArgumentTuple& u,
									 const JacobianTuple& jac, 
									 RangeType& s )
		{ 
			
			s = 0;
		//s[0]=dF/drho  s[1]=-dF/dphi
      const double dtStiff = model_.thetaSource( en, time, x, u[uVar], u[sigmaVar], s );
			
			return dtStiff;
		} 

		template <class QuadratureImp, class ArgumentTupleVector > 
		void initializeIntersection(const Intersection& it,
																const double time,
																const QuadratureImp& quadInner, 
																const QuadratureImp& quadOuter,
																const ArgumentTupleVector& uLeftVec,
																const ArgumentTupleVector& uRightVec) 
		{
		}

		template <class QuadratureImp, class ArgumentTupleVector > 
		void initializeBoundary(const Intersection& it,
														const double time,
														const QuadratureImp& quadInner, 
														const ArgumentTupleVector& uLeftVec)
		{
		}

		template <class QuadratureImp,class ArgumentTuple, class JacobianTuple>
		double boundaryFlux(const Intersection& it,
												double time, 
												const QuadratureImp& faceQuadInner,
												const int quadPoint,
												const ArgumentTuple& uLeft,
												const JacobianTuple& jacLeft,
												RangeType& gLeft,
												JacobianRangeType& gDiffLeft ) const   /*@LST0E@*/
		{
			const FaceDomainType& x = faceQuadInner.localPoint( quadPoint );
			const DomainType normal = it.integrationOuterNormal(x);
			
			JacobianRangeType diffmatL(0.),diffmatR(0.);
 	
		 	model_.boundaryallenCahnDiffusion(uLeft[uVar],uLeft[sigmaVar] , diffmatL);
			
			diffmatL.mv(normal, gLeft);
			gDiffLeft = 0;
			

			// add penalty term ( enVolume() is available since we derive from
			//    DiscreteModelDefaultWithInsideOutside)
			const double factor = acpenalty_  ;
			
			


			double diffTimeStep(0.);
			return diffTimeStep;
	
		}

		template <class ArgumentTuple, class JacobianTuple>    /*@LST1S@*/
		void analyticalFlux(const EntityType& en,
												const double time, const DomainType& x,
												const ArgumentTuple& u, 
												const JacobianTuple& jac,
												JacobianRangeType& f)
		{
				model_.allenCahnDiffusion(u[uVar],u[sigmaVar], f);
		}                     

           
		template< class QuadratureImp,
							class ArgumentTuple, 
							class JacobianTuple >          /*@LST0S@*/
		double numericalFlux(const Intersection& it,
												 const double time,
												 const QuadratureImp& faceQuadInner,
												 const QuadratureImp& faceQuadOuter,
												 const int quadPoint, 
												 const ArgumentTuple& uLeft,
												 const ArgumentTuple& uRight,
												 const JacobianTuple& jacLeft,
												 const JacobianTuple& jacRight,
												 RangeType& gLeft,
												 RangeType& gRight,
												 JacobianRangeType& gDiffLeft,
												 JacobianRangeType& gDiffRight )
		{
			

			const FaceDomainType& x = faceQuadInner.localPoint( quadPoint );
			const DomainType normal = it.integrationOuterNormal(x);
			
			JacobianRangeType diffmatL(0.),diffmatR(0.);
 	
		 	model_.allenCahnDiffusion(uLeft[uVar],uLeft[sigmaVar] , diffmatL);
 			model_.allenCahnDiffusion(uRight[uVar],uLeft[sigmaVar], diffmatR);
			
			diffmatL.mv(normal, gLeft);
			diffmatR.umv(normal, gLeft);
			gLeft*=0.5;
			gRight=gLeft;
			gDiffLeft = 0;
      gDiffRight = 0;     
			

			// add penalty term ( enVolume() is available since we derive from
			//    DiscreteModelDefaultWithInsideOutside)
			const double factor = acpenalty_  ;
			
	 		double jmp( uLeft[uVar][dimDomain+1] );
			jmp -= uRight[uVar][dimDomain+1];
			RangeType jump(0.);
			jump[1]=jmp;
			gLeft.axpy(factor, jump);
			


			double diffTimeStep(0.);
			return diffTimeStep;
		}

		template <class ArgumentTuple, class JacobianTuple >
    void analyticalFlux( const EntityType& en,
                         const double time,
                         const DomainType& x,
                         const ArgumentTuple& u,
                         const JacobianTuple& jac,
                         JacobianRangeType& f ) const
    {
			JacobianRangeType diffmatrix;
		
			model_.allenCahndiffusion(u[sigmaVar], diffmatrix);
			
			f = diffmatrix;
    }




	private:
		const double acpenalty_;
		const Model& model_;
	};
}//end namespce
#endif

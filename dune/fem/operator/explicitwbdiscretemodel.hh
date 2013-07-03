#ifndef DUNE_PHASEFIELD_EXPLICITWBDISCRETEMODEL_HH
#define DUNE_PHASEFIELD_EXPLICITWBDISCRETEMODEL_HH

// Dune includes
#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>

// Dune-Fem includes
#include <dune/fem/space/discontinuousgalerkin.hh>
#include <dune/fem/pass/localdg/discretemodel.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/misc/boundaryidentifier.hh>
#include <dune/fem/misc/fmatrixconverter.hh>

// local includes
#include "projectiondiscretemodelcommon.hh"
#include <dune/fem/fluxes/ldgfluxwellbalanced.hh>

//*************************************************************
namespace Dune {


  // for the 3rd Pass
	template< class Model, 
						class NumFlux, 
						int polOrd, int passUId,int passProjId,int passGradId,
						bool advectionPartExists, bool diffusionPartExists >
	class ExplicitPhasefieldLDGModel;

  
	template <class Model, class NumFlux,
						int polOrd, int passUId, int passProjId, int passGradId, 
						bool advectionPartExists, bool diffusionPartExists >
	struct ExplicitPhasefieldLDGTraits
		: public AdvectionProjTraits< Model, NumFlux, polOrd, passUId,passProjId, passGradId, advectionPartExists>
          
	{
		typedef AdvectionDiffusionLDGModel
		< Model, NumFlux, polOrd, passUId,passProjId, passGradId, 
			advectionPartExists, diffusionPartExists >                   DGDiscreteModelType;
	};


	// AdvectionDiffusionLDGModel
	//---------------------------
	template< class Model, 
						class NumFlux, 
						int polOrd, int passUId, int passProjId,int passGradId,
						bool advectionPartExists, bool diffusionPartExists >
	class ExplicitPhasefieldLDGModel :
		public AdvectionProjModel< Model, NumFlux, polOrd, passUId, passProjId, passGradId, advectionPartExists > 
	{
	public:
		typedef ExplicitPhasefieldLDGTraits< Model, 
																				 NumFlux, 
																				 polOrd,
																				 passUId, 
																				 passProjId,
																				 passGradId,
																				 advectionPartExists, 
																				 diffusionPartExists >  Traits; 
		
		typedef AdvectionProjModel< Model, 
																NumFlux, 
																polOrd, 
																passUId,
																passProjId,
																passGradId, 
																advectionPartExists>    BaseType;
		
		using BaseType :: inside ;
		using BaseType :: outside ;
		using BaseType :: model_;
		using BaseType :: uBnd_;

		using BaseType :: uVar;
    
		// defined in AdvectionModel 
		using BaseType :: maxAdvTimeStep_ ;
		using BaseType :: maxDiffTimeStep_ ;


		integral_constant< int, passGradId> sigmaVar;

		integral_constant< int, passProjId> thetaVar;
    
	public:
		enum { dimDomain = Traits :: dimDomain };
		enum { dimRange  = Traits :: dimRange };

		enum { advection = advectionPartExists };
		enum { diffusion = diffusionPartExists }<F2>;

		typedef FieldVector< double, dimDomain >               DomainType;
		typedef FieldVector< double, dimDomain-1 >             FaceDomainType;

#if defined TESTOPERATOR
		static const bool ApplyInverseMassOperator = false;
#else
		static const bool ApplyInverseMassOperator = true;
#endif

		typedef typename Traits :: GridPartType                            GridPartType;
		typedef typename Traits :: GridType                                GridType;
		typedef typename GridPartType :: IntersectionIteratorType          IntersectionIterator;
		typedef typename IntersectionIterator :: Intersection              Intersection;
		typedef typename GridType :: template Codim< 0 > :: Entity         EntityType;
		typedef typename GridType :: template Codim< 0 > :: EntityPointer  EntityPointerType;
		typedef typename Traits :: RangeFieldType                          RangeFieldType;
		typedef typename Traits :: DomainFieldType                         DomainFieldType;
		typedef typename Traits :: RangeType                               RangeType;
		typedef typename Traits :: JacobianRangeType                       JacobianRangeType;

	
		typedef typename Traits :: DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

		typedef LDGDiffusionFlux< DiscreteFunctionSpaceType, Model> DiffusionFluxType;
		enum { evaluateJacobian = true };

  public:

		ExplicitPhasefieldLDGModel(const Model& mod,
															 const NumFlux& numflux)
			: BaseType( mod, numflux)
		{
    }

    bool hasSource() const
		{                
			return ( model_.hasNonStiffSource(); );
		} 

		bool hasFlux() const { return advection;  };      
	
    template <class ArgumentTuple, class JacobianTuple >
		double source( const EntityType& en,
									 const double time, 
									 const DomainType& x,
									 const ArgumentTuple& u, 
									 const JacobianTuple& jac, 
									 RangeType& s ) const
		{
			s = 0;

      //we need \sigma_phi= u[sigmaVar][dimDomain+1]
			//        \theta_1=u[thetaVar][1]
			//        \nablan\theta_2 jac[thetatVar]

      double dtEst = std::numeric_limits< double > :: max();
      //nonStiffSource should give the source term for explicit treatment
			const double dtStiff = model_.nonStiffSource( en, time, x, u[uVar],u[sigmaVar],u[thetaVar],jac[thetaVar],jac[uVar], s );
		  
      dtEst = ( dtStiff > 0 ) ? dtStiff : dtEst;
			maxDiffTimeStep_ = std::max( dtStiff, maxDiffTimeStep_ );
	
			
			// return the fastest wave from source terms
			return dtEst;
		}

		void switchUpwind() const 
		{ 
			BaseType :: switchUpwind(); 
			diffFlux_.switchUpwind(); 
		}

	public:
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
			// advection
		  const FaceDomainType& x = faceQuadInner.localPoint( quadPoint );
     
      const DomainType normal = it.integrationOuterNormal( x );
  
      double wave = BaseType :: numericalFlux( it, time, faceQuadInner, faceQuadOuter,
																										 quadPoint, uLeft, uRight, jacLeft, jacRight, 
																										 gLeft, gRight, gDiffLeft, gDiffRight );
      
		  RangeType nonCons(0.);
 
      RangeType average(0.);
      double phiLeft,phiRight;
      double rhoLeft,rhoRight;
      rhoLeft  = uLeft[uVar][0];
      rhoRight = uRight[uVar][0];
      phiLeft  = uLeft[uVar][dimDomain+1];
      phiRight = uRight[uVar][dimDomain+1];
      phiLeft/=rhoLeft;
      phiRight/=rhoRight;

      
      // {{rho}}
      average[1]=uLeft[uVar][0]+uRight[uVar][0];
      average*=0.5;
   
      //[[\mu]]
      for(int i=0;i<dimDomain;i++)
      {  
        nonCons[i+1]=normal[i];
        nonCons[i+1]*=uLeft[thetaVar][0]-uRight[thetaVar][0];
        nonCons[i+1]*=average[1];
      }   
 
      average[1]=0.5*(uLeft[thetaVar][1]+uRight[thetaVar][1])*(phiLeft-phiRight);
   
      //nonCon
#if USEJACOBIAN
      nonCons[1]=average[1];
#else
      nonCons[1]=0;
#endif

    
      //factor comes from the meanvalue of the testfunctions
      nonCons*=0.5;
          
      //{{\theta}}[[phi]]
      gLeft +=nonCons;
      gRight-=nonCons;



			gDiffLeft  = 0;
			gDiffRight = 0;

			maxAdvTimeStep_  = std::max( wave, maxAdvTimeStep_ );
			maxDiffTimeStep_ = std::max( diffTimeStep, maxDiffTimeStep_ );
			return std::max( wave, diffTimeStep );
		}


		/**
		 * @brief same as numericalFlux() but for fluxes over boundary interfaces
		 */
		template <class QuadratureImp, 
							class ArgumentTuple, class JacobianTuple>
		double boundaryFlux(const Intersection& it,
												const double time, 
												const QuadratureImp& faceQuadInner,
												const int quadPoint,
												const ArgumentTuple& uLeft,
												const JacobianTuple& jacLeft,
												RangeType& gLeft,
												JacobianRangeType& gDiffLeft ) const   /*@LST0E@*/
		{

			// advection

			const double wave = BaseType :: 
				boundaryFlux( it, time, faceQuadInner, quadPoint,
											uLeft, jacLeft, gLeft, gDiffLeft );
                                  
			// diffusion
      
			double diffTimeStep = 0.0;

			bool hasBoundaryValue = 
				model_.hasBoundaryValue( it, time, faceQuadInner.localPoint(0) );

			gDiffLeft = 0;


			maxAdvTimeStep_  = std::max( wave, maxAdvTimeStep_ );
			maxDiffTimeStep_ = std::max( diffTimeStep, maxDiffTimeStep_ );

			return std::max( wave, diffTimeStep );

		}
		template <class ArgumentTuple, class JacobianTuple >
		void analyticalFlux( const EntityType& en,
												 const double time, 
												 const DomainType& x,
												 const ArgumentTuple& u, 
												 const JacobianTuple& jac, 
												 JacobianRangeType& f ) const
		{
      // advection
			BaseType :: analyticalFlux( en, time, x, u, jac, f );

  	}
	};                                              /*@LST0E@*/

}
#include "thetadiscretemodel.hh"
#endif

#ifndef DUNE_PHASEFIELD_IMPLICITWBDISCRETEMODEL_HH
#define DUNE_PHASEFIELD_IMPLICITWBDISCRETEMODEL_HH

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

 
	// AdvectionDiffusionLDGModel
	//---------------------------
  // for the 3rd Pass
	template< class Model, 
						class NumFlux, 
						int polOrd, int passUId,int passProjId,int passGradId,
						bool advectionPartExists, bool diffusionPartExists >
	class ImplicitPhasefieldLDGModel;


	// AdvectionDiffusionLDGTraits
	//----------------------------
  
	template <class Model, class NumFlux,
						int polOrd, int passUId, int passProjId, int passGradId, 
						bool advectionPartExists, bool diffusionPartExists >
	struct ImplicitPhasefieldLDGTraits 
		: public AdvectionProjTraits< Model, NumFlux, polOrd, passUId,passProjId, passGradId, advectionPartExists>
          
	{
		typedef ImplicitPhasefieldLDGModel
		< Model, NumFlux, polOrd, passUId,passProjId, passGradId, 
			advectionPartExists, diffusionPartExists >                   DGDiscreteModelType;
	};


	//ImplicitDGModel
	//---------------------------
	template< class Model, 
						class NumFlux, 
						int polOrd, int passUId, int passProjId,int passGradId,
						bool advectionPartExists, bool diffusionPartExists >
	class ImplicitPhasefieldLDGModel :
		public AdvectionProjModel< Model, NumFlux, polOrd, passUId, passProjId, passGradId, advectionPartExists > 
	{
	public:
		typedef ImplicitPhasefieldLDGTraits< Model, 
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
		enum { diffusion = diffusionPartExists };

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
		/**
		 * @brief constructor
		 */
		ImplicitPhasefieldLDGModel(const Model& mod,
															 const NumFlux& numf,
															 DiffusionFluxType& diffflux)
			: BaseType( mod, numf ),
				diffFlux_( diffflux ),
				penalty_( 1.0 ),
        switch_(Fem::Parameter::getValue<double>("phasefield.thetaswitch")),
				cflDiffinv_( 8.0 * ( polOrd + 1) )
		{
    }
		bool hasSource() const
		{                
			return ( model_.hasNonStiffSource() || model_.hasStiffSource() );
		} 

		bool hasFlux() const { return advection || diffusion; };      
		/**
		 * @brief analytical flux function$
		 */
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
      //stiffSource should be the monotone part of df/dphi
			const double dtStiff = model_.stiffSource( en, time, x, u[uVar],u[sigmaVar],u[thetaVar],jac[thetaVar],jac[uVar], s );
		  
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
		  const FaceDomainType& x = faceQuadInner.localPoint( quadPoint );
     
      const DomainType normal = it.integrationOuterNormal( x );
  
     
			// diffusion
			double diffTimeStep = 0.0;
			if( diffusion ) 
				{ 
					RangeType dLeft, dRight;
					diffTimeStep = diffFlux_.numericalFlux(it, *this,
																								 time, faceQuadInner, faceQuadOuter, quadPoint,
																								 uLeft[ uVar ], uRight[ uVar ], 
																								 uLeft[ sigmaVar ], uRight[ sigmaVar ], 
																								 dLeft, dRight,
																								 gDiffLeft, gDiffRight);

					gLeft  = dLeft;
					gRight = dRight;
				}
//should the noncon part be here?
#if 0
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
 //     average[1]=(phiLeft-phiRight);
   
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
#endif

			gDiffLeft  = 0;
			gDiffRight = 0;

			maxDiffTimeStep_ = std::max( diffTimeStep, maxDiffTimeStep_ );
			return  diffTimeStep;
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

			if( diffusion && hasBoundaryValue ) 
				{
					// diffusion boundary flux for Dirichlet boundaries 
					RangeType dLeft;
					diffTimeStep = diffFlux_.boundaryFlux(it, 
																								*this, 
																								time, faceQuadInner, quadPoint,
																								uLeft[ uVar ], uBnd_, // is set during call of  BaseType::boundaryFlux
																								uLeft[ sigmaVar ], 
																								dLeft,
																								gDiffLeft);
					gLeft += dLeft;
				}
			else if ( diffusion )
				{
					RangeType diffBndFlux;
					model_.diffusionBoundaryFlux( it, time, faceQuadInner.localPoint(quadPoint),
																				uLeft[uVar], jacLeft[uVar], diffBndFlux );
					gLeft += diffBndFlux;
				}
			else
				gDiffLeft = 0;


			maxAdvTimeStep_  = std::max( wave, maxAdvTimeStep_ );
			maxDiffTimeStep_ = std::max( diffTimeStep, maxDiffTimeStep_ );

			return std::max( wave, diffTimeStep );

		}
		/*@LST0S@*/
		/**
		 * @brief analytical flux function$
		 */
		template <class ArgumentTuple, class JacobianTuple >
		void analyticalFlux( const EntityType& en,
												 const double time, 
												 const DomainType& x,
												 const ArgumentTuple& u, 
												 const JacobianTuple& jac, 
												 JacobianRangeType& f ) const
		{
				JacobianRangeType diffmatrix;
				model_.diffusion(en, time, x, u[ uVar ],u[sigmaVar], diffmatrix);
				// ldg case 
				f = diffmatrix;
		}
	protected:
		mutable DiffusionFluxType& diffFlux_;
		const double penalty_;
    const double switch_; 
    const double cflDiffinv_;
	};                                              /*@LST0E@*/

}
#include "thetadiscretemodel.hh"
#endif

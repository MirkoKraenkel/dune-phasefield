#ifndef DUNE_FEM_DG_LDGFLUXTHETA_HH
#define DUNE_FEM_DG_LDGFLUXTHETA_HH

// local includes
#include <dune/fem-dg/operator/fluxes/diffusionflux.hh>

namespace Dune {


  // LDGDiffusionFlux
  //-----------------

  template <class DiscreteFunctionSpaceImp, 
            class Model>
  class LDGDiffusionFlux :
    public DGDiffusionFluxBase< DiscreteFunctionSpaceImp, Model>
  {
    typedef DGDiffusionFluxBase< DiscreteFunctionSpaceImp, Model> BaseType;
  public:
    typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;

    enum { dimDomain = DiscreteFunctionSpaceType :: dimDomain };
    enum { dimRange  = DiscreteFunctionSpaceType :: dimRange };

    typedef typename DiscreteFunctionSpaceType :: DomainType           DomainType;
    typedef typename DiscreteFunctionSpaceType :: RangeFieldType       RangeFieldType;
    typedef typename DiscreteFunctionSpaceType :: DomainFieldType      DomainFieldType;
    typedef typename DiscreteFunctionSpaceType :: RangeType            RangeType;
    typedef typename DiscreteFunctionSpaceType :: JacobianRangeType    JacobianRangeType;

    typedef FieldVector< DomainFieldType, dimDomain-1 > FaceDomainType;

    typedef typename DiscreteFunctionSpaceType :: GridPartType         GridPartType;
    typedef typename GridPartType :: IntersectionIteratorType          IntersectionIterator;
    typedef typename IntersectionIterator :: Intersection              Intersection;
    typedef typename GridPartType :: GridType                          GridType;
    typedef typename DiscreteFunctionSpaceType :: EntityType           EntityType;
    enum { dimGradRange = dimDomain * dimRange };
    enum { polOrd = DiscreteFunctionSpaceType :: polynomialOrder };

    // type of gradient space 
    typedef typename BaseType :: DiscreteGradientSpaceType  DiscreteGradientSpaceType;

    typedef typename DiscreteGradientSpaceType :: RangeType GradientRangeType;
    typedef typename DiscreteGradientSpaceType :: JacobianRangeType GradientJacobianType;

    // jacobians of the functions do not have to be evaluated for this flux 
    enum { evaluateJacobian = false };

  private:
    // no copying 
    LDGDiffusionFlux(const LDGDiffusionFlux& other);
  protected:
    using BaseType :: determineDirection;
    using BaseType :: model_;
    using BaseType :: cflDiffinv_;
    using BaseType :: dimensionFactor_;
    using BaseType :: nonconformingFactor_;
    using BaseType :: numericalFlux ;


  public:
    /**
     * @brief constructor
     */
    LDGDiffusionFlux(GridPartType& gridPart,
                     const Model& model ) :
      BaseType( model, true ),
      penalty_(Parameter::getValue<double>("dgdiffusionflux.penalty")),
      // Set CFL number for penalty term (compare diffusion in first pass)
      penaltyTerm_( std::abs(  penalty_ ) > 0 )
    {
      if( Parameter :: verbose () ) 
				{
					std::cout << "LDGDiffusionFlux: penalty = " << penalty_ << std::endl;
				}
    }

    //! returns true if lifting has to be calculated 
    const bool hasLifting () const { return false; }

    void diffusionFluxPenalty( std::ostream& out ) const
    {
      out <<penalty_;
    }

  public:
    /**
     * @brief flux function on interfaces between cells
     *
     * @param intersection intersection
     * @param time current time given by TimeProvider
     * @param x coordinate of required evaluation local to \c intersection
     * @param uLeft DOF evaluation on this side of \c intersection
     * @param uRight DOF evaluation on the other side of \c intersection
     * @param gLeft result for this side of \c intersection
     * @param gRight result for the other side of \c intersection
     * @return wave speed estimate (multiplied with the integration element of the intersection).
     *         To estimate the time step |T|/wave is used
     */
    template <class QuadratureImp>
    double gradientNumericalFlux(
																 const Intersection& intersection,
																 const EntityType& inside,
																 const EntityType& outside,
																 const double time,
																 const QuadratureImp& faceQuadInner,
																 const QuadratureImp& faceQuadOuter,
																 const int quadPoint, 
																 const RangeType& uLeft,
																 const RangeType& uRight,
																 GradientRangeType& gLeft,
																 GradientRangeType& gRight,
																 GradientJacobianType& gDiffLeft,
																 GradientJacobianType& gDiffRight) const
    {
      const FaceDomainType& x = faceQuadInner.localPoint( quadPoint );
      const DomainType normal = intersection.integrationOuterNormal( x );

      // determine side 
      //      const bool useInterior = determineDirection( normal );
      const bool useInterior = determineDirection(false,0.,0.,intersection );
	
      GradientJacobianType diffmatrix; 

      if( useInterior ) 
				{
					model_.jacobian(inside,         /* inside entity */
													time,           /* for time dependent diffusion */
													faceQuadInner.point( quadPoint ),      /* inside point on intersection */
													uLeft,          /* { u_(x^-) } */
													diffmatrix      /* return diffusion tensor */
                          );
				}
      else 
				{
					model_.jacobian(outside,       /* outside entity */
													time,          /* for time dependent diffusion */
													faceQuadOuter.point( quadPoint ),  /* outside point on intersection */
													uRight,        /* { u_(x^+) } */
													diffmatrix     /* return diffusion tensor */
													);
				}

      // mutliply with normal 
      diffmatrix.mv(normal, gLeft);

      // copy flux 
      gRight = gLeft;

      gDiffLeft = 0;
      gDiffRight = 0;

      // time step is set in 2nd pass 
      return 0.0;
    }

    /*
     * @brief numerical flux for u given as u_h
		 */
    template <class QuadratureImp> 
    double gradientBoundaryFlux(const Intersection& intersection,
                                const EntityType& inside,
                                const double time, 
                                const QuadratureImp& faceQuadInner,
                                const int quadPoint,
                                const RangeType& uLeft,
                                const RangeType& uBnd,
                                GradientRangeType& gLeft,
                                GradientJacobianType& gDiffLeft) const 
    {
      const FaceDomainType& x = faceQuadInner.localPoint( quadPoint );
      const DomainType normal = intersection.integrationOuterNormal(x);
      const DomainType& xglInside  = faceQuadInner.point( quadPoint );
      GradientJacobianType diffmatrix;

      // determine side 
      //const bool useInterior = determineDirection( normal );

      // get apropriate value 
      //const RangeType&  uVal = ( useInterior ) ? uLeft : uBnd ;

      // for the numerical diffusion flux \tilde u
      // one uses \tilde u = g_D where g_D is Dirichlet boundary data
      model_.jacobian( inside, 
                       time,
                       xglInside,
                       uBnd,
                       diffmatrix
											 );

      // apply normal 
      diffmatrix.mv(normal, gLeft); 

      // gDiffLeft is not needed for dual formalation (LDG)
      gDiffLeft = 0;

      // time step is set in 2nd pass 
      return 0.0;
    }


    /**
     * @brief flux function on interfaces between cells
     *
     * @param intersection intersection
     * @param time current time given by TimeProvider
     * @param x coordinate of required evaluation local to \c intersection
     * @param uLeft DOF evaluation on this side of \c intersection
     * @param uRight DOF evaluation on the other side of \c intersection
     * @param gLeft result for this side of \c intersection
     * @param gRight result for the other side of \c intersection
     * @return wave speed estimate (multiplied with the integration element of the intersection).
     *         To estimate the time step |T|/wave is used
     */
    template <class DiscreteModelImp, 
              class QuadratureImp>
    inline double numericalFlux(const Intersection& intersection,
																const DiscreteModelImp& discreteModel,
																const double time,
																const QuadratureImp& faceQuadInner,
																const QuadratureImp& faceQuadOuter,
																const int quadPoint, 
																const RangeType& uLeft,
																const RangeType& uRight,
																const GradientRangeType& sigmaLeft,
																const GradientRangeType& sigmaRight,
																const GradientRangeType& thetaLeft,
																const GradientRangeType& thetaRight,
																RangeType& gLeft,
																RangeType& gRight,
																JacobianRangeType& gDiffLeft, // not used here (only for primal passes)
																JacobianRangeType& gDiffRight ) const
    {
      FieldMatrixConverter< GradientRangeType, JacobianRangeType> jacLeft( sigmaLeft );
      FieldMatrixConverter< GradientRangeType, JacobianRangeType> jacRight( sigmaRight );
      FieldMatrixConverter< GradientRangeType, JacobianRangeType> tensLeft( thetaLeft );
      FieldMatrixConverter< GradientRangeType, JacobianRangeType> tensRight( thetaRight );
   
      return numericalFlux( intersection, discreteModel, time,
                            faceQuadInner, faceQuadOuter, quadPoint,
                            uLeft, uRight, jacLeft, jacRight,
														tensLeft,tensRight,
                            gLeft, gRight, gDiffLeft, gDiffRight );
    }


    template< class DiscreteModelImp, 
              class QuadratureImp, class JacobianRangeTypeImp >
    double numericalFlux(const Intersection& intersection,
                         const DiscreteModelImp& discreteModel,
                         const double time,
                         const QuadratureImp& faceQuadInner,
                         const QuadratureImp& faceQuadOuter,
                         const int quadPoint, 
                         const RangeType& uLeft,
                         const RangeType& uRight,
                         const JacobianRangeTypeImp& sigmaLeft,
                         const JacobianRangeTypeImp& sigmaRight,
												 const JacobianRangeTypeImp& thetaLeft,
                         const JacobianRangeTypeImp& thetaRight,
												 RangeType& gLeft,
                         RangeType& gRight,
                         JacobianRangeType& gDiffLeft, // not used here (only for primal passes)
                         JacobianRangeType& gDiffRight ) const
    {
      const FaceDomainType& x = faceQuadInner.localPoint( quadPoint );
      const DomainType normal = intersection.integrationOuterNormal(x);

      const EntityType& inside  = discreteModel.inside();
      const EntityType& outside = discreteModel.outside();

			/**********************************
			 * Diffusion sigma Flux (Pass 2)  *
			 **********************************/
      JacobianRangeType diffmatrix;
      JacobianRangeType tension;
      tension=thetaLeft;


      tension+=thetaRight;
      
      tension*=-0.5;
       

      // determine side (needs to be opposite of above)
      const bool useExterior = ! determineDirection(false,0.,0.,intersection );

      if( useExterior ) 
				{
					model_.diffusion( inside, time, 
														faceQuadInner.point( quadPoint ),
														uLeft, sigmaLeft, diffmatrix);
	

				}
      else 
				{
					model_.diffusion( outside, time, 
														faceQuadOuter.point( quadPoint ),
														uRight, sigmaRight, diffmatrix);
				}

      

      // apply normal 
      diffmatrix.mv(normal, gLeft);
      tension.mv(normal,gLeft);

      //////////////////////////////////////////////////////////
      //
      //  --Time step calculation 
      //
      //////////////////////////////////////////////////////////
      const double faceLengthSqr = normal.two_norm2();

      const double faceVolumeEstimate = dimensionFactor_ *
        ( intersection.conforming() ) ? faceLengthSqr : ( nonconformingFactor_ * faceLengthSqr );

      const double diffTimeLeft =
        model_.diffusionTimeStep( intersection,
																	discreteModel.enVolume(),
																	faceVolumeEstimate,
																	time, x, uLeft );

      const double diffTimeRight =
        model_.diffusionTimeStep( intersection, 
																	discreteModel.nbVolume(),
																	faceVolumeEstimate,
																	time, x, uRight );

      //////////////////////////////////////////////////////////
      //
      //  --Penalty Term
      //

      //////////////////////////////////////////////////////////

      // take minimum to proceed 
      const double diffTimeStep = std::max( diffTimeLeft, diffTimeRight );
      if( penaltyTerm_ )
				{
					// add penalty term ( enVolume() is available since we derive from
					//    DiscreteModelDefaultWithInsideOutside)
					const double factor = penalty_ * diffTimeStep ;

					RangeType jump( uLeft );
					jump -= uRight;
					gLeft.axpy(factor, jump);
				}

      gRight = gLeft ;

      // gDiffLeft should be 0 in case of LDG
      gDiffLeft = 0;
      gDiffRight = 0;

      // timestep restict to diffusion timestep
      // WARNING: reconsider this
      return diffTimeStep * cflDiffinv_;
    }


    /**
     * @brief same as numericalFlux() but for fluxes over boundary interfaces
     */
    template <class DiscreteModelImp, 
              class QuadratureImp> 
    inline double boundaryFlux(const Intersection& intersection,
															 const DiscreteModelImp& discreteModel,  
															 const double time, 
															 const QuadratureImp& faceQuadInner,
															 const int quadPoint,
															 const RangeType& uLeft,
															 const RangeType& uRight,
															 const GradientRangeType& sigmaLeft,
															 const GradientRangeType& thetaLeft,
															 RangeType& gLeft,
															 JacobianRangeType& gDiffLeft) const   /*@LST0E@*/
    {
      FieldMatrixConverter< GradientRangeType, JacobianRangeType> jacLeft( sigmaLeft );
      FieldMatrixConverter< GradientRangeType, JacobianRangeType> tensLeft(thetaLeft );
      return boundaryFlux( intersection, discreteModel, time,
                           faceQuadInner, quadPoint, uLeft, uRight,
                           jacLeft, tensLeft,gLeft, gDiffLeft );
    }

    /**
     * @brief same as numericalFlux() but for fluxes over boundary interfaces
     */
    template< class DiscreteModelImp, 
              class QuadratureImp, class JacobianRangeTypeImp > 
    double boundaryFlux(const Intersection& intersection,
                        const DiscreteModelImp& discreteModel,  
                        const double time, 
                        const QuadratureImp& faceQuadInner,
                        const int quadPoint,
                        const RangeType& uLeft,
                        const RangeType& uRight,
                        const JacobianRangeTypeImp& sigmaLeft,
											  const JacobianRangeTypeImp& thetaLeft,
												RangeType& gLeft,
                        JacobianRangeType& gDiffLeft) const   /*@LST0E@*/
    {
      // get local point 
      const FaceDomainType& x = faceQuadInner.localPoint( quadPoint );
      const DomainType normal = intersection.integrationOuterNormal(x);

      const EntityType& inside  = discreteModel.inside();

      /****************************/
      /* Diffusion (Pass 2)       */
      /****************************/
      JacobianRangeType diffmatrix;
     

      diffmatrix.mv(normal, gLeft);
      

      //////////////////////////////////////////////////////////
      //
      //  --Time step calculation 
      //
      //////////////////////////////////////////////////////////
      const double faceVolumeEstimate = normal.two_norm2();

      const double diffTimeStep =
        model_.diffusionTimeStep( intersection,
																	discreteModel.enVolume(),
																	faceVolumeEstimate,
																	time, x, uLeft );

      //////////////////////////////////////////////////////////
      //
      //  --Penalty Term
      //
      //////////////////////////////////////////////////////////
      if( penaltyTerm_ ) 
				{
					// add penalty term
					const double factor = penalty_ * diffTimeStep;

					RangeType jump( uLeft );
					jump -= uRight;
					gLeft.axpy(factor, jump);
				}

      // gDiffLeft should be 0 in case of LDG
      gDiffLeft = 0;

      return diffTimeStep * cflDiffinv_; 
    }

  protected:
    const double penalty_;
    const bool penaltyTerm_;

  }; // end LDGDiffusionFlux                        

} // end namespace 
#endif

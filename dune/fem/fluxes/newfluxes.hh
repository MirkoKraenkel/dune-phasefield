#ifndef DUNE_FEM_DG_EULERFLUXES_HH
#define DUNE_FEM_DG_EULERFLUXES_HH
  
// system includes
#include <string>
#include <cmath>

#include <dune/common/fvector.hh>





////////////////////////////////////////////////////////
//
// Dune interface for Euler fluxes 
//
////////////////////////////////////////////////////////


// Wellbalanced Flux
//--------

template< class Model >
class WBFlux
{
public:
  typedef Model                                       ModelType;
  enum { dimDomain = Model::dimDomain };
  enum { dimRange = Model::dimRange };
  typedef typename Model::Traits                      Traits;
  typedef typename Traits::GridType                   GridType;
  typedef typename GridType::ctype                    ctype;
  typedef typename Traits::EntityType                 EntityType;
  typedef typename Traits::EntityPointerType          EntityPointerType;

  typedef typename Traits::DomainType                 DomainType;
  typedef typename Traits::FaceDomainType             FaceDomainType;
  typedef typename Traits::RangeType                  RangeType;
  typedef typename Traits::GradientRangeType          GradientType;
 
	typedef typename Traits::FluxRangeType              FluxRangeType;
	typedef typename Traits::ThetaRangeType             ThetaRangeType; 
  
	WBFlux( const Model& mod )
    : model_(mod),
      visc_(Dune::Fem::Parameter::getValue<double>("phasefield.addvisc",1)),
			numViscMu_(Dune::Fem::Parameter::getValue<double>("phasefield.nuvisc",0.))
  {
  }

  static std::string name () { return "nonconservativtranssport"; }

  const Model& model() const { return model_; }

  // Return value: maximum wavespeed*length of integrationOuterNormal
  // gLeft,gRight are fluxed * length of integrationOuterNormal
  template< class Intersection, class QuadratureImp >
  inline double numericalFlux(const Intersection& intersection,
                              const EntityType& inside,
                              const EntityType& outside,
                              const double time, 
                              const QuadratureImp& faceQuadInner,
                              const QuadratureImp& faceQuadOuter,
                              const int quadPoint,
                              const RangeType& uLeft, 
                              const RangeType& uRight,
								              const ThetaRangeType& thetaLeft,
								              const ThetaRangeType& thetaRight,
								              RangeType& gLeft,
                              RangeType& gRight) const 
  {
    const FaceDomainType& x = faceQuadInner.localPoint( quadPoint );
    DomainType normal = intersection.integrationOuterNormal(x);  
    const double len = normal.two_norm();
    normal *= 1./len;
        
    double rhoLeft  = uLeft[0];
    double rhoRight = uRight[0];
    double phiLeft  = uLeft[dimDomain+1];
    double phiRight = uRight[dimDomain+1];

    double vLeft[dimDomain],vRight[dimDomain];

    for(int i=0; i<dimDomain; i++)
    {
      vLeft[i]  = uLeft[1+i]/rhoLeft;
      vRight[i] = uRight[1+i]/rhoRight;
    }

    

    RangeType visc;
		ThetaRangeType newvisc,thetaFluxLeft,thetaFluxRight;
   	FluxRangeType anaflux;
    model_.advection( inside, time, faceQuadInner.point( quadPoint ),
                      uLeft, anaflux );
     // set gLeft 
    anaflux.mv( normal, gLeft );
      
    
    model_.advection( outside, time, faceQuadOuter.point( quadPoint ),
                      uRight, anaflux );
    //add F(uleft) 
    anaflux.umv( normal, gLeft );
#if  0
    model_.thetaSource( inside, time, faceQuadInner.point( quadPoint ),
                      uLeft, thetaFluxLeft );

    model_.thetaSource( outside, time, faceQuadOuter.point( quadPoint ),
                 uRight,thetaFluxRight );

#endif

    double maxspeedl, maxspeedr, maxspeed;
    double viscparal, viscparar, viscpara;
    
    const DomainType xGlobal = intersection.geometry().global(x);
    
    model_.maxSpeed( normal, time, xGlobal, 
                     uLeft, viscparal, maxspeedl );
    model_.maxSpeed( normal, time, xGlobal,
                     uRight, viscparar, maxspeedr );


    maxspeed = (maxspeedl > maxspeedr) ? maxspeedl : maxspeedr;
    viscpara = (viscparal > viscparar) ? viscparal : viscparar;
    viscpara*=visc_;
    visc = uRight;
    visc -= uLeft;

    visc *= viscpara;
    for(int i=1; i<dimDomain+1;i++)
 		{
     		gLeft[i] -= visc[i];
    }

    
     // \delta\mu  consider sign!!!!!!!!
    newvisc=thetaRight[0];
    newvisc-=thetaLeft[0];
    newvisc*=viscpara*numViscMu_; 
     
   gLeft[0]-=newvisc[0];
   gLeft *= 0.5*len; 
   gRight = gLeft;
 
   RangeType nonConLeft(0.),nonConRight(0.);
  
   nonConFlux( normal,
               len,        
               rhoLeft,
               rhoRight, 
               vLeft,
               vRight,
               thetaLeft, 
               thetaRight,
               phiLeft,
               phiRight,
               nonConLeft,
               nonConRight);
    
   gLeft  -= nonConLeft;
   gRight += nonConRight;
   
   return maxspeed * len;
  }


  inline void nonConFlux( const DomainType& normal,
                          const double length,  
                          const double rhoLeft,
                          const double rhoRight,
                          const double* vLeft,
                          const double* vRight,
                          const ThetaRangeType& thetaLeft,
                          const ThetaRangeType& thetaRight,
                          const double phiLeft,
                          const double phiRight,
                          RangeType& nonConLeft,
                          RangeType& nonConRight) const
  {
      double tauLeft,tauRight;
      tauLeft  = thetaLeft[1];
      tauRight = thetaRight[1];
      
      //[[\mu]]
      double jumpMu=thetaLeft[0]-thetaRight[0];
     
      //[\phi]
      double jumpPhi=phiLeft-phiRight;
      
      double vLeftNormal{0.},vRightNormal{0.};

      nonConLeft=nonConRight=0.;
      for(int i=0;i<dimDomain;i++)
      {  
        nonConLeft[i+1]=normal[i];
        vLeftNormal+=normal[i]*vLeft[i];
        vRightNormal+=normal[i]*vRight[i];
      }     
    
      nonConRight=nonConLeft;
      
      nonConLeft *=(jumpMu*rhoLeft -jumpPhi*tauLeft);
      nonConRight*=(jumpMu*rhoRight-jumpPhi*tauRight);
      
      nonConLeft[dimDomain+1]=vLeftNormal*jumpPhi; 
      nonConRight[dimDomain+1]=vRightNormal*jumpPhi;   
      nonConLeft*=0.5*length;
      nonConRight*=0.5*length;
  }
                           
                           



 protected:
  const Model& model_;
  const double visc_;
	const double numViscMu_;
};
#endif // file declaration

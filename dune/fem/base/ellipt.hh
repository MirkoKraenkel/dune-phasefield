#ifndef ELLIPTIC_ALGORITHM_HH
#define ELLIPTIC_ALGORITHM_HH
#include <config.h>

#ifdef LOCALDEBUG
static double sum_ = 0.;
static double sum2_ = 0.;
static double localMaxRatio_ = 0.;
static double localMinRatio_ = 1e+100;
static double maxRatioOfSums = 0.;
static double minRatioOfSums = 1e+100;
#endif

#ifndef NDEBUG 
// enable fvector and fmatrix checking
#define DUNE_ISTL_WITH_CHECKING
#endif

// include std libs
#include <iostream>
#include <string>

// dune-fem includes
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/operator/projection/l2projection.hh>
#include <dune/fem/gridpart/gridpart.hh>
#include <dune/fem/solver/odesolver.hh>

// dune-fem-dg includes
#include <dune/fem-dg/operator/dg/dgoperatorchoice.hh>

// include local header files
#include "../base/base.hh"
#include "../base/baseevolution.hh"


using namespace Dune;                                        


template <class GridImp,
          class ProblemTraits, 
          int polOrd>             
struct ElliptTraits 
{
  enum { polynomialOrder = polOrd };

  // type of Grid
  typedef GridImp                                   GridType;

  // Choose a suitable GridView
  typedef DGAdaptiveLeafGridPart< GridType >       GridPartType;

  // problem dependent types 
  typedef typename ProblemTraits :: template Traits< GridPartType > :: InitialDataType  InitialDataType;
  typedef typename ProblemTraits :: template Traits< GridPartType > :: ModelType        ModelType;
  typedef typename ProblemTraits :: template Traits< GridPartType > :: FluxType         FluxType;
  static const Dune :: DGDiffusionFluxIdentifier DiffusionFluxId = 
    ProblemTraits :: template Traits< GridPartType > :: PrimalDiffusionFluxId ;

  // The DG Operator (using 2 Passes)
  // The first operator is sum of the other two
  // The other two are needed for semi-implicit time discretization
  typedef DGAdvectionDiffusionOperator< ModelType, FluxType,
                                        DiffusionFluxId,  
                                        polynomialOrder >            DgType; /*@LST1E@*/
  typedef DGAdvectionOperator< ModelType, FluxType,
                               DiffusionFluxId,  
                               polynomialOrder >                     DgAdvectionType; /*@LST1E@*/
  typedef DGDiffusionOperator< ModelType, FluxType,
                               DiffusionFluxId,  
                               polynomialOrder >                     DgDiffusionType; /*@LST1E@*/

  // The discrete function for the unknown solution is defined in the DgOperator
  typedef typename DgType :: DestinationType                         DiscreteFunctionType;
  // ... as well as the Space type
  typedef typename DgType :: SpaceType                               DiscreteSpaceType;

  // The ODE Solvers                                                         /*@LST1S@*/
  typedef DuneODE :: OdeSolverInterface< DiscreteFunctionType >      OdeSolverInterfaceType;
  typedef DuneODE :: ImplicitOdeSolver< DiscreteFunctionType >       ImplicitOdeSolverType; /*@\label{dg:}@*/
  typedef DuneODE :: ExplicitOdeSolver< DiscreteFunctionType >       ExplicitOdeSolverType; /*@\label{}@*/
  typedef DuneODE :: SemiImplicitOdeSolver< DiscreteFunctionType >   SemiImplicitOdeSolverType; /*@\label{}@*/


  // type of restriction/prolongation projection for adaptive simulations 
  typedef RestrictProlongDefault< DiscreteFunctionType > RestrictionProlongationType;
};


template <class GridImp,
          class ProblemTraits, 
          int polynomialOrder>             
struct EllipticAlgorithm
{
  // my traits class 
  typedef ElliptTraits< GridImp, ProblemTraits, polynomialOrder> Traits ;

  // type of Grid
  typedef typename Traits :: GridType                 GridType;

  // Choose a suitable GridView
  typedef typename Traits :: GridPartType             GridPartType;

  // initial data type 
  typedef typename Traits :: InitialDataType          InitialDataType;

  // An analytical version of our model
  typedef typename Traits :: ModelType                 ModelType;

  // The flux for the discretization of advection terms
  typedef typename Traits :: FluxType                  FluxType;

  // The DG space operator
  // The first operator is sum of the other two
  // The other two are needed for semi-implicit time discretization
  typedef typename Traits :: DgType                    DgType;
  typedef typename Traits :: DgAdvectionType           DgAdvectionType;
  typedef typename Traits :: DgDiffusionType           DgDiffusionType;

  // The discrete function for the unknown solution is defined in the DgOperator
  typedef typename Traits :: DiscreteFunctionType      DiscreteFunctionType;

  // ... as well as the Space type
  typedef typename Traits :: DiscreteSpaceType         DiscreteSpaceType;

  typedef tuple< const DiscreteFunctionType* >           IOTupleType;
  typedef Dune::DataWriter<GridType,IOTupleType>            DataWriterType;

  EllipticAlgorithm(GridType& grid) :
    grid_( grid ),
    gridPart_( grid_ ),
    space_( gridPart_ ),
    solution_("solution", space_ ),
    problem_( ProblemTraits::problem() ),
    model_( new ModelType( problem() ) ),
    convectionFlux_( *model_ ),
    dgOperator_(grid_, convectionFlux_),
    dgAdvectionOperator_(grid_, convectionFlux_),
    dgDiffusionOperator_(grid_, convectionFlux_),
    eocId_( FemEoc::addEntry(std::string("$L^2$-error")) )
  {
    std::string filename( Parameter::commonOutputPath() );
    filename += "/run.gnu";
    runFile_.open( filename.c_str() );
    if( ! runFile_.is_open() ) 
    {
      std::cerr << filename << "runfile not open" << std::endl;
      abort();
    }
    runFile_ << "# h | elements | CPU time | iter | l_min | l_max | cond  | L2 error" << std::endl;
    // initializeFemEoc( description() );
  }                                                                        /*@LST1E@*/

  //! return reference to discrete space 
  DiscreteSpaceType & space() { return space_; }                    /*@LST0E@*/

  //! returns data prefix for EOC loops ( default is loop )
  virtual std::string dataPrefix() const 
  {
    return problem_->dataPrefix();
  }

  // gather information from the space operator, the time integratior
  // and the problem to output before each table in tex file
  std::string description() const
  {
    std::string latexInfo;

    latexInfo = dgAdvectionOperator_.description()
                + dgDiffusionOperator_.description();

    std::stringstream odeInfo; 

    latexInfo += odeInfo.str() 
                  + "\n"
                  + problem_->description()
                  + "\n\n";

    return latexInfo;
  }


  //! default time loop implementation, overload for changes 
  void operator()(double& averagedt,
                  double& mindt, double& maxdt, size_t& counter,
                  int& total_newton_iterations, int& total_ils_iterations,
                  int& max_newton_iterations, int& max_ils_iterations)
  {
    numbers_.resize( 0 );

    // calculate grid width
    const double h = GridWidth::calcGridWidth(gridPart_);
    numbers_.push_back( h );

    const double size = grid_.size(0);
    numbers_.push_back( size );
    
    //assert( solution_.space().size() > 0 );
    // solveOperator( dgDiffusionOperator_ , solution_ , numbers_ );
    solveOperator( dgOperator_ , solution_ , numbers_ );

    //DGL2ProjectionImpl :: project( problem(), solution_ );
    return ;
  }


  //! finalize computation by calculating errors and EOCs 
  void finalize( const int eocloop )
  {
    //---- Adapter for exact solution ------------------------------------------
    typedef Dune::GridFunctionAdapter< InitialDataType, GridPartType >
        GridExactSolutionType;

    // create grid function adapter 
    GridExactSolutionType ugrid( "exact solution", problem(), gridPart_,
                                 DiscreteSpaceType :: polynomialOrder + 1 );

    // create L2 - Norm 
    Dune::L2Norm< GridPartType > l2norm( gridPart_ );

    // calculate L2 - Norm 
    const double l2error = l2norm.distance( ugrid, solution_ );

    numbers_.push_back( l2error );

    for(size_t i=0; i<numbers_.size(); ++i)
      runFile_ << numbers_[ i ] << " ";

    runFile_ << std::endl;

    // store values 
    std::vector<double> errors;
    errors.push_back( l2error );

    // submit error to the FEM EOC calculator 
    Dune::FemEoc :: setErrors(eocId_, errors);
  }

  const InitialDataType& problem() const { assert( problem_ ); return *problem_; }
  const ModelType& model() const { assert( model_ ); return *model_; }
  const DiscreteFunctionType& solution() const { return solution_; }

protected:
  GridType& grid_;
  GridPartType gridPart_;      // reference to grid part, i.e. the leaf grid 
  DiscreteSpaceType space_;    // the discrete function space
  // the solution 
  DiscreteFunctionType   solution_;
  // InitialDataType is a Dune::Operator that evaluates to $u_0$ and also has a
  // method that gives you the exact solution.
  const InitialDataType*  problem_;
  ModelType*              model_;
  // Initial flux for advection discretization (UpwindFlux)
  FluxType                convectionFlux_;
  DgType                  dgOperator_;
  DgAdvectionType         dgAdvectionOperator_;
  DgDiffusionType         dgDiffusionOperator_;
  const unsigned int      eocId_;
  std::vector<double> numbers_;
  std::ofstream runFile_;
};
#endif // FEMHOWTO_STEPPER_HH

// include host specific macros set by configure script   /*@LST0S@*/
#include <config.h>

// include std libs
#include <iostream>
#include <string>

// Dune includes
#include <dune/fem/misc/threads/threadmanager.hh>
#include <dune/fem/operator/projection/l2projection.hh>
#include <dune/fem/gridpart/common/gridpart.hh>
// local includes
// definition of simulate() 
#if AC
#include "allencahn.hh"
#else
#include "phasefield.hh"
#endif




/**
 * @return 0 we don't program bugs. :)
 */
int main(int argc, char ** argv, char ** envp) {        

  /* Initialize MPI (always do this even if you are not using MPI) */
  Dune::Fem::MPIManager :: initialize( argc, argv );
  try {

        // *** Initialization
        Dune::Fem::Parameter::append(argc,argv);

        if (argc == 2)
        {
          Dune::Fem::Parameter::append(argv[1]);
        }
        else
        {
          Dune::Fem::Parameter::append("parameter");
        }

    // get number of desired threads (default is 1)
    int numThreads = Dune::Fem::Parameter::getValue< int >("fem.parallel.numberofthreads", 1);
    Dune :: Fem :: ThreadManager :: setMaxNumberThreads( numThreads );


    simulation :: simulate();  
	  
    // write parameters used 
    Dune::Fem::Parameter::write("parameter.log");
  }
  catch (Dune::Exception &e) 
  {                           
    std::cerr << e << std::endl;
    return 1;
  } 
  catch (...) 
  {
    std::cerr << "Generic exception!" << std::endl;
    return 2;
  }                                                     

  return 0;
}


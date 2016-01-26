
#include <iomanip>
#include <chrono>

#include "comm/global.hh"
#include "EigenSystem.hh"
#include "EigenSystemFactory.hh"
#include "EigenSolverFactory.hh"
#include "AleaTypedefs.hh"
#include "DeviceTraits.hh"

#include "Kokkos_Core.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "MatrixMarket_Tpetra.hpp"

//#include "cuda.h"
//#include "cuda_runtime_api.h"

using namespace alea;
using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::milliseconds;

int main( int argc, char *argv[] )
{
    profugus::initialize(argc,argv);
    Kokkos::initialize(argc,argv);

    high_resolution_clock::time_point time_start =
        high_resolution_clock::now();

    // Read ParameterList from file
    Teuchos::RCP<Teuchos::ParameterList> pl;
    VALIDATE( argc > 1, "USAGE: xalea input_file.xml");
    pl = Teuchos::getParametersFromXmlFile(argv[1]);
    CHECK( pl != Teuchos::null );

    // Initialize Kokkos device
    //DeviceTraits<DEVICE>::initialize(pl);

    Teuchos::RCP<EigenSystem> system =
        alea::EigenSystemFactory::buildEigenSystem(pl);
    Teuchos::RCP<const MATRIX> myA = system->getMatrix();
    Teuchos::RCP<const MV> initial_guess = system->getInitialGuess();

    Teuchos::RCP<MV> x( new MV(myA->getDomainMap(),1) );

    high_resolution_clock::time_point time_end = high_resolution_clock::now();
    if( profugus::node() == 0 )
    {
        std::cout << "Matrix loading took "
            << duration_cast<milliseconds>(time_end-time_start).count()
            << " milliseconds" << std::endl;
    }

    time_start = high_resolution_clock::now();

    std::string solver_type = pl->get("solver_type","mcsa");
    Teuchos::RCP<alea::AleaSolver> solver =
        alea::EigenSolverFactory::buildSolver(solver_type,myA,pl);

    time_end = high_resolution_clock::now();
    if( profugus::node()== 0 )
    {
        std::cout << "Solver setup took "
            << duration_cast<milliseconds>(time_end-time_start).count()
            << " milliseconds" << std::endl;
    }

    time_start = high_resolution_clock::now();
    solver->apply(*initial_guess,*x);
    time_end = high_resolution_clock::now();
    if( profugus::node() == 0 )
    {
        std::cout << "Solve took "
            << duration_cast<milliseconds>(time_end-time_start).count()
            << " milliseconds" << std::endl;
    }


    FILE *pFile;
    pFile = fopen("eigenvector.mtx", "w");

    size_t N = myA -> getGlobalNumRows();
    
    Teuchos::ArrayRCP<SCALAR> xp = x->getDataNonConst(0);

    fprintf(pFile, "%%%MatrixMarket matrix array real general \n");
    fprintf(pFile, "%d %d \n", N,1);
    for (unsigned int index = 0; index < N; ++index)
         fprintf (pFile, "%14.15f \n", xp[index]);
    
    fclose(pFile);

/*  int nDevices;

  cudaGetDeviceCount(&nDevices);
  printf("The number of GPU devices is: %d\n", nDevices);
  for (int i = 0; i < nDevices; i++) {
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, i);
    printf("Device Number: %d\n", i);
    printf("  Device name: %s\n", prop.name);
    printf("  Memory Clock Rate (KHz): %d\n",
           prop.memoryClockRate);
    printf("  Memory Bus Width (bits): %d\n",
           prop.memoryBusWidth);
    printf("  Peak Memory Bandwidth (GB/s): %f\n\n",
           2.0*prop.memoryClockRate*(prop.memoryBusWidth/8)/1.0e6);
    printf("  globalL1CacheSupported:  %d\n", prop.globalL1CacheSupported);
    printf("  localL1CacheSupported:  %d\n", prop.localL1CacheSupported);
  }
*/

    Kokkos::finalize();
    profugus::finalize();
    return 0;
}



#include <iomanip>
#include <chrono>

#include "comm/global.hh"
#include "LinearSystem.hh"
#include "LinearSystemFactory.hh"
#include "LinearSolverFactory.hh"
#include "AleaTypedefs.hh"
#include "DeviceTraits.hh"

#include "Kokkos_Core.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "MatrixMarket_Tpetra.hpp"

#include "cuda.h"
#include "cuda_runtime_api.h"

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

    Teuchos::RCP<LinearSystem> system =
        alea::LinearSystemFactory::buildLinearSystem(pl);
    Teuchos::RCP<const MATRIX> myA = system->getMatrix();
    Teuchos::RCP<const MV> b = system->getRhs();

    Teuchos::RCP<MV> x( new MV(myA->getDomainMap(),1) );

    high_resolution_clock::time_point time_end = high_resolution_clock::now();
    if( profugus::node() == 0 )
    {
        std::cout << "Matrix loading took "
            << duration_cast<milliseconds>(time_end-time_start).count()
            << " milliseconds" << std::endl;
    }

    time_start = high_resolution_clock::now();

    unsigned int num_blocks   = solver_pl->get("num_blocks",10);
    VALIDATE( num_blocks >= 2, "Minimal number of partitions is 2" );
    
    unsigned int overlap      = solver_pl->get("overlap",0.0);
    VALIDATE( overlap>= 0.0 && overlap<=1.0, 
            "The percentage of overlapping must be a number between 0 and 1");
    
    std::string inner_solver = solver_pl->get("inner_solver","richardson");
    VALIDATE( inner_solver == "richardson" || inner_solver =="monte_carlo", 
             "The type of inner solver provided is not valid for Multi-Splitting" );
             
             
    //measure the size of the problem 
    unsigned int N = myA -> getNodeMaxNumRowEntries();
             
    //determine the number of rows that overlaps between adjacent subdomains         
    unsigned int overlapping = overlap * N;        
    
    //determine the number of rows for each subdomain
    unsigned int block_size = N / num_blocks;        
     
    Teuchos::ArrayRCP< ENDPOINTS > partitions(num_blocks);
       
    unsigned int p = 0;
    
    partitions[p][0]=0;
    partitions[p][1]=block_size - 1;
    p=1;
    
    while(p!=num_blocks - 1)
    {
    	partitions[p][0] = partitions[p-1][1] + 1 - overlapping;
    	partitions[p][1] = partitions[p-1][1] + 1 - overlapping + block_size;
    }
             
    partitions[p][0] = partitions[p-1][1] + 1 - overlapping;
    partitions[p][1] = N-1;

    Kokkos::finalize();
    profugus::finalize();
    return 0;
}


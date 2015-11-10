
#include <iomanip>
#include <chrono>

#include "comm/global.hh"
#include "LinearSystem.hh"
#include "LinearSystemFactory.hh"
#include "LinearSolverFactory.hh"
#include "LinearSystem_MultiSplitting.hh"
#include "MultiSplitting.hh"
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
    VALIDATE( argc > 1, "USAGE: xalea_multisplit input_file.xml");
    pl = Teuchos::getParametersFromXmlFile(argv[1]);
    CHECK( pl != Teuchos::null );

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

/*    LinearSystem_MultiSplitting L_MS(pl);
 *
 *        splitting split = L_MS.buildSplitting(pl, 0);*/

    MultiSplitting MS(pl);
    MS.solve(x);

/*    for(unsigned int i=0; i!=b_data.size(); ++i)
       std::cout<<b_data[i]<<std::endl;
*/

    Kokkos::finalize();
    profugus::finalize();
    return 0;
}


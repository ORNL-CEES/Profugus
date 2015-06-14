
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
#include "MatrixMarket_Tpetra.hpp"

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

    std::string solver_type = pl->get("solver_type","mcsa");
    Teuchos::RCP<alea::AleaSolver> solver =
        alea::LinearSolverFactory::buildSolver(solver_type,myA,pl);

    time_end = high_resolution_clock::now();
    if( profugus::node()== 0 )
    {
        std::cout << "Solver setup took "
            << duration_cast<milliseconds>(time_end-time_start).count()
            << " milliseconds" << std::endl;
    }

    time_start = high_resolution_clock::now();
    solver->apply(*b,*x);
    time_end = high_resolution_clock::now();
    if( profugus::node() == 0 )
    {
        std::cout << "Solve took "
            << duration_cast<milliseconds>(time_end-time_start).count()
            << " milliseconds" << std::endl;
    }

    // Compute final residual
    Teuchos::RCP<MV> r( new MV(myA->getDomainMap(),1) );
    myA->apply(*x,*r);
    r->update(1.0,*b,-1.0);
    Teuchos::ArrayRCP<SCALAR> res_norm(1), b_norm(1);
    r->norm2(res_norm());
    b->norm2(b_norm());
    if( profugus::node() == 0 )
    {
        std::cout << "Final relative residual norm: "
                  << res_norm[0]/b_norm[0] << std::endl;
    }

    //std::ofstream sol_file;
    //sol_file.open("solution.txt");
    std::string file_sol("solution.mtx");
    std::string file_rhs("rhs.mtx");
    Tpetra::MatrixMarket::Writer<MV>::writeDenseFile (file_sol,x);
    Tpetra::MatrixMarket::Writer<MV>::writeDenseFile (file_rhs,b);
    std::cout<<"I have printed somewhere"<<std::endl;
    //sol_file.close();

    // Finalize Kokkos device
    //DeviceTraits<DEVICE>::finalize();

    Kokkos::finalize();
    profugus::finalize();
    return 0;
}


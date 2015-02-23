
#include <iomanip>
#include <chrono>

#include "comm/global.hh"
#include "LinearSystem.hh"
#include "LinearSystemFactory.hh"
#include "LinearSolverFactory.hh"
#include "AleaTypedefs.hh"
#include "DeviceTraits.hh"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

using namespace alea;
using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::milliseconds;

int main( int argc, char *argv[] )
{
    profugus::initialize(argc,argv);

    high_resolution_clock::time_point time_start =
        high_resolution_clock::now();

    // Read ParameterList from file
    Teuchos::RCP<Teuchos::ParameterList> pl;
    TEUCHOS_ASSERT( argc > 1 );
    pl = Teuchos::getParametersFromXmlFile(argv[1]);
    TEUCHOS_ASSERT( pl != Teuchos::null );

    // Initialize Kokkos device
    DeviceTraits<DEVICE>::initialize(pl);

    Teuchos::RCP<LinearSystem> system =
        alea::LinearSystemFactory::buildLinearSystem(pl);
    Teuchos::RCP<const MATRIX> myA = system->getMatrix();
    Teuchos::RCP<const MV> b = system->getRhs();

    Teuchos::RCP<MV> x( new MV(myA->getDomainMap(),1) );

    high_resolution_clock::time_point time_end = high_resolution_clock::now();
    if( profugus::node() == 0 )
    {
        std::cout << "Matrix loading and eigenvalue calc took "
            << duration_cast<milliseconds>(time_end-time_start).count()
            << " milliseconds" << std::endl;
    }

    // Disable all output
    Teuchos::sublist(pl,"Polynomial")->set("verbosity","none");
    Teuchos::sublist(pl,"Belos")->set("verbosity","none");
    Teuchos::sublist(pl,"Synthetic Acceleration")->set("verbosity","none");
    Teuchos::sublist(pl,"Richardson")->set("verbosity","none");
    Teuchos::sublist(pl,"Chebyshev")->set("verbosity","none");
    Teuchos::sublist(pl,"Monte Carlo")->set("verbosity","none");

    Teuchos::RCP<Teuchos::ParameterList> multirun_pl =
        Teuchos::sublist(pl,"Multirun");

    int num_history_levels  = multirun_pl->get("Num History Levels",1);
    int start_num_histories = multirun_pl->get("Starting Histories",100);
    int history_multiplier  = multirun_pl->get("History Multiplier",2);
    std::vector<int> num_histories(num_history_levels);
    num_histories[0] = start_num_histories;
    for( int i=1; i<num_history_levels; ++i )
        num_histories[i] = num_histories[i-1] * history_multiplier;

    int num_length_levels = multirun_pl->get("Num Length Levels",1);
    int start_length      = multirun_pl->get("Starting Length",1);
    std::string length_type = multirun_pl->get("Length Increase Type",
                                               "Multiplicative");
    int length_factor = multirun_pl->get("Length Increase Factor",2);
    std::string output_format = multirun_pl->get("Output Format","standard");
    std::vector<int> max_length(num_length_levels);
    max_length[0] = start_length;
    for( int i=1; i<num_length_levels; ++i )
    {
        if( length_type == "Multiplicative" )
        {
            TEUCHOS_ASSERT( length_factor > 1 );
            max_length[i] = max_length[i-1] * length_factor;
        }
        else if( length_type == "Additive" )
        {
            TEUCHOS_ASSERT( length_factor > 0 );
            max_length[i] = max_length[i-1] + length_factor;
        }
    }

    std::string solver_type = pl->get("solver_type","mcsa");
    pl->set("verbosity","none");

    if( profugus::node() == 0 )
    {
        std::cout << "Num histories: ";
        for( int j=0; j<num_history_levels; ++j )
            std::cout << num_histories[j] << " ";
        std::cout << std::endl;

        std::cout << "Max history length: ";
        for( int i=0; i<num_length_levels; ++i )
            std::cout << max_length[i] << " ";
        std::cout << std::endl;
    }

    Teuchos::RCP<Teuchos::ParameterList> poly_pl =
        Teuchos::sublist(pl,"Polynomial");
    Teuchos::RCP<Teuchos::ParameterList> mc_pl =
        Teuchos::sublist(pl,"Monte Carlo");

    // Flush output stream on every write
    std::cout << std::unitbuf;

    std::string itemsep = " ";
    std::string lineend = "";
    if( output_format == "latex" )
    {
        itemsep = " & ";
        lineend = " \\\\";
    }

    for( int j=0; j<num_history_levels; ++j )
    {
        mc_pl->set("num_histories",num_histories[j]);
        std::cout << num_histories[j];
        for( int i=0; i<num_length_levels; ++i )
        {
            int iters=0;
            int ms=0;
            // Catch exceptions to allow multirun to proceed
            // In particular, certain polynomial objects may fail at
            // construction time and throw exceptions.
            // An example is the GMRES polynomial constructed by the
            // normal equations which may be too ill conditioned to
            // solve at high polynomial orders.
            try
            {
                poly_pl->set("polynomial_order",max_length[i]);

                Teuchos::RCP<alea::AleaSolver> solver =
                    alea::LinearSolverFactory::buildSolver(solver_type,myA,pl);

                time_start = high_resolution_clock::now();
                solver->apply(*b,*x);
                time_end = high_resolution_clock::now();

                iters = solver->getNumIters();
                ms = duration_cast<milliseconds>(time_end-time_start).count();
            }
            catch(...)
            {
            }

            if( profugus::node() == 0 )
            {
                std::cout << itemsep;
                if( iters > 0 )
                {
                    std::cout << iters << "(" << ms << ")";
                }
                else
                {
                    std::cout << "-";
                }
            }
        }
        if( profugus::node() == 0 )
            std::cout << lineend << std::endl;
    }

    // Initialize Kokkos device
    DeviceTraits<DEVICE>::finalize();

    profugus::finalize();
    return 0;
}


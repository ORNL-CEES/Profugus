//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/mc/Anderson_Solver.t.hh
 * \author Thomas M. Evans
 * \date   Tue Apr 07 20:53:11 2015
 * \brief  Anderson_Solver template method definitions.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_mc_Anderson_Solver_t_hh
#define MC_mc_Anderson_Solver_t_hh

#include <iostream>
#include <iomanip>

#include "Teuchos_Array.hpp"

#include "harness/DBC.hh"
#include "comm/Timing.hh"
#include "spn/MatrixTraits.hh"
#include "Anderson_Solver.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template<class T>
Anderson_Solver<T>::Anderson_Solver(RCP_Std_DB db)
    : d_db(db)
    , d_nodes(profugus::nodes())
    , d_node(profugus::node())
    , d_k_anderson(0.0)
{
    REQUIRE(!d_db.is_null());

    // Process the solver options
    VALIDATE(db->isSublist("anderson_db"), "No anderson_db nested "
             << "database defined for Anderson solver.");

    // Get the Anderson database
    auto &adb = db->sublist("anderson_db");

    // validate that there is a fission source mesh defined
    VALIDATE(adb.isParameter("x_bounds"),
             "Failed to define x-boundaries for anderson mesh");
    VALIDATE(adb.isParameter("y_bounds"),
             "Failed to define y-boundaries for anderson mesh");
    VALIDATE(adb.isParameter("z_bounds"),
             "Failed to define z-boundaries for anderson mesh");

    // Set Anderson solver defaults
    adb.template get<int>("Storage Depth", 2);
    adb.template get<double>("Mixing Parameter", 0.8);
    adb.template get<int>("Acceleration Start Iteration", 0);
    adb.template get<double>("tolerance", 1.0e-3);

    // Make a set-constant communicator
    profugus::split(d_node, 0, d_set_comm);

#ifdef ENSURE_ON
    profugus::set_internal_comm(d_set_comm);
    ENSURE(profugus::nodes() == 1);
    ENSURE(profugus::node() == 0);
    profugus::reset_internal_comm();
#endif
    ENSURE(profugus::nodes() == d_nodes);
    ENSURE(profugus::node()  == d_node);
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Set the underlying fixed-source transporter and fission source.
 */
template<class T>
void Anderson_Solver<T>::set(SP_Source_Transporter transporter,
                             SP_Fission_Source     source)
{
    typedef Teuchos::Array<double> OneDArray_dbl;

    // get a reference to the tallier
    b_tallier = transporter->tallier();
    INSIST(b_tallier, "The tallier has not been assigned.");
    CHECK(b_tallier->geometry() && b_tallier->physics());
    CHECK(!b_tallier->is_built());

    // get initial k and build keff tally
    double init_keff = d_db->get("keff_init", 1.0);
    VALIDATE(init_keff >= 0., "Initial keff guess (keff_init="
             << init_keff << ") must be nonnegative.");
    d_keff_tally = std::make_shared<Keff_Tally>(
        init_keff, b_tallier->physics());

    // Get the anderson mesh boundaries
    auto adb = Teuchos::sublist(d_db, "anderson_db");
    auto xb  = adb->template get<OneDArray_dbl>("x_bounds").toVector();
    auto yb  = adb->template get<OneDArray_dbl>("y_bounds").toVector();
    auto zb  = adb->template get<OneDArray_dbl>("z_bounds").toVector();

    // Make the Anderson mesh
    auto mesh = std::make_shared<profugus::Cartesian_Mesh>(xb, yb, zb);
    CHECK(mesh);

    // build the map for the solver
    profugus::set_internal_comm(d_set_comm);
    RCP_MAP map = MatrixTraits<T>::build_map(
        mesh->num_cells() + 1, mesh->num_cells() + 1);
    profugus::reset_internal_comm();
    CHECK(!map.is_null());

    // Build the operator
    d_operator = Teuchos::rcp(
        new Operator(transporter, source, mesh, map, d_set_comm));

    // Make the Anderson nonlinear solver
    d_anderson = std::make_shared<Anderson_t>(d_operator, adb);

    // store the integrated cycle weight (const number of requested particles
    // per cycle)
    d_Np = static_cast<double>(source->Np());
    CHECK(d_Np > 0.0);

    ENSURE(!d_operator.is_null());
    ENSURE(b_tallier->is_built());
}

//---------------------------------------------------------------------------//
// DERIVED INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief Solve the problem..
 */
template<class T>
void Anderson_Solver<T>::solve()
{
    using std::cout;

    REQUIRE(!d_operator.is_null());

    SCOPED_TIMER("MC::Anderson_Solver.solve");

    // Set output parameters
    cout.precision(6);
    cout.setf(std::ios::internal);

    // >>> INITIALIZE THE PROBLEM
    initialize();

    // >>> RUN A SINGLE INACTIVE CYCLE
    run_inactive();

    // >>> SOLVE THE EIGENVALUE PROBLEM WITH ANDERSON
    anderson_solve();

    // >>> RUN ACTIVE CYCLES
    run_active();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Call to reset the solver and tallies for another calculation.
 */
template<class T>
void Anderson_Solver<T>::reset()
{
    NOT_IMPLEMENTED("Reset not available.");
}

//---------------------------------------------------------------------------//
// PRIVATE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Initialize the k-eigenvalue problem.
 *
 * The initialize method adds the keffective tally to the tallier and builds
 * the initial/fission source.
 */
template<class T>
void Anderson_Solver<T>::initialize()
{
    INSIST(!b_tallier->is_built(), "The given tallier can only be built once.");
    INSIST(!b_tallier->is_finalized(), "The given tallier was used in "
           "a prior transport solve without reset() being called.");

    // Add the keff tally to the ACTIVE cycle tallier and build it
    b_tallier->add_pathlength_tally(d_keff_tally);
    b_tallier->build();
    CHECK(b_tallier->is_built());
    CHECK(b_tallier->num_pathlength_tallies() >= 1);

    // build the source
    auto source = d_operator->source();
    CHECK(source);
    if (source->is_initial_source())
    {
        // build the initial fission source
        source->build_initial_source();
    }

    ENSURE(b_tallier->num_pathlength_tallies() >= 1);
    ENSURE(d_keff_tally->cycle_count() == 0);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Run inactive cycle(s).
 */
template<class T>
void Anderson_Solver<T>::run_inactive()
{
    using std::endl; using std::cout;
    using std::setw; using std::fixed; using std::scientific;

    // Timers (some required, some optional)
    Timer cycle_timer;

    // First set a null tallier
    auto first_tallier = std::make_shared<Tallier_t>();
    first_tallier->set(b_tallier->geometry(), b_tallier->physics());
    first_tallier->add_pathlength_tally(d_keff_tally);
    first_tallier->build();
    CHECK(first_tallier->is_built());
    CHECK(first_tallier->num_tallies() == 1);

    // Set the tallier in the operator
    d_operator->set_tallier(first_tallier);

    // Run a cycle
    cout << ">>> Beginning inactive cycles." << endl;
    cout << endl;
    cout << " Cycle  k_cycle   Time (s) " << endl;

    for (int cycle = 0; cycle < 1; ++cycle)
    {
        // Iterate and time
        cycle_timer.start();
        d_operator->iterate(d_keff_tally->latest());
        cycle_timer.stop();

        if (d_node == 0)
        {
            cout << fixed      << setw(4) << cycle << "   "
                 << fixed      << setw(8) << d_keff_tally->latest() << "  "
                 << scientific << setw(9) << cycle_timer.TIMER_CLOCK()
                 << endl;
        }
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Do Anderson solve.
 */
template<class T>
void Anderson_Solver<T>::anderson_solve()
{
    using std::endl; using std::cout;
    using std::setw; using std::fixed; using std::scientific;

    REQUIRE(d_anderson);

    // First set a null tallier
    auto null_tallier = std::make_shared<Tallier_t>();
    null_tallier->set(b_tallier->geometry(), b_tallier->physics());
    null_tallier->build();
    CHECK(null_tallier->is_built());
    CHECK(null_tallier->num_tallies() == 0);

    // Set the null tallier
    d_operator->set_tallier(null_tallier);

    // Initialize the Anderson operator and get a solution vector
    auto v = d_operator->initialize_Anderson();

    // Solve using Anderson
    d_anderson->solve(v);

    // Finalize the operator after Anderson solve in order to resume active
    // cycles
    d_k_anderson = d_operator->finalize_Anderson(*v);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Run active cycles.
 */
template<class T>
void Anderson_Solver<T>::run_active()
{
    using std::endl; using std::cout;
    using std::setw; using std::fixed; using std::scientific;

    REQUIRE(b_tallier);
    REQUIRE(b_tallier->is_built());

    // set cycle numbers
    const int num_total = d_db->get("num_cycles", 50);

    // Timers (some required, some optional)
    Timer cycle_timer;

    // Set the tallier to run active cycles
    b_tallier->begin_active_cycles();
    CHECK(d_keff_tally->cycle_count() == 0);

    // Set the keff tally to start active cycles with k from Anderson
    d_keff_tally->set_keff(d_k_anderson);

    // Set the active tallier in the operator
    d_operator->set_tallier(b_tallier);

    cout << endl;
    cout << ">>> Beginning active cycles." << endl;
    cout << endl;
    cout << " Cycle  k_cycle    k_avg    k_variance    Time (s) "
         << endl;

    // Solve active cycles
    for (int cycle = 0; cycle < num_total; ++cycle)
    {
        // Iterate and time
        cycle_timer.start();
        d_operator->iterate(d_keff_tally->latest());
        d_operator->update_source();
        cycle_timer.stop();

        // Print diagnostic output
        if (d_node == 0)
        {
            cout << fixed      << setw(4)  << cycle << "   "
                 << fixed      << setw(8)  << d_keff_tally->latest() << "  "
                 << fixed      << setw(8)  << d_keff_tally->mean() << "  "
                 << scientific << setw(11) << d_keff_tally->variance() << "  "
                 << scientific << setw(11) << cycle_timer.TIMER_CLOCK()
                 << endl;
        }
    }
    cout << endl;
    CHECK(num_total == d_keff_tally->cycle_count());

    profugus::global_barrier();

    // >>> FINALIZE THE TALLIES
    CHECK(!b_tallier->is_finalized());
    b_tallier->finalize(num_total * d_Np);

    ENSURE(b_tallier->is_finalized());
}

} // end namespace profugus

#endif // MC_mc_Anderson_Solver_t_hh

//---------------------------------------------------------------------------//
// end of MC/mc/Anderson_Solver.t.hh
//---------------------------------------------------------------------------//

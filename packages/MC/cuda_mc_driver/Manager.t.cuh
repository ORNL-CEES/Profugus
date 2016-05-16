//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc_driver/Manager.t.cuh
 * \author Thomas M. Evans
 * \date   Wed Jun 18 11:21:16 2014
 * \brief  Manager member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_driver_Manager_t_hh
#define cuda_mc_driver_Manager_t_hh

#include <fstream>
#include <iomanip>

#include "Teuchos_XMLParameterListHelpers.hpp"

#include "harness/DBC.hh"
#include "comm/Timing.hh"
#include "comm/global.hh"
#include "utils/Serial_HDF5_Writer.hh"
#include "utils/Parallel_HDF5_Writer.hh"
#include "solvers/LinAlgTypedefs.hh"
#include "cuda_mc/Fission_Source.hh"
#include "cuda_mc/Uniform_Source.hh"
#include "cuda_mc/KCode_Solver.hh"
#include "Manager.hh"

namespace cuda_mc
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template <class Geometry>
Manager<Geometry>::Manager()
    : d_node(profugus::node())
    , d_nodes(profugus::nodes())
{
}

//---------------------------------------------------------------------------//
/*!
 * \brief Setup the problem.
 * \param master Problem parameters.
 */
template <class Geometry>
void Manager<Geometry>::setup(RCP_ParameterList master)
{
    SCOPED_TIMER("Manager.setup");

    SCREEN_MSG("Building and initializing geometry, physics, "
               << "variance reduction, and tallies");

    // use the problem builder to setup the problem
    Prob_Builder builder;
    builder.setup(master);

    // get the problem database from the problem-builder
    d_db = builder.problem_db();
    CHECK(!d_db.is_null());

    // store the problem name
    d_problem_name = d_db->get("problem_name", std::string("MC"));
    d_db->setName(d_problem_name + "-PROBLEM");

    // get the geometry and physics
    d_geometry = builder.get_geometry();
    d_physics  = builder.get_physics();
    CHECK(d_geometry);
    CHECK(d_physics);

    // get the variance reduction
    auto var_reduction = builder.get_var_reduction();
    CHECK(var_reduction);

    // get the external source shape (it could be null)
    auto shape = builder.get_source_shape();

    // problem type
    std::string prob_type = shape ? "fixed" : "eigenvalue";

    // set the problem type in the final db
    d_db->set("problem_type", prob_type);

    // build the random controller
    d_rnd_control = std::make_shared<RNG_Control_t>(
        d_db->template get<int>("seed", 32442));

    SCREEN_MSG("Building " << prob_type << " solver");

    // get the tallier
    SP_Tallier tallier = builder.get_tallier();

    // make the transporter
    SP_Transporter transporter(std::make_shared<Transporter_t>(
                                   d_db, d_geometry, d_physics));
    transporter->set(tallier);
    transporter->set(var_reduction);

    // build the appropriate solver (default is eigenvalue)
    if (prob_type == "eigenvalue")
    {
        // make the fission source
        SP_Fission_Source source(
            std::make_shared<Fission_Source_t>(
                d_db, d_geometry, d_physics, d_rnd_control));

        // >>> determine eigensolver

        // Standard K-Code
        {
            // make the solver
            auto kcode_solver =
                std::make_shared<profugus::KCode_Solver<Geom_t> >(d_db);

            // set the solver
            kcode_solver->set(transporter, source);

            // set hybrid acceleration
            kcode_solver->set(builder.get_acceleration());

            // assign the base solver
            d_keff_solver = kcode_solver;
        }

        // assign the base solver
        d_solver = d_keff_solver;
    }
    else if (prob_type == "fixed")
    {
        // make the uniform source
        std::shared_ptr<profugus::Uniform_Source<Geometry> > source(
            std::make_shared<profugus::Uniform_Source<Geometry> >(
                d_db, d_geometry, d_physics, d_rnd_control));
        source->build_source(shape);

        // make the solver
        d_fixed_solver = std::make_shared<Fixed_Source_Solver_t>();

        // set it
        d_fixed_solver->set(transporter, source);

        // assign the base solver
        d_solver = d_fixed_solver;
    }
    else
    {
        throw profugus::assertion(
            "Undefined problem type; choose eigenvalue or fixed");
    }

    ENSURE(d_geometry);
    ENSURE(d_physics);
    ENSURE(d_solver);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Solve the problem.
 */
template <class Geometry>
void Manager<Geometry>::solve()
{
    if (d_db->template get<bool>("do_transport", true))
    {
        SCOPED_TIMER("Manager.solve");

        SCREEN_MSG("Executing solver");

        // run the appropriate solver
        if (d_keff_solver)
        {
            CHECK(!d_fixed_solver);
            d_keff_solver->solve();
        }
        else if (d_fixed_solver)
        {
            CHECK(!d_keff_solver);
            d_fixed_solver->solve();
        }
        else
        {
            throw profugus::assertion("No legitimate solver built");
        }
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Do output.
 */
template <class Geometry>
void Manager<Geometry>::output()
{
    using std::string;

    SCOPED_TIMER("Manager.output");

    SCREEN_MSG("Outputting data");

    // >>> OUTPUT FINAL DATABASE

    // output the final database
    if (d_node == 0)
    {
        std::ostringstream m;
        m << d_problem_name << "_db.xml";
        Teuchos::writeParameterListToXmlFile(*d_db, m.str());
    }

    profugus::global_barrier();

    // >>> OUTPUT SOLUTION (only available if HDF5 is on)
#ifdef USE_HDF5

    // Output filename
    std::ostringstream m;
    m << d_problem_name << "_output.h5";
    string outfile = m.str();

    // scalar output for kcode
    if (d_keff_solver)
    {
        // get the kcode tally
        auto keff = d_keff_solver->keff_tally();
        CHECK(keff);

        // make the hdf5 file
        profugus::Serial_HDF5_Writer writer;
        writer.open(outfile);

        // output scalar quantities
        writer.begin_group("keff");

        writer.write(string("mean"), keff->mean());
        writer.write(string("variance"), keff->variance());
        writer.write(string("num_active_cycles"),
                     static_cast<int>(keff->cycle_count()));
        writer.write(string("cycle_estimates"), keff->all_keff());

        // do diagnostics on acceleration if it exists
        if (d_keff_solver->acceleration())
        {
            d_keff_solver->acceleration()->diagnostics(writer);
        }

        writer.end_group();

        writer.close();
    }

#endif // USE_HDF5
}

} // end namespace cuda_mc

#endif // cuda_mc_driver_Manager_t_hh

//---------------------------------------------------------------------------//
//                 end of Manager.t.cuh
//---------------------------------------------------------------------------//

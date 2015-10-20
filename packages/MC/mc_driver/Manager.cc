//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc_driver/Manager.cc
 * \author Thomas M. Evans
 * \date   Wed Jun 18 11:21:16 2014
 * \brief  Manager member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <fstream>
#include <iomanip>

#include "Teuchos_XMLParameterListHelpers.hpp"

#include "harness/DBC.hh"
#include "comm/Timing.hh"
#include "comm/global.hh"
#include "utils/Serial_HDF5_Writer.hh"
#include "utils/Parallel_HDF5_Writer.hh"
#include "physics/Uniform_Source.hh"
#include "Manager.hh"

namespace mc
{
//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
Manager::Manager()
    : d_node(profugus::node())
    , d_nodes(profugus::nodes())
{
}

//---------------------------------------------------------------------------//
/*!
 * \brief Setup the problem.
 * \param xml_file
 */
void Manager::setup(const std::string &xml_file)
{
    SCOPED_TIMER("Manager.setup");

    SCREEN_MSG("Reading xml file -> " << xml_file);
    SCREEN_MSG("Building and initializing geometry, physics, "
               << "variance reduction, and tallies");

    // use the problem builder to setup the problem
    Problem_Builder builder;
    builder.setup(xml_file);

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

    // output the geometry
    if (d_db->get<bool>("output_geometry", false))
    {
        std::ostringstream g;
        g << d_db->get<std::string>("problem_name") << "_geo.out";

        if (d_node == 0)
        {
            // geo file
            std::ofstream gfile(g.str().c_str(), std::ofstream::out);

            // output the geometry
            d_geometry->array().output(gfile);
        }
    }

    // get the variance reduction
    auto var_reduction = builder.get_var_reduction();
    CHECK(var_reduction);

    // get the external source shape (it could be null)
    auto shape = builder.get_source_shape();

    // problem type
    std::string prob_type = "fixed";

    // set the problem type in the final db
    d_db->set("problem_type", prob_type);

    // get the tallier
    d_tallier = builder.get_tallier();

    // make the transporter
    SP_Transporter transporter(std::make_shared<Transporter_t>(
                                   d_db, d_geometry, d_physics));
    transporter->set(d_tallier);
    transporter->set(var_reduction);

    REQUIRE(prob_type == "fixed");

    // make the uniform source
    std::shared_ptr<profugus::Uniform_Source> source(
	std::make_shared<profugus::Uniform_Source>(
	    d_db, d_geometry, d_physics));
    source->build_source(shape);

    // make the solver
    d_fixed_solver = std::make_shared<Fixed_Source_Solver_t>();

    // set it
    d_fixed_solver->set(transporter, source);

    // assign the base solver
    d_solver = d_fixed_solver;


#ifdef USE_HDF5
    // Output filename
    std::ostringstream m;
    m << d_problem_name << "_output.h5";
    d_output_name = m.str();

    SCREEN_MSG("Initializing " << d_output_name << " output file");

    // make the hdf5 file
    profugus::Serial_HDF5_Writer writer;
    writer.open(d_output_name);
    writer.close();
#else
    ADD_WARNING("HDF5 not available in this build, tally output will be off.");
#endif

    ENSURE(d_geometry);
    ENSURE(d_physics);
    ENSURE(d_solver);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Solve the problem.
 */
void Manager::solve()
{
    if (d_db->get<bool>("do_transport", true))
    {
        SCOPED_TIMER("Manager.solve");

        SCREEN_MSG("Executing solver");

        REQUIRE(d_fixed_solver);
	d_fixed_solver->solve();
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Do output.
 */
void Manager::output()
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

    // Only do output if we did transport
    if (!d_db->get<bool>("do_transport"))
    {
        return;
    }

    // Output all other end-of-problem tallies
    d_tallier->output(d_output_name);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Write timing to output.
 */
void Manager::timing()
{
    SCREEN_MSG("Outputting timing data");

#ifdef USE_HDF5

    // make the hdf5 file
    profugus::Serial_HDF5_Writer writer;
    writer.open(d_output_name, profugus::HDF5_IO::APPEND);

    // make the timing group
    writer.begin_group("timing");

    const auto &keys = profugus::Timing_Diagnostics::timer_keys();
    for (auto k : keys)
    {
        writer.write(k, profugus::Timing_Diagnostics::timer_value(k));
    }

    writer.end_group();

    writer.close();

#endif // USE_HDF5

}

} // end namespace mc

//---------------------------------------------------------------------------//
//                 end of Manager.cc
//---------------------------------------------------------------------------//

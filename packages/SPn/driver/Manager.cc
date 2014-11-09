//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   driver/Manager.cc
 * \author Thomas M. Evans
 * \date   Fri Mar 14 11:32:36 2014
 * \brief  Manager member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "EpetraExt_RowMatrixOut.h"
#include "Epetra_RowMatrix.h"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "harness/DBC.hh"
#include "comm/Timing.hh"
#include "comm/global.hh"
#include "utils/String_Functions.hh"
#include "utils/Parallel_HDF5_Writer.hh"
#include "utils/Definitions.hh"
#include "spn/Dimensions.hh"
#include "spn/SpnSolverBuilder.hh"
#include "Manager.hh"

namespace spn
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
    SCREEN_MSG("Building and initializing mesh");

    // use the problem builder to setup the problem
    Problem_Builder builder;
    builder.setup(xml_file);

    // get the problem database from the problem-builder
    d_db = builder.problem_db();

    // store the problem name
    d_problem_name = d_db->get("problem_name", std::string("SPn"));
    d_db->setName(d_problem_name + "-PROBLEM");

    // get the mesh objects from the builder
    d_mesh    = builder.mesh();
    d_indexer = builder.indexer();
    d_gdata   = builder.global_data();

    // get the material database from the problem builder
    d_mat = builder.mat_db();

    // get the source (it will be null for eigenvalue problems)
    d_external_source = builder.source();

    // build the problem dimensions
    d_dim = Teuchos::rcp(new profugus::Dimensions(d_db->get("SPn_order", 1)));

    // problem type
    std::string prob_type = d_external_source.is_null() ? "eigenvalue" :
                            "fixed";

    // if a timestep control is defined then this is a time-dependent problem
    if (d_db->isSublist("timestep control"))
    {
        INSIST(prob_type == "fixed",
                "Cannot do time-dependent eigenproblems.");
        prob_type = "fixed_tdep";
    }

    // set the problem type in the final db
    d_db->set("problem_type", prob_type);

    SCREEN_MSG("Building " << prob_type << " solver");

    // default linear solver type (stratimikios)
    d_db->get("solver_type", std::string("stratimikos"));

    // build the solver
    d_solver_base =
        profugus::SpnSolverBuilder::build(prob_type,d_db);

    // set default as a forward (not adjoint) run
    bool adjoint = d_db->get("adjoint", false);

    // setup the solver
    d_solver_base->setup(d_dim, d_mat, d_mesh, d_indexer, d_gdata, adjoint);

    // make the state
    d_state = Teuchos::rcp(
        new profugus::State(d_mesh, d_mat->xs().num_groups()));

    ENSURE(!d_mesh.is_null());
    ENSURE(!d_indexer.is_null());
    ENSURE(!d_gdata.is_null());
    ENSURE(!d_mat.is_null());
    ENSURE(!d_dim.is_null());
    ENSURE(!d_state.is_null());
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

        // solve problem (source will be null for eigenvalue problems)
        d_solver_base->solve(d_external_source);

        // write the solution vector into the state
        d_solver_base->write_state(*d_state);
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Do output.
 */
void Manager::output()
{
    using def::I; using def::J; using def::K;

    REQUIRE(!d_state.is_null());
    REQUIRE(!d_db.is_null());

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

    // >>> OUTPUT SOLUTION

    // Output filename
    std::ostringstream m;
    m << d_problem_name << "_output.h5";
    std::string outfile = m.str();

    // get a constant reference to the state
    const profugus::State &state = *d_state;

    // group offset (if doing a truncated range)
    auto g_first = d_db->get<int>("g_first");
    CHECK(1 + d_db->get<int>("g_last") - g_first == state.num_groups());

    // output the fluxes if Parallel HDF5 is available
#ifdef H5_HAVE_PARALLEL
    {
        // make the parallel hdf5 writer
        profugus::Parallel_HDF5_Writer writer;
        writer.open(outfile);

        // make the decomposition for parallel output (state is ordered
        // i->j->k)
        profugus::HDF5_IO::Decomp d(3, profugus::HDF5_IO::COLUMN_MAJOR);
        d.ndims     = 3;
        d.global[I] = d_gdata->num_cells(I);
        d.global[J] = d_gdata->num_cells(J);
        d.global[K] = d_gdata->num_cells(K);
        d.local[I]  = d_mesh->num_cells_dim(I);
        d.local[J]  = d_mesh->num_cells_dim(J);
        d.local[K]  = d_mesh->num_cells_dim(K);
        d.offset[I] = d_indexer->offset(I);
        d.offset[J] = d_indexer->offset(J);
        d.offset[K] = 0; // currently do not partition in K

        // make a group for the fluxes
        writer.begin_group("fluxes");

        // loop over groups and write the fluxes
        for (int g = 0, Ng = d_mat->xs().num_groups(); g < Ng; ++g)
        {
            // dataset name
            std::ostringstream f;
            f << "group_" << g + g_first;

            // get the group fluxes
            profugus::State::const_View_Field flux = state.flux(g, g);
            CHECK(!flux.is_null());
            CHECK(flux.size() == d_mesh->num_cells());

            // write that data in parallel
            writer.write(f.str(), d, flux.getRawPtr());

            profugus::global_barrier();
        }

        writer.end_group();
        writer.close();
    }
#endif

    profugus::global_barrier();

    // >>> OUTPUT MATRICES
    if (d_db->get<bool>("output_matrices", false))
    {
        d_solver_base->write_problem_to_file();
    }
}

} // end namespace spn

//---------------------------------------------------------------------------//
//                 end of Manager.cc
//---------------------------------------------------------------------------//

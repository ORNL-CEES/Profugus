//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   driver/Manager.cc
 * \author Thomas M. Evans
 * \date   Fri Mar 14 11:32:36 2014
 * \brief  Manager member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Teuchos_XMLParameterListHelpers.hpp"

#include "harness/DBC.hh"
#include "comm/Timing.hh"
#include "comm/global.hh"
#include "utils/String_Functions.hh"
#include "spn/Dimensions.hh"
#include "Manager.hh"

namespace profugus
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
    d_dim = Teuchos::rcp(new Dimensions(d_db->get("SPn_order", 1)));

    // problem type
    std::string prob_type = d_external_source.is_null() ? "eigenvalue" :
                            "fixed";

    // set the problem type in the final db
    d_db->set("problem_type", prob_type);

    SCREEN_MSG("Building " << prob_type << " solver");

    // default linear solver type (stratimikios)
    d_db->get("solver_type", std::string("stratimikos"));

    // build the appropriate solver (default is eigenvalue)
    if (prob_type == "eigenvalue")
    {
        d_eigen_solver = Teuchos::rcp(new Eigenvalue_Solver(d_db));
        d_solver_base  = d_eigen_solver;
    }
    else if (prob_type == "fixed")
    {
        d_fixed_solver = Teuchos::rcp(new Fixed_Source_Solver(d_db));
        d_solver_base  = d_fixed_solver;
    }
    else
    {
        throw profugus::assertion(
            "Undefined problem type; choose eigenvalue or fixed");
    }

    // setup the solver
    d_solver_base->setup(d_dim, d_mat, d_mesh, d_indexer, d_gdata);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Solve the problem.
 */
void Manager::solve()
{
    SCOPED_TIMER("Manager.solve");

    SCREEN_MSG("Executing solver");

    // run the appropriate solver
    if (!d_eigen_solver.is_null())
    {
        Check (d_fixed_solver.is_null());
        d_eigen_solver->solve();
    }
    else
    {
        Check (!d_fixed_solver.is_null());
        Check (d_eigen_solver.is_null());
        Check (!d_external_source.is_null());
        d_fixed_solver->solve(*d_external_source);
    }

    // write the solution vector into the state
    d_solver_base->write_state(*d_state);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Do output.
 */
void Manager::output()
{
    SCOPED_TIMER("Manager.output");

    SCREEN_MSG("Outputting data");

    // output the final database
    if (d_node == 0)
    {
        std::ostringstream m;
        m << d_problem_name << "_db.xml";
        Teuchos::writeParameterListToXmlFile(*d_db, m.str());
    }
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Manager.cc
//---------------------------------------------------------------------------//

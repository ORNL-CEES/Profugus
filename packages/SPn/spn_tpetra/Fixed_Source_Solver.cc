//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   spn_tpetra/Fixed_Source_Solver.cc
 * \author Thomas M. Evans
 * \date   Mon Feb 17 21:00:33 2014
 * \brief  Fixed_Source_Solver member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <string>

#include "harness/DBC.hh"
#include "utils/String_Functions.hh"
#include "Linear_System_FV.hh"
#include "Fixed_Source_Solver.hh"

namespace profugus
{
namespace tpetra
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
Fixed_Source_Solver::Fixed_Source_Solver(RCP_ParameterList db)
    : Base(db)
    , d_solver(b_db)
{
    ENSURE(!b_db.is_null());
}

//---------------------------------------------------------------------------//
// PUBLIC INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief Setup the solver.
 *
 * Calls to this function builds the linear SPN system.
 */
void Fixed_Source_Solver::setup(RCP_Dimensions  dim,
                                RCP_Mat_DB      mat,
                                RCP_Mesh        mesh,
                                RCP_Indexer     indexer,
                                RCP_Global_Data data)
{
    REQUIRE(!b_db.is_null());
    REQUIRE(!dim.is_null());
    REQUIRE(!mat.is_null());
    REQUIRE(!mesh.is_null());
    REQUIRE(!indexer.is_null());
    REQUIRE(!data.is_null());

    // build the linear system (we only provide finite volume for now)
    std::string &eqn_type = b_db->get("eqn_type", std::string("fv"));

    if (profugus::to_lower(eqn_type) == "fv")
    {
        b_system = Teuchos::rcp(
            new Linear_System_FV(b_db, dim, mat, mesh, indexer, data));
    }
    else
    {
        std::string msg = "Undefined equation type: " + eqn_type;
        throw profugus::assertion(msg);
    }
    CHECK(!b_system.is_null());

    // build the matrix
    b_system->build_Matrix();

    // register the operator with the solver
    d_solver.set_operator(b_system->get_Operator());

    // allocate the left-hand side solution vector
    d_lhs = Teuchos::rcp(new Vector_t(b_system->get_Map()));

    ENSURE(b_system->get_Map()->getNodeNumElements() ==
            d_lhs->getLocalLength());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Solve the SPN equations for a given external source.
 */
void Fixed_Source_Solver::solve(const External_Source &q)
{
    REQUIRE(!b_system.is_null());
    REQUIRE(!d_lhs.is_null());

    // null lhs vector
    d_lhs->putScalar(0.0);

    // make the right-hand side vector based on the source
    b_system->build_RHS(q);
    CHECK(b_system->get_RHS()->getLocalLength() == d_lhs->getLocalLength());

    // solve the problem
    d_solver.solve(d_lhs, b_system->get_RHS());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Write scalar flux into the state.
 */
void Fixed_Source_Solver::write_state(State_t &state)
{
    REQUIRE(state.mesh().num_cells() *
             b_system->get_dims()->num_equations() * state.num_groups()
             <= d_lhs->getLocalLength());

    Base::write_u_into_state(*d_lhs, state);
}

} // end namespace profugus
} // end namespace tpetra

//---------------------------------------------------------------------------//
//                 end of Fixed_Source_Solver.cc
//---------------------------------------------------------------------------//

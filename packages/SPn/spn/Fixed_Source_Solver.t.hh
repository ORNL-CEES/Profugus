//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   spn/Fixed_Source_Solver.t.hh
 * \author Thomas M. Evans
 * \date   Mon Feb 17 21:00:33 2014
 * \brief  Fixed_Source_Solver template member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef spn_Fixed_Source_Solver_t_hh
#define spn_Fixed_Source_Solver_t_hh

#include <string>

#include "harness/DBC.hh"
#include "utils/String_Functions.hh"
#include "solvers/LinAlgTypedefs.hh"
#include "Linear_System_FV.hh"
#include "Fixed_Source_Solver.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template <class T>
Fixed_Source_Solver<T>::Fixed_Source_Solver(RCP_ParameterList db)
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
template <class T>
void Fixed_Source_Solver<T>::setup(RCP_Dimensions  dim,
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
            new Linear_System_FV<EpetraTypes>(
                b_db, dim, mat, mesh, indexer, data));
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
    d_lhs = VectorTraits<T>::build_vector(b_system->get_Map());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Solve the SPN equations for a given external source.
 */
template <class T>
void Fixed_Source_Solver<T>::solve(const External_Source &q)
{
    REQUIRE(!b_system.is_null());
    REQUIRE(!d_lhs.is_null());

    // null lhs vector
    VectorTraits<T>::put_scalar(d_lhs,0.0);

    // make the right-hand side vector based on the source
    b_system->build_RHS(q);
    CHECK( VectorTraits<T>::local_length(b_system->get_RHS()) ==
           VectorTraits<T>::local_length(d_lhs) );

    // solve the problem
    d_solver.solve(d_lhs, b_system->get_RHS());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Write scalar flux into the state.
 */
template <class T>
void Fixed_Source_Solver<T>::write_state(State &state)
{
    REQUIRE(state.mesh().num_cells() *
             b_system->get_dims()->num_equations() * state.num_groups()
             <= VectorTraits<T>::local_length(d_lhs) );

    Base::write_u_into_state(d_lhs, state);
}

} // end namespace profugus

#endif // spn_Fixed_Source_Solver_t_hh

//---------------------------------------------------------------------------//
//                 end of Fixed_Source_Solver.cc
//---------------------------------------------------------------------------//

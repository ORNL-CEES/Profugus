//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Fission_Matrix_Acceleration.t.hh
 * \author Thomas M. Evans
 * \date   Thu Nov 13 19:53:21 2014
 * \brief  Fission_Matrix_Acceleration template member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef mc_Fission_Matrix_Acceleration_t_hh
#define mc_Fission_Matrix_Acceleration_t_hh

#include "spn/Dimensions.hh"
#include "spn/SpnSolverBuilder.hh"
#include "spn/Linear_System_FV.hh"
#include "spn/Eigenvalue_Solver.hh"
#include "Fission_Matrix_Acceleration.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template<class T>
Fission_Matrix_Acceleration_Impl<T>::Fission_Matrix_Acceleration_Impl()
{
}

//---------------------------------------------------------------------------//
// PUBLIC INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief Build the SPN problem.
 *
 * \param builder
 */
template<class T>
void Fission_Matrix_Acceleration_Impl<T>::build_problem(
    const Problem_Builder_t &builder)
{
    // get the problem database from the problem-builder
    b_db = builder.problem_db();

    // get the mesh objects from the builder
    b_mesh    = builder.mesh();
    b_indexer = builder.indexer();
    b_gdata   = builder.global_data();

    // get the material database from the problem builder
    b_mat = builder.mat_db();

    // build the problem dimensions
    d_dim = Teuchos::rcp(new profugus::Dimensions(b_db->get("SPn_order", 1)));

    // default linear solver type (stratimikios)
    b_db->get("solver_type", std::string("stratimikos"));

    // make the forward and adjoint states
    b_forward = Teuchos::rcp(new profugus::State(
                                 b_mesh, b_mat->xs().num_groups()));
    b_adjoint = Teuchos::rcp(new profugus::State(
                                 b_mesh, b_mat->xs().num_groups()));

    // make the linear system for this problem
    d_system = Teuchos::rcp(
        new Linear_System_FV<T>(b_db, d_dim, b_mat, b_mesh, b_indexer, b_gdata));

    // make the matrices (A,B) for the SPN problem, Ap = (1/k)Bp
    d_system->build_Matrix();
    d_system->build_fission_matrix();

    ENSURE(!b_mesh.is_null());
    ENSURE(!b_indexer.is_null());
    ENSURE(!b_gdata.is_null());
    ENSURE(!b_mat.is_null());
    ENSURE(!d_dim.is_null());
    ENSURE(!b_adjoint.is_null());
    ENSURE(!b_forward.is_null());
    ENSURE(!d_system.is_null());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Initialize the acceleration at the beginning of the K-code solve.
 *
 * This call invokes two eigenvalue solves:
 * \f[
   \mathbf{A}\phi = \frac{1}{k}\mathbf{B}\phi\:,
 * \f]
 * and
 * \f[
   \mathbf{A}^{\dagger}\phi^{\dagger} =
   \frac{1}{k}\mathbf{B}^{\dagger}\phi^{\dagger}\:,
 * \f]
 */
template<class T>
void Fission_Matrix_Acceleration_Impl<T>::initialize()
{
    typedef Eigenvalue_Solver<T>               Solver_t;
    typedef typename Solver_t::External_Source Source_t;

    // make a "null" external source to pass to the solver
    Teuchos::RCP<const Source_t> null_source;
    CHECK(null_source.is_null());

    // make an eigenvalue solver
    Teuchos::RCP<Solver_t> eigensolver =
        Teuchos::rcp_dynamic_cast<Solver_t>(
            SpnSolverBuilder::build("eigenvalue", b_db));

    // do the adjoint solve first so that the operators in the linear system
    // are set back to forward (not transposed) when the solves are complete

    // >>> ADJOINT SOLVE

    // setup the solver for the adjoint solve
    eigensolver->setup(b_mat, b_mesh, b_indexer, b_gdata, d_system, true);

    // solve the adjoint problem
    eigensolver->solve(null_source);

    // write the results into the adjoint state
    eigensolver->write_state(*b_adjoint);

    // >>> FORWARD SOLVE

    // setup the solver for the forward solve
    eigensolver->setup(b_mat, b_mesh, b_indexer, b_gdata, d_system, false);

    // solve the forward problem
    eigensolver->solve(null_source);

    // write the results into the forward state
    eigensolver->write_state(*b_forward);
}

} // end namespace profugus

#endif // mc_Fission_Matrix_Acceleration_t_hh

//---------------------------------------------------------------------------//
//                 end of Fission_Matrix_Acceleration.t.hh
//---------------------------------------------------------------------------//

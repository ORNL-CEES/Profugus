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

#include <string>

#include "harness/Soft_Equivalence.hh"
#include "spn/Dimensions.hh"
#include "spn/SpnSolverBuilder.hh"
#include "solvers/LinearSolverBuilder.hh"
#include "solvers/PreconditionerBuilder.hh"
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

    // make the linear system for this problem
    d_system = Teuchos::rcp(
        new Linear_System_FV<T>(b_db, d_dim, b_mat, b_mesh, b_indexer, b_gdata));

    // make the adjoint state
    d_adjoint = VectorTraits<T>::build_vector(d_system->get_Map());

    // make the matrices (A,B) for the SPN problem, Ap = (1/k)Bp
    d_system->build_Matrix();
    d_system->build_fission_matrix();

    ENSURE(!b_mesh.is_null());
    ENSURE(!b_indexer.is_null());
    ENSURE(!b_gdata.is_null());
    ENSURE(!b_mat.is_null());
    ENSURE(!d_dim.is_null());
    ENSURE(!d_adjoint.is_null());
    ENSURE(!d_system.is_null());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Initialize the acceleration at the beginning of the K-code solve.
 *
 * This call invokes one eigenvalue solve:
 * \f[
   \mathbf{A}^{\dagger}\phi^{\dagger} =
   \frac{1}{k}\mathbf{B}^{\dagger}\phi^{\dagger}\:,
 * \f]
 */
template<class T>
void Fission_Matrix_Acceleration_Impl<T>::initialize(RCP_ParameterList mc_db)
{
    typedef Eigenvalue_Solver<T>               Solver_t;
    typedef typename Solver_t::External_Source Source_t;

    REQUIRE(mc_db->isSublist("fission_matrix_db"));

    // setup linear solver settings
    solver_db(Teuchos::sublist(mc_db, "fission_matrix_db"));

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

    // store the eigenvalue
    d_keff = eigensolver->get_eigenvalue();
    CHECK(d_keff > 0.0);

    // get the eigenvector
    auto eigenvector =
        VectorTraits<T>::get_data(eigensolver->get_eigenvector());
    auto adjoint =
        VectorTraits<T>::get_data_nonconst(d_adjoint);
    CHECK(eigenvector.size() == adjoint.size());

    // copy local storage
    adjoint.assign(eigenvector);

    // set the system back to forward
    d_system->set_adjoint(false);

    // make the ShiftedOperator
    d_operator = Teuchos::rcp(new ShiftedOperator_t);
    d_operator->set_operator(d_system->get_Operator());
    d_operator->set_rhs_operator(d_system->get_fission_matrix());
    d_operator->set_shift(1.0 / d_keff);

    // build the linear solver
    d_solver = LinearSolverBuilder<T>::build_solver(mc_db);
    d_solver->set_operator(d_operator);

    // build the preconditioner
    auto preconditioner = PreconditionerBuilder<T>::build_preconditioner(
        d_system->get_Operator(), mc_db);

    // set the preconditioner
    d_solver->set_preconditioner(preconditioner);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Start cycle initialization with fission container at \f$l\f$.
 *
 * Calculate beginning of cycle fission density for use in acceleration at end
 * of cycle.
 */
template<class T>
void Fission_Matrix_Acceleration_Impl<T>::start_cycle(
    const Fission_Site_Container &f)
{
}

//---------------------------------------------------------------------------//
/*!
 * \brief End cycle acceleration to create fission container at \f$l+1\f$.
 *
 * Build new fission source based on SPN acceleration.
 */
template<class T>
void Fission_Matrix_Acceleration_Impl<T>::end_cycle(
    Fission_Site_Container &f)
{
}

//---------------------------------------------------------------------------//
// PRIVATE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Default solver functions.
 *
 * The default settings for the solver is
 *
 * \verbatim

   solver_type: "stratimikos"
   Preconditioner: "ml"
   Stratimikos:
        Linear Solver Type: "Belos"
        Preconditioner Type: "None"

   \endverbatim
 *
 * The Stratimikos solver should \b not define a preconditioner.
 */
template<class T>
void Fission_Matrix_Acceleration_Impl<T>::solver_db(
    RCP_ParameterList mc_db)
{
    using std::string;

    // set the default solver type
    mc_db->get("solver_type", string("stratimikos"));

    // set the default preconditioner
    mc_db->get("Preconditioner", string("ml"));

    // setup the stratimikos sublist

}

} // end namespace profugus

#endif // mc_Fission_Matrix_Acceleration_t_hh

//---------------------------------------------------------------------------//
//                 end of Fission_Matrix_Acceleration.t.hh
//---------------------------------------------------------------------------//

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
#include <sstream>
#include <algorithm>
#include <vector>
#include <cmath>

#include "Teuchos_XMLParameterListHelpers.hpp"
#include "AnasaziOperatorTraits.hpp"

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

    // make the linear system for this problem
    d_system = Teuchos::rcp(
        new Linear_System_FV<T>(b_db, d_dim, b_mat, b_mesh, b_indexer, b_gdata));

    // make the adjoint state
    d_adjoint = VTraits::build_vector(d_system->get_Map());

    // build the work vector and solution vectors
    d_work = VTraits::build_vector(d_system->get_Map());
    d_g    = VTraits::build_vector(d_system->get_Map());

    // make the matrices (A,B) for the SPN problem, Ap = (1/k)Bp
    d_system->build_Matrix();
    d_system->build_fission_matrix();

    // make the external source that will be used to calculate B\phi for the
    // acceleration equation
    d_q = std::make_shared<Isotropic_Source>(b_mesh->num_cells());

    // allocate space for the external source shapes, ids, and fields
    d_q_field.resize(b_mesh->num_cells());
    d_q_shapes.resize(b_mesh->num_cells());

    // number of groups
    int num_groups = b_mat->xs().num_groups();

    // build the source shapes (chi) [The SPN problem is setup so that the
    // material ids go from [0,N)]
    const auto &xs = b_mat->xs();
    for (int mid = 0; mid < xs.num_mat(); ++mid)
    {
        CHECK(xs.has(mid));

        // get chi (they can be zero)
        const auto &chi = xs.vector(mid, XS::CHI);
        CHECK(chi.length() == num_groups);
        CHECK(d_q_shapes[mid].empty());

        // get the pointer to the underlying data
        const auto *chi_p = chi.values();

        // store chi
        d_q_shapes[mid].insert(d_q_shapes[mid].end(),
                               chi_p, chi_p + num_groups);
    }

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

    // get the fission matrix db
    RCP_ParameterList fmdb = Teuchos::sublist(mc_db, "fission_matrix_db");

    // setup linear solver settings
    solver_db(fmdb);

    // make a "null" external source to pass to the solver
    Teuchos::RCP<const Source_t> null_source;
    CHECK(null_source.is_null());

    // make an eigenvalue solver
    Teuchos::RCP<Solver_t> eigensolver =
        Teuchos::rcp_dynamic_cast<Solver_t>(
            SpnSolverBuilder::build("eigenvalue", b_db));

    // >>> ADJOINT SOLVE

    // setup the solver for the adjoint solve
    eigensolver->setup(b_mat, b_mesh, b_indexer, b_gdata, d_system, true);

    // solve the adjoint problem
    eigensolver->solve(null_source);

    // store the eigenvalue
    d_keff = eigensolver->get_eigenvalue();
    CHECK(d_keff > 0.0);

    // get the eigenvector
    auto eigenvector = VTraits::get_data(eigensolver->get_eigenvector());
    auto adjoint     = VTraits::get_data_nonconst(d_adjoint);
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
    d_solver = LinearSolverBuilder<T>::build_solver(fmdb);
    d_solver->set_operator(d_operator);

    // build the preconditioner
    auto preconditioner = PreconditionerBuilder<T>::build_preconditioner(
        d_system->get_Operator(), fmdb);

    // set the preconditioner
    if (!preconditioner.is_null())
    {
        d_solver->set_preconditioner(preconditioner);
    }
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
    double                        k_l,
    const Fission_Site_Container &f)
{
    REQUIRE(k_l > 0.0);

    // store the beginning-of-cycle eigenvalue
    d_k_l = k_l;

    // store B\phi^l at the beginning of the cycle
    auto Bphi_l = build_Bphi(f);
    CHECK(VTraits::local_length(Bphi_l) == VTraits::local_length(d_work));

    // now deep copy into local space
    ATraits::MvAddMv(1.0, *Bphi_l, 0.0, *Bphi_l, *d_work);
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
    typedef Anasazi::OperatorTraits<double, MultiVector_t, Operator_t> OT;

    // calculate B\phi^l at l+1/2 (end of MC cycle)
    auto Bphi = build_Bphi(f);
    CHECK(VTraits::local_length(Bphi) == VTraits::local_length(d_work));

    // make some useful name references
    RCP_Vector Bphi_l = d_work;
    RCP_Vector rhs    = Bphi;

    // now calculate k for the solvability condition to the correction
    // equation
    std::vector<double> numerator(1, 0.0), denominator(1, 0.0);
    ATraits::MvDot(*Bphi,   *d_adjoint, numerator);
    ATraits::MvDot(*Bphi_l, *d_adjoint, denominator);
    CHECK(denominator[0] != 0.0);

    // calculate the solvability k
    double k = d_k_l * numerator[0] / denominator[0];
    CHECK(k != 0.0);

    // build the RHS vector
    ATraits::MvAddMv(1.0/k, *Bphi, -1.0/d_k_l, *Bphi_l, *rhs);

    // solve the system
    d_solver->solve(d_g, rhs);

    // ensure orthagonality with adjoint vector, (x*, g) = 0 by applying a
    // correction (the preconditioned solve may not preserve this)
    ATraits::MvDot(*d_adjoint, *d_g, numerator);
    ATraits::MvDot(*d_adjoint, *d_adjoint, denominator);
    CHECK(denominator[0] != 0.0);

    // do the correction
    // ATraits::MvAddMv(1.0, *d_g, -numerator[0]/denominator[0], *d_adjoint, *d_g);

    // calculate the residual
    OT::Apply(*d_operator, *d_g, *d_work);
    ATraits::MvAddMv(1.0, *rhs, -1.0, *d_work, *d_work);
    ATraits::MvNorm(*d_work, numerator);
    std::cout << numerator[0] << std::endl;

#ifdef REMEMBER_ON
    std::vector<double> ortho(1, 0.0);
    ATraits::MvDot(*d_adjoint, *d_g, ortho);
    // VALIDATE(std::fabs(ortho[0] < 1.0e-8),
    //         "Orthogonality of correction fails versus adjoint eigenvector, "
    //         "(x*,g) = " << ortho[0] << " which is > 1.0e-8");
    std::cout << ortho[0] << std::endl;
#endif
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

    // get user-specified tolerances and iterations that we will over-ride
    // later
    double tol     = mc_db->get("tolerance", 1.0e-6);
    double max_itr = mc_db->get("max_itr", 500);

    // set the defaults for the linear solver
    std::ostringstream m;
    m << "<ParameterList name='fission_matrix_db'>\n"
      << "<Parameter name='Preconditioner' type='string' value='ml'/>\n"
      << "<Parameter name='solver_type' type='string' value='stratimikos'/>\n"
      << "<ParameterList name='Stratimikos'>\n"
      << " <Parameter name='Linear Solver Type' type='string' value='Belos'/>\n"
      << " <Parameter name='Preconditioner Type' type='string' value='None'/>\n"
      << "</ParameterList>\n"
      << "</ParameterList>";
    const string pldefault(m.str());

    // Convert string to a Teuchos PL
    RCP_ParameterList default_pl =
        Teuchos::getParametersFromXmlString(pldefault);

    // Insert defaults into pl, leaving existing values in tact
    mc_db->setParametersNotAlreadySet(*default_pl);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build the Bphi mat-vec product from the fission sites.
 */
template<class T>
auto Fission_Matrix_Acceleration_Impl<T>::build_Bphi(
    const Fission_Site_Container &f) -> RCP_Vector
{
    REQUIRE(d_q_field.size() == b_mesh->num_cells());

    // initialize the source field
    std::fill(d_q_field.begin(), d_q_field.end(), 0.0);

    // dimension vector for each cell
    Mesh::Dim_Vector ijk;

    // loop through the fission sites and add them up in each mesh cell
    for (const auto &site : f)
    {
        CHECK(b_mesh->find_cell(site.r, ijk));

        // find the cell containing the site
        b_mesh->find_cell(site.r, ijk);
        CHECK(b_mesh->convert(ijk[0], ijk[1], ijk[2]) < d_q_field.size());

        // add the site to the field
        ++d_q_field[b_mesh->convert(ijk[0], ijk[1], ijk[2])];
    }

    // normalize by volume
    for (int cell = 0; cell < b_mesh->num_cells(); ++cell)
    {
        d_q_field[cell] /= b_mesh->volume(cell);
    }

    // setup the source
    d_q->set(b_mat->matids(), d_q_shapes, d_q_field);

    // use the linear system to calculate a right-hand side vector from this
    // source, which is B\phi in SPN space
    d_system->build_RHS(*d_q);

    // return the RHS vector
    ENSURE(!d_system->get_RHS().is_null());
    return d_system->get_RHS();
}

} // end namespace profugus

#endif // mc_Fission_Matrix_Acceleration_t_hh

//---------------------------------------------------------------------------//
//                 end of Fission_Matrix_Acceleration.t.hh
//---------------------------------------------------------------------------//

//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   MC/mc/Fission_Matrix_Solver.t.hh
 * \author Thomas M. Evans
 * \date   Mon Dec 08 17:18:34 2014
 * \brief  Fission_Matrix_Solver template member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_mc_Fission_Matrix_Solver_t_hh
#define MC_mc_Fission_Matrix_Solver_t_hh

#include <algorithm>
#include <string>
#include <vector>
#include <cmath>
#include <sstream>

#include "Teuchos_XMLParameterListHelpers.hpp"
#include "AnasaziOperatorTraits.hpp"

#include "harness/DBC.hh"
#include "comm/global.hh"
#include "utils/Definitions.hh"
#include "solvers/LinearSolverBuilder.hh"
#include "solvers/PreconditionerBuilder.hh"
#include "Fission_Matrix_Solver.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Build the solver.
 */
template<class Geometry, class T>
Fission_Matrix_Solver<Geometry,T>::Fission_Matrix_Solver(
        RCP_ParameterList fm_db,
        RCP_Mesh          mesh,
        RCP_Indexer       indexer,
        RCP_Mat_DB        mat,
        RCP_Linear_System system,
        SP_Cart_Mesh      global_mesh,
        double            keff)
    : d_mesh(mesh)
    , d_indexer(indexer)
    , d_mat(mat)
    , d_system(system)
    , d_global_mesh(global_mesh)
{
    REQUIRE(!d_mesh.is_null());
    REQUIRE(!d_indexer.is_null());
    REQUIRE(!d_mat.is_null());
    REQUIRE(!d_system.is_null());
    REQUIRE(d_global_mesh);
    REQUIRE(keff > 0.0);
    REQUIRE(d_global_mesh->num_cells() >= d_mesh->num_cells());

    // setup linear solver settings and defaults
    solver_db(fm_db);

    // build the work vector and solution vectors
    d_work = VTraits::build_vector(d_system->get_Map());
    d_g    = VTraits::build_vector(d_system->get_Map());

    // make the external source that will be used to calculate B\phi for the
    // acceleration equation
    d_q = std::make_shared<Isotropic_Source>(d_mesh->num_cells());

    // allocate space for the external source shapes, ids, and fields
    d_q_field.resize(d_mesh->num_cells());
    d_q_shapes.resize(d_mesh->num_cells());

    // make a global q field for calculating the fission density on each
    // domain
    d_global_q_field.resize(d_global_mesh->num_cells());

    // number of groups
    int num_groups = d_mat->xs().num_groups();

    // build the source shapes (chi) [The SPN problem is setup so that the
    // material ids go from [0,N)]
    const auto &xs = d_mat->xs();
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

    // make the ShiftedOperator
    d_operator = Teuchos::rcp(new ShiftedOperator_t);
    d_operator->set_operator(d_system->get_Operator());
    d_operator->set_rhs_operator(d_system->get_fission_matrix());
    d_operator->set_shift(1.0 / keff);

    // build the linear solver
    d_solver = LinearSolverBuilder<T>::build_solver(fm_db);
    d_solver->set_operator(d_operator);

    // build the preconditioner
    auto preconditioner = PreconditionerBuilder<T>::build_preconditioner(
        d_system->get_Operator(), fm_db);

    // set the preconditioner
    if (!preconditioner.is_null())
    {
        d_solver->set_preconditioner(preconditioner);
    }
}

//---------------------------------------------------------------------------//
// PUBLIC MEMBERS
//---------------------------------------------------------------------------//
/*!
 * \brief Set the eigenvectors of the SPN system.
 */
template<class Geometry, class T>
void Fission_Matrix_Solver<Geometry,T>::set_eigenvectors(RCP_Vector forward,
                                                         RCP_Vector adjoint)
{
    REQUIRE(!forward.is_null());
    REQUIRE(!adjoint.is_null());

    d_forward = forward;
    d_adjoint = adjoint;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Set the fission source at the beginning of the cycle.
 */
template<class Geometry, class T>
void Fission_Matrix_Solver<Geometry,T>::set_u_begin(
        const Fission_Site_Container &f,
        double                        k_l)
{
    REQUIRE(k_l > 0.0);

    // store the beginning-of-cycle eigenvalue
    d_k_l = k_l;

    // store B\phi^l at the beginning of the cycle
    auto Bphi_l = build_Bphi(f);
    CHECK(VTraits::local_length(Bphi_l) == VTraits::local_length(d_work));

    // now deep copy into local space
    std::vector<int> index(1, 0);
    ATraits::SetBlock(*Bphi_l, index, *d_work);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Solve for the correction.
 */
template<class Geometry, class T>
void Fission_Matrix_Solver<Geometry,T>::solve(const Fission_Site_Container &f)
{
    typedef Anasazi::OperatorTraits<double, MultiVector_t, Operator_t> OT;

    REQUIRE(!d_forward.is_null());
    REQUIRE(!d_adjoint.is_null());

    // calculate B\phi^l at l+1/2 (end of MC cycle)
    auto Bphi = build_Bphi(f);
    CHECK(VTraits::local_length(Bphi) == VTraits::local_length(d_work));

    // make some useful name references
    RCP_Vector Bphi_l = d_work;
    RCP_Vector rhs    = Bphi;

    // now calculate k for the solvability condition to the correction
    // equation
    std::vector<double> numerator(1, 0.0), denominator(1, 0.0);
    ATraits::MvDot(*Bphi, *d_adjoint, numerator);
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
    ATraits::MvDot(*d_adjoint, *d_forward, denominator);
    CHECK(denominator[0] != 0.0);

    // orthogonalize the correction vector
    ATraits::MvAddMv(1.0, *d_g, -numerator[0]/denominator[0], *d_forward, *d_g);

#ifdef REMEMBER_ON
    // Check orthogonality
    std::vector<double> ortho(1, 0.0);
    ATraits::MvDot(*d_adjoint, *d_g, ortho);
    VALIDATE(std::fabs(ortho[0] < 1.0e-8),
             "Orthogonality of correction fails versus adjoint eigenvector, "
             << "(x*,g) = " << ortho[0] << " which is > 1.0e-8");

    // check the residual
    ATraits::MvNorm(*rhs, denominator);
    OT::Apply(*d_operator, *d_g, *d_work);
    ATraits::MvAddMv(1.0, *rhs, -1.0, *d_work, *d_work);
    ATraits::MvNorm(*d_work, numerator);
    VALIDATE(numerator[0]/denominator[0] < 1.0e-3,
             "Residual does not preserve the original operator solution, "
             << " r = " << numerator[0]);
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
template<class Geometry, class T>
void Fission_Matrix_Solver<Geometry,T>::solver_db(RCP_ParameterList fm_db)
{
    using std::string;

    // get user-specified tolerances and iterations that we will over-ride
    // later
    double tol     = fm_db->get("tolerance", 1.0e-6);
    double max_itr = fm_db->get("max_itr", 500);

    // set the defaults for the linear solver
    std::ostringstream m;
    m << "<ParameterList name='acceleration'>\n"
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
    fm_db->setParametersNotAlreadySet(*default_pl);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build the Bphi mat-vec product from the fission sites.
 */
template<class Geometry, class T>
auto Fission_Matrix_Solver<Geometry,T>::build_Bphi(
    const Fission_Site_Container &f) -> RCP_Vector
{
    using def::I; using def::J; using def::K;

    REQUIRE(d_q_field.size() == d_mesh->num_cells());
    REQUIRE(d_global_q_field.size() == d_global_mesh->num_cells());

    // initialize the source fields
    std::fill(d_global_q_field.begin(), d_global_q_field.end(), 0.0);

    // dimension vector for each cell
    Mesh::Dim_Vector ijk;

    // loop through the global fission sites and add them up in each global
    // mesh cell
    for (const auto &site : f)
    {
        CHECK(d_global_mesh->find(site.r, ijk));

        // find the cell containing the site
        d_global_mesh->find(site.r, ijk);
        CHECK(d_global_mesh->index(ijk[0], ijk[1], ijk[2]) <
              d_global_q_field.size());

        // add the site to the global field
        ++d_global_q_field[d_global_mesh->index(ijk[0], ijk[1], ijk[2])];
    }

    // reduce the global field
    profugus::global_sum(d_global_q_field.data(), d_global_mesh->num_cells());

    // map the global field to the local field and normalize by volume
    for (int k = 0; k < d_mesh->num_cells_dim(K); ++k)
    {
        for (int j = 0; j < d_mesh->num_cells_dim(J); ++j)
        {
            for (int i = 0; i < d_mesh->num_cells_dim(I); ++i)
            {
                // get the global index for this cell
                int global = d_indexer->l2g(i, j, k);
                int local  = d_indexer->l2l(i, j, k);
                CHECK(local == d_mesh->convert(i, j, k));
                CHECK(global < d_global_q_field.size());

                // get the volume
                double vol = d_mesh->volume(local);
                CHECK(vol == d_global_mesh->volume(global));
                CHECK(vol > 0.0);

                // make the local fission source density
                d_q_field[local] = d_global_q_field[global] / vol;
            }
        }
    }

    // setup the source
    d_q->set(d_mat->matids(), d_q_shapes, d_q_field);

    // use the linear system to calculate a right-hand side vector from this
    // source, which is B\phi in SPN space
    d_system->build_RHS(*d_q);

    // return the RHS vector
    ENSURE(!d_system->get_RHS().is_null());
    return d_system->get_RHS();
}

} // end namespace profugus

#endif // MC_mc_Fission_Matrix_Solver_t_hh

//---------------------------------------------------------------------------//
//                 end of Fission_Matrix_Solver.t.hh
//---------------------------------------------------------------------------//

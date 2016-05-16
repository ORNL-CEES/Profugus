//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   MC/mc/Fission_Matrix_Acceleration.t.hh
 * \author Thomas M. Evans
 * \date   Thu Nov 13 19:53:21 2014
 * \brief  Fission_Matrix_Acceleration template member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_mc_Fission_Matrix_Acceleration_t_hh
#define MC_mc_Fission_Matrix_Acceleration_t_hh

#include <cmath>
#include <algorithm>

#include "harness/DBC.hh"
#include "harness/Soft_Equivalence.hh"
#include "comm/global.hh"
#include "spn/Dimensions.hh"
#include "spn/SpnSolverBuilder.hh"
#include "spn/Linear_System_FV.hh"
#include "spn/Moment_Coefficients.hh"
#include "spn/Eigenvalue_Solver.hh"
#include "Global_RNG.hh"
#include "Fission_Matrix_Acceleration.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template<class Geometry, class T>
Fission_Matrix_Acceleration_Impl<Geometry,T>::Fission_Matrix_Acceleration_Impl()
    : d_cycle_ctr(0)
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
template<class Geometry, class T>
void Fission_Matrix_Acceleration_Impl<Geometry,T>::build_problem(
    const Problem_Builder_t &builder)
{
    // get the problem database from the problem-builder
    b_db = builder.problem_db();

    // get the mesh objects from the builder
    b_mesh    = builder.mesh();
    b_indexer = builder.indexer();
    b_gdata   = builder.global_data();

    // build a global SPN mesh (for mapping fission sites to partitioned SPN
    // fields)
    b_global_mesh = std::make_shared<Cartesian_Mesh>(
        b_gdata->edges(0), b_gdata->edges(1), b_gdata->edges(2));
    CHECK(b_global_mesh->num_cells() == b_gdata->num_cells());

    // get the material database from the problem builder
    b_mat = builder.mat_db();

    // build the problem dimensions
    d_dim = Teuchos::rcp(new profugus::Dimensions(b_db->get("SPn_order", 1)));

    // make the linear system for this problem
    d_system = Teuchos::rcp(
        new Linear_System_FV<T>(b_db, d_dim, b_mat, b_mesh, b_indexer, b_gdata));

    // make the forward and adjoint eigenvectors
    d_forward = VTraits::build_vector(d_system->get_Map());
    d_adjoint = VTraits::build_vector(d_system->get_Map());

    // make the matrices (A,B) for the SPN problem, Ap = (1/k)Bp
    d_system->build_Matrix();
    d_system->build_fission_matrix();

    // build the weight corrections in the FM mesh
    d_nu.resize(b_global_mesh->num_cells());

    ENSURE(!b_mesh.is_null());
    ENSURE(!b_indexer.is_null());
    ENSURE(!b_gdata.is_null());
    ENSURE(!b_mat.is_null());
    ENSURE(!d_dim.is_null());
    ENSURE(!d_adjoint.is_null());
    ENSURE(!d_forward.is_null());
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
template<class Geometry, class T>
void Fission_Matrix_Acceleration_Impl<Geometry,T>::initialize(
        RCP_ParameterList mc_db)
{
    typedef Eigenvalue_Solver<T>               Solver_t;
    typedef typename Solver_t::External_Source Source_t;

    REQUIRE(mc_db->isSublist("fission_matrix_db"));
    REQUIRE(mc_db->sublist("fission_matrix_db").isSublist("acceleration"));

    // get the fission matrix db->acceleration
    RCP_ParameterList fmdb = Teuchos::sublist(
        Teuchos::sublist(mc_db, "fission_matrix_db"), "acceleration");

    // set the damping (defaults to 1.0)
    d_beta = fmdb->get("damping", 1.0);

    // set the begin/end cycle
    d_cycle_ctr   = 0;
    d_cycle_begin = fmdb->get("begin_cycle", 1);
    d_cycle_end   = fmdb->get(
        "end_cycle", mc_db->template get<int>("num_inactive_cycles", 10));
    d_accelerate  = false;

    // total number (inactive + active) of cycles
    int num_cycles = mc_db->template get<int>("num_cycles", 50);

    // setup to generate an initial fission source distribution
    d_initial_source = fmdb->get("initial_source", false);
    VALIDATE(d_cycle_begin < num_cycles || d_initial_source,
             "No valid acceleration, the beginning cycle of "
             << d_cycle_begin << " is >= the number of cycles, "
             << num_cycles << " and no initial source distribution requested.");

    // determine if we should update the particle population
    b_update_Np = fmdb->get("update_Np", false);

    // make a "null" external source to pass to the solver
    Teuchos::RCP<const Source_t> null_source;
    CHECK(null_source.is_null());

    // make an eigenvalue solver, use the database from the SPN problem
    Teuchos::RCP<Solver_t> eigensolver =
        Teuchos::rcp_dynamic_cast<Solver_t>(
            SpnSolverBuilder::build("eigenvalue", b_db));

    // >>> FORWARD SOLVE

    // setup the solver for the forward solve (puts the operators back to
    // forward)
    eigensolver->setup(b_mat, b_mesh, b_indexer, b_gdata, d_system, false);

    // solve the adjoint problem
    eigensolver->solve(null_source);

    // store the eigenvalue
    d_keff = eigensolver->get_eigenvalue();
    CHECK(d_keff > 0.0);

    // assign the forward eigenvector
    std::vector<int> index(1, 0);
    ATraits::SetBlock(*eigensolver->get_eigenvector(), index, *d_forward);

    // >>> ADJOINT SOLVE and BUILD FISSION_MATRIX_SOLVER

    // only do this if we eventually plan to do fission matrix acceleration
    if (d_cycle_begin < num_cycles)
    {
        // setup the solver for the adjoint solve
        eigensolver->setup(b_mat, b_mesh, b_indexer, b_gdata, d_system, true);

        // solve the adjoint problem
        eigensolver->solve(null_source);
        CHECK(profugus::soft_equiv(
                  eigensolver->get_eigenvalue(), d_keff, 1.0e-3));

        // get the eigenvector and assign it
        ATraits::SetBlock(*eigensolver->get_eigenvector(), index, *d_adjoint);

        // >>> BUILD THE FISSION_MATRIX_SOLVER

        // set the system back to forward
        d_system->set_adjoint(false);

        // build and set the fission matrix solver
        d_fm_solver = Teuchos::rcp(
            new FM_Solver_t(
                fmdb, b_mesh, b_indexer, b_mat, d_system, b_global_mesh,
                d_keff));
        d_fm_solver->set_eigenvectors(d_forward, d_adjoint);
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Start cycle initialization with fission container at \f$l\f$.
 *
 * Calculate beginning of cycle fission density for use in acceleration at end
 * of cycle.
 */
template<class Geometry, class T>
void Fission_Matrix_Acceleration_Impl<Geometry,T>::start_cycle(
    double                        k_l,
    const Fission_Site_Container &f)
{
    REQUIRE(!d_fm_solver.is_null());

    d_accelerate = false;

   if (d_cycle_ctr >= d_cycle_begin && d_cycle_ctr < d_cycle_end)
    {
        REQUIRE(!d_fm_solver.is_null());
        d_fm_solver->set_u_begin(f, k_l);
        d_accelerate = true;
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief End cycle acceleration to create fission container at \f$l+1\f$.
 *
 * Build new fission source based on SPN acceleration.
 */
template<class Geometry, class T>
void Fission_Matrix_Acceleration_Impl<Geometry,T>::end_cycle(
    Fission_Site_Container &f)
{
    // increment
    ++d_cycle_ctr;

    // return if not accelerating
    if (!d_accelerate)
    {
        // calculate and store the number of fission sites
        b_current_Np = f.size();
        profugus::global_sum(&b_current_Np, 1);
        return;
    }

    REQUIRE(!d_fm_solver.is_null());
    REQUIRE(Global_RNG::d_rng.assigned());
    REQUIRE(d_nu.size() == b_global_mesh->num_cells());

    // initialize d_nu multiplication and norm
    std::fill(d_nu.begin(), d_nu.end(), 0.0);
    double l2_norm = 0.0;
    double max_c   = 0.0;

    // solve the acceleration equation
    d_fm_solver->solve(f);
    CHECK(VTraits::local_length(d_fm_solver->get_g()) ==
          VTraits::local_length(d_forward));

    // push back the latest solve iteration count
    d_iterations.push_back(d_fm_solver->iterations());

    // convert the correction vector from SPN space to fission sites
    convert_g(d_fm_solver->get_g());

    // get the local fission density at l+1/2
    auto fis_den = d_fm_solver->current_f();
    CHECK(fis_den.size() == b_mesh->num_cells());

    // number of local cells
    int Nc[3] = {b_mesh->num_cells_dim(0),
                 b_mesh->num_cells_dim(1),
                 b_mesh->num_cells_dim(2)};
    CHECK(Nc[0]*Nc[1]*Nc[2] == b_mesh->num_cells());

    // build the multiplicative correction
    for (int k = 0; k < Nc[2]; ++k)
    {
        for (int j = 0; j < Nc[1]; ++j)
        {
            for (int i = 0; i < Nc[0]; ++i)
            {
                // get the local and global cell indices
                int global = b_indexer->l2g(i, j, k);
                int local  = b_indexer->l2l(i, j, k);
                CHECK(global < d_nu.size());
                CHECK(local == b_mesh->convert(i, j, k));
                CHECK(fis_den[local] >= 0.0);

                // get the l2 norm of the correction (locally)
                l2_norm += d_nu[global] * d_nu[global];

                // only make multiplicative correction in cells with fission
                // density

                // NOTE: because the correction is multiplicative, there could
                // be no fission sites in a cell due to sampling, but a
                // nonzero correction, there should only be a zero correction
                // in regions with no fissionable material, of course there
                // will be no MC fission sites there.  For now, we do not
                // correct in cells that have no MC fission sites, even if
                // there is a correction [we will look at this later]
                if (fis_den[local] > 0.0)
                {
                    // d_nu is currently in neutrons/cc so we don't need to
                    // multiply the fission density by volume->the volumes
                    // cancel out
                    d_nu[global] = 1.0 + d_beta * d_nu[global] / fis_den[local];

                    // calculate the maximum correction
                    max_c = std::max(max_c, d_nu[global]);
                }
            }
        }
    }

    // add up the norm over every domain
    profugus::global_sum(&l2_norm, 1);

    // get the maximum correction across all domains
    profugus::global_max(&max_c, 1);

    // calculate the l2 norm and push it back onto the stack
    l2_norm = std::sqrt(l2_norm);
    d_norms.push_back(l2_norm);

    // push the maxium correction onto the stack
    d_correction.push_back(max_c);

    // reduce d_nu so that every domain has the correction
    profugus::global_sum(d_nu.data(), d_nu.size());

    // make a new fission site container
    Fission_Site_Container nf;
    CHECK(nf.empty());

    // dimension vector for each cell
    Mesh::Dim_Vector ijk;

    // get the global RNG (by reference)
    auto rng = Global_RNG::d_rng;

    // correct the fission sites
    while (!f.empty())
    {
        // get the fission site of the back
        const auto &site = f.back();
        CHECK(b_global_mesh->find(site.r, ijk));

        // find the cell containing the site
        b_global_mesh->find(site.r, ijk);
        int cell = b_global_mesh->index(ijk[0], ijk[1], ijk[2]);
        CHECK(cell < d_nu.size());
        CHECK(cell == b_indexer->g2g(ijk[0], ijk[1], ijk[2]));

        // set d_nu[cell] to 0.0 if it is negative (we should count these)!
        if (d_nu[cell] < 0.0) d_nu[cell] = 0.0;

        // sample to determine the number of sites at this location
        int n = d_nu[cell];
        CHECK(static_cast<double>(n) + 1.0 - d_nu[cell] >= 0.0 &&
              static_cast<double>(n) + 1.0 - d_nu[cell] <= 1.0);

        // with probability n+1-nu there will be n sites; with propability
        // nu-n there will be n+1 sites
        if (rng.ran() < d_nu[cell] - static_cast<double>(n))
        {
            ++n;
        }

        // add n sites to the new fission bank
        for (int m = 0; m < n; ++m)
        {
            nf.push_back(site);
        }

        // pop the fission site from the container
        f.pop_back();
    }

    // swap the new and old containers
    std::swap(nf, f);

    // calculate and store the number of fission sites
    b_current_Np = f.size();
    profugus::global_sum(&b_current_Np, 1);
    d_Np.push_back(b_current_Np);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build the initial fission source.
 */
template<class Geometry, class T>
void Fission_Matrix_Acceleration_Impl<Geometry,T>::build_initial_source(
    Fission_Source_t &source)
{
    // if we aren't using an initial source, then do "standard" fission source
    // initialization and return
    if (!d_initial_source)
    {
        source.build_initial_source();
        return;
    }

    REQUIRE(b_global_mesh);

    // otherwise use the forward eigenvector to sample the inital fission
    // distribution

    // temporarily use d_nu to hold the fission distribution from the SPN
    // solution
    std::fill(d_nu.begin(), d_nu.end(), 0.0);
    convert_eigenvector(d_forward);

    // reduce d_nu so that every domain has the fission density
    profugus::global_sum(d_nu.data(), d_nu.size());

    // build the initial source
    source.build_initial_source(b_global_mesh,
                                Teuchos::ArrayView<const double>(d_nu));
}

//---------------------------------------------------------------------------//
#ifdef USE_HDF5
/*!
 * \brief Write diagnostics to HDF5 file.
 */
template<class Geometry, class T>
void Fission_Matrix_Acceleration_Impl<Geometry,T>::diagnostics(
    Serial_HDF5_Writer &writer) const
{
    writer.begin_group("fission_mat_acceleration");

    // write the eigenvalue of the SPN problem
    writer.write("SPn_k", d_keff);

    // only output data if acceleration has been performed
    if (!d_iterations.empty())
    {
        // write the iteration history
        writer.write("iterations", d_iterations);

        // write the correction l2 norm
        writer.write("correction_L2", d_norms);

        // write the L_inf norm of the multiplicative correction
        writer.write("mult_correction_Linf", d_correction);

        // write the history of post-corrected fission site numbers
        writer.write("num_fission_sites", d_Np);
    }

    writer.end_group();
}
#endif

//---------------------------------------------------------------------------//
// PRIVATE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Convert g to fission density.
 */
template<class Geometry, class T>
void Fission_Matrix_Acceleration_Impl<Geometry,T>::convert_g(
        RCP_Const_Vector gc)
{
    REQUIRE(b_mesh->num_cells() * d_system->get_dims()->num_equations()
            * b_mat->xs().num_groups() <=
            VectorTraits<T>::local_length(gc));
    REQUIRE(d_nu.size() == b_global_mesh->num_cells());

    // SPN moments
    double u_m[4] = {0.0, 0.0, 0.0, 0.0};

    // get cross sections from the database
    const auto &xs = b_mat->xs();

    // number of equations (moments) in the SPN solution
    int N = d_system->get_dims()->num_equations();

    // number of groups
    int Ng = xs.num_groups();

    // number of local cells
    int Nc[3] = {b_mesh->num_cells_dim(0),
                 b_mesh->num_cells_dim(1),
                 b_mesh->num_cells_dim(2)};
    CHECK(Nc[0]*Nc[1]*Nc[2] == b_mesh->num_cells());

    // the u vector is ordered cell->moments->group
    Teuchos::ArrayRCP<const double> gcv = VectorTraits<T>::get_data(gc);

    // loop over cells on this domain
    for (int k = 0; k < Nc[2]; ++k)
    {
        for (int j = 0; j < Nc[1]; ++j)
        {
            for (int i = 0; i < Nc[0]; ++i)
            {
                // get the local and global cell indices
                int global = b_indexer->l2g(i, j, k);
                int local  = b_indexer->l2l(i, j, k);
                CHECK(global < d_nu.size());
                CHECK(local == b_mesh->convert(i, j, k));

                // get nu-sigma_f for this cell
                const auto &nu_sigf = xs.vector(b_mat->matid(local),
                                                XS::NU_SIG_F);
                CHECK(nu_sigf.length() == Ng);

                // loop over groups
                for (int g = 0; g < Ng; ++g)
                {
                    // loop over moments
                    for (int n = 0; n < N; ++n)
                    {
                        CHECK(d_system->index(g, n, local) < gcv.size());

                        // get the moment from the solution vector
                        u_m[n] = gcv[d_system->index(g, n, local)];
                    }

                    // do f^T * phi_0
                    d_nu[global] +=  Moment_Coefficients::u_to_phi(
                        u_m[0], u_m[1], u_m[2], u_m[3]) * nu_sigf[g];
                }
            }
        }
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Convert eigenvector to fission density.
 */
template<class Geometry, class T>
void Fission_Matrix_Acceleration_Impl<Geometry,T>::convert_eigenvector(
    RCP_Const_Vector ev)
{
    REQUIRE(b_mesh->num_cells() * d_system->get_dims()->num_equations()
            * b_mat->xs().num_groups() <=
            VectorTraits<T>::local_length(ev));
    REQUIRE(d_nu.size() == b_global_mesh->num_cells());

    // SPN moments
    double u_m[4] = {0.0, 0.0, 0.0, 0.0};

    // get cross sections from the database
    const auto &xs = b_mat->xs();

    // number of equations (moments) in the SPN solution
    int N = d_system->get_dims()->num_equations();

    // number of groups
    int Ng = xs.num_groups();

    // number of local cells
    int Nc[3] = {b_mesh->num_cells_dim(0),
                 b_mesh->num_cells_dim(1),
                 b_mesh->num_cells_dim(2)};
    CHECK(Nc[0]*Nc[1]*Nc[2] == b_mesh->num_cells());

    // the u vector is ordered cell->moments->group
    auto v = VTraits::get_data(ev);

    // >>> Determine normalization of eigenvector

    double norm[] = {0.0, 0.0};
    double flux   = 0.0;

    // loop over cells on this domain
    for (int cell = 0; cell < b_mesh->num_cells(); ++cell)
    {
        // loop over groups
        for (int g = 0; g < Ng; ++g)
        {
            // loop over moments
            for (int n = 0; n < N; ++n)
            {
                CHECK(d_system->index(g, n, cell) < v.size());

                // get the moment from the solution vector
                u_m[n] = v[d_system->index(g, n, cell)];
            }

            // get the flux for this cell/group, store in d_nu
            flux =  Moment_Coefficients::u_to_phi(
                u_m[0], u_m[1], u_m[2], u_m[3]);

            // calculate the norms
            norm[0] += flux * flux; // for normalization
            norm[1] += flux;        // for sign
        }
    }

    // do a global reduction on the normalization
    profugus::global_sum(norm, 2);
    CHECK(norm[0] > 0.0);

    // calculate the normalization
    double norm_f = (1.0 / std::sqrt(norm[0])) * (std::fabs(norm[1]) / norm[1]);

    // >>> Calculate the fission density

    // loop over cells on this domain
    for (int k = 0; k < Nc[2]; ++k)
    {
        for (int j = 0; j < Nc[1]; ++j)
        {
            for (int i = 0; i < Nc[0]; ++i)
            {
                // get the local and global cell indices
                int global = b_indexer->l2g(i, j, k);
                int local  = b_indexer->l2l(i, j, k);
                CHECK(global < d_nu.size());
                CHECK(local == b_mesh->convert(i, j, k));

                // get nu-sigma_f for this cell
                const auto &nu_sigf = xs.vector(b_mat->matid(local),
                                                XS::NU_SIG_F);
                CHECK(nu_sigf.length() == Ng);

                // loop over groups
                for (int g = 0; g < Ng; ++g)
                {
                    // loop over moments
                    for (int n = 0; n < N; ++n)
                    {
                        CHECK(d_system->index(g, n, local) < v.size());

                        // get the moment from the solution vector
                        u_m[n] = v[d_system->index(g, n, local)];
                    }

                    // do f^T * phi_0
                    d_nu[global] +=  Moment_Coefficients::u_to_phi(
                        u_m[0], u_m[1], u_m[2], u_m[3]) * norm_f * nu_sigf[g];
                }
            }
        }
    }
}

} // end namespace profugus

#endif // MC_mc_Fission_Matrix_Acceleration_t_hh

//---------------------------------------------------------------------------//
//                 end of Fission_Matrix_Acceleration.t.hh
//---------------------------------------------------------------------------//

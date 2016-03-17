//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/mc/Anderson_Operator.t.hh
 * \author Thomas M. Evans
 * \date   Tue Apr 07 21:05:56 2015
 * \brief  Anderson_Operator template method definitions.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_mc_Anderson_Operator_t_hh
#define MC_mc_Anderson_Operator_t_hh

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>

#include "spn/VectorTraits.hh"
#include "spn/MatrixTraits.hh"
#include "Keff_Tally.hh"
#include "Anderson_Operator.hh"
#include "Global_RNG.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template<class Geometry, class T>
Anderson_Operator<Geometry,T>::Anderson_Operator(
        RCP_PL                   pl,
        SP_Source_Transporter    transporter,
        SP_Fission_Source        fission_source,
        SP_Cart_Mesh             mesh,
        RCP_MAP                  map,
        profugus::Communicator_t set_comm)
    : Base(map)
    , d_pl(pl)
    , d_transporter(transporter)
    , d_source(fission_source)
    , d_mesh(mesh)
    , d_Np(static_cast<double>(fission_source->Np()))
    , d_nodes(profugus::nodes())
    , d_node(profugus::node())
    , d_set_comm(set_comm)
{
    REQUIRE(d_transporter);
    REQUIRE(d_source);
    REQUIRE(d_mesh);

    // Fission site container
    d_fission_sites = d_source->create_fission_site_container();
    CHECK(d_fission_sites);
    CHECK(d_fission_sites->empty());

    // Build the g' vector stored for each iteration (dimensioned only over
    // the mesh, not over the full solution vector)
    profugus::set_internal_comm(d_set_comm);
    RCP_MAP mesh_map = MatrixTraits<T>::build_map(d_mesh->num_cells(),
                                                  d_mesh->num_cells());
    profugus::reset_internal_comm();
    d_gp = VectorTraits<T>::build_vector(mesh_map);
    ENSURE(!d_gp.is_null());
    ENSURE(profugus::nodes() == d_nodes);

    d_use_tally = d_pl->get("use_tally",true);

    d_tallies_built = false;
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Do a transport iteration.
 */
template<class Geometry, class T>
void Anderson_Operator<Geometry,T>::iterate(double k) const
{
    REQUIRE( d_tallier );
    REQUIRE( d_tallier->is_built() && !d_tallier->is_finalized() );
    REQUIRE(d_fission_sites && d_fission_sites->empty());

    d_transporter->set(d_tallier);

    // assign the current source state in the transporter (this will generally
    // be a pass through, but it gives the solver a chance to update
    // quantities if the source changes from cycle-to-cycle)
    d_transporter->assign_source(d_source);

    // set the solver to sample fission sites
    d_transporter->sample_fission_sites(d_fission_sites, k);

    // initialize keff tally to the beginning of the cycle
    REQUIRE( d_tallier );
    d_tallier->begin_cycle();

    // solve the fixed source problem using the transporter
    d_transporter->solve();

    // do end-of-cycle tally processing including global sum Note: this is the
    // total *requested* number of particles, not the actual number of
    // histories. Each processor adjusts the particle weights so that the
    // total weight emitted, summed over all processors, is d_Np.
    d_tallier->end_cycle(d_Np);

}

//---------------------------------------------------------------------------//
/*!
 * \brief Update the fission source.
 */
template<class Geometry, class T>
void Anderson_Operator<Geometry,T>::update_source() const
{
    REQUIRE(!d_fission_sites->empty());

    // build a new source from the fission site distribution
    d_source->build_source(d_fission_sites);

    ENSURE(d_fission_sites);
    ENSURE(d_fission_sites->empty());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build tallies necessary for operator apply
 */
template<class Geometry, class T>
void Anderson_Operator<Geometry,T>::build_tallies()
{
    // Need Tallier with Keff and Fission tallies
    d_tallier = std::make_shared<Tallier_t>();
    d_tallier->set( d_source->geometry(), d_source->physics() );
    auto src_tallier = d_transporter->tallier();
    for( const auto &tally : *src_tallier )
    {
        if( tally->name() == "keff" )
        {
            d_keff_tally = std::dynamic_pointer_cast<Keff_Tally_t>(tally);
            ENSURE( d_keff_tally );
            d_tallier->add_pathlength_tally(d_keff_tally);
        }
    }
    ENSURE( d_tallier->num_pathlength_tallies() == 1 );

    // Build Fission_Tally
    d_fisn_tally = std::make_shared<Fission_Tally_t>(
        d_source->physics());
    d_fisn_tally->set_mesh(d_mesh);
    d_tallier->add_pathlength_tally(d_fisn_tally);
    d_tallier->build();
    ENSURE( d_tallier->is_built() );
    ENSURE( d_tallier->num_tallies() == 2 );
    ENSURE( d_tallier->num_pathlength_tallies() == 2 );

    d_tallies_built = true;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Initialize Anderson solve.
 */
template<class Geometry, class T>
auto Anderson_Operator<Geometry,T>::initialize_Anderson() -> RCP_MV
{
    REQUIRE( d_tallies_built );

    // Make the solution vector
    auto sol_vec = VectorTraits<T>::build_vector(Base::d_domain_map);

    // Get the data from sol_vec (this is the total vector g + k)
    auto v = VectorTraits<T>::get_data_nonconst(sol_vec);

    // Get an ArrayRCP to gp (just over g)
    auto gp = VectorTraits<T>::get_data_nonconst(d_gp);
    CHECK(gp.size() + 1 == v.size());

    // Make gp from the existing cycle fission sites
    restrict(*d_fission_sites, gp());

    // Write gp into v to initialize the solution vector
    std::copy(gp.begin(), gp.end(), v.begin());

    // Store the latest k iterate in the initial solution vector
    v[d_mesh->num_cells()] = d_keff_tally->latest();
    CHECK(v[d_mesh->num_cells()] > 0.0);

    // return v
    return sol_vec;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Finalize the Anderson solver.
 *
 * This sets up the fission source for the continued running of active cycles
 * after Anderson solution.
 *
 * \return k from Anderson solve
 */
template<class Geometry, class T>
double Anderson_Operator<Geometry,T>::finalize_Anderson(const MV &v)
{
    // Number of cells in the grid
    int nc = d_mesh->num_cells();

    // Get ArrayRCP's to the data
    auto in  = VectorTraits<T>::get_data(Teuchos::rcpFromRef(v));
    CHECK(in.size() == nc + 1);

    // Make view of data
    auto g = Teuchos::arrayView(in.get(), nc);
    CHECK(g.size() == nc);

    // get last k from v
    double k = in[nc];
    CHECK(k > 0.0);

    // Apply prolongation P: g->f, to make fission sites for active cycles
    prolongate(g, *d_fission_sites);

    // update the source
    update_source();

    // return the last k iterate
    return k;
}

//---------------------------------------------------------------------------//
// OPERATOR INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief Operator apply for solves.
 */
template<class Geometry, class T>
void Anderson_Operator<Geometry,T>::ApplyImpl(const MV &x, MV &y) const
{
    REQUIRE( d_tallies_built );

    // Number of cells in the grid
    int nc = d_mesh->num_cells();

    // Get ArrayRCP's to the data
    auto in  = VectorTraits<T>::get_data(Teuchos::rcpFromRef(x));
    auto out = VectorTraits<T>::get_data_nonconst(Teuchos::rcpFromRef(y));
    CHECK(in.size() == nc + 1);
    CHECK(out.size() == in.size());

    // Make views of data
    auto g = Teuchos::arrayView(in.get(), nc);
    CHECK(g.size() == nc);

    std::cout << "Anderson in vector: ";
    for( auto x : g )
        std::cout << x << " ";
    std::cout << std::endl;

    // This is an ArrayRCP
    auto gp = VectorTraits<T>::get_data_nonconst(d_gp);
    CHECK(gp.size() == nc);

    // get k from x
    double k = in[nc];
    if (d_node == 0)
    {
        std::cout << "Anderson k iterate: " << std::fixed << std::setw(8)
                  << k << std::endl;
    }
    CHECK(k > 0.0);

    // Apply prolongation P: g->f
    prolongate(g, *d_fission_sites);

    // Update the fission source for the next Monte Carlo transport
    update_source();

    // Reset tallier
    d_fisn_tally->reset();

    // Do a Monte Carlo iteration
    iterate(k);

    // Restrict to get new g = Rf
    restrict(*d_fission_sites, gp());

    std::cout << "New fission sites: ";
    for( auto x : gp() )
        std::cout << x << " ";
    std::cout << std::endl;

    std::cout << "Ran " << d_source->num_run() << " histories" << std::endl;
    d_fisn_tally->finalize(d_source->num_run());
    std::cout << "Fission tally result: ";
    for( auto x : d_fisn_tally->results() )
        std::cout << x.first << " ";
    std::cout << std::endl;

    // >>> Update F(g,k)

    double transport_eig = d_keff_tally->latest();
    REQUIRE( transport_eig > 0.0 );

    // F(g)
    for (int cell = 0; cell < nc; ++cell)
    {
        out[cell] = gp[cell] - g[cell];
    }

    // F(k)
    // Take norms of g and gp
    double norm_gp = d_blas.NRM2( nc, gp.getRawPtr(), 1 );
    double norm_g  = d_blas.NRM2( nc, g.getRawPtr(),  1 );

    CHECK(norm_g > 0.0);
    out[nc] = transport_eig - k;

    // Fix screwy formatting from NOX
    std::cout << std::scientific;
}

//---------------------------------------------------------------------------//
// PRIVATE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Prolongate, \f$\mathbf{P}: g\rightarrow f\f$.
 *
 * The prolongation operation is defined
 * \f[
   \mathbf{P}g = f' \frac{g}{g'}\:.
 * \f]
 */
template<class Geometry, class T>
void Anderson_Operator<Geometry,T>::prolongate(const_View              g,
                                               Fission_Site_Container &f) const
{
    REQUIRE(g.size() == d_mesh->num_cells());

    // Make views of g' (the old vector).
    auto gp = VectorTraits<T>::get_data(d_gp);
    CHECK(gp.size() == g.size());

    // make a new fission site container
    Fission_Site_Container nf;
    CHECK(nf.empty());

    // dimension vector for each cell
    Cartesian_Mesh::Dim_Vector ijk;

    // get the global RNG (by reference)
    auto rng = Global_RNG::d_rng;

    // Multiplicative correction
    double nu = 0.0;

    // correct the fission sites
    while (!f.empty())
    {
        // get the fission site of the back
        const auto &site = f.back();
        CHECK(d_mesh->find(site.r, ijk));

        // find the cell containing the site
        d_mesh->find(site.r, ijk);
        int cell = d_mesh->index(ijk[0], ijk[1], ijk[2]);
        CHECK(cell < gp.size());

        // calculate the multiplicative correction
        if (gp[cell] > 0.0)
            nu = g[cell] / gp[cell];

        // sample to determine the number of sites at this location
        int n = nu;
        CHECK(static_cast<double>(n) + 1.0 - nu >= 0.0 &&
              static_cast<double>(n) + 1.0 - nu <= 1.0);

        // with probability n+1-nu there will be n sites; with propability
        // nu-n there will be n+1 sites
        if (rng.ran() < nu - static_cast<double>(n))
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

    ENSURE(nf.empty());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Restrict, \f$\mathbf{R}f = g\f$.
 *
 * The restriction operation is applied by adding all of the fission sites in
 * each cell and dividing by volume to get a fission density.
 */
template<class Geometry, class T>
void Anderson_Operator<Geometry,T>::restrict(
        const Fission_Site_Container &f,
        View                          g) const
{
    REQUIRE(g.size() == d_mesh->num_cells());

    // Initialize g to zero
    std::fill(g.begin(), g.end(), 0.0);

    if( d_use_tally )
    {
        // Get results of fission tally
        auto fisn_result = d_fisn_tally->results();

        // Update g
        double g_sum = 0.0;
        for( int cell = 0; cell < d_mesh->num_cells(); ++cell )
        {
            g[cell] = fisn_result[cell].first;
            g_sum += g[cell];
        }

        // Normalize the distribution
        // This prevents wild swings in number of fission sites
        // in early iterations
        for( int cell = 0; cell < d_mesh->num_cells(); ++cell )
        {
            g[cell] *= static_cast<double>(d_mesh->num_cells()) / g_sum;
        }
    }
    else
    {
        // dimension vector for each cell
        Cartesian_Mesh::Dim_Vector ijk;

        // loop through the global fission sites and add them up in each global
        // mesh cell
        for (const auto &site : f)
        {
            CHECK(d_mesh->find(site.r, ijk));

            // find the cell containing the site
            d_mesh->find(site.r, ijk);
            CHECK(d_mesh->index(ijk[0], ijk[1], ijk[2]) < g.size());

            // add the site to the global field
            g[d_mesh->index(ijk[0], ijk[1], ijk[2])] += 1.0;

        }

        double num_cells = static_cast<double>(d_mesh->num_cells());
        double num_p     = static_cast<double>(d_source->Np());

        // Normalize
        for (int cell = 0; cell < d_mesh->num_cells(); ++cell)
        {
            g[cell] *= num_cells / num_p;
        }
    }
}

} // end namespace profugus

#endif // MC_mc_Anderson_Operator_t_hh

//---------------------------------------------------------------------------//
// end of MC/mc/Anderson_Operator.t.hh
//---------------------------------------------------------------------------//

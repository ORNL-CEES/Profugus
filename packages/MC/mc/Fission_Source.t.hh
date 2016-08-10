//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   MC/mc/Fission_Source.t.hh
 * \author Thomas M. Evans
 * \date   Mon May 05 14:22:46 2014
 * \brief  Fission_Source template member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_mc_Fission_Source_t_hh
#define MC_mc_Fission_Source_t_hh

#include <algorithm>
#include <numeric>

#include "Teuchos_Array.hpp"

#include "Fission_Source.hh"

#include "harness/DBC.hh"
#include "harness/Diagnostics.hh"
#include "comm/global.hh"
#include "comm/Timing.hh"
#include "utils/Constants.hh"
#include "mc/Global_RNG.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template <class Geometry>
Fission_Source<Geometry>::Fission_Source(RCP_Std_DB     db,
                                         SP_Geometry    geometry,
                                         SP_Physics     physics,
                                         SP_RNG_Control rng_control)
    : Base(geometry, physics, rng_control)
    , d_fission_rebalance(std::make_shared<Fission_Rebalance_t>())
    , d_np_requested(0)
    , d_np_total(0)
    , d_np_domain(0)
    , d_wt(0.0)
    , d_num_left(0)
    , d_num_run(0)
{
    using def::I; using def::J; using def::K;

    REQUIRE(b_geometry);
    REQUIRE(b_physics);
    REQUIRE(b_rng_control);

    // Boundaries in -X, +X, -Y, +Y, -Z, +Z
    Teuchos::Array<double> extents(6, 0.);

    // Assign large extents that will be trimmed to the geometry by default
    extents[0] = extents[2] = extents[4] = -std::numeric_limits<double>::max();
    extents[1] = extents[3] = extents[5] =  std::numeric_limits<double>::max();

    extents = db->get("init_fission_src", extents);

    VALIDATE(extents.size() == 6,
             "Fission source must have 6 entries, but it has "
             << extents.size());

    // get the low and upper bounds of the geometry
    auto box = geometry->get_extents();
    const Space_Vector &low_edge  = box.lower();
    const Space_Vector &high_edge = box.upper();

    for (int i = 0; i < 3; ++i)
    {
        double lower = extents[2 * i];
        double upper = extents[2 * i + 1];

        lower = std::max(lower, low_edge[i]);
        upper = std::min(upper, high_edge[i]);

        d_lower[i] = lower;
        d_width[i] = upper - lower;

        VALIDATE(d_width[i] > 0., "Fission source width for axis " << i
                 << " has non-positive width " << d_width[i]
                 << " (lower=" << lower << ", upper=" << upper << ")");
    }

    // store the total number of requested particles per cycle
    d_np_requested = static_cast<size_type>(db->get("Np", 1000));
    VALIDATE(d_np_requested > 0, "Number of source particles ("
            << d_np_requested << ") must be positive");

    // initialize the total for the first cycle
    d_np_total = d_np_requested;
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Build the initial fission source.
 */
template <class Geometry>
void Fission_Source<Geometry>::build_initial_source()
{
    // send an empty mesh and view
    SP_Cart_Mesh     mesh;
    Const_Array_View view;
    build_initial_source(mesh, view);

    ENSURE(d_wt >= 1.0);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build the initial source from a mesh distribution.
 */
template <class Geometry>
void Fission_Source<Geometry>::build_initial_source(SP_Cart_Mesh     mesh,
                                                    Const_Array_View fis_dens)
{
    REQUIRE(d_np_total > 0);

    SCOPED_TIMER("MC::Fission_Source.build_initial_source");

    // set the fission site container to an unassigned state
    d_fission_sites = SP_Fission_Sites();
    CHECK(!d_fission_sites);

    // make the RNG for this cycle
    make_RNG();

    // build the domain-replicated fission source
    build_DR(mesh, fis_dens);

    // set counters
    d_num_left = d_np_domain;
    d_num_run  = 0;

    // weight per particle
    d_wt = static_cast<double>(d_np_requested) /
           static_cast<double>(d_np_total);

    // initialize starting cell
    d_current_cell = 0;

    profugus::global_barrier();

    ENSURE(d_wt > 0.0);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build a source from a fission site container.
 *
 * \param fission_sites
 */
template <class Geometry>
void Fission_Source<Geometry>::build_source(SP_Fission_Sites &fission_sites)
{
    REQUIRE(fission_sites);

    // build an empty fission site container if one doesn't exist (meaning
    // that we haven't yet initialized from an existing fission source)
    if (!d_fission_sites)
    {
        d_fission_sites = std::make_shared<Fission_Site_Container>();
    }

    // the internal fission source should be empty
    REQUIRE(d_fission_sites->empty());

    SCOPED_TIMER("MC::Fission_Source.build_source");

    // swap the input fission sites with the internal storage fissino sites
    d_fission_sites.swap(fission_sites);

    // rebalance across sets (when number of blocks per set > 1; the
    // set-rebalance may try to do some load-balancing when it can, that is
    // why this call should comm after the gather; otherwise the
    // load-balancing could provide poor results)
    d_fission_rebalance->rebalance(*d_fission_sites);

    // get the number of fission sites on this domain, on this set, and
    // globally from the fission-rebalance
    d_np_domain = d_fission_rebalance->num_fissions();
    d_np_total  = d_fission_rebalance->num_global_fissions();
    CHECK(d_np_domain >= d_fission_sites->size()); // there could be multiple
                                                    // fissions at a single
                                                    // site

    // make the RNG for this cycle
    make_RNG();

    // set counters
    d_num_left = d_np_domain;
    d_num_run  = 0;

    // weight per particle
    d_wt = static_cast<double>(d_np_requested) /
           static_cast<double>(d_np_total);

    profugus::global_barrier();

    ENSURE(d_wt > 0.0);
    ENSURE(fission_sites);
    ENSURE(fission_sites->empty());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Create a fission site container.
 */
template <class Geometry>
auto Fission_Source<Geometry>::create_fission_site_container() const
    -> SP_Fission_Sites
{
    SP_Fission_Sites fs(std::make_shared<Fission_Site_Container>());
    ENSURE(fs);
    ENSURE(fs->empty());
    return fs;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get a particle from the source.
*/
template <class Geometry>
auto Fission_Source<Geometry>::get_particle() -> SP_Particle
{
    using def::I; using def::J; using def::K;

    REQUIRE(d_wt > 0.0);
    REQUIRE(profugus::Global_RNG::d_rng.assigned());

    // particle
    SP_Particle p;
    CHECK(!p);

    // return a null particle if no source
    if (!d_num_left)
    {
        ENSURE(d_fission_sites ? d_fission_sites->empty() : true);
        return p;
    }

    SCOPED_TIMER_2("MC::Fission_Source.get_particle");

    // make a particle
    p = std::make_shared<Particle_t>();

    // use the global rng on this domain for the random number generator
    p->set_rng(profugus::Global_RNG::d_rng);
    RNG rng = p->rng();

    // material id
    int matid = 0;

    // particle position and isotropic direction
    Space_Vector r, omega;

    // sample the angle isotropically
    Base::sample_angle(omega, rng);

    // sample flag
    bool sampled;

    // if there is a fission site container than get the particle from there;
    // otherwise assume this is an initial source
    if (!is_initial_source())
    {
        CHECK(!d_fission_sites->empty());

        // get the last element in the site container
        Fission_Site &fs = d_fission_sites->back();

        // get the location of the physics site
        r = b_physics->fission_site(fs);

        // intialize the geometry state
        b_geometry->initialize(r, omega, p->geo_state());

        // get the material id
        matid = b_geometry->matid(p->geo_state());

        // initialize the physics state at the fission site
        sampled = b_physics->initialize_fission(fs, *p);
        CHECK(sampled);

        // pop this fission site from the list
        d_fission_sites->pop_back();
    }
    else
    {
        matid = sample_geometry(r, omega, *p, rng);
    }

    // set the material id in the particle
    p->set_matid(matid);

    // set particle weight
    p->set_wt(d_wt);

    // make particle alive
    p->live();

    // update counters
    d_num_left--;
    d_num_run++;

    ENSURE(p->matid() == matid);
    return p;
}

//---------------------------------------------------------------------------//
// PRIVATE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Build the initial source for Domain Replicated decompositions.
 *
 * In DR decompositions the number of sets is equal to the number of domains
 * (1 block per set).  Thus, the number of particles per set is equal to the
 * number of particles per domain.
 */
template <class Geometry>
void Fission_Source<Geometry>::build_DR(SP_Cart_Mesh     mesh,
                                        Const_Array_View fis_dens)
{
    // calculate the number of particles per domain and set (equivalent)
    d_np_domain = d_np_total / b_nodes;

    // recalculate the total number of particles (we want the same number of
    // particles in each domain, so the total may change slightly from the
    // requested value)
    d_np_total = d_np_domain * b_nodes;

    // if there is a mesh then do stratified sampling to calculate the initial
    // fission distribution
    if (mesh)
    {
        // Check extents of Cartesian mesh vs. geometry
        auto geom_extents = b_geometry->get_extents();
        double tol = 1e-10;

        // X bounds
        double mesh_bnd = mesh->low_corner(def::I);
        double geom_bnd = geom_extents.lower()[def::I];
        VALIDATE(mesh_bnd > geom_bnd - tol,
                 "Lower x bound of mesh (" << mesh_bnd << ") is below MC "
                 "geometry lower x bound (" << geom_bnd << ")");
        mesh_bnd = mesh->high_corner(def::I);
        geom_bnd = geom_extents.upper()[def::I];
        VALIDATE(mesh_bnd < geom_bnd + tol,
                 "Upper x bound of mesh (" << mesh_bnd << ") is above MC "
                 "geometry upper x bound (" << geom_bnd << ")");

        // Y bounds
        mesh_bnd = mesh->low_corner(def::J);
        geom_bnd = geom_extents.lower()[def::J];
        VALIDATE(mesh_bnd > geom_bnd - tol,
                 "Lower y bound of mesh (" << mesh_bnd << ") is below MC "
                 "geometry lower y bound (" << geom_bnd << ")");
        mesh_bnd = mesh->high_corner(def::J);
        geom_bnd = geom_extents.upper()[def::J];
        VALIDATE(mesh_bnd < geom_bnd + tol,
                 "Upper y bound of mesh (" << mesh_bnd << ") is above MC "
                 "geometry upper y bound (" << geom_bnd << ")");

        // Z bounds
        mesh_bnd = mesh->low_corner(def::K);
        geom_bnd = geom_extents.lower()[def::K];
        VALIDATE(mesh_bnd > geom_bnd - tol,
                 "Lower z bound of mesh (" << mesh_bnd << ") is below MC "
                 "geometry lower z bound (" << geom_bnd << ")");
        mesh_bnd = mesh->high_corner(def::K);
        geom_bnd = geom_extents.upper()[def::K];
        VALIDATE(mesh_bnd < geom_bnd + tol,
                 "Upper z bound of mesh (" << mesh_bnd << ") is above MC "
                 "geometry upper z bound (" << geom_bnd << ")");

        REQUIRE(mesh->num_cells() == fis_dens.size());

        // number of cells in the mesh
        int num_cells = mesh->num_cells();

        // determine the total number of fissions
        double fissions = 0.0;
        for (int cell = 0; cell < num_cells; ++cell)
        {
            fissions += fis_dens[cell] * mesh->volume(cell);
        }
        CHECK(fissions > 0.0);

        // allocate fission distribution
        Vec_Int n(num_cells, 0);

        // pre-sample sites on this domain
        double nu            = 0.0;
        int    new_np_domain = 0;
        for (int cell = 0; cell < num_cells; ++cell)
        {
            // calculate the expected number of sites in this cell
            nu = fis_dens[cell] * mesh->volume(cell) / fissions * d_np_domain;

            // there can be n or n+1 sites; with probability n+1-nu there will
            // be n sites, with probability nu-n there will be n+1 sites
            n[cell] = nu;
            if (Global_RNG::d_rng.ran() < nu - static_cast<double>(n[cell]))
            {
                ++n[cell];
            }

            // add up the number of particles on this domain
            new_np_domain += n[cell];
        }

        // store the distributions persistently
        std::swap(n, d_fis_dist);
        CHECK(d_fis_dist.size() == num_cells);

        // update the number of particles globally and on the domain
        d_np_domain = new_np_domain;
        d_np_total  = d_np_domain;
        profugus::global_sum(&d_np_total, 1);

        // store the mesh for later use
        d_fis_mesh = mesh;
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Sample geometry to get a particle.
 */
template <class Geometry>
int Fission_Source<Geometry>::sample_geometry(Space_Vector       &r,
                                              const Space_Vector &omega,
                                              Particle_t         &p,
                                              RNG_t               rng)
{
    using def::I; using def::J; using def::K;

    // sampled complete flag
    bool sampled = false;

    // material id
    int matid = 0;

    // >>> Sample the full geometry
    if (!d_fis_mesh)
    {
        CHECK(d_fis_dist.empty());

        // Sanity check that some fissionable material exists
        std::vector<int> matids;
        b_physics->xs().get_matids(matids);
        CHECK( matids.size() > 0 );
        bool fission_found = false;
        for( auto matid : matids )
        {
            fission_found = fission_found | b_physics->is_fissionable(matid);
        }
        INSIST(fission_found,"No fissionable material found in problem.");

        // sample the geometry until a fission site is found (if there is no
        // fission in a given domain the number of particles on that domain is
        // zero, and we never get here) --> so, fission sampling should always
        // be successful
        int attempts = 0;
        while (!sampled)
        {
            attempts++;
            INSIST(attempts < 1000,
                   "Failed to locate viable fission site after 1000 samples.");

            // sample a point in the geometry
            r[I] = d_width[I] * rng.ran() + d_lower[I];
            r[J] = d_width[J] * rng.ran() + d_lower[J];
            r[K] = d_width[K] * rng.ran() + d_lower[K];

            // intialize the geometry state
            b_geometry->initialize(r, omega, p.geo_state());

            // get the material id
            matid = b_geometry->matid(p.geo_state());

            // try initializing fission here, if it is successful we are
            // finished
            if (b_physics->initialize_fission(matid, p))
            {
                sampled = true;
            }

            DIAGNOSTICS_TWO(integers["fission_src_samples"]++);
        }
    }
    // >>> Sample the mesh source
    else
    {
        CHECK(d_fis_dist.size() == d_fis_mesh->num_cells());
        CHECK(d_current_cell < d_fis_dist.size());
        CHECK(d_fis_dist[d_current_cell] >= 0);

        // determine the particle birth cell
        while (d_fis_dist[d_current_cell] == 0)
        {
            ++d_current_cell;
            CHECK(d_current_cell < d_fis_dist.size());
        }

        // get the logical cell indices in the mesh
        auto ijk = d_fis_mesh->cardinal(d_current_cell);

        // get the low/high bounds for the cell
        double xdims[] = {d_fis_mesh->edges(I)[ijk[I]],
                          d_fis_mesh->edges(I)[ijk[I] + 1]};
        double ydims[] = {d_fis_mesh->edges(J)[ijk[J]],
                          d_fis_mesh->edges(J)[ijk[J] + 1]};
        double zdims[] = {d_fis_mesh->edges(K)[ijk[K]],
                          d_fis_mesh->edges(K)[ijk[K] + 1]};

        // sample the appropriate cell until a fissionable region is found
        while (!sampled)
        {
            // sample a point in the geometry
            r[I] = (xdims[1] - xdims[0]) * rng.ran() + xdims[0];
            r[J] = (ydims[1] - ydims[0]) * rng.ran() + ydims[0];
            r[K] = (zdims[1] - zdims[0]) * rng.ran() + zdims[0];

            // intialize the geometry state
            b_geometry->initialize(r, omega, p.geo_state());

            // get the material id
            matid = b_geometry->matid(p.geo_state());

            // try initializing fission here, if it is successful we are
            // finished
            if (b_physics->initialize_fission(matid, p))
            {
                sampled = true;
            }

            DIAGNOSTICS_TWO(integers["fission_src_samples"]++);
        }

        // subtract a particle from the current cell
        --d_fis_dist[d_current_cell];
        CHECK(d_fis_dist[d_current_cell] >= 0);
    }

    return matid;
}

} // end namespace profugus

#endif // MC_mc_Fission_Source_t_hh

//---------------------------------------------------------------------------//
//                 end of Fission_Source.t.hh
//---------------------------------------------------------------------------//

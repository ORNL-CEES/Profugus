//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Fission_Source.cc
 * \author Thomas M. Evans
 * \date   Mon May 05 14:22:46 2014
 * \brief  Fission_Source member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

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
Fission_Source::Fission_Source(RCP_Std_DB     db,
                               SP_Geometry    geometry,
                               SP_Physics     physics,
                               SP_RNG_Control rng_control)
    : Base(geometry, physics, rng_control)
    , d_fission_rebalance(std::make_shared<Fission_Rebalance>())
    , d_np_requested(0)
    , d_np_total(0)
    , d_np_domain(0)
    , d_wt(0.0)
    , d_num_left(0)
    , d_num_run(0)
{
    using def::I; using def::J; using def::K;

    Require (b_geometry);
    Require (b_physics);
    Require (b_rng_control);

    // Boundaries in -X, +X, -Y, +Y, -Z, +Z
    Teuchos::Array<double> extents(6, 0.);

    // Assign large extents that will be trimmed to the geometry by default
    extents[0] = extents[2] = extents[4] = -std::numeric_limits<double>::max();
    extents[1] = extents[3] = extents[5] =  std::numeric_limits<double>::max();

    extents = db->get("init_fission_src", extents);

    Validate(extents.size() == 6,
             "Fission source must have 6 entries, but it has "
             << extents.size());

    // get the underlying geometry array
    const auto &array = b_geometry->array();

    // get the low and upper bounds of the geometry
    const Space_Vector &low_edge = array.corner();
    Space_Vector high_edge       = array.corner();
    high_edge[I] += array.pitch(I);
    high_edge[J] += array.pitch(J);
    high_edge[K] += array.height();

    for (int i = 0; i < 3; ++i)
    {
        double lower = extents[2 * i];
        double upper = extents[2 * i + 1];

        lower = std::max(lower, low_edge[i]);
        upper = std::min(upper, high_edge[i]);

        d_lower[i] = lower;
        d_width[i] = upper - lower;

        Validate(d_width[i] > 0., "Fission source width for axis " << i
                 << " has non-positive width " << d_width[i]
                 << " (lower=" << lower << ", upper=" << upper << ")");

    }

    // store the total number of requested particles per cycle
    d_np_requested = static_cast<size_type>(db->get("Np", 1000));
    Validate(d_np_requested > 0, "Number of source particles ("
            << d_np_requested << ") must be positive");

    // initialize the total for the first cycle
    d_np_total = d_np_requested;
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Build the initial source.
 */
void Fission_Source::build_initial_source()
{
    Require (d_np_total > 0);

    SCOPED_TIMER("MC::Fission_Source.build_initial_source");

    // set the fission site container to an unassigned state
    d_fission_sites = SP_Fission_Sites();
    Check (!d_fission_sites);

    // make the RNG for this cycle
    make_RNG();

    // build the domain-replicated fission source
    build_DR();

    // set counters
    d_num_left = d_np_domain;
    d_num_run  = 0;

    // weight per particle
    d_wt = static_cast<double>(d_np_requested) /
           static_cast<double>(d_np_total);

    profugus::global_barrier();

    Ensure (d_wt >= 1.0);
    Ensure (d_wt > 0.0);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build a source from a fission site container.
 *
 * \param fission_sites
 */
void Fission_Source::build_source(SP_Fission_Sites &fission_sites)
{
    Require (fission_sites);

    // build an empty fission site container if one doesn't exist (meaning
    // that we haven't yet initialized from an existing fission source)
    if (!d_fission_sites)
    {
        d_fission_sites = std::make_shared<Fission_Site_Container>();
    }

    // the internal fission source should be empty
    Require (d_fission_sites->empty());

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
    Check (d_np_domain >= d_fission_sites->size()); // there could be multiple
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

    Ensure (d_wt > 0.0);
    Ensure (fission_sites);
    Ensure (fission_sites->empty());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Create a fission site container.
 */
Fission_Source::SP_Fission_Sites
Fission_Source::create_fission_site_container() const
{
    SP_Fission_Sites fs(std::make_shared<Fission_Site_Container>());
    Ensure (fs);
    Ensure (fs->empty());
    return fs;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get a particle from the source.
*/
Fission_Source::SP_Particle Fission_Source::get_particle()
{
    using def::I; using def::J; using def::K;

    Require (d_wt > 0.0);
    Require (profugus::Global_RNG::d_rng.assigned());

    // particle
    SP_Particle p;
    Check (!p);

    // return a null particle if no source
    if (!d_num_left)
    {
        Ensure (d_fission_sites ? d_fission_sites->empty() : true);
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
        Check (!d_fission_sites->empty());

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
        Check (sampled);

        // pop this fission site from the list
        d_fission_sites->pop_back();
    }
    else
    {
        // sample the geometry until a fission site is found (if there is no
        // fission in a given domain the number of particles on that domain is
        // zero, and we never get here) --> so, fission sampling should always
        // be successful
        sampled = false;
        while (!sampled)
        {
            // sample a point in the geometry
            sample_r(r, rng);

            // intialize the geometry state
            b_geometry->initialize(r, omega, p->geo_state());

            // get the material id
            matid = b_geometry->matid(p->geo_state());

            // try initializing fission here, if it is successful we are
            // finished
            if (b_physics->initialize_fission(matid, *p))
            {
                sampled = true;
            }

            DIAGNOSTICS_TWO(integers["fission_src_samples"]++);
        }
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

    Ensure (p->matid() == matid);
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
void Fission_Source::build_DR()
{
    // calculate the number of particles per domain and set (equivalent)
    d_np_domain = d_np_total / b_nodes;

    // recalculate the total number of particles (we want the same number of
    // particles in each domain, so the total may change slightly from the
    // requested value)
    d_np_total = d_np_domain * b_nodes;
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Fission_Source.cc
//---------------------------------------------------------------------------//
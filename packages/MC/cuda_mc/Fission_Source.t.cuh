//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Fission_Source.t.cuh
 * \author Steven Hamilton
 * \date   Mon May 05 14:22:46 2014
 * \brief  Fission_Source template member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Fission_Source_t_cuh
#define cuda_mc_Fission_Source_t_cuh

#include <algorithm>
#include <numeric>

#include "Teuchos_Array.hpp"

#include "Fission_Source.cuh"
#include "Sampler.cuh"

#include "harness/DBC.hh"
#include "harness/Diagnostics.hh"
#include "comm/global.hh"
#include "comm/Timing.hh"
#include "utils/Constants.hh"

namespace cuda_mc
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template <class Geometry>
Fission_Source<Geometry>::Fission_Source(RCP_Std_DB     db,
                                         SDP_Geometry   geometry,
                                         SDP_Physics    physics)
    : Base(geometry)
    , d_fission_rebalance(std::make_shared<Fission_Rebalance>())
    , d_physics_host(physics)
    , d_nodes(profugus::nodes())
    , d_have_sites(false)
    , d_wt(0.0)
{
    using def::I; using def::J; using def::K;

    REQUIRE(b_geometry);

    REQUIRE(d_physics_host.get_host_ptr());
    REQUIRE(d_physics_host.get_device_ptr());
    d_physics = physics.get_device_ptr();

    // Boundaries in -X, +X, -Y, +Y, -Z, +Z
    Teuchos::Array<double> extents(6, 0.);

    // Assign large extents that will be trimmed to the geometry by default
    extents[0] = extents[2] = extents[4] = -std::numeric_limits<double>::max();
    extents[1] = extents[3] = extents[5] =  std::numeric_limits<double>::max();

    extents = db->get("init_fission_src", extents);

    INSIST(extents.size() == 6, "Fission source must have 6 entries");

    // get the low and upper bounds of the geometry
    auto low_edge  = geometry.get_host_ptr()->lower();
    auto high_edge = geometry.get_host_ptr()->upper();

    // X Bounds
    double lower = extents[0];
    double upper = extents[1];

    lower = std::max(lower, low_edge[I]);
    upper = std::min(upper, high_edge[I]);

    d_lower[I] = lower;
    d_width[I] = upper - lower;

    // Y Bounds
    lower = extents[2];
    upper = extents[3];

    lower = std::max(lower, low_edge[J]);
    upper = std::min(upper, high_edge[J]);

    d_lower[J] = lower;
    d_width[J] = upper - lower;

    // Z Bounds
    lower = extents[4];
    upper = extents[5];

    lower = std::max(lower, low_edge[K]);
    upper = std::min(upper, high_edge[K]);

    d_lower[K] = lower;
    d_width[K] = upper - lower;

    INSIST(d_width[I] > 0., "Fission source x width is non-positive");
    INSIST(d_width[J] > 0., "Fission source y width is non-positive");
    INSIST(d_width[K] > 0., "Fission source z width is non-positive");

    // store the total number of requested particles per cycle
    d_np_requested = 1000;
    if (db->isType<int>("Np"))
        d_np_requested = db->get<int>("Np");
    else if (db->isType<size_type>("Np"))
        d_np_requested = db->get<size_type>("Np");
    else if (db->isParameter("Np"))
        VALIDATE(false,"Unrecognized type for parameter Np.");
    INSIST(d_np_requested > 0., "Number of source particles must be positive");

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
    REQUIRE(d_np_total > 0);

    // build the domain-replicated fission source
    build_DR();

    d_np_left = d_np_domain;

    // weight per particle
    d_wt = static_cast<double>(d_np_requested) /
           static_cast<double>(d_np_total);

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
void Fission_Source<Geometry>::build_source(SP_Fission_Site_Vec &fission_sites)
{
    REQUIRE(fission_sites);

    // Make fission site vector if it hasn't been built yet
    if( !d_fission_site_vec )
        d_fission_site_vec = std::make_shared<Fission_Site_Vector>();

    SCOPED_TIMER("MC::Fission_Source.build_source");

    // swap the input fission sites with the internal storage fission sites
    d_fission_site_vec.swap(fission_sites);

    // rebalance across sets (when number of blocks per set > 1; the
    // set-rebalance may try to do some load-balancing when it can, that is
    // why this call should comm after the gather; otherwise the
    // load-balancing could provide poor results)
    std::vector<Fission_Site> host_sites(d_fission_site_vec->size());
    thrust::copy(d_fission_site_vec->begin(),
                 d_fission_site_vec->end(),
                 host_sites.begin());
    d_fission_rebalance->rebalance(host_sites);
    *d_fission_site_vec = host_sites;

    // get the number of fission sites on this domain, on this set, and
    // globally from the fission-rebalance
    d_np_domain = d_fission_rebalance->num_fissions();
    d_np_total  = d_fission_rebalance->num_global_fissions();
    CHECK(d_np_domain >= d_fission_site_vec->size()); // there could be multiple
                                                      // fissions at a single
                                                      // site

    d_np_left = d_np_domain;

    // weight per particle
    d_wt = static_cast<double>(d_np_requested) /
           static_cast<double>(d_np_total);

    profugus::global_barrier();

    // Assign raw pointer for on-device access
    d_fission_sites = d_fission_site_vec->data().get();
    d_have_sites = true;

    ENSURE(d_wt > 0.0);
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
void Fission_Source<Geometry>::build_DR()
{
    // calculate the number of particles per domain and set (equivalent)
    d_np_domain = d_np_total / d_nodes;

    // recalculate the total number of particles (we want the same number of
    // particles in each domain, so the total may change slightly from the
    // requested value)
    d_np_total = d_np_domain * d_nodes;
}

} // end namespace cuda_mc

#endif // cuda_mc_Fission_Source_t_cuh

//---------------------------------------------------------------------------//
//                 end of Fission_Source.t.cuh
//---------------------------------------------------------------------------//

//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Uniform_Source.cc
 * \author Thomas M. Evans
 * \date   Tue May 06 16:43:26 2014
 * \brief  Uniform_Source member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <numeric>

#include "Teuchos_Array.hpp"

#include "harness/Soft_Equivalence.hh"
#include "harness/DBC.hh"
#include "comm/Timing.hh"
#include "comm/global.hh"
#include "Sampler.hh"
#include "Uniform_Source.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 *
 * \param db
 * \param geometry
 * \param physics
 * \param rng_control
 */
Uniform_Source::Uniform_Source(RCP_Std_DB     db,
                               SP_Geometry    geometry,
                               SP_Physics     physics,
                               SP_RNG_Control rng_control)
    : Base(geometry, physics, rng_control)
    , d_erg_cdf(b_physics->num_groups(), 0.0)
    , d_np_requested(0)
    , d_np_total(0)
    , d_np_domain(0)
    , d_wt(1.0)
    , d_np_run(0)
{
    REQUIRE(!db.is_null());

    // store the total number of requested particles
    d_np_requested = static_cast<size_type>(db->get("Np", 1000));
    VALIDATE(d_np_requested > 0., "Number of source particles ("
            << d_np_requested << ") must be positive");

    // initialize the total
    d_np_total = d_np_requested;

    // get the spectral shape
    const auto &shape = db->get(
        "spectral_shape", Teuchos::Array<double>(b_physics->num_groups(), 1.0));
    CHECK(shape.size() == d_erg_cdf.size());

    // calculate the normalization
    double norm = std::accumulate(shape.begin(), shape.end(), 0.0);
    CHECK(norm > 0.0);

    // assign to the shape cdf
    norm  = 1.0 / norm;
    d_erg_cdf[0] = shape[0] * norm;
    for ( int n = 1; n < d_erg_cdf.size(); ++n )
    {
        d_erg_cdf[n] = d_erg_cdf[n-1] + shape[n] * norm;
    }
    CHECK(profugus::soft_equiv(1.0, d_erg_cdf.back(), 1.0e-6));

    // initialize timers in this class, which may be necessary because domains
    // with no source will not make this timer otherwise
#if UTILS_TIMING > 0
    profugus::Timing_Diagnostics::update_timer(
        "profugus::Uniform_Source.get_particle", 0.0);
#endif
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Build the initial source.
 * \param geometric_shape
 */
void Uniform_Source::build_source(SP_Shape geometric_shape)
{
    REQUIRE(geometric_shape);

    SCOPED_TIMER("profugus::Uniform_Source.build_source");

    // store the spatial shape
    d_geo_shape = geometric_shape;

    // build the source based on domain replication
    build_DR();

    // set counters
    d_np_run.store(0);

    profugus::global_barrier();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get a particle from the source.
 */
Uniform_Source::SP_Particle Uniform_Source::get_particle()
{
    using def::I; using def::J; using def::K;

    REQUIRE(d_wt > 0.0);
    REQUIRE(d_geo_shape);

    // Get how many particles we have run and add one for this one.
    size_type np_run = d_np_run.fetch_add( 1 );

    // unassigned particle
    SP_Particle p;

    // return a null particle if no source. this may be more than one if
    // several threads have recently modified the run count atomic
    if (np_run >= d_np_domain)
    {
        return p;
    }

    // make a particle
    p = std::make_shared<Particle_t>();

    // create a unique rng for the particle
    int stream_id = np_run*b_nodes + b_node;
    p->set_rng( this->b_rng_control->rng(stream_id) );
    auto rng = p->rng();

    // material id
    int matid = 0;

    // particle position and direction
    Space_Vector r, omega;

    // sample the angle isotropically
    Base::sample_angle(omega, rng);

    // sample the geometry shape-->we should not get here if there are no
    // particles on this domain
    r = d_geo_shape->sample(rng);

    // intialize the geometry state
    b_geometry->initialize(r, omega, p->geo_state());

    // get the material id
    matid = b_geometry->matid(p->geo_state());

    // initialize the physics state by manually sampling the group
    int group = sampler::sample_discrete_CDF(
        d_erg_cdf.size(), &d_erg_cdf[0], rng.ran());
    CHECK(group < b_physics->num_groups());
    p->set_group(group);

    // set the material id in the particle
    p->set_matid(matid);

    // set particle weight
    p->set_wt(d_wt);

    // make particle alive
    p->live();

    ENSURE(p->matid() == matid);
    return p;
}

//---------------------------------------------------------------------------//
// PRIVATE IMPLEMENTATION
//---------------------------------------------------------------------------//
/*!
 * \brief Build domain replicated source.
 */
void Uniform_Source::build_DR()
{
    // calculate the number of particles per domain and set (equivalent)
    d_np_domain = d_np_total / b_nodes;

    // recalculate the total number of particles (we want the same number of
    // particles in each domain, so the total may change slightly from the
    // requested value)
    d_np_total  = d_np_domain * b_nodes;
}

} // end namespace profugus

//---------------------------------------------------------------------------//
//                 end of Uniform_Source.cc
//---------------------------------------------------------------------------//
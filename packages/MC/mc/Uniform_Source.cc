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
#include "Global_RNG.hh"
#include "Sampler.hh"
#include "Uniform_Source.hh"

namespace profugus
{

//---------------------------------------------------------------------------//
// STATIC VARIABLES
//---------------------------------------------------------------------------//

Source::size_type Uniform_Source::d_np_left = 0;
Source::size_type Uniform_Source::d_np_run  = 0;

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
    REMEMBER(double sum = 0.0);
    norm  = 1.0 / norm;
    int n = 0;
    for (double &c : d_erg_cdf)
    {
        c = shape[n] * norm;
        REMEMBER(sum += c);
        ++n;
    }
    ENSURE(soft_equiv(sum, 1.0));

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
    REQUIRE(!in_thread_parallel_region());

    SCOPED_TIMER("profugus::Uniform_Source.build_source");

    // store the spatial shape
    d_geo_shape = geometric_shape;

    // make the RNG for this cycle
    Base::make_RNG();

    // build the source based on domain replication
    build_DR();

    // set counters
    d_np_left = d_np_domain;
    d_np_run  = 0;

#pragma omp parallel copyin(d_np_left, d_np_run)
    {
        // store number of threads in the team and thread id
        int nt = num_current_threads();
        int id = thread_id();

        // extra particles
        int extra = d_np_left % nt;

        // base number of particles per thread
        d_np_left = d_np_left / nt;

        // add extra to threads
        if (id < extra)
            ++d_np_left;
    }

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
    REQUIRE(profugus::Global_RNG::d_rng.assigned());
    REQUIRE(d_geo_shape);

    // unassigned particle
    SP_Particle p;

    // return a null particle if no source
    if (!d_np_left)
    {
        return p;
    }

    SCOPED_TIMER_2("MC::Uniform_Source.get_particle");

    // make a particle
    p = std::make_shared<Particle_t>();

    // use the global rng on this domain for the random number generator
    p->set_rng(profugus::Global_RNG::d_rng);
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

    // update counters
    --d_np_left;
    ++d_np_run;

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

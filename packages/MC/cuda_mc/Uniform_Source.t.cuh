//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Uniform_Source.t.cuh
 * \author Steven Hamilton
 * \date   Tue May 06 16:43:26 2014
 * \brief  Uniform_Source member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Uniform_Source_t_cuh
#define cuda_mc_Uniform_Source_t_cuh

#include <numeric>

#include "Teuchos_Array.hpp"

#include "utils/View_Field.hh"
#include "cuda_utils/CudaDBC.hh"
#include "cuda_utils/Utility_Functions.hh"
#include "comm/Timing.hh"
#include "comm/global.hh"
#include "Sampler.cuh"
#include "Uniform_Source.cuh"

namespace cuda_mc
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
 */
template <class Geometry>
Uniform_Source<Geometry>::Uniform_Source(RCP_Std_DB     db,
                                         SDP_Geometry   geometry)
    : Base(geometry)
    , d_erg_cdf(db->get<int>("num_groups"))
    , d_np_requested(0)
    , d_np_total(0)
    , d_np_domain(0)
    , d_wt(1.0)
    , d_np_left(0)
    , d_np_run(0)
    , d_nodes(profugus::nodes())
{
    REQUIRE(!db.is_null());

    // store the total number of requested particles
    d_np_requested = static_cast<size_type>(db->get("Np", 1000));
    INSIST(d_np_requested > 0., "Number of source particles must be positive");

    // get the spectral shape
    d_num_groups = db->get<int>("num_groups");
    const auto &shape = db->get(
        "spectral_shape", Teuchos::Array<double>(d_num_groups, 1.0));
    CHECK(shape.size() == d_num_groups);

    // calculate the normalization
    double norm = std::accumulate(shape.begin(), shape.end(), 0.0);
    CHECK(norm > 0.0);

    // assign to the shape cdf
    double sum = 0.0;
    norm  = 1.0 / norm;
    int n = 0;
    Vec_Dbl erg_cdf(d_num_groups,0.0);
    for (double &c : erg_cdf)
    {
        double val = shape[n] * norm;
        sum += val;
        c = sum;
        ++n;
    }
    ENSURE(cuda::utility::soft_equiv(sum, 1.0));

    // Copy to device
    d_erg_cdf.assign(profugus::make_view(erg_cdf));

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
template <class Geometry>
void Uniform_Source<Geometry>::build_source(SP_Shape geometric_shape)
{
    REQUIRE(geometric_shape);

    SCOPED_TIMER("profugus::Uniform_Source.build_source");

    // store the spatial shape
    d_geo_shape_host = SDP_Shape(geometric_shape);
    d_geo_shape = d_geo_shape_host.get_device_ptr();

    // make the RNG for this cycle
    Base::make_RNG();

    // build the source based on domain replication
    build_DR();

    // set counters
    d_np_left = d_np_domain;
    d_np_run  = 0;

    profugus::global_barrier();
}

//---------------------------------------------------------------------------//
// PRIVATE IMPLEMENTATION
//---------------------------------------------------------------------------//
/*!
 * \brief Build domain replicated source.
 */
template <class Geometry>
void Uniform_Source<Geometry>::build_DR()
{
    // calculate the number of particles per domain and set (equivalent)
    d_np_domain = d_np_total / d_nodes;

    // recalculate the total number of particles (we want the same number of
    // particles in each domain, so the total may change slightly from the
    // requested value)
    d_np_total  = d_np_domain * d_nodes;
}

} // end namespace cuda_mc

#endif // cuda_mc_Uniform_Source_t_cuh

//---------------------------------------------------------------------------//
//                 end of Uniform_Source.t.cuh
//---------------------------------------------------------------------------//

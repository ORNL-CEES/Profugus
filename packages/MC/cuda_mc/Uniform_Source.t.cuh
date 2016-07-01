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
#include "cuda_utils/Launch_Args.t.cuh"
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
{
    REQUIRE(!db.is_null());

    // store the total number of requested particles
    d_np_total = 1000;
    if (db->isType<int>("Np"))
        d_np_total = db->get<int>("Np");
    else if (db->isType<size_type>("Np"))
        d_np_total = db->get<size_type>("Np");
    else if (db->isParameter("Np"))
        VALIDATE(false,"Unrecognized type for parameter Np.");
    INSIST(d_np_total > 0., "Number of source particles must be positive");
    d_np_domain = d_np_total / profugus::nodes();

    // get the spectral shape
    d_num_groups = db->get<int>("num_groups");
    const auto &shape = db->get(
        "spectral_shape", Teuchos::Array<double>(d_num_groups, 1.0));
    CHECK(shape.size() == d_num_groups);

    d_erg_cdf_vec.resize(d_num_groups);

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
    d_erg_cdf_vec = erg_cdf;
    d_erg_cdf = d_erg_cdf_vec.data().get();

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
void Uniform_Source<Geometry>::build_source(SDP_Shape geometric_shape)
{
    REQUIRE(geometric_shape.get_host_ptr());
    REQUIRE(geometric_shape.get_device_ptr());

    SCOPED_TIMER("profugus::Uniform_Source.build_source");

    // store the spatial shape
    d_geo_shape = geometric_shape.get_device_ptr();

    d_np_left = d_np_domain;
    profugus::global_barrier();
}

} // end namespace cuda_mc

#endif // cuda_mc_Uniform_Source_t_cuh

//---------------------------------------------------------------------------//
//                 end of Uniform_Source.t.cuh
//---------------------------------------------------------------------------//

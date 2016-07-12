//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Keff_Tally.t.cuh
 * \author Steven Hamilton
 * \date   Wed May 14 13:29:40 2014
 * \brief  Keff_Tally member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Keff_Tally_t_cuh
#define cuda_mc_Keff_Tally_t_cuh

#include "Keff_Tally.cuh"

#include "cuda_utils/CudaDBC.hh"
#include "comm/global.hh"

namespace cuda_mc
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Kcode solver should construct this with initial keff estimate.
 */
template <class Geometry>
Keff_Tally<Geometry>::Keff_Tally(double         keff_init,
                                 SDP_Physics    physics)
    : d_keff_cycle(keff_init)
{
    d_physics = physics.get_device_ptr();
    REQUIRE( d_physics );

    // reset tally
    reset();
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Calculate average keff over active cycles.
 *
 * This should really only be called after begin_active_cycles; but either way
 * it will return the average keff since either initialization of this tally or
 * the last call to begin_active_cycles.
 *
 * This is only meaningful after one complete cycle. If called before then,
 * we'll return an arbitrary value.
 */
template <class Geometry>
double Keff_Tally<Geometry>::mean() const
{
    if (d_cycle < 1)
        return -1.;

    double keff = d_keff_sum / d_cycle;

    ENSURE(keff >= 0.);
    return keff;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Calculate variance of average keff estimate over active cycles.
 *
 * This is only meaningful after two complete cycles. If called before then,
 * we'll return an arbitrary value.
 */
template <class Geometry>
double Keff_Tally<Geometry>::variance() const
{
    if (d_cycle < 2)
        return d_keff_sum * d_keff_sum;

    double var = (d_keff_sum_sq - d_keff_sum * d_keff_sum / d_cycle)
                    / (d_cycle * (d_cycle - 1));

    ENSURE(var >= 0.);
    return var;
}

//---------------------------------------------------------------------------//
// PUBLIC INTERFACE FUNCTIONS
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/*!
 * \brief Begin active cycles in a kcode calculation.
 *
 * This resets the accumulated keff statistics.
 */
template <class Geometry>
void Keff_Tally<Geometry>::begin_active_cycles()
{
    d_cycle       = 0;
    d_keff_sum    = 0.;
    d_keff_sum_sq = 0.;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Begin a new keff cycle
 *
 * This clears the current accumulated path lengths.
 */
template <class Geometry>
void Keff_Tally<Geometry>::begin_cycle(size_type num_particles)
{
    REQUIRE(num_particles>0);
    d_keff_cycle = 0.;
    d_keff_data.resize(num_particles);
    thrust::fill(d_keff_data.begin(),d_keff_data.end(),0.0);
    d_thread_keff = d_keff_data.data().get();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Mark the end of a keff cycle
 *
 * This performs a global sum across processors, because the provided
 * num_particles is the total number over all domains (blocks plus sets).
 *
 * We accumulate the sum and sum-of-squares of the keff so that we can calculate
 * averages and variances.
 */
template <class Geometry>
void Keff_Tally<Geometry>::end_cycle(double num_particles)
{
    REQUIRE(num_particles > 0.);

    // Keff estimate is total nu-sigma-f reaction rate / num particles
    d_keff_cycle = thrust::reduce( d_keff_data.begin(), d_keff_data.end(),
                                   0.0, thrust::plus<double>() );
    d_keff_cycle /= num_particles;

    // Do a global sum (since num_particles is global)
    profugus::global_sum(d_keff_cycle);

    // Accumulate first and second moments of cycles, as well as counter
    ++d_cycle;
    d_keff_sum    += d_keff_cycle;
    d_keff_sum_sq += d_keff_cycle * d_keff_cycle;

    // Store keff estimate
    d_all_keff.push_back(d_keff_cycle);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Clear accumulated keff values
 *
 * This will still preserve the latest keff value so that existing fission sites
 * will be sampled correctly.
 */
template <class Geometry>
void Keff_Tally<Geometry>::reset()
{
    d_cycle       = 0;
    d_keff_sum    = 0.;
    d_keff_sum_sq = 0.;

    d_all_keff.clear();
}

} // end namespace cuda_mc

#endif // cuda_mc_Keff_Tally_t_cuh

//---------------------------------------------------------------------------//
//                 end of Keff_Tally.t.cuh
//---------------------------------------------------------------------------//

//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Tallier.i.cuh
 * \author Thomas M. Evans
 * \date   Mon May 12 12:15:30 2014
 * \brief  Tallier member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Tallier_i_cuh
#define cuda_mc_Tallier_i_cuh

#include "cuda_utils/CudaDBC.hh"
#include "Tallier.cuh"

namespace cuda_mc
{

//---------------------------------------------------------------------------//
/*!
 * \brief Process path-length tally events.
 *
 * \param step step-length
 * \param p particle
 */
template <class Geometry>
__device__ void Tallier<Geometry>::path_length(double            step,
                                               const Particle_t &p)
{
    DEVICE_REQUIRE(step >= 0.0);

    // accumulate results for all pathlength tallies
    if (d_cell_tally)
        d_cell_tally->accumulate(step, p);
    if (d_keff_tally)
        d_keff_tally->accumulate(step, p);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Tally any source events.
 *
 * \param p particle
 */
template <class Geometry>
__device__ void Tallier<Geometry>::source(const Particle_t &p)
{
}

//---------------------------------------------------------------------------//
/*!
 * \brief Perform all end-history tally tasks.
 */
template <class Geometry>
__device__ void Tallier<Geometry>::end_history()
{
}

} // end namespace cuda_mc

#endif // cuda_mc_Tallier_t_cuh

//---------------------------------------------------------------------------//
//                 end of Tallier.t.cuh
//---------------------------------------------------------------------------//

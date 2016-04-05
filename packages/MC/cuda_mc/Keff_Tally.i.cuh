//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Keff_Tally.i.cuh
 * \author Steven Hamilton
 * \date   Wed May 14 13:29:40 2014
 * \brief  Keff_Tally inline member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Keff_Tally_i_cuh
#define cuda_mc_Keff_Tally_i_cuh

#include "Keff_Tally.cuh"

#include "cuda_utils/CudaDBC.hh"
#include "cuda_utils/Utility_Functions.hh"

namespace cuda_mc
{

//---------------------------------------------------------------------------//
// PUBLIC INTERFACE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Track particle and do tallying.
 *
 * This only uses the particle's weight. The accumulated tally is
 * \f[
   k_{\mbox{eff}} = wl\nu\sigma_{\mbox{f}}
 * \f]
 * where \f$l\f$ is the step-length.
 */
template <class Geometry>
__device__ void Keff_Tally<Geometry>::accumulate(double            step,
                                                 const Particle_t &p)
{
    REQUIRE(d_physics);

    cuda::utility::atomic_add_double( &d_keff_cycle,
        p.wt() * step * d_physics->total(profugus::physics::NU_FISSION, p) );
}

} // end namespace cuda_mc

#endif // cuda_mc_Keff_Tally_i_cuh

//---------------------------------------------------------------------------//
//                 end of Keff_Tally.i.cuh
//---------------------------------------------------------------------------//

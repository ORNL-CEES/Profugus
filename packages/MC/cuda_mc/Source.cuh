//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Source.cuh
 * \author Steven Hamilton
 * \date   Mon May 05 14:22:41 2014
 * \brief  Source class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Source_cuh
#define cuda_mc_Source_cuh

#include <memory>
#include <curand_kernel.h>

#include "utils/Definitions.hh"
#include "cuda_utils/Definitions.hh"
#include "cuda_utils/Constants.hh"
#include "cuda_utils/Shared_Device_Ptr.hh"

namespace cuda_mc
{

//===========================================================================//
/*!
 * \class Source
 * \brief Base class definition for Monte Carlo sources.
 */
//===========================================================================//

template <class Geometry>
class Source
{
  public:
    //@{
    //! Useful typedefs.
    typedef Geometry                            Geometry_t;
    typedef typename Geometry_t::Geo_State_t    Geo_State_t;
    typedef typename Geometry_t::Space_Vector   Space_Vector;
    typedef curandState_t                       RNG_t;
    typedef def::size_type                      size_type;
    //@}

    //@{
    //! Smart pointers
    typedef cuda::Shared_Device_Ptr<Geometry_t> SDP_Geometry;
    //@}

  protected:
    // >>> DATA

    // Geometry and physics.
    SDP_Geometry b_geometry_host;
    Geometry_t  *b_geometry;

    // Sample isotropic angle.
    __device__ void sample_angle(cuda::Space_Vector &omega, RNG_t rng)
    {
        omega.z        = 1.0 - 2.0 * curand_uniform_double(&rng);
        double phi      = cuda::constants::two_pi *
                          curand_uniform_double(&rng);
        double sintheta = sqrt(1.0 - omega.z * omega.z);

        omega.x = sintheta * cos(phi);
        omega.y = sintheta * sin(phi);
    }

    // Calculate random number offsets.
    void make_RNG();

    // Node ids.
    int b_node, b_nodes;

  public:
    // Constructor.
    Source(SDP_Geometry    geometry);

    // Derived classes should implement the following functions,
    // but because this is an on-device class there can be NO virtual
    // functions.  Toggling between derived classes must be handled
    // through templates.

    //! Get a particle from the source.
    //Particle get_particle();

    //! Whether the source has finished emitting all its particles.
    //bool empty() const;

    //! Number of particles to transport in the source on the current domain.
    //size_type num_to_transport() const;

    //! Total number of particles to transport in the entire problem/cycle.
    //size_type total_num_to_transport() const;

    // >>> INHERITED INTERFACE

    //! Get the geometry.
    const Geometry_t& geometry() const { return *b_geometry; }

    //! Number of random number streams generated so far (inclusive).
    int num_streams() const { return d_rng_stream; }

  private:
    // >>> DATA

    // Offsets used for random number generator selection.
    int d_rng_stream;
};

} // end namespace cuda_mc

#endif // cuda_mc_Source_cuh

//---------------------------------------------------------------------------//
//                 end of Source.cuh
//---------------------------------------------------------------------------//

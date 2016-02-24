//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Particle.cuh
 * \author Thomas M. Evans
 * \date   Fri Apr 25 11:26:16 2014
 * \brief  Particle class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Particle_cuh
#define cuda_mc_Particle_cuh

#include "CudaUtils/cuda_utils/Definitions.hh"
#include "MC/mc/Definitions.hh"
#include "CudaUtils/cuda_utils/CudaDBC.hh"

#include <curand_kernel.h>

namespace cuda_mc
{

//===========================================================================//
/*!
 * \class Particle
 * \brief Particle class for MC transport.
 */
/*!
 * \example mc/test/tstParticle.cc
 *
 * Test of Particle.
 */
//===========================================================================//

template <class Geometry>
class Particle
{
  public:
    //@{
    //! Typedefs.
    typedef cuda::Space_Vector              Space_Vector;
    typedef profugus::events::Event         Event_Type;
    typedef typename Geometry::Geo_State_t  Geo_State_t;
    typedef curandState_t                   RNG_State;
    //@}

  private:
    // >>> DATA

    // Material id in current region.
    int d_matid;

    // Particle group index.
    int d_group;

    // Particle weight.
    double d_wt;

    // Curand state
    RNG_State d_rng;

    // Alive/dead status.
    bool d_alive;

    // Latest particle event.
    Event_Type d_event;

    // Particle geometric state.
    Geo_State_t d_geo_state;

  public:
    // Constructor
    __device__ Particle(){};

    // >>> PARTICLE FUNCTIONS

    //! Set a new weight.
    __device__ void set_wt(double wt) { d_wt = wt; }

    //! Set particle group.
    __device__ void set_group(int g) { d_group = g; }

    //! Multiply weight.
    __device__ void multiply_wt(double wt) { d_wt *= wt; }

    //! Set a new random number generator.
    __device__ void set_rng(const RNG_State &rng) { d_rng = rng; }

    //! Set the particle event flag.
    __device__ void set_event(Event_Type event) { d_event = event; }

    //! Set the material id of the region occupied by the particle.
    __device__ void set_matid(int matid) { d_matid = matid; }

    //! Kill the particle.
    __device__ void kill() { d_alive = false; }

    //! Set particle status to alive.
    __device__ void live() { d_alive = true; }

    //@{
    //! Get a handle to the geometric state of the particle.
    __device__ Geo_State_t& geo_state() { return d_geo_state; }
    __device__ const Geo_State_t& geo_state() const { return d_geo_state; }
    //@}

    //@{
    //! Access particle data.
    __device__ bool alive() const { return d_alive; }
    __device__ double wt() const { return d_wt; }
    __device__ const RNG_State& rng() const { return d_rng; }
    __device__ double ran() const
    {
        RNG_State rng_state = d_rng;
        return curand_uniform_double(&rng_state);
    }
    __device__ Event_Type event() const { return d_event; }
    __device__ int matid() const { return d_matid; }
    __device__ int group() const { return d_group; }
    //@}
};

} // end namespace cuda_mc

#endif // cuda_mc_Particle_cuh

//---------------------------------------------------------------------------//
//                 end of Particle.cuh
//---------------------------------------------------------------------------//

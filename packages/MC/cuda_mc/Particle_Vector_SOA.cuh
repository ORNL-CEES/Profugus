//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/cuda_mc/Particle_Vector_SOA.cuh
 * \author Steven Hamilton
 * \date   Tue Mar 07 09:28:41 2017
 * \brief  Particle_Vector_SOA class declaration.
 * \note   Copyright (c) 2017 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_cuda_mc_Particle_Vector_SOA_cuh
#define MC_cuda_mc_Particle_Vector_SOA_cuh

#include <thrust/device_vector.h>
#include "MC/mc/Definitions.hh"
#include "Utils/comm/Timing.hh"
#include "CudaUtils/cuda_utils/Device_View_Field.hh"
#include "CudaUtils/cuda_utils/Device_Memory_Manager.hh"
#include "CudaUtils/cuda_utils/Utility_Functions.hh"
#include "Definitions.hh"
#include "RNG_Control.cuh"

namespace cuda_mc
{

//===========================================================================//
/*!
 * \class Particle_Vector_SOA
 * \brief Class for managing collection of particle objects.
 */
//===========================================================================//

template <class Geometry>
class Particle_Vector_SOA
{
  public:
    //@{
    //! Typedefs
    typedef cuda_utils::Space_Vector                Space_Vector;
    typedef profugus::events::Event                 Event_Type;
    typedef typename Geometry::Geo_State_Vector_t   Geo_State_Vector_t;
    //@}

  private:

    // >>> DATA
    cuda::Device_View_Field<int>          d_matids;
    cuda::Device_View_Field<int>          d_groups;
    cuda::Device_View_Field<double>       d_wts;
    cuda::Device_View_Field<bool>         d_alive;
    cuda::Device_View_Field<Event_Type>   d_events;
    Geo_State_Vector_t                    d_geo_state_vec;
    RNG_Control                           d_rng_control;

  public:

    // Constructor
    Particle_Vector_SOA(cuda::Device_View_Field<int> matids,
                        cuda::Device_View_Field<int> groups,
                        cuda::Device_View_Field<double> wts,
                        cuda::Device_View_Field<bool> alive,
                        cuda::Device_View_Field<Event_Type> events,
                        Geo_State_Vector_t geo_state_vec,
                        RNG_Control rng_control)
      : d_matids(matids)
      , d_groups(groups)
      , d_wts(wts)
      , d_alive(alive)
      , d_events(events)
      , d_geo_state_vec(geo_state_vec)
      , d_rng_control(rng_control)
    {
    }

    // Number of particles currently allocated
    int size() const {return d_matids.size();}

    // >>> PARTICLE FUNCTIONS

    //! Set a new weight.
    __device__ void set_wt(int ind, double wt)
    {
        DEVICE_REQUIRE(ind < d_wts.size());
        d_wts[ind] = wt;
    }

    //! Set particle group.
    __device__ void set_group(int ind, int g)
    {
        DEVICE_REQUIRE(ind < d_groups.size());
        d_groups[ind] = g;
    }

    //! Multiply weight.
    __device__ void multiply_wt(int ind, double wt)
    {
        DEVICE_REQUIRE(ind < d_wts.size());
        d_wts[ind] *= wt;
    }

    //! Set the particle event flag.
    __device__ void set_event(int ind, Event_Type event)
    {
        DEVICE_REQUIRE(ind < d_events.size());
        d_events[ind] = event;
    }

    //! Set the material id of the region occupied by the particle.
    __device__ void set_matid(int ind, int matid)
    {
        DEVICE_REQUIRE(ind < d_matids.size());
        d_matids[ind] = matid;
    }

    //! Kill the particle.
    __device__ void kill(int ind)
    {
        DEVICE_REQUIRE(ind < d_alive.size());
        d_alive[ind] = false;
    }

    //! Set particle status to alive.
    __device__ void live(int ind)
    {
        DEVICE_REQUIRE(ind < d_alive.size());
        d_alive[ind] = true;
    }

    //@{
    //! Get a handle to the geometric state of the particle.
    __device__ Geo_State_Vector_t& geo_states()
    {
        return d_geo_state_vec;
    }
    __device__ const Geo_State_Vector_t& geo_states() const
    {
        return d_geo_state_vec;
    }
    //@}

    //@{
    //! Access particle data.
    __device__ bool alive(int ind) const
    {
        DEVICE_REQUIRE(ind < d_alive.size());
        return d_alive[ind];
    }
    __device__ double wt(int ind) const
    {
        DEVICE_REQUIRE(ind < d_wts.size());
        return d_wts[ind];
    }
    __device__ RNG_State_t * rng(int ind)
    {
        // Ignore index provided, use thread id
        return d_rng_control.get_state(cuda::utility::thread_id());
    }
    __device__ double ran(int ind)
    {
        return curand_uniform_double(rng(ind));
    }
    __device__ Event_Type event(int ind) const
    {
        DEVICE_REQUIRE(ind < d_events.size());
        return d_events[ind];
    }
    __device__ int matid(int ind) const
    {
        DEVICE_REQUIRE(ind < d_matids.size());
        return d_matids[ind];
    }
    __device__ int group(int ind) const
    {
        DEVICE_REQUIRE(ind < d_groups.size());
        return d_groups[ind];
    }
    //@}

};

//===========================================================================//
/*!
 * \class Particle_Vector_SOA_DMM
 * \brief Memory manager for Particle_Vector_SOA.
 */
//===========================================================================//

template <class Geometry>
class Particle_Vector_SOA_DMM :
    public cuda::Device_Memory_Manager<Particle_Vector_SOA<Geometry>>
{
  public:

    typedef Particle_Vector_SOA<Geometry>             Particle_Vector_SOA_t;
    typedef profugus::events::Event                   Event_Type;
    typedef typename Geometry::Geo_State_Vector_t     Geo_State_Vector_t;
    typedef typename Geometry::Geo_State_Vector_DMM_t Geo_State_Vector_DMM_t;

    // Constructor
    Particle_Vector_SOA_DMM(int rng_seed)
      : d_rng_control(rng_seed)
    {}

    // Number of particles currently allocated
    int size() const {return d_matids.size();}

    // Memory manager interface
    Particle_Vector_SOA_t device_instance()
    {
        return Particle_Vector_SOA_t(cuda::make_view(d_matids),
                                     cuda::make_view(d_groups),
                                     cuda::make_view(d_wts),
                                     cuda::make_view(d_alive),
                                     cuda::make_view(d_events),
                                     d_geo_states_dmm.device_instance(),
                                     d_rng_control.device_instance());
    }

    // Initialize vector for specified number of states
    void initialize(int num_states)
    {
        SCOPED_TIMER("MC::Particle_Vector::initialize");
        REQUIRE(num_states > 0);
        d_matids.resize(    num_states);
        d_groups.resize(    num_states);
        d_wts.resize(       num_states);
        d_alive.resize(     num_states);
        d_events.resize(    num_states);
        d_geo_states_dmm.initialize(num_states);
        d_rng_control.initialize(num_states);
    }

  private:

    thrust::device_vector<int>          d_matids;
    thrust::device_vector<int>          d_groups;
    thrust::device_vector<double>       d_wts;
    thrust::device_vector<bool>         d_alive;
    thrust::device_vector<Event_Type>   d_events;
    Geo_State_Vector_DMM_t              d_geo_states_dmm;
    RNG_Control_DMM                     d_rng_control;
};

//---------------------------------------------------------------------------//
} // end namespace cuda_mc

//---------------------------------------------------------------------------//
#endif // MC_cuda_mc_Particle_Vector_SOA_cuh

//---------------------------------------------------------------------------//
// end of MC/cuda_mc/Particle_Vector_SOA.cuh
//---------------------------------------------------------------------------//

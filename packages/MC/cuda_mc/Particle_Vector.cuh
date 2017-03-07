//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/cuda_mc/Particle_Vector.cuh
 * \author Steven Hamilton
 * \date   Tue Mar 07 09:28:41 2017
 * \brief  Particle_Vector class declaration.
 * \note   Copyright (c) 2017 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_cuda_mc_Particle_Vector_cuh
#define MC_cuda_mc_Particle_Vector_cuh

#include <thrust/device_vector.h>
#include "Particle.cuh"
#include "CudaUtils/cuda_utils/Device_View_Field.hh"
#include "CudaUtils/cuda_utils/Device_Memory_Manager.hh"

namespace cuda_mc
{

//===========================================================================//
/*!
 * \class Particle_Vector
 * \brief Class for managing collection of particle objects.
 */
/*!
 * \example cuda_mc/test/tstParticle_Vector.cc
 *
 * Test of Particle_Vector.
 */
//===========================================================================//

template <class Geometry>
class Particle_Vector
{
  public:
    //@{
    //! Typedefs
    typedef cuda_utils::Space_Vector            Space_Vector;
    typedef profugus::events::Event             Event_Type;
    typedef typename Geometry::Geo_State_t      Geo_State_t;
    typedef Particle<Geometry>                  Particle_t;
    typedef cuda::Device_View_Field<Particle_t> Particle_View;
    //@}

  private:

    // >>> DATA
    Particle_View d_particles;

  public:

    // Constructor
    Particle_Vector(Particle_View particles)
      : d_particles(particles)
    {
    }

    // Number of particles currently allocated
    int size() const {return d_particles.size();}

    // >>> PARTICLE FUNCTIONS

    //! Set a new weight.
    __device__ void set_wt(int ind, double wt)
    {
        DEVICE_REQUIRE(ind < d_particles.size());
        d_particles[ind].set_wt(wt);
    }

    //! Set particle group.
    __device__ void set_group(int ind, int g)
    {
        DEVICE_REQUIRE(ind < d_particles.size());
        d_particles[ind].set_group(g);
    }

    //! Multiply weight.
    __device__ void multiply_wt(int ind, double wt)
    {
        DEVICE_REQUIRE(ind < d_particles.size());
        d_particles[ind].multiply_wt(wt);
    }

    //! Set a new random number generator.
    __device__ void set_rng(int ind, RNG_State_t *rng)
    {
        DEVICE_REQUIRE(ind < d_particles.size());
        d_particles[ind].set_rng(rng);
    }

    //! Set the particle event flag.
    __device__ void set_event(int ind, Event_Type event)
    {
        DEVICE_REQUIRE(ind < d_particles.size());
        d_particles[ind].set_event(event);
    }

    //! Set the material id of the region occupied by the particle.
    __device__ void set_matid(int ind, int matid)
    {
        DEVICE_REQUIRE(ind < d_particles.size());
        d_particles[ind].set_matid(matid);
    }

    //! Kill the particle.
    __device__ void kill(int ind)
    {
        DEVICE_REQUIRE(ind < d_particles.size());
        d_particles[ind].kill();
    }

    //! Set particle status to alive.
    __device__ void live(int ind)
    {
        DEVICE_REQUIRE(ind < d_particles.size());
        d_particles[ind].live();
    }

    //@{
    //! Get a handle to the geometric state of the particle.
    __device__ Geo_State_t& geo_state(int ind)
    {
        DEVICE_REQUIRE(ind < d_particles.size());
        return d_particles[ind].geo_state();
    }
    __device__ const Geo_State_t& geo_state(int ind) const
    {
        DEVICE_REQUIRE(ind < d_particles.size());
        return d_particles[ind].geo_state();
    }
    //@}

    //@{
    //! Access particle data.
    __device__ bool alive(int ind) const
    {
        DEVICE_REQUIRE(ind < d_particles.size());
        return d_particles[ind].alive();
    }
    __device__ double wt(int ind) const
    {
        DEVICE_REQUIRE(ind < d_particles.size());
        return d_particles[ind].wt();
    }
    __device__ RNG_State_t * rng(int ind) const
    {
        DEVICE_REQUIRE(ind < d_particles.size());
        return d_particles[ind].rng();
    }
    __device__ double ran(int ind) const
    {
        DEVICE_REQUIRE(ind < d_particles.size());
        return d_particles[ind].ran();
    }
    __device__ Event_Type event(int ind) const
    {
        DEVICE_REQUIRE(ind < d_particles.size());
        return d_particles[ind].event();
    }
    __device__ int matid(int ind) const
    {
        DEVICE_REQUIRE(ind < d_particles.size());
        return d_particles[ind].matid();
    }
    __device__ int group(int ind) const
    {
        DEVICE_REQUIRE(ind < d_particles.size());
        return d_particles[ind].group();
    }
    //@}

};

//===========================================================================//
/*!
 * \class Particle_Vector_DMM
 * \brief Memory manager for Particle_Vector.
 */
//===========================================================================//

template <class Geometry>
class Particle_Vector_DMM :
    public cuda::Device_Memory_Manager<Particle_Vector<Geometry>>
{
  public:

    typedef Particle_Vector<Geometry> Particle_Vector_t;
    typedef Particle<Geometry>        Particle_t;

    // Constructor
    Particle_Vector_DMM(){}

    // Number of particles currently allocated
    int size() const {return d_particles.size();}

    // Memory manager interface
    Particle_Vector_t device_instance()
    {
        return Particle_Vector_t(cuda::make_view(d_particles));
    }

    // Initialize vector for specified number of states
    void initialize(int num_states)
    {
        REQUIRE(num_states > 0);
        d_particles.resize(num_states);
    }

  private:

    thrust::device_vector<Particle_t> d_particles;

};

//---------------------------------------------------------------------------//
} // end namespace cuda_mc

//---------------------------------------------------------------------------//
#endif // MC_cuda_mc_Particle_Vector_cuh

//---------------------------------------------------------------------------//
// end of MC/cuda_mc/Particle_Vector.cuh
//---------------------------------------------------------------------------//

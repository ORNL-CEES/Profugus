//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Particle_Vector.hh
 * \author Stuart Slattery
 * \brief  Particle vector class definition.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Particle_Vector_hh
#define cuda_mc_Particle_Vector_hh

#include "utils/Definitions.hh"

#include "cuda_utils/CudaMacros.hh"
#include "cuda_utils/CudaDBC.hh"

#include "rng/RNG.hh"

#include "Definitions.hh"

#include <Teuchos_Array.hpp>

#include <curand_kernel.h>

#include <thrust/device_vector.h>

namespace cuda_profugus
{
//===========================================================================//
/*!
 * \class Particle_Vector
 * \brief Particle_Vector class for MC transport.
 */
/*!
 * \example mc/test/tstParticle_Vector.cc
 *
 * Test of Particle.
 */
//===========================================================================//
template <class Geometry>
class Particle_Vector
{
  public:
    //@{
    //! Typedefs.
    typedef typename Geometry::Geo_State_t Geo_State_t;
    typedef events::Event Event_t;
    //@}

  private:
    // >>> DATA

    // Number of particles in the vector.
    int d_size;

    // Material id in current region.
    int* d_matid;

    // Particle group index.
    int* d_group;

    // Particle weight.
    double* d_wt;

    // Random number generator.
    curandState* d_rng;

    // Alive/dead status.
    bool* d_alive;

    // Particle geometric state.
    Geo_State_t* d_geo_state;

    // Latest particle event.
    Event_t* d_event;

    // Particle statistical batch id.
    int* d_batch;

    // Distance of last particle step.
    double* d_step;

    // Distance to next collision in mean-free-paths.
    double* d_dist_mfp;

    // Event indices.
    int* d_event_indices;

    // Number of particles with a given event. Host only.
    Teuchos::Array<int> d_event_sizes;

    // Local id vector for sorting.
    thrust::device_vector<int> d_event_lid;    

    // Stencil vector for sorting.
    thrust::device_vector<int> d_event_stencil;

    // Steering vector for sorting.
    thrust::device_vector<int> d_event_steering;

  public:

    // >>> HOST API

    // Constructor
    Particle_Vector( const int num_particle, const profugus::RNG& rng );
    
    // Destructor.
    ~Particle_Vector();

    // Sort the local indices by event key, effectively sorting the vector
    // using the default stream.
    void sort_by_event( const int sort_size );

    // Get the number of particles with a given event on the host.
    int get_event_size( const events::Event event ) const;

    // >>> DEVICE API

    //! Get the number of particles in the vector.
    PROFUGUS_HOST_DEVICE_FUNCTION
    int size() const { return d_size; }

    //! Get the particle indices with a given event.
    PROFUGUS_DEVICE_FUNCTION
    int* event_indices( const events::Event event ) const
    {
        REQUIRE( event < events::END_EVENT );
        return d_event_indices + event*d_size;
    }

    //! Get a uniform random number on [0,1] from a particle's RNG.
    PROFUGUS_DEVICE_FUNCTION
    double ran( const std::size_t i )
    {
	REQUIRE( i < d_size );
	return curand_uniform( &d_rng[i] );
    }

    //! Set the weight of a particle.
    PROFUGUS_DEVICE_FUNCTION
    void set_wt( const std::size_t i, const double wt ) const 
    { 
	REQUIRE( i < d_size );
	d_wt[i] = wt; 
    }

    //! Multiply the weight of a particle.
    PROFUGUS_DEVICE_FUNCTION
    void multiply_wt( const std::size_t i, const double wt ) const 
    { 
	REQUIRE( i < d_size );
	d_wt[i] *= wt; 
    }

    //! Get the weight of a particle.
    PROFUGUS_DEVICE_FUNCTION
    double wt( const std::size_t i ) const 
    { 
	REQUIRE( i < d_size );
	return d_wt[i]; 
    }

    //! Get the group of a particle.
    PROFUGUS_DEVICE_FUNCTION
    int group( const std::size_t i ) const 
    { 
	REQUIRE( i < d_size );
	return d_group[i]; 
    }

    //! Set the group of a particle.
    PROFUGUS_DEVICE_FUNCTION
    void set_group( const std::size_t i, const int group )
    { 
	REQUIRE( i < d_size );
	d_group[i] = group; 
    }

    //! Get the matid of a particle.
    PROFUGUS_DEVICE_FUNCTION
    int matid( const std::size_t i ) const 
    { 
	REQUIRE( i < d_size );
	return d_matid[i]; 
    }

    //! Set the matid of a particle.
    PROFUGUS_DEVICE_FUNCTION
    void set_matid( const std::size_t i, const int matid )
    { 
	REQUIRE( i < d_size );
	d_matid[i] = matid;
    }

    //! Get the alive status of a particle.
    PROFUGUS_DEVICE_FUNCTION
    bool alive( const std::size_t i ) const 
    { 
	REQUIRE( i < d_size );
	return d_alive[i]; 
    }

    //! Set the particle status to alive.
    PROFUGUS_DEVICE_FUNCTION
    void live( const std::size_t i )
    { 
	REQUIRE( i < d_size );
	d_alive[i] = true; 
    }

    //! Kill a particle.
    PROFUGUS_DEVICE_FUNCTION
    void kill( const std::size_t i )
    { 
	REQUIRE( i < d_size );
	d_alive[i] = false; 
    }

    //! Get the event of a particle.
    PROFUGUS_DEVICE_FUNCTION
    Event_t event( const std::size_t i ) const 
    { 
	REQUIRE( i < d_size );
	return d_event[i];
    }

    //! Set the event of a particle.
    PROFUGUS_DEVICE_FUNCTION
    void set_event( const std::size_t i, const Event_t event )
    { 
	REQUIRE( i < d_size );
	d_event[i] = event; 
    }

    //! Get the geometry state of a particle.
    PROFUGUS_DEVICE_FUNCTION
    const Geo_State_t& geo_state( const std::size_t i ) const 
    { 
	REQUIRE( i < d_size );
	return d_geo_state[i]; 
    }
    PROFUGUS_DEVICE_FUNCTION
    Geo_State_t& geo_state( const std::size_t i )
    { 
	REQUIRE( i < d_size );
	return d_geo_state[i]; 
    }
    //@}

    //! Get the particle batch.
    PROFUGUS_DEVICE_FUNCTION
    int batch( const std::size_t i ) const
    {
	REQUIRE( i < d_size );
	return d_batch[i];
    }

    //! Set the particle batch.
    PROFUGUS_DEVICE_FUNCTION
    void set_batch( const std::size_t i, const int batch )
    {
	REQUIRE( i < d_size );
	d_batch[i] = batch;
    }

    //! Get the particle step.
    PROFUGUS_DEVICE_FUNCTION
    double step( const std::size_t i ) const
    {
	REQUIRE( i < d_size );
	return d_step[i];
    }

    //! Set the particle step.
    PROFUGUS_DEVICE_FUNCTION
    void set_step( const std::size_t i, const double step )
    {
	REQUIRE( i < d_size );
	d_step[i] = step;
    }

    //! Get the particle distance to collision.
    PROFUGUS_DEVICE_FUNCTION
    double dist_mfp( const std::size_t i ) const
    {
	REQUIRE( i < d_size );
	return d_dist_mfp[i];
    }

    //! Set the particle distance to collision.
    PROFUGUS_DEVICE_FUNCTION
    void set_dist_mfp( const std::size_t i, const double dist_mfp )
    {
	REQUIRE( i < d_size );
	d_dist_mfp[i] = dist_mfp;
    }

  private:

    // Event stencil functor.
    struct Stencil_Functor
    {
        events::Event d_event;

        PROFUGUS_HOST_DEVICE_FUNCTION
        int operator()(event::Event e )
        { return (e == d_event) 1 : 0; }
    };
};

//---------------------------------------------------------------------------//

} // end namespace cuda_profugus

#endif // cuda_mc_Particle_Vector_hh

//---------------------------------------------------------------------------//
//                 end of Particle_Vector.hh
//---------------------------------------------------------------------------//

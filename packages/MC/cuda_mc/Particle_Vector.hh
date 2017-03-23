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

    // Particle local index linked to the events.
    int* d_lid;

    // Particle statistical batch id.
    int* d_batch;

    // Distance of last particle step.
    double* d_step;

    // Distance to next collision in mean-free-paths.
    double* d_dist_mfp;

    // Event bounds. Updated at every sort.
    int* d_event_bounds;

    // Number of particles with a given event.
    int* d_num_event;

    // Event bins.
    int* d_event_bins;

    // Number of particles with a given event. Host only.
    Teuchos::Array<int> d_event_sizes;

  public:

    // >>> HOST API

    // Constructor
    Particle_Vector( const int num_particle, const profugus::RNG& rng );

    // Destructor.
    ~Particle_Vector();

    // Return if the vector is empty. If the vector has no events of any kind
    // - not event dead particles, it is empty. We initialize everything to
    // dead so it is not empty.
    bool empty() const;

    // Get the number of particles that are not dead.
    int num_alive() const;

    // Sort the local indices by event key, effectively sorting the vector
    // using the default stream.
    void sort_by_event( const int sort_size );

    // Get the number of particles with a given event on the host.
    int get_event_size( const events::Event event ) const;

    // >>> DEVICE API

    //! Get the number of particles in the vector.
    PROFUGUS_HOST_DEVICE_FUNCTION
    int size() const { return d_size; }

    //! Get the lower bound of an event in the vector.
    PROFUGUS_DEVICE_FUNCTION
    int event_lower_bound( const events::Event event ) const
    {
        DEVICE_REQUIRE( event < events::END_EVENT );
        return d_event_bounds[ event ];
    }

    //! Get the particle indices with a given event.
    PROFUGUS_DEVICE_FUNCTION
    int* event_indices( const events::Event event ) const
    {
        DEVICE_REQUIRE( event < events::END_EVENT );
        return d_lid + event_lower_bound(event);
    }

    //! Get a uniform random number on [0,1] from a particle's RNG.
    PROFUGUS_DEVICE_FUNCTION
    double ran( const int i )
    {
	DEVICE_REQUIRE( i < d_size );
	return curand_uniform( &d_rng[i] );
    }

    //! Set the weight of a particle.
    PROFUGUS_DEVICE_FUNCTION
    void set_wt( const int i, const double wt ) const
    {
	DEVICE_REQUIRE( i < d_size );
	d_wt[i] = wt;
    }

    //! Multiply the weight of a particle.
    PROFUGUS_DEVICE_FUNCTION
    void multiply_wt( const int i, const double wt ) const
    {
	DEVICE_REQUIRE( i < d_size );
	d_wt[i] *= wt;
    }

    //! Get the weight of a particle.
    PROFUGUS_DEVICE_FUNCTION
    double wt( const int i ) const
    {
	DEVICE_REQUIRE( i < d_size );
	return d_wt[i];
    }

    //! Get the group of a particle.
    PROFUGUS_DEVICE_FUNCTION
    int group( const int i ) const
    {
	DEVICE_REQUIRE( i < d_size );
	return d_group[i];
    }

    //! Set the group of a particle.
    PROFUGUS_DEVICE_FUNCTION
    void set_group( const int i, const int group )
    {
	DEVICE_REQUIRE( i < d_size );
	d_group[i] = group;
    }

    //! Get the matid of a particle.
    PROFUGUS_DEVICE_FUNCTION
    int matid( const int i ) const
    {
	DEVICE_REQUIRE( i < d_size );
	return d_matid[i];
    }

    //! Set the matid of a particle.
    PROFUGUS_DEVICE_FUNCTION
    void set_matid( const int i, const int matid )
    {
	DEVICE_REQUIRE( i < d_size );
	d_matid[i] = matid;
    }

    //! Get the alive status of a particle.
    PROFUGUS_DEVICE_FUNCTION
    bool alive( const int i ) const
    {
	DEVICE_REQUIRE( i < d_size );
	return d_alive[i];
    }

    //! Set the particle status to alive.
    PROFUGUS_DEVICE_FUNCTION
    void live( const int i )
    {
	DEVICE_REQUIRE( i < d_size );
	d_alive[i] = true;
    }

    //! Kill a particle.
    PROFUGUS_DEVICE_FUNCTION
    void kill( const int i )
    {
	DEVICE_REQUIRE( i < d_size );
	d_alive[i] = false;
    }

    //! Get the event of a particle.
    PROFUGUS_DEVICE_FUNCTION
    Event_t event( const int i ) const
    {
	DEVICE_REQUIRE( i < d_size );
	return d_event[i];
    }

    //! Set the event of a particle.
    PROFUGUS_DEVICE_FUNCTION
    void set_event( const int i, const Event_t event )
    {
#ifdef __NVCC__
	DEVICE_REQUIRE( i < d_size );

        // Set the event.
	d_event[i] = event;

        // Bin this particle with its event.
        int* address = &d_num_event[event];
        int bin_size = atomicAdd( address, 1 );
        d_event_bins[ event*d_size + bin_size ] = i;
#endif
    }

    //! Get the geometry state of a particle.
    PROFUGUS_DEVICE_FUNCTION
    const Geo_State_t& geo_state( const int i ) const
    {
	DEVICE_REQUIRE( i < d_size );
	return d_geo_state[i];
    }
    PROFUGUS_DEVICE_FUNCTION
    Geo_State_t& geo_state( const int i )
    {
	DEVICE_REQUIRE( i < d_size );
	return d_geo_state[i];
    }
    //@}

    //! Get the particle batch.
    PROFUGUS_DEVICE_FUNCTION
    int batch( const int i ) const
    {
	DEVICE_REQUIRE( i < d_size );
	return d_batch[i];
    }

    //! Set the particle batch.
    PROFUGUS_DEVICE_FUNCTION
    void set_batch( const int i, const int batch )
    {
	DEVICE_REQUIRE( i < d_size );
	d_batch[i] = batch;
    }

    //! Get the particle step.
    PROFUGUS_DEVICE_FUNCTION
    double step( const int i ) const
    {
	DEVICE_REQUIRE( i < d_size );
	return d_step[i];
    }

    //! Set the particle step.
    PROFUGUS_DEVICE_FUNCTION
    void set_step( const int i, const double step )
    {
	DEVICE_REQUIRE( i < d_size );
	d_step[i] = step;
    }

    //! Get the particle distance to collision.
    PROFUGUS_DEVICE_FUNCTION
    double dist_mfp( const int i ) const
    {
	DEVICE_REQUIRE( i < d_size );
	return d_dist_mfp[i];
    }

    //! Set the particle distance to collision.
    PROFUGUS_DEVICE_FUNCTION
    void set_dist_mfp( const int i, const double dist_mfp )
    {
	DEVICE_REQUIRE( i < d_size );
	d_dist_mfp[i] = dist_mfp;
    }
};

//---------------------------------------------------------------------------//

} // end namespace cuda_profugus

#endif // cuda_mc_Particle_Vector_hh

//---------------------------------------------------------------------------//
//                 end of Particle_Vector.hh
//---------------------------------------------------------------------------//

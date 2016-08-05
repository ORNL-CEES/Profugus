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

#include "Particle.hh"

#include "utils/Definitions.hh"

#include "cuda_utils/CudaMacros.hh"
#include "cuda_utils/CudaDBC.hh"

#include "rng/RNG.hh"

#include "Definitions.hh"

#include <Teuchos_Array.hpp>

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
    typedef Particle<Geometry> Particle_t;
    //@}

  private:
    // >>> DATA

    // Number of particles in the vector.
    int d_size;

    // Particles.
    Particle_t* d_particles;

    // Particle local index linked to the events.
    int* d_lid;

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
    // - not even dead particles, it is empty. We initialize everything to
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
        REQUIRE( event < events::END_EVENT );
        return d_event_bounds[ event ];
    }

    //! Get a particle.
    PROFUGUS_DEVICE_FUNCTION
    Particle_t& particle( const int i )
    {
	REQUIRE( i < d_size );
	REQUIRE( d_lid[i] < d_size );
        return d_particles[ d_lid[i] ];
    }

    //! Get a particle.
    PROFUGUS_DEVICE_FUNCTION
    const Particle_t& particle( const int i ) const
    {
	REQUIRE( i < d_size );
	REQUIRE( d_lid[i] < d_size );
        return d_particles[ d_lid[i] ];
    }
    
    //! Get a uniform random number on [0,1] from a particle's RNG.
    PROFUGUS_DEVICE_FUNCTION
    double ran( const int i )
    {
	return particle(i).ran();
    }

    //! Set the weight of a particle.
    PROFUGUS_DEVICE_FUNCTION
    void set_wt( const int i, const double wt ) const 
    { 
        particle(i).set_wt( wt );
    }

    //! Multiply the weight of a particle.
    PROFUGUS_DEVICE_FUNCTION
    void multiply_wt( const int i, const double wt ) const 
    { 
        particle(i).multiply_wt( wt );
    }

    //! Get the weight of a particle.
    PROFUGUS_DEVICE_FUNCTION
    double wt( const int i ) const 
    { 
        return particle(i).wt();
    }

    //! Get the group of a particle.
    PROFUGUS_DEVICE_FUNCTION
    int group( const int i ) const 
    { 
	return particle(i).group();
    }

    //! Set the group of a particle.
    PROFUGUS_DEVICE_FUNCTION
    void set_group( const int i, const int group )
    { 
        particle(i).set_group( group );
    }

    //! Get the matid of a particle.
    PROFUGUS_DEVICE_FUNCTION
    int matid( const int i ) const 
    { 
        return particle(i).matid();
    }

    //! Set the matid of a particle.
    PROFUGUS_DEVICE_FUNCTION
    void set_matid( const int i, const int matid )
    { 
        particle(i).set_matid( matid );
    }

    //! Get the alive status of a particle.
    PROFUGUS_DEVICE_FUNCTION
    bool alive( const int i ) const 
    { 
        return particle(i).alive();
    }

    //! Set the particle status to alive.
    PROFUGUS_DEVICE_FUNCTION
    void live( const int i )
    { 
        particle(i).live();
    }

    //! Kill a particle.
    PROFUGUS_DEVICE_FUNCTION
    void kill( const int i )
    { 
        particle(i).kill();
    }

    //! Get the event of a particle.
    PROFUGUS_DEVICE_FUNCTION
    Event_t event( const int i ) const 
    { 
        return particle(i).event();
    }

    //! Set the event of a particle.
    PROFUGUS_DEVICE_FUNCTION
    void set_event( const int i, const Event_t event )
    { 
#ifdef __NVCC__
	REQUIRE( i < d_size );

        // Set the event.
        particle(i).set_event( event );

        // Bin this particle with its event.
        int* address = &d_num_event[event];
        int bin_size = atomicAdd( address, 1 );
        d_event_bins[ event*d_size + bin_size ] = d_lid[i];
#endif
    }

    //! Get the geometry state of a particle.
    PROFUGUS_DEVICE_FUNCTION
    const Geo_State_t& geo_state( const int i ) const 
    { 
	return particle(i).geo_state();
    }
    PROFUGUS_DEVICE_FUNCTION
    Geo_State_t& geo_state( const int i )
    { 
	return particle(i).geo_state();        
    }
    //@}

    //! Get the particle batch.
    PROFUGUS_DEVICE_FUNCTION
    int batch( const int i ) const
    {
        return particle(i).batch();
    }

    //! Set the particle batch.
    PROFUGUS_DEVICE_FUNCTION
    void set_batch( const int i, const int batch )
    {
        particle(i).set_batch( batch );
    }

    //! Get the particle step.
    PROFUGUS_DEVICE_FUNCTION
    double step( const int i ) const
    {
        return particle(i).step();
    }

    //! Set the particle step.
    PROFUGUS_DEVICE_FUNCTION
    void set_step( const int i, const double step )
    {
        particle(i).set_step( step );
    }

    //! Get the particle distance to collision.
    PROFUGUS_DEVICE_FUNCTION
    double dist_mfp( const int i ) const
    {
        return particle(i).dist_mfp();
    }

    //! Set the particle distance to collision.
    PROFUGUS_DEVICE_FUNCTION
    void set_dist_mfp( const int i, const double dist_mfp )
    {
        particle(i).set_dist_mfp( dist_mfp );
    }
};

//---------------------------------------------------------------------------//

} // end namespace cuda_profugus

#endif // cuda_mc_Particle_Vector_hh

//---------------------------------------------------------------------------//
//                 end of Particle_Vector.hh
//---------------------------------------------------------------------------//

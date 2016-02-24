//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   cuda_mc/test/Particle_Vector_Tester.hh
 * \author Stuart Slattery
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_test_Particle_Vector_Tester_hh
#define cuda_mc_test_Particle_Vector_Tester_hh

#include "cuda_utils/Shared_Device_Ptr.hh"

#include "rng/RNG.hh"

#include "cuda_geometry/Mesh_Geometry.hh"

#include "../Particle_Vector.hh"

#include <Teuchos_Array.hpp>

//---------------------------------------------------------------------------//
class Particle_Vector_Tester
{
  public:

    typedef cuda_profugus::Particle_Vector<cuda_profugus::Mesh_Geometry> Particle_Vector;
    typedef typename Particle_Vector::Event_t Event_t;
    typedef typename Particle_Vector::Geo_State_t Geo_State_t;

    Particle_Vector_Tester( const int num_particle, const profugus::RNG& rng );

    // get the underlying vector.
    cuda::Shared_Device_Ptr<Particle_Vector>& get_vector()
    { return d_vector; }

    // get the size of the vector
    int size() const { return d_size; }

    // get a vector of random numbers for the vector.
    Teuchos::Array<double> ran();

    // set the entire vector to the same weight.
    void set_wt( const double wt );

    // mulitply the entire vector by a weight.
    void multiply_wt( const Teuchos::Array<double>& wt );

    // get a vector of weights.
    Teuchos::Array<double> wt();

    // get a vector of groups.
    Teuchos::Array<int> group();

    // Set the entire vector to a group.
    void set_group( const int group );

    // get a vector of matids.
    Teuchos::Array<int> matid();

    // Set the entire vector to a matid.
    void set_matid( const int matid );

    // get a vector of alive status.
    Teuchos::Array<int> alive();

    // set the whole vector to live.
    void live();

    // kill the whole vector.
    void kill();

    // get a vector of events
    Teuchos::Array<Event_t> event();

    // set a vector of events
    void set_event( const Teuchos::Array<Event_t>& events );

    // sort the vector by event.
    void sort_by_event() { d_vector.get_host_ptr()->sort_by_event(); }

    // get the particles with an event.
    void get_event_particles( const Event_t event, 
			      std::size_t& start_index,
			      std::size_t& num_particle ) const
    { d_vector.get_host_ptr()->get_event_particles(event,start_index,num_particle); }

    // set the geometry state for the whole vector
    void set_geo_state( const Geo_State_t& geo_state );
    
    // get the geometry state for the whole vector
    Teuchos::Array<Geo_State_t> geo_state();

    // get a vector of batches.
    Teuchos::Array<int> batch();

    // Set the entire vector to a batch.
    void set_batch( const int batch );

  private:
    
    int d_size;

    cuda::Shared_Device_Ptr<Particle_Vector> d_vector;
};

//---------------------------------------------------------------------------//

#endif // cuda_mc_test_Particle_Vector_Tester_hh

//---------------------------------------------------------------------------//
//                 end of cuda_utils/Particle_Vector_Tester.hh
//---------------------------------------------------------------------------//

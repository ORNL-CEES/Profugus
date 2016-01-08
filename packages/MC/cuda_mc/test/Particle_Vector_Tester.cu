//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   cuda_mc/test/Particle_Vector_Tester.cu
 * \author Stuart Slattery
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Particle_Vector_Tester.hh"

#include "cuda_utils/Hardware.hh"

#include <cuda_runtime.h>

//---------------------------------------------------------------------------//
// CUDA Kernels
//---------------------------------------------------------------------------//
__global__ void ran_kernel( Particle_Vector_Tester::Particle_Vector* vector, 
			    double* ran )
{
    int i = threadIdx.x;
    ran[i] = vector->ran( i );
}

//---------------------------------------------------------------------------//
__global__ void set_wt_kernel( Particle_Vector_Tester::Particle_Vector* vector, 
			       double wt )
{
    int i = threadIdx.x;
    vector->set_wt( i, wt );
}

//---------------------------------------------------------------------------//
__global__ void multiply_wt_kernel( Particle_Vector_Tester::Particle_Vector* vector, 
				    double* wt )
{
    int i = threadIdx.x;
    vector->multiply_wt( i, wt[i] );
}

//---------------------------------------------------------------------------//
__global__ void wt_kernel( Particle_Vector_Tester::Particle_Vector* vector, 
			   double* wt )
{
    int i = threadIdx.x;
    wt[i] = vector->wt( i );
}

//---------------------------------------------------------------------------//
__global__ void group_kernel( Particle_Vector_Tester::Particle_Vector* vector, 
			      int* group )
{
    int i = threadIdx.x;
    group[i] = vector->group( i );
}

//---------------------------------------------------------------------------//
__global__ void set_group_kernel( Particle_Vector_Tester::Particle_Vector* vector, 
				  int group )
{
    int i = threadIdx.x;
    vector->set_group( i, group );
}

//---------------------------------------------------------------------------//
__global__ void matid_kernel( Particle_Vector_Tester::Particle_Vector* vector, 
			      int* matid )
{
    int i = threadIdx.x;
    matid[i] = vector->matid( i );
}

//---------------------------------------------------------------------------//
__global__ void set_matid_kernel( Particle_Vector_Tester::Particle_Vector* vector, 
				  int matid )
{
    int i = threadIdx.x;
    vector->set_matid( i, matid );
}

//---------------------------------------------------------------------------//
__global__ void alive_kernel( Particle_Vector_Tester::Particle_Vector* vector, 
			      int* alive )
{
    int i = threadIdx.x;
    alive[i] = vector->alive( i );
}

//---------------------------------------------------------------------------//
__global__ void live_kernel( Particle_Vector_Tester::Particle_Vector* vector )
{
    int i = threadIdx.x;
    vector->live( i );
}

//---------------------------------------------------------------------------//
__global__ void kill_kernel( Particle_Vector_Tester::Particle_Vector* vector )
{
    int i = threadIdx.x;
    vector->kill( i );
}

//---------------------------------------------------------------------------//
__global__ void set_event_kernel( Particle_Vector_Tester::Particle_Vector* vector, 
				  typename Particle_Vector_Tester::Event_t* event )
{
    typename std::size_t i = threadIdx.x;
    vector->set_event( i, event[i] );
}

//---------------------------------------------------------------------------//
__global__ void event_kernel( Particle_Vector_Tester::Particle_Vector* vector, 
			      typename Particle_Vector_Tester::Event_t* event )
{
    typename std::size_t i = threadIdx.x;
    event[i] = vector->event( i );
}

//---------------------------------------------------------------------------//
__global__ void set_geo_state_kernel( Particle_Vector_Tester::Particle_Vector* vector, 
				      typename Particle_Vector_Tester::Geo_State_t geo_state )
{
    typename std::size_t i = threadIdx.x;
    vector->geo_state( i ) = geo_state;
}

//---------------------------------------------------------------------------//
__global__ void geo_state_kernel( Particle_Vector_Tester::Particle_Vector* vector, 
				  typename Particle_Vector_Tester::Geo_State_t* geo_state )
{
    typename std::size_t i = threadIdx.x;
    geo_state[i] = vector->geo_state( i );
}

//---------------------------------------------------------------------------//
// Particle_Vector_Tester
//---------------------------------------------------------------------------//
Particle_Vector_Tester::Particle_Vector_Tester( const int num_particle, 
						const profugus::RNG& rng )
    : d_size( num_particle )
{
    // Acquire hardware for the test.
    cuda::Hardware<cuda::arch::Device>::acquire();

    // Create the vector after hardware acquisition.
    d_vector = cuda::Shared_Device_Ptr<Particle_Vector>( num_particle, rng );
}

//---------------------------------------------------------------------------//
// get a vector of random numbers for the vector.
Teuchos::Array<double> Particle_Vector_Tester::ran()
{
    double* device_ran;
    cudaMalloc( (void**) &device_ran, d_size * sizeof(double) );
    ran_kernel<<<1,d_size>>>( d_vector.get_device_ptr(), device_ran );

    Teuchos::Array<double> host_ran( d_size );
    cudaMemcpy( host_ran.getRawPtr(), device_ran, d_size * sizeof(double),
		cudaMemcpyDeviceToHost );

    cudaFree( device_ran );
    return host_ran;
}

//---------------------------------------------------------------------------//
// set the entire vector to the same weight.
void Particle_Vector_Tester::set_wt( const double wt )
{
    set_wt_kernel<<<1,d_size>>>( d_vector.get_device_ptr(), wt );
}

//---------------------------------------------------------------------------//
// mulitply each particle weight by a different value.
void Particle_Vector_Tester::multiply_wt( const Teuchos::Array<double>& wt )
{
    double* device_wt;
    cudaMalloc( (void**) &device_wt, d_size * sizeof(double) );
    cudaMemcpy( device_wt, wt.getRawPtr(), d_size * sizeof(double),
		cudaMemcpyHostToDevice );
    
    multiply_wt_kernel<<<1,d_size>>>( d_vector.get_device_ptr(), device_wt );

    cudaFree( device_wt );
}

//---------------------------------------------------------------------------//
// get a vector of weights.
Teuchos::Array<double> Particle_Vector_Tester::wt()
{
    double* device_wt;
    cudaMalloc( (void**) &device_wt, d_size * sizeof(double) );
    wt_kernel<<<1,d_size>>>( d_vector.get_device_ptr(), device_wt );

    Teuchos::Array<double> host_wt( d_size );
    cudaMemcpy( host_wt.getRawPtr(), device_wt, d_size * sizeof(double),
		cudaMemcpyDeviceToHost );

    cudaFree( device_wt );
    return host_wt;
}

//---------------------------------------------------------------------------//
// get a vector of groups.
Teuchos::Array<int> Particle_Vector_Tester::group()
{
    int* device_group;
    cudaMalloc( (void**) &device_group, d_size * sizeof(int) );
    group_kernel<<<1,d_size>>>( d_vector.get_device_ptr(), device_group );

    Teuchos::Array<int> host_group( d_size );
    cudaMemcpy( host_group.getRawPtr(), device_group, d_size * sizeof(int),
		cudaMemcpyDeviceToHost );

    cudaFree( device_group );
    return host_group;
}

//---------------------------------------------------------------------------//
// Set the entire vector to a group.
void Particle_Vector_Tester::set_group( const int group )
{
    set_group_kernel<<<1,d_size>>>( d_vector.get_device_ptr(), group );
}

//---------------------------------------------------------------------------//
// get a vector of matids.
Teuchos::Array<int> Particle_Vector_Tester::matid()
{
    int* device_matid;
    cudaMalloc( (void**) &device_matid, d_size * sizeof(int) );
    matid_kernel<<<1,d_size>>>( d_vector.get_device_ptr(), device_matid );

    Teuchos::Array<int> host_matid( d_size );
    cudaMemcpy( host_matid.getRawPtr(), device_matid, d_size * sizeof(int),
		cudaMemcpyDeviceToHost );

    cudaFree( device_matid );
    return host_matid;
}

//---------------------------------------------------------------------------//
// Set the entire vector to a matid.
void Particle_Vector_Tester::set_matid( const int matid )
{
    set_matid_kernel<<<1,d_size>>>( d_vector.get_device_ptr(), matid );
}

//---------------------------------------------------------------------------//
// get a vector of alive status.
Teuchos::Array<int> Particle_Vector_Tester::alive()
{
    int* device_alive;
    cudaMalloc( (void**) &device_alive, d_size * sizeof(int) );
    alive_kernel<<<1,d_size>>>( d_vector.get_device_ptr(), device_alive );

    Teuchos::Array<int> host_alive( d_size );
    cudaMemcpy( host_alive.getRawPtr(), device_alive, d_size * sizeof(int),
		cudaMemcpyDeviceToHost );

    cudaFree( device_alive );
    return host_alive;
}

//---------------------------------------------------------------------------//
// set the whole vector to live.
void Particle_Vector_Tester::live()
{
    live_kernel<<<1,d_size>>>( d_vector.get_device_ptr() );
}

//---------------------------------------------------------------------------//
// kill the whole vector.
void Particle_Vector_Tester::kill()
{
    kill_kernel<<<1,d_size>>>( d_vector.get_device_ptr() );
}

//---------------------------------------------------------------------------//
// get a vector of events
Teuchos::Array<typename Particle_Vector_Tester::Event_t> 
Particle_Vector_Tester::event()
{
    Event_t* device_event;
    cudaMalloc( (void**) &device_event, d_size * sizeof(Event_t) );
    event_kernel<<<1,d_size>>>( d_vector.get_device_ptr(), device_event );

    Teuchos::Array<Event_t> host_event( d_size );
    cudaMemcpy( host_event.getRawPtr(), device_event, d_size * sizeof(Event_t),
		cudaMemcpyDeviceToHost );

    cudaFree( device_event );
    return host_event;
}

//---------------------------------------------------------------------------//
// set a vector of events
void Particle_Vector_Tester::set_event( const Teuchos::Array<Event_t>& events )
{
    Event_t* device_event;
    cudaMalloc( (void**) &device_event, d_size * sizeof(Event_t) );
    cudaMemcpy( device_event, events.getRawPtr(), d_size * sizeof(Event_t),
		cudaMemcpyHostToDevice );
    
    set_event_kernel<<<1,d_size>>>( d_vector.get_device_ptr(), device_event );

    cudaFree( device_event );
}

//---------------------------------------------------------------------------//
// set the geometry state for the whole vector
void Particle_Vector_Tester::set_geo_state( const Geo_State_t& geo_state )
{
    set_geo_state_kernel<<<1,d_size>>>( d_vector.get_device_ptr(), geo_state );
}

//---------------------------------------------------------------------------//
// get the geometry state for the whole vector
Teuchos::Array<typename Particle_Vector_Tester::Geo_State_t> 
Particle_Vector_Tester::geo_state()
{
    Geo_State_t* device_geo_state;
    cudaMalloc( (void**) &device_geo_state, d_size * sizeof(Geo_State_t) );
    geo_state_kernel<<<1,d_size>>>( d_vector.get_device_ptr(), device_geo_state );

    Teuchos::Array<Geo_State_t> host_geo_state( d_size );
    cudaMemcpy( host_geo_state.getRawPtr(), device_geo_state, d_size * sizeof(Geo_State_t),
		cudaMemcpyDeviceToHost );

    cudaFree( device_geo_state );
    return host_geo_state;
}

//---------------------------------------------------------------------------//
//                 end of cuda_mc/Particle_Vector_Tester.cu
//---------------------------------------------------------------------------//

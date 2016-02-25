//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   cuda_mc/test/Particle_Vector_Tester.cu
 * \author Stuart Slattery
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Particle_Vector_Tester.hh"

#include "cuda_utils/Hardware.hh"
#include "cuda_utils/CudaDBC.hh"

#include <cuda_runtime.h>

//---------------------------------------------------------------------------//
// CUDA Kernels
//---------------------------------------------------------------------------//
__global__ void ran_kernel( Particle_Vector_Tester::Particle_Vector* vector, 
			    double* ran,
			    const int num_particle )
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if ( i < num_particle ) ran[i] = vector->ran( i );
}

//---------------------------------------------------------------------------//
__global__ void set_wt_kernel( Particle_Vector_Tester::Particle_Vector* vector, 
			       double wt,
			       const int num_particle )
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if ( i < num_particle ) vector->set_wt( i, wt );
}

//---------------------------------------------------------------------------//
__global__ void multiply_wt_kernel( Particle_Vector_Tester::Particle_Vector* vector, 
				    double* wt,
				    const int num_particle )
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if ( i < num_particle ) vector->multiply_wt( i, wt[i] );
}

//---------------------------------------------------------------------------//
__global__ void wt_kernel( Particle_Vector_Tester::Particle_Vector* vector, 
			   double* wt,
			   const int num_particle )
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if ( i < num_particle ) wt[i] = vector->wt( i );
}

//---------------------------------------------------------------------------//
__global__ void group_kernel( Particle_Vector_Tester::Particle_Vector* vector, 
			      int* group,
			      const int num_particle )
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if ( i < num_particle ) group[i] = vector->group( i );
}

//---------------------------------------------------------------------------//
__global__ void set_group_kernel( Particle_Vector_Tester::Particle_Vector* vector, 
				  int group,
				  const int num_particle )
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if ( i < num_particle ) vector->set_group( i, group );
}

//---------------------------------------------------------------------------//
__global__ void matid_kernel( Particle_Vector_Tester::Particle_Vector* vector, 
			      int* matid,
			      const int num_particle )
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if ( i < num_particle ) matid[i] = vector->matid( i );
}

//---------------------------------------------------------------------------//
__global__ void set_matid_kernel( Particle_Vector_Tester::Particle_Vector* vector, 
				  int matid,
				  const int num_particle )
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if ( i < num_particle ) vector->set_matid( i, matid );
}

//---------------------------------------------------------------------------//
__global__ void alive_kernel( Particle_Vector_Tester::Particle_Vector* vector, 
			      int* alive,
			      const int num_particle )
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if ( i < num_particle ) alive[i] = vector->alive( i );
}

//---------------------------------------------------------------------------//
__global__ void live_kernel( Particle_Vector_Tester::Particle_Vector* vector,
			     const int num_particle )
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if ( i < num_particle ) vector->live( i );
}

//---------------------------------------------------------------------------//
__global__ void kill_kernel( Particle_Vector_Tester::Particle_Vector* vector,
			     const int num_particle )
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if ( i < num_particle ) vector->kill( i );
}

//---------------------------------------------------------------------------//
__global__ void set_event_kernel( Particle_Vector_Tester::Particle_Vector* vector, 
				  typename Particle_Vector_Tester::Event_t* event,
				  const int num_particle )
{
    typename std::size_t i = threadIdx.x + blockIdx.x * blockDim.x;
    if ( i < num_particle ) vector->set_event( i, event[i] );
}

//---------------------------------------------------------------------------//
__global__ void event_kernel( Particle_Vector_Tester::Particle_Vector* vector, 
			      typename Particle_Vector_Tester::Event_t* event,
			      const int num_particle )
{
    typename std::size_t i = threadIdx.x + blockIdx.x * blockDim.x;
    if ( i < num_particle ) event[i] = vector->event( i );
}

//---------------------------------------------------------------------------//
__global__ void set_geo_state_kernel( Particle_Vector_Tester::Particle_Vector* vector, 
				      typename Particle_Vector_Tester::Geo_State_t geo_state,
				      const int num_particle )
{
    typename std::size_t i = threadIdx.x + blockIdx.x * blockDim.x;
    if ( i < num_particle ) vector->geo_state( i ) = geo_state;
}

//---------------------------------------------------------------------------//
__global__ void geo_state_kernel( Particle_Vector_Tester::Particle_Vector* vector, 
				  typename Particle_Vector_Tester::Geo_State_t* geo_state,
				  const int num_particle )
{
    typename std::size_t i = threadIdx.x + blockIdx.x * blockDim.x;
    if ( i < num_particle ) geo_state[i] = vector->geo_state( i );
}

//---------------------------------------------------------------------------//
__global__ void batch_kernel( Particle_Vector_Tester::Particle_Vector* vector, 
			      int* batch,
			      const int num_particle )
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if ( i < num_particle ) batch[i] = vector->batch( i );
}

//---------------------------------------------------------------------------//
__global__ void set_batch_kernel( Particle_Vector_Tester::Particle_Vector* vector, 
				  int batch,
				  const int num_particle )
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if ( i < num_particle ) vector->set_batch( i, batch );
}

//---------------------------------------------------------------------------//
__global__ void step_kernel( Particle_Vector_Tester::Particle_Vector* vector, 
			     double* step,
			     const int num_particle )
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if ( i < num_particle ) step[i] = vector->step( i );
}

//---------------------------------------------------------------------------//
__global__ void set_step_kernel( Particle_Vector_Tester::Particle_Vector* vector, 
				 double step,
				 const int num_particle )
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if ( i < num_particle ) vector->set_step( i, step );
}

//---------------------------------------------------------------------------//
// Particle_Vector_Tester
//---------------------------------------------------------------------------//
Particle_Vector_Tester::Particle_Vector_Tester( const int num_particle, 
						const profugus::RNG& rng )
    : d_size( num_particle )
{
    // Acquire hardware for the test.
    if ( !cuda::Hardware<cuda::arch::Device>::have_acquired() )
	cuda::Hardware<cuda::arch::Device>::acquire();

    // Create the vector after hardware acquisition.
    d_vector = cuda::shared_device_ptr<Particle_Vector>( num_particle, rng );
}

//---------------------------------------------------------------------------//
// get a vector of random numbers for the vector.
Teuchos::Array<double> Particle_Vector_Tester::ran()
{
    double* device_ran;
    cudaMalloc( (void**) &device_ran, d_size * sizeof(double) );

    unsigned int threads_per_block = 
	cuda::Hardware<cuda::arch::Device>::num_cores_per_mp();
    unsigned int num_blocks = d_size / threads_per_block;
    if ( d_size % threads_per_block > 0 ) ++num_blocks;

    ran_kernel<<<num_blocks,threads_per_block>>>( 
	d_vector.get_device_ptr(), device_ran, d_size );

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
    unsigned int threads_per_block = 
	cuda::Hardware<cuda::arch::Device>::num_cores_per_mp();
    unsigned int num_blocks = d_size / threads_per_block;
    if ( d_size % threads_per_block > 0 ) ++num_blocks;

    set_wt_kernel<<<num_blocks,threads_per_block>>>( d_vector.get_device_ptr(), wt, d_size );
}

//---------------------------------------------------------------------------//
// mulitply each particle weight by a different value.
void Particle_Vector_Tester::multiply_wt( const Teuchos::Array<double>& wt )
{
    double* device_wt;
    cudaMalloc( (void**) &device_wt, d_size * sizeof(double) );
    cudaMemcpy( device_wt, wt.getRawPtr(), d_size * sizeof(double),
		cudaMemcpyHostToDevice );
    
    unsigned int threads_per_block = 
	cuda::Hardware<cuda::arch::Device>::num_cores_per_mp();
    unsigned int num_blocks = d_size / threads_per_block;
    if ( d_size % threads_per_block > 0 ) ++num_blocks;

    multiply_wt_kernel<<<num_blocks,threads_per_block>>>( 
	d_vector.get_device_ptr(), device_wt, d_size );

    cudaFree( device_wt );
}

//---------------------------------------------------------------------------//
// get a vector of weights.
Teuchos::Array<double> Particle_Vector_Tester::wt()
{
    double* device_wt;
    cudaMalloc( (void**) &device_wt, d_size * sizeof(double) );

    unsigned int threads_per_block = 
	cuda::Hardware<cuda::arch::Device>::num_cores_per_mp();
    unsigned int num_blocks = d_size / threads_per_block;
    if ( d_size % threads_per_block > 0 ) ++num_blocks;

    wt_kernel<<<num_blocks,threads_per_block>>>( 
	d_vector.get_device_ptr(), device_wt, d_size );

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
 
    unsigned int threads_per_block = 
	cuda::Hardware<cuda::arch::Device>::num_cores_per_mp();
    unsigned int num_blocks = d_size / threads_per_block;
    if ( d_size % threads_per_block > 0 ) ++num_blocks;

    group_kernel<<<num_blocks,threads_per_block>>>( 
	d_vector.get_device_ptr(), device_group, d_size );

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
    unsigned int threads_per_block = 
	cuda::Hardware<cuda::arch::Device>::num_cores_per_mp();
    unsigned int num_blocks = d_size / threads_per_block;
    if ( d_size % threads_per_block > 0 ) ++num_blocks;

    set_group_kernel<<<num_blocks,threads_per_block>>>( 
	d_vector.get_device_ptr(), group, d_size );
}

//---------------------------------------------------------------------------//
// get a vector of matids.
Teuchos::Array<int> Particle_Vector_Tester::matid()
{
    int* device_matid;
    cudaMalloc( (void**) &device_matid, d_size * sizeof(int) );

    unsigned int threads_per_block = 
	cuda::Hardware<cuda::arch::Device>::num_cores_per_mp();
    unsigned int num_blocks = d_size / threads_per_block;
    if ( d_size % threads_per_block > 0 ) ++num_blocks;

    matid_kernel<<<num_blocks,threads_per_block>>>( 
	d_vector.get_device_ptr(), device_matid, d_size );

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
    unsigned int threads_per_block = 
	cuda::Hardware<cuda::arch::Device>::num_cores_per_mp();
    unsigned int num_blocks = d_size / threads_per_block;
    if ( d_size % threads_per_block > 0 ) ++num_blocks;

    set_matid_kernel<<<num_blocks,threads_per_block>>>( 
	d_vector.get_device_ptr(), matid, d_size );
}

//---------------------------------------------------------------------------//
// get a vector of alive status.
Teuchos::Array<int> Particle_Vector_Tester::alive()
{
    int* device_alive;
    cudaMalloc( (void**) &device_alive, d_size * sizeof(int) );

    unsigned int threads_per_block = 
	cuda::Hardware<cuda::arch::Device>::num_cores_per_mp();
    unsigned int num_blocks = d_size / threads_per_block;
    if ( d_size % threads_per_block > 0 ) ++num_blocks;

    alive_kernel<<<num_blocks,threads_per_block>>>( 
	d_vector.get_device_ptr(), device_alive, d_size );

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
    unsigned int threads_per_block = 
	cuda::Hardware<cuda::arch::Device>::num_cores_per_mp();
    unsigned int num_blocks = d_size / threads_per_block;
    if ( d_size % threads_per_block > 0 ) ++num_blocks;

    live_kernel<<<num_blocks,threads_per_block>>>( d_vector.get_device_ptr(), d_size );
}

//---------------------------------------------------------------------------//
// kill the whole vector.
void Particle_Vector_Tester::kill()
{
    unsigned int threads_per_block = 
	cuda::Hardware<cuda::arch::Device>::num_cores_per_mp();
    unsigned int num_blocks = d_size / threads_per_block;
    if ( d_size % threads_per_block > 0 ) ++num_blocks;

    kill_kernel<<<num_blocks,threads_per_block>>>( d_vector.get_device_ptr(), d_size );
}

//---------------------------------------------------------------------------//
// get a vector of events
Teuchos::Array<typename Particle_Vector_Tester::Event_t> 
Particle_Vector_Tester::event()
{
    Event_t* device_event;
    cudaMalloc( (void**) &device_event, d_size * sizeof(Event_t) );

    unsigned int threads_per_block = 
	cuda::Hardware<cuda::arch::Device>::num_cores_per_mp();
    unsigned int num_blocks = d_size / threads_per_block;
    if ( d_size % threads_per_block > 0 ) ++num_blocks;

    event_kernel<<<num_blocks,threads_per_block>>>(
	d_vector.get_device_ptr(), device_event, d_size );

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
    
    unsigned int threads_per_block = 
	cuda::Hardware<cuda::arch::Device>::num_cores_per_mp();
    unsigned int num_blocks = d_size / threads_per_block;
    if ( d_size % threads_per_block > 0 ) ++num_blocks;

    set_event_kernel<<<num_blocks,threads_per_block>>>( 
	d_vector.get_device_ptr(), device_event, d_size );

    cudaFree( device_event );
}

//---------------------------------------------------------------------------//
// set the geometry state for the whole vector
void Particle_Vector_Tester::set_geo_state( const Geo_State_t& geo_state )
{
    unsigned int threads_per_block = 
	cuda::Hardware<cuda::arch::Device>::num_cores_per_mp();
    unsigned int num_blocks = d_size / threads_per_block;
    if ( d_size % threads_per_block > 0 ) ++num_blocks;

    set_geo_state_kernel<<<num_blocks,threads_per_block>>>( 
	d_vector.get_device_ptr(), geo_state, d_size );
}

//---------------------------------------------------------------------------//
// get the geometry state for the whole vector
Teuchos::Array<typename Particle_Vector_Tester::Geo_State_t> 
Particle_Vector_Tester::geo_state()
{
    Geo_State_t* device_geo_state;
    cudaMalloc( (void**) &device_geo_state, d_size * sizeof(Geo_State_t) );

    unsigned int threads_per_block = 
	cuda::Hardware<cuda::arch::Device>::num_cores_per_mp();
    unsigned int num_blocks = d_size / threads_per_block;
    if ( d_size % threads_per_block > 0 ) ++num_blocks;

    geo_state_kernel<<<num_blocks,threads_per_block>>>( 
	d_vector.get_device_ptr(), device_geo_state, d_size );

    Teuchos::Array<Geo_State_t> host_geo_state( d_size );
    cudaMemcpy( host_geo_state.getRawPtr(), device_geo_state, d_size * sizeof(Geo_State_t),
		cudaMemcpyDeviceToHost );

    cudaFree( device_geo_state );
    return host_geo_state;
}

//---------------------------------------------------------------------------//
// get a vector of batches.
Teuchos::Array<int> Particle_Vector_Tester::batch()
{
    int* device_batch;
    cudaMalloc( (void**) &device_batch, d_size * sizeof(int) );

    unsigned int threads_per_block = 
	cuda::Hardware<cuda::arch::Device>::num_cores_per_mp();
    unsigned int num_blocks = d_size / threads_per_block;
    if ( d_size % threads_per_block > 0 ) ++num_blocks;

    batch_kernel<<<num_blocks,threads_per_block>>>( 
	d_vector.get_device_ptr(), device_batch, d_size );

    Teuchos::Array<int> host_batch( d_size );
    cudaMemcpy( host_batch.getRawPtr(), device_batch, d_size * sizeof(int),
		cudaMemcpyDeviceToHost );

    cudaFree( device_batch );
    return host_batch;
}

//---------------------------------------------------------------------------//
// Set the entire vector to a batch.
void Particle_Vector_Tester::set_batch( const int batch )
{
    unsigned int threads_per_block = 
	cuda::Hardware<cuda::arch::Device>::num_cores_per_mp();
    unsigned int num_blocks = d_size / threads_per_block;
    if ( d_size % threads_per_block > 0 ) ++num_blocks;

    set_batch_kernel<<<num_blocks,threads_per_block>>>( 
	d_vector.get_device_ptr(), batch, d_size );
}

//---------------------------------------------------------------------------//
// get a vector of stepes.
Teuchos::Array<int> Particle_Vector_Tester::step()
{
    double* device_step;
    cudaMalloc( (void**) &device_step, d_size * sizeof(double) );

    unsigned int threads_per_block = 
	cuda::Hardware<cuda::arch::Device>::num_cores_per_mp();
    unsigned int num_blocks = d_size / threads_per_block;
    if ( d_size % threads_per_block > 0 ) ++num_blocks;

    step_kernel<<<num_blocks,threads_per_block>>>( 
	d_vector.get_device_ptr(), device_step, d_size );

    Teuchos::Array<double> host_step( d_size );
    cudaMemcpy( host_step.getRawPtr(), device_step, d_size * sizeof(double),
		cudaMemcpyDeviceToHost );

    cudaFree( device_step );
    return host_step;
}

//---------------------------------------------------------------------------//
// Set the entire vector to a step.
void Particle_Vector_Tester::set_step( const double step )
{
    unsigned int threads_per_block = 
	cuda::Hardware<cuda::arch::Device>::num_cores_per_mp();
    unsigned int num_blocks = d_size / threads_per_block;
    if ( d_size % threads_per_block > 0 ) ++num_blocks;

    set_step_kernel<<<num_blocks,threads_per_block>>>( 
	d_vector.get_device_ptr(), step, d_size );
}

//---------------------------------------------------------------------------//
//                 end of cuda_mc/Particle_Vector_Tester.cu
//---------------------------------------------------------------------------//

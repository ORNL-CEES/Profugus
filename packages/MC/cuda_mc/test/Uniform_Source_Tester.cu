//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   cuda_mc/test/Uniform_Source_Tester.cu
 * \author Stuart Slattery
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Uniform_Source_Tester.hh"

#include "cuda_utils/Hardware.hh"
#include "cuda_utils/CudaDBC.hh"

#include <Teuchos_Array.hpp>

#include <cuda_runtime.h>

//---------------------------------------------------------------------------//
// CUDA Kernels
//---------------------------------------------------------------------------//
__global__ void wt_kernel( Uniform_Source_Tester::Particle_Vector* vector, 
			   double* wt )
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    wt[i] = vector->wt( i );
}

//---------------------------------------------------------------------------//
__global__ void matid_kernel( Uniform_Source_Tester::Particle_Vector* vector, 
			      int* matid )
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    matid[i] = vector->matid( i );
}

//---------------------------------------------------------------------------//
__global__ void alive_kernel( Uniform_Source_Tester::Particle_Vector* vector, 
			      int* alive )
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    alive[i] = vector->alive( i );
}

//---------------------------------------------------------------------------//
__global__ void group_kernel( Uniform_Source_Tester::Particle_Vector* vector, 
			      int* group )
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    group[i] = vector->group( i );
}

//---------------------------------------------------------------------------//
__global__ void event_kernel( Uniform_Source_Tester::Particle_Vector* vector, 
			      typename Uniform_Source_Tester::Event_t* event )
{
    typename std::size_t i = threadIdx.x + blockIdx.x * blockDim.x;
    event[i] = vector->event( i );
}

//---------------------------------------------------------------------------//
__global__ void batch_kernel( Uniform_Source_Tester::Particle_Vector* vector, 
			      int* batch )
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    batch[i] = vector->batch( i );
}

//---------------------------------------------------------------------------//
__global__ void kill_kernel( Uniform_Source_Tester::Particle_Vector* vector )
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    vector->set_event( i, cuda_profugus::events::DEAD );
}

//---------------------------------------------------------------------------//
// Uniform_Source_Tester
//---------------------------------------------------------------------------//
Uniform_Source_Tester::Uniform_Source_Tester( 
    const std::vector<double>& x_edges,
    const std::vector<double>& y_edges,
    const std::vector<double>& z_edges,
    const int matid,
    const int vector_size,
    const profugus::RNG& rng,
    const int num_group )
    : d_size( vector_size )
{
    // Acquire hardware for the test.
    cuda_utils::Hardware<cuda_utils::arch::Device>::acquire();

    // Create the geometry host .
    std::shared_ptr<Geometry_DMM> host_geom = 
	std::make_shared<Geometry_DMM>( x_edges, y_edges, z_edges );
    int num_cells = host_geom->num_cells();

    // Set matids with the geometry on the host.
    std::vector<int> matids( num_cells, matid );
    host_geom->set_matids( matids );

    // Create a device copy of the geometry.
    d_geometry = cuda_utils::shared_device_ptr<Geometry>(host_geom->device_instance() );

    // Create the particle vector.
    d_particles = cuda_utils::shared_device_ptr<Particle_Vector>( vector_size, rng );

    // Create the source shape.
    d_shape = cuda_utils::shared_device_ptr<Shape>( 
	*std::min_element(x_edges.begin(),x_edges.end()),
	*std::max_element(x_edges.begin(),x_edges.end()),
	*std::min_element(y_edges.begin(),y_edges.end()),
	*std::max_element(y_edges.begin(),y_edges.end()),
	*std::min_element(z_edges.begin(),z_edges.end()),
	*std::max_element(z_edges.begin(),z_edges.end()) );
}

//---------------------------------------------------------------------------//
// get a vector of weights.
Teuchos::Array<double> Uniform_Source_Tester::wt()
{
    double* device_wt;
    cudaMalloc( (void**) &device_wt, d_size * sizeof(double) );
    int num_threads = 64;
    int num_block = d_size / num_threads;
    wt_kernel<<<num_block,num_threads>>>( 
	d_particles.get_device_ptr(), device_wt );

    Teuchos::Array<double> host_wt( d_size );
    cudaMemcpy( host_wt.getRawPtr(), device_wt, d_size * sizeof(double),
		cudaMemcpyDeviceToHost );

    cudaFree( device_wt );
    return host_wt;
}

//---------------------------------------------------------------------------//
// get a vector of matids.
Teuchos::Array<int> Uniform_Source_Tester::matid()
{
    int* device_matid;
    cudaMalloc( (void**) &device_matid, d_size * sizeof(int) );
    int num_threads = 64;
    int num_block = d_size / num_threads;

    matid_kernel<<<num_block,num_threads>>>( 
	d_particles.get_device_ptr(), device_matid );

    Teuchos::Array<int> host_matid( d_size );
    cudaMemcpy( host_matid.getRawPtr(), device_matid, d_size * sizeof(int),
		cudaMemcpyDeviceToHost );

    cudaFree( device_matid );
    return host_matid;
}

//---------------------------------------------------------------------------//
// get a vector of groups.
Teuchos::Array<int> Uniform_Source_Tester::group()
{
    int* device_group;
    cudaMalloc( (void**) &device_group, d_size * sizeof(int) );
    int num_threads = 64;
    int num_block = d_size / num_threads;

    group_kernel<<<num_block,num_threads>>>( 
	d_particles.get_device_ptr(), device_group );

    Teuchos::Array<int> host_group( d_size );
    cudaMemcpy( host_group.getRawPtr(), device_group, d_size * sizeof(int),
		cudaMemcpyDeviceToHost );

    cudaFree( device_group );
    return host_group;
}

//---------------------------------------------------------------------------//
// get a vector of alive status.
Teuchos::Array<int> Uniform_Source_Tester::alive()
{
    int* device_alive;
    cudaMalloc( (void**) &device_alive, d_size * sizeof(int) );
    int num_threads = 64;
    int num_block = d_size / num_threads;

    alive_kernel<<<num_block,num_threads>>>( 
	d_particles.get_device_ptr(), device_alive );

    Teuchos::Array<int> host_alive( d_size );
    cudaMemcpy( host_alive.getRawPtr(), device_alive, d_size * sizeof(int),
		cudaMemcpyDeviceToHost );

    cudaFree( device_alive );
    return host_alive;
}

//---------------------------------------------------------------------------//
// get a vector of events
Teuchos::Array<typename Uniform_Source_Tester::Event_t> 
Uniform_Source_Tester::event()
{
    Event_t* device_event;
    cudaMalloc( (void**) &device_event, d_size * sizeof(Event_t) );
    int num_threads = 64;
    int num_block = d_size / num_threads;

    event_kernel<<<num_block,num_threads>>>(
	d_particles.get_device_ptr(), device_event );

    Teuchos::Array<Event_t> host_event( d_size );
    cudaMemcpy( host_event.getRawPtr(), device_event, d_size * sizeof(Event_t),
		cudaMemcpyDeviceToHost );

    cudaFree( device_event );
    return host_event;
}

//---------------------------------------------------------------------------//
// get a vector of batches.
Teuchos::Array<int> Uniform_Source_Tester::batch()
{
    int* device_batch;
    cudaMalloc( (void**) &device_batch, d_size * sizeof(int) );

    int num_threads = 64;
    int num_block = d_size / num_threads;

    batch_kernel<<<num_block,num_threads>>>( 
	d_particles.get_device_ptr(), device_batch );

    Teuchos::Array<int> host_batch( d_size );
    cudaMemcpy( host_batch.getRawPtr(), device_batch, d_size * sizeof(int),
		cudaMemcpyDeviceToHost );

    cudaFree( device_batch );
    return host_batch;
}

//---------------------------------------------------------------------------//
// kill all the particles.
void Uniform_Source_Tester::kill_particles()
{
    int num_threads = 64;
    int num_block = d_size / num_threads;
    kill_kernel<<<num_block,num_threads>>>( d_particles.get_device_ptr() );
}

//---------------------------------------------------------------------------//
//                 end of cuda_mc/Uniform_Source_Tester.cu
//---------------------------------------------------------------------------//

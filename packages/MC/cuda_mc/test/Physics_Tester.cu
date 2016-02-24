//---------------------------------*-CUDA-*----------------------------------//
/*!
 * \file   cuda_mc/test/Physics_Tester.cu
 * \author Stuart Slattery
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Physics_Tester.hh"

#include "cuda_utils/Hardware.hh"
#include "cuda_utils/CudaDBC.hh"
#include "cuda_utils/Utility_Functions.hh"
#include "cuda_utils/Memory.cuh"

#include <Teuchos_Array.hpp>

#include <cuda_runtime.h>

//---------------------------------------------------------------------------//
// CUDA Kernels
//---------------------------------------------------------------------------//
__global__ void geometry_initialize_kernel( 
    Physics_Tester::Particle_Vector* vector, 
    const Physics_Tester::Geometry* geometry,
    const Physics_Tester::Space_Vector r,
    const Physics_Tester::Space_Vector d,
    const int matid,
    const int num_particle )
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if ( i < num_particle ) 
    {
	geometry->initialize( r, d, vector->geo_state(i) );
	vector->set_matid( i, matid );
	vector->set_event( i, cuda_profugus::events::COLLISION );
	vector->set_wt( i, 0.9 );
    }
}

//---------------------------------------------------------------------------//
__global__ void sample_group_kernel( 
    Physics_Tester::Particle_Vector* vector,
    const double* cdf,
    const int cdf_size,
    const int num_particle )
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if ( i < num_particle ) 
    {
	int group =
	    cuda::utility::sample_discrete_CDF( cdf_size, cdf, vector->ran(i) );
	vector->set_group( i, group );
    }
}

//---------------------------------------------------------------------------//
// Physics_Tester
//---------------------------------------------------------------------------//
Physics_Tester::Physics_Tester( 
    const std::vector<double>& x_edges,
    const std::vector<double>& y_edges,
    const std::vector<double>& z_edges,
    const int vector_size,
    const profugus::RNG& rng,
    Teuchos::ParameterList& db,
    const profugus::XS& xs )
    : d_size( vector_size )
    , d_particle_tester( vector_size, rng )
{
    // Acquire hardware for the test.
    if ( !cuda::Hardware<cuda::arch::Device>::have_acquired() )
	cuda::Hardware<cuda::arch::Device>::acquire();

    // Create the geometry on the host.
    std::shared_ptr<Geometry> host_geom = 
	std::make_shared<Geometry>( x_edges, y_edges, z_edges );
    int num_cells = host_geom->num_cells();

    // Set matids with the geometry on the host.
    std::vector<typename Geometry::matid_type> matids( num_cells, 0 );
    host_geom->set_matids( matids );

    // Create a device copy of the geometry.
    d_geometry = cuda::Shared_Device_Ptr<Geometry>( host_geom );

    // Create the physics.
    d_physics = cuda::shared_device_ptr<Physics>( db, xs );

    // Set the geometry with the physics.
    d_physics.get_host_ptr()->set_geometry( d_geometry );
}

//---------------------------------------------------------------------------//
// Initialize particles with the geometry and set to collide.
void Physics_Tester::geometry_initialize( 
    const Space_Vector r, const Space_Vector d, const int matid )
{
    unsigned int threads_per_block = 
	cuda::Hardware<cuda::arch::Device>::num_cores_per_mp();
    unsigned int num_blocks = d_size / threads_per_block;
    if ( d_size % threads_per_block > 0 ) ++num_blocks;

    geometry_initialize_kernel<<<num_blocks,threads_per_block>>>(
	particles().get_device_ptr(), d_geometry.get_device_ptr(),
	r, d, matid, d_size );
}

//---------------------------------------------------------------------------//
// Sample a cdf and set the particle group.
void Physics_Tester::sample_group( const std::vector<double>& cdf )
{
    // copy the cdf to the device.
    double* device_cdf;
    cuda::memory::Malloc( device_cdf, cdf.size() );
    cuda::memory::Copy_To_Device( device_cdf, cdf.data(), cdf.size() );

    // Sample the cdf and set the particle groups.
    unsigned int threads_per_block = 
	cuda::Hardware<cuda::arch::Device>::num_cores_per_mp();
    unsigned int num_blocks = d_size / threads_per_block;
    if ( d_size % threads_per_block > 0 ) ++num_blocks;

    sample_group_kernel<<<num_blocks,threads_per_block>>>(
	particles().get_device_ptr(), device_cdf, cdf.size(), d_size );

    // free allocated data.
    cuda::memory::Free( device_cdf );
}

//---------------------------------------------------------------------------//
//                 end of cuda_mc/Physics_Tester.cu
//---------------------------------------------------------------------------//

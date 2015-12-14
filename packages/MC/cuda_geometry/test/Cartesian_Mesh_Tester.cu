//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/cuda_geometry/test/Cartesian_Mesh_Tester.cu
 * \author Steven Hamilton
 * \date   Mon Dec 14 13:28:26 2015
 * \brief  Cartesian_Mesh_Tester class definitions.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "utils/View_Field.hh"

#include "../Cartesian_Mesh.hh"
#include "Cartesian_Mesh_Tester.hh"

namespace cuda_profugus
{

__global__ void compute_volumes_kernel(Cartesian_Mesh mesh,
                                       int            num_points,
                                       const int     *cells,
                                       double        *volumes)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if( tid < num_points )
    {
        volumes[tid] = mesh.volume(cells[tid]);
    }
}

//---------------------------------------------------------------------------//
// Constructor
//---------------------------------------------------------------------------//
Cartesian_Mesh_Tester::Cartesian_Mesh_Tester( const Vec_Dbl &x_edges,
                                              const Vec_Dbl &y_edges,
                                              const Vec_Dbl &z_edges )
{
    d_mesh = std::make_shared<Cartesian_Mesh>(x_edges,y_edges,z_edges);
}

//---------------------------------------------------------------------------//
// Compute volumes of specified cells
//---------------------------------------------------------------------------//
void Cartesian_Mesh_Tester::compute_volumes( const Vec_Int &host_cells,
                                                   Vec_Dbl &host_volumes) const
{
    REQUIRE( host_cells.size() == host_volumes.size() );

    int num_points = host_cells.size();

    // Create memroy on device
    typedef cuda::arch::Device Arch;
    cuda::Device_Vector<Arch,int>    device_cells(profugus::make_view(host_cells));
    cuda::Device_Vector<Arch,double> device_volumes(num_points);

    // Execute kernel
    compute_volumes_kernel<<<1,num_points>>>(
            *d_mesh,
             num_points,
             device_cells.data(),
             device_volumes.data());

    REQUIRE( cudaGetLastError() == cudaSuccess );

    // Copy volumes back to host
    device_volumes.to_host(profugus::make_view(host_volumes));
}

//---------------------------------------------------------------------------//
} // end namespace cuda_profugus

//---------------------------------------------------------------------------//
// end of MC/cuda_geometry/test/Cartesian_Mesh_Tester.cu
//---------------------------------------------------------------------------//

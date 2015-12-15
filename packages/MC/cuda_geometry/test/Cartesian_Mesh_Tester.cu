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

__global__ void compute_indices_kernel(Cartesian_Mesh mesh,
                                       int            num_vals,
                                       const int     *ii,
                                       const int     *jj,
                                       const int     *kk,
                                       int           *cells)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if( tid < num_vals )
    {
        size_t cell;
        mesh.index(ii[tid], jj[tid], kk[tid], cell);
        cells[tid] = cell;
    }
}

__global__ void compute_cardinals_kernel(Cartesian_Mesh mesh,
                                        int            num_vals,
                                        const int     *cells,
                                        int           *ii,
                                        int           *jj,
                                        int           *kk)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if( tid < num_vals )
    {
        mesh.cardinal(cells[tid], ii[tid], jj[tid], kk[tid]);
    }
}

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
// Compute indices of specified cells
//---------------------------------------------------------------------------//
void Cartesian_Mesh_Tester::compute_indices( const Vec_Int &host_ii,
                                             const Vec_Int &host_jj,
                                             const Vec_Int &host_kk,
                                                   Vec_Int &host_cells ) const
{
    int num_points = host_ii.size();
    REQUIRE( num_points == host_jj.size() );
    REQUIRE( num_points == host_kk.size() );
    REQUIRE( num_points == host_cells.size() );

    // Create memroy on device
    typedef cuda::arch::Device Arch;
    cuda::Device_Vector<Arch,int> device_ii(profugus::make_view(host_ii));
    cuda::Device_Vector<Arch,int> device_jj(profugus::make_view(host_jj));
    cuda::Device_Vector<Arch,int> device_kk(profugus::make_view(host_kk));
    cuda::Device_Vector<Arch,int> device_cells(num_points);

    // Execute kernel
    compute_indices_kernel<<<1,num_points>>>( *d_mesh,
                                               num_points,
                                               device_ii.data(),
                                               device_jj.data(),
                                               device_kk.data(),
                                               device_cells.data() );

    REQUIRE( cudaGetLastError() == cudaSuccess );

    // Copy data back to host
    device_cells.to_host(profugus::make_view(host_cells));
}

//---------------------------------------------------------------------------//
// Compute cardinal indices of specified cells
//---------------------------------------------------------------------------//
void Cartesian_Mesh_Tester::compute_cardinals( const Vec_Int &host_cells,
                                                     Vec_Int &host_ii,
                                                     Vec_Int &host_jj,
                                                     Vec_Int &host_kk ) const
{
    int num_points = host_cells.size();
    REQUIRE( num_points == host_ii.size() );
    REQUIRE( num_points == host_jj.size() );
    REQUIRE( num_points == host_kk.size() );

    // Create memroy on device
    typedef cuda::arch::Device Arch;
    cuda::Device_Vector<Arch,int> device_cells(profugus::make_view(host_cells));
    cuda::Device_Vector<Arch,int> device_ii(num_points);
    cuda::Device_Vector<Arch,int> device_jj(num_points);
    cuda::Device_Vector<Arch,int> device_kk(num_points);

    // Execute kernel
    compute_cardinals_kernel<<<1,num_points>>>( *d_mesh,
                                                 num_points,
                                                 device_cells.data(),
                                                 device_ii.data(),
                                                 device_jj.data(),
                                                 device_kk.data());

    REQUIRE( cudaGetLastError() == cudaSuccess );

    // Copy data back to host
    device_ii.to_host(profugus::make_view(host_ii));
    device_jj.to_host(profugus::make_view(host_jj));
    device_kk.to_host(profugus::make_view(host_kk));
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

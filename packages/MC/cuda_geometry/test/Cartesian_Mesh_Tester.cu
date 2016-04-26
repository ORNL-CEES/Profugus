//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/cuda_geometry/test/Cartesian_Mesh_Tester.cu
 * \author Steven Hamilton
 * \date   Mon Dec 14 13:28:26 2015
 * \brief  Cartesian_Mesh_Tester class definitions.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include <cuda.h>
#include <cuda_runtime.h>

#include "utils/View_Field.hh"
#include "gtest/Gtest_Functions.hh"
#include "cuda_utils/CudaDBC.hh"
#include "cuda_utils/Device_Vector_device.t.hh"

#include "../Cartesian_Mesh.hh"
#include "Cartesian_Mesh_Tester.hh"

typedef profugus::geometry::cell_type cell_type;
typedef cuda::Space_Vector            Point;
typedef cuda::Coordinates             Coords;
typedef cuda_profugus::Cartesian_Mesh Cartesian_Mesh;

__global__ void compute_indices_kernel(Cartesian_Mesh mesh,
                                       int            num_vals,
                                       const int     *ii,
                                       const int     *jj,
                                       const int     *kk,
                                       cell_type     *cells)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if( tid < num_vals )
    {
        cell_type cell;
        mesh.index(ii[tid], jj[tid], kk[tid], cell);
        cells[tid] = cell;
    }
}

__global__ void compute_cardinals_kernel(Cartesian_Mesh mesh,
                                        int              num_vals,
                                        const cell_type *cells,
                                        int             *ii,
                                        int             *jj,
                                        int             *kk)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if( tid < num_vals )
    {
        mesh.cardinal(cells[tid], ii[tid], jj[tid], kk[tid]);
    }
}

__global__ void compute_volumes_kernel(Cartesian_Mesh mesh,
                                       int              num_points,
                                       const cell_type *cells,
                                       double          *volumes)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if( tid < num_points )
    {
        volumes[tid] = mesh.volume(cells[tid]);
    }
}

__global__ void compute_coords_kernel(Cartesian_Mesh   mesh,
                                      int              num_points,
                                      const Point     *points,
                                      Coords          *coords)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if( tid < num_points )
    {
        mesh.find_upper(points[tid],coords[tid]);
    }
}

namespace
{

// Build Cartesian_Mesh
std::shared_ptr<Cartesian_Mesh> get_mesh()
{
    std::vector<double> x_edges = {0.0, 0.1, 0.6, 0.9, 1.0};
    std::vector<double> y_edges = {-1.0, -0.6, 0.0};
    std::vector<double> z_edges = {2.0, 2.6, 3.4, 4.0};
    
    auto mesh = std::make_shared<Cartesian_Mesh>(
        x_edges,y_edges,z_edges);
    return mesh;
}

}

//---------------------------------------------------------------------------//
// Compute indices of specified cells
//---------------------------------------------------------------------------//
void Cartesian_Mesh_Tester::test_index()
{
    auto mesh = get_mesh();

    std::vector<int> host_ii = {0, 1, 2, 3};
    std::vector<int> host_jj = {1, 0, 1, 0};
    std::vector<int> host_kk = {0, 0, 2, 1};

    int num_points = host_ii.size();

    // Create memroy on device
    typedef cuda::arch::Device Arch;
    cuda::Device_Vector<Arch,int> device_ii(profugus::make_view(host_ii));
    cuda::Device_Vector<Arch,int> device_jj(profugus::make_view(host_jj));
    cuda::Device_Vector<Arch,int> device_kk(profugus::make_view(host_kk));
    cuda::Device_Vector<Arch,cell_type> device_cells(num_points);

    // Execute kernel
    compute_indices_kernel<<<1,num_points>>>( *mesh,
                                               num_points,
                                               device_ii.data(),
                                               device_jj.data(),
                                               device_kk.data(),
                                               device_cells.data() );

    REQUIRE( cudaGetLastError() == cudaSuccess );

    // Copy data back to host
    std::vector<cell_type> host_cells(host_ii.size());
    device_cells.to_host(profugus::make_view(host_cells));

    std::vector<cell_type> expected_cells = {4, 1, 22, 11};
    EXPECT_VEC_EQ(expected_cells, host_cells);
}

//---------------------------------------------------------------------------//
// Compute cardinal indices of specified cells
//---------------------------------------------------------------------------//
void Cartesian_Mesh_Tester::test_cardinal()
{
    auto mesh = get_mesh();

    std::vector<cell_type> host_cells = {4, 1, 22, 11};

    int num_points = host_cells.size();

    // Create memroy on device
    typedef cuda::arch::Device Arch;
    cuda::Device_Vector<Arch,cell_type> device_cells(
        profugus::make_view(host_cells));
    cuda::Device_Vector<Arch,int>       device_ii(num_points);
    cuda::Device_Vector<Arch,int>       device_jj(num_points);
    cuda::Device_Vector<Arch,int>       device_kk(num_points);

    // Execute kernel
    compute_cardinals_kernel<<<1,num_points>>>( *mesh,
                                                 num_points,
                                                 device_cells.data(),
                                                 device_ii.data(),
                                                 device_jj.data(),
                                                 device_kk.data() );

    REQUIRE( cudaGetLastError() == cudaSuccess );

    // Copy data back to host
    std::vector<int> host_ii(num_points);
    std::vector<int> host_jj(num_points);
    std::vector<int> host_kk(num_points);
    device_ii.to_host(profugus::make_view(host_ii));
    device_jj.to_host(profugus::make_view(host_jj));
    device_kk.to_host(profugus::make_view(host_kk));

    std::vector<int> expected_ii = {0, 1, 2, 3};
    std::vector<int> expected_jj = {1, 0, 1, 0};
    std::vector<int> expected_kk = {0, 0, 2, 1};
    EXPECT_VEC_EQ(expected_ii, host_ii);
    EXPECT_VEC_EQ(expected_jj, host_jj);
    EXPECT_VEC_EQ(expected_kk, host_kk);
}

//---------------------------------------------------------------------------//
// Compute volumes of specified cells
//---------------------------------------------------------------------------//
void Cartesian_Mesh_Tester::test_volume()
{
    auto mesh = get_mesh();

    std::vector<cell_type> host_cells = {4, 1, 22, 11};

    int num_points = host_cells.size();

    // Create memroy on device
    typedef cuda::arch::Device Arch;
    cuda::Device_Vector<Arch,cell_type> device_cells(
        profugus::make_view(host_cells));
    cuda::Device_Vector<Arch,double> device_volumes(num_points);

    // Execute kernel
    compute_volumes_kernel<<<1,num_points>>>( *mesh,
                                               num_points,
                                               device_cells.data(),
                                               device_volumes.data() );

    REQUIRE( cudaGetLastError() == cudaSuccess );

    // Copy data back to host
    std::vector<double> host_volumes(num_points);
    device_volumes.to_host(profugus::make_view(host_volumes));

    std::vector<double> expected_volumes = {0.1 * 0.6 * 0.6,
                                            0.5 * 0.4 * 0.6,
                                            0.3 * 0.6 * 0.6,
                                            0.1 * 0.4 * 0.8};
    
    EXPECT_VEC_SOFT_EQ(expected_volumes,host_volumes);
}

//---------------------------------------------------------------------------//
// Test locating points in mesh
//---------------------------------------------------------------------------//
void Cartesian_Mesh_Tester::test_find_upper()
{
    auto x_edges = {0.0, 0.10, 0.25, 0.30, 0.42};
    auto y_edges = {0.0, 0.20, 0.40, 0.50};
    auto z_edges = {-0.1, 0.0, 0.15, 0.50};

    auto mesh = std::make_shared<Cartesian_Mesh>(x_edges,y_edges,z_edges);

    std::vector<Point> host_points = {{-1.0, 1.0, 0.1}, {0.2, 0.45, 0.4}};
    int num_points = host_points.size();

    // Create memroy on device
    typedef cuda::arch::Device Arch;
    cuda::Device_Vector<Arch,Point> device_points(
        profugus::make_view(host_points));
    cuda::Device_Vector<Arch,Coords> device_coords(num_points);

    // Execute kernel
    compute_coords_kernel<<<1,num_points>>>( *mesh,
                                              num_points,
                                              device_points.data(),
                                              device_coords.data() );

    REQUIRE( cudaGetLastError() == cudaSuccess );

    // Copy data back to host
    std::vector<Coords> host_coords(num_points);
    device_coords.to_host(profugus::make_view(host_coords));

    std::vector<Coords> expected_coords = {{-1, 3, 1}, {1, 2, 2}};
    
    EXPECT_EQ(expected_coords[0].i, host_coords[0].i);
    EXPECT_EQ(expected_coords[0].j, host_coords[0].j);
    EXPECT_EQ(expected_coords[0].k, host_coords[0].k);
    EXPECT_EQ(expected_coords[1].i, host_coords[1].i);
    EXPECT_EQ(expected_coords[1].j, host_coords[1].j);
    EXPECT_EQ(expected_coords[1].k, host_coords[1].k);
}

//---------------------------------------------------------------------------//
// end of MC/cuda_geometry/test/Cartesian_Mesh_Tester.cu
//---------------------------------------------------------------------------//

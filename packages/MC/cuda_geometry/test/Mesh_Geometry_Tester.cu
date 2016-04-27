//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/cuda_geometry/test/Mesh_Geometry_Tester.cu
 * \author Steven Hamilton
 * \date   Mon Dec 14 13:28:26 2015
 * \brief  Mesh_Geometry_Tester class definitions.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "gtest/Gtest_Functions.hh"
#include "utils/View_Field.hh"
#include "cuda_utils/Definitions.hh"
#include "cuda_utils/Device_Vector_device.t.hh"
#include "cuda_utils/CudaDBC.hh"

#include "../Mesh_Geometry.hh"
#include "Mesh_Geometry_Tester.hh"

typedef profugus::geometry::cell_type  cell_type;
typedef cuda_utils::Space_Vector       Point;
typedef cuda_utils::Coordinates        Coords;
typedef cuda_profugus::Mesh_Geometry   Mesh_Geometry;

// Get volume from the mesh for each specified cell
__global__ void compute_volumes_kernel(Mesh_Geometry    mesh,
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

// Compute the matid for particle at each specified spatial location
__global__ void compute_matids_kernel(Mesh_Geometry   mesh,
                                      int             num_points,
                                      const Point    *points,
                                      int            *matids)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if( tid < num_points )
    {
        // Create and initialize state on each thread
        // We're only testing matids so direction doesn't matter
        cuda_profugus::Mesh_State state;
        Point dir = {1.0, 0.0, 0.0};
        mesh.initialize(points[tid],dir,state);

        // Get matid
        matids[tid] = mesh.matid(state);
    }
}

// Compute the distance to boundary
__global__ void compute_distance_kernel(Mesh_Geometry   mesh,
                                        int             num_points,
                                        const Point    *points,
                                        const Point    *dirs,
                                        double         *distances,
                                        Coords         *next_ijk)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if( tid < num_points )
    {
        // Create and initialize state on each thread
        cuda_profugus::Mesh_State state;
        mesh.initialize(points[tid],dirs[tid],state);

        distances[tid] = mesh.distance_to_boundary(state);
        next_ijk[tid] = state.next_ijk;
    }
}

// Compute the distance to boundary
__global__ void compute_move_to_surf_kernel(Mesh_Geometry   mesh,
                                            int             num_points,
                                            const Point    *points,
                                            const Point    *dirs,
                                            Coords         *ijk)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if( tid < num_points )
    {
        // Create and initialize state on each thread
        cuda_profugus::Mesh_State state;
        mesh.initialize(points[tid],dirs[tid],state);

        mesh.distance_to_boundary(state);
        mesh.move_to_surface(state);
        ijk[tid] = state.ijk;
    }
}

namespace
{

// Build Mesh_Geometry
std::shared_ptr<Mesh_Geometry> get_mesh()
{
    std::vector<double> x_edges = {0.0, 0.1, 0.6, 0.9, 1.0};
    std::vector<double> y_edges = {-1.0, -0.6, 0.0};
    std::vector<double> z_edges = {2.0, 2.6, 3.4, 4.0};
    
    auto mesh = std::make_shared<Mesh_Geometry>(x_edges,y_edges,z_edges);
    return mesh;
}

}

//---------------------------------------------------------------------------//
// Compute volumes of specified cells
//---------------------------------------------------------------------------//
void Mesh_Geometry_Tester::test_volume()
{
    auto mesh = get_mesh();

    std::vector<cell_type> host_cells = {4, 1, 22, 11};
    int num_points = host_cells.size();

    // Create memory on device
    typedef cuda::arch::Device Arch;
    cuda::Device_Vector<Arch,cell_type> device_cells(
        profugus::make_view(host_cells));
    cuda::Device_Vector<Arch,double> device_volumes(num_points);

    // Execute kernel
    compute_volumes_kernel<<<1,num_points>>>( *mesh,
                                               num_points,
                                               device_cells.data(),
                                               device_volumes.data());

    REQUIRE( cudaGetLastError() == cudaSuccess );

    // Copy volumes back to host
    std::vector<double> host_volumes(num_points);
    device_volumes.to_host(profugus::make_view(host_volumes));

    std::vector<double> expected_volumes = {0.1 * 0.6 * 0.6,
                                            0.5 * 0.4 * 0.6,
                                            0.3 * 0.6 * 0.6,
                                            0.1 * 0.4 * 0.8};

    EXPECT_VEC_SOFT_EQ( expected_volumes, host_volumes );
}

//---------------------------------------------------------------------------//
// Compute matids of specified cells
//---------------------------------------------------------------------------//
void Mesh_Geometry_Tester::test_matid()
{
    auto mesh = get_mesh();

    std::vector<int> all_matids = {1, 3, 2, 0,
                                   3, 1, 4, 1,
                                   2, 5, 2, 1,
                                   0, 1, 2, 3,
                                   1, 2, 3, 4,
                                   2, 3, 4, 5};

    std::vector<Point> host_points = {{0.7,  -0.9,  2.1},
                                      {0.5,  -0.5,  2.5},
                                      {0.99, -0.01, 3.99},
                                      {0.05, -0.8,  2.4}};

    int num_points = host_points.size();
        
    // Set matids with mesh (host call)
    mesh->set_matids(all_matids);

    // Create memroy on device
    typedef cuda::arch::Device Arch;
    cuda::Device_Vector<Arch,Point> device_points(
        profugus::make_view(host_points));
    cuda::Device_Vector<Arch,matid_type> device_cell_matids(num_points);

    // Execute kernel
    compute_matids_kernel<<<1,num_points>>>(
            *mesh,
             num_points,
             device_points.data(),
             device_cell_matids.data());

    REQUIRE( cudaGetLastError() == cudaSuccess );

    // Copy matids back to host
    std::vector<int> host_cell_matids(num_points);
    device_cell_matids.to_host(profugus::make_view(host_cell_matids));

    std::vector<int> expected_matids = {2, 1, 5, 1};

    EXPECT_VEC_EQ( expected_matids, host_cell_matids);
}

//---------------------------------------------------------------------------//
// Test distance to boundary
//---------------------------------------------------------------------------//
void Mesh_Geometry_Tester::test_dist_to_bdry()
{
    // Build mesh
    auto x_edges = {0.0, 0.10, 0.25, 0.30, 0.42};
    auto y_edges = {0.0, 0.20, 0.40, 0.50};
    auto z_edges = {-0.1, 0.0, 0.15, 0.50};
    auto mesh = std::make_shared<Mesh_Geometry>(x_edges,y_edges,z_edges);


    double sqrt_half = sqrt(0.5);
    std::vector<Point> host_points = {{0.01, 0.01, -0.01},
                                      {0.26, 0.35, -0.01}};
    std::vector<Point> host_dirs   = {{1.0,  0.0,   0.0},
                                      {sqrt_half, sqrt_half, 0.0}};

    int num_points = host_points.size();

    // Create memroy on device
    typedef cuda::arch::Device Arch;
    cuda::Device_Vector<Arch,Point> device_points(
        profugus::make_view(host_points));
    cuda::Device_Vector<Arch,Point> device_dirs(
        profugus::make_view(host_dirs));
    cuda::Device_Vector<Arch,double> device_distances(num_points);
    cuda::Device_Vector<Arch,Coords> device_coords(num_points);

    // Execute kernel
    compute_distance_kernel<<<1,num_points>>>(
            *mesh,
             num_points,
             device_points.data(),
             device_dirs.data(),
             device_distances.data(),
             device_coords.data());

    REQUIRE( cudaGetLastError() == cudaSuccess );

    // Copy matids back to host
    std::vector<double> host_distances(num_points);
    device_distances.to_host(profugus::make_view(host_distances));
    std::vector<Coords> host_coords(num_points);
    device_coords.to_host(profugus::make_view(host_coords));

    std::vector<double> expected_distances = {0.1 - 0.01, 0.04 / sqrt_half};
    std::vector<Coords> expected_coords = {{1, 0, 0}, {3, 1, 0}};

    EXPECT_VEC_SOFT_EQ(expected_distances, host_distances);
    EXPECT_EQ(expected_coords[0].i, host_coords[0].i);
    EXPECT_EQ(expected_coords[0].j, host_coords[0].j);
    EXPECT_EQ(expected_coords[0].k, host_coords[0].k);
    EXPECT_EQ(expected_coords[1].i, host_coords[1].i);
    EXPECT_EQ(expected_coords[1].j, host_coords[1].j);
    EXPECT_EQ(expected_coords[1].k, host_coords[1].k);
}

//---------------------------------------------------------------------------//
// Test move to boundary
//---------------------------------------------------------------------------//
void Mesh_Geometry_Tester::test_move_to_surf()
{
    // Build mesh
    auto x_edges = {0.0, 0.10, 0.25, 0.30, 0.42};
    auto y_edges = {0.0, 0.20, 0.40, 0.50};
    auto z_edges = {-0.1, 0.0, 0.15, 0.50};
    auto mesh = std::make_shared<Mesh_Geometry>(x_edges,y_edges,z_edges);


    double sqrt_half = sqrt(0.5);
    std::vector<Point> host_points = {{0.01, 0.01, -0.01},
                                      {0.26, 0.35, -0.01}};
    std::vector<Point> host_dirs   = {{1.0,  0.0,   0.0},
                                      {sqrt_half, sqrt_half, 0.0}};

    int num_points = host_points.size();

    // Create memroy on device
    typedef cuda::arch::Device Arch;
    cuda::Device_Vector<Arch,Point> device_points(
        profugus::make_view(host_points));
    cuda::Device_Vector<Arch,Point> device_dirs(
        profugus::make_view(host_dirs));
    cuda::Device_Vector<Arch,Coords> device_coords(num_points);

    // Execute kernel
    compute_move_to_surf_kernel<<<1,num_points>>>(
            *mesh,
             num_points,
             device_points.data(),
             device_dirs.data(),
             device_coords.data());

    REQUIRE( cudaGetLastError() == cudaSuccess );

    // Copy matids back to host
    std::vector<Coords> host_coords(num_points);
    device_coords.to_host(profugus::make_view(host_coords));

    std::vector<Coords> expected_coords = {{1, 0, 0}, {3, 1, 0}};

    EXPECT_EQ(expected_coords[0].i, host_coords[0].i);
    EXPECT_EQ(expected_coords[0].j, host_coords[0].j);
    EXPECT_EQ(expected_coords[0].k, host_coords[0].k);
    EXPECT_EQ(expected_coords[1].i, host_coords[1].i);
    EXPECT_EQ(expected_coords[1].j, host_coords[1].j);
    EXPECT_EQ(expected_coords[1].k, host_coords[1].k);
}

//---------------------------------------------------------------------------//
// end of MC/cuda_geometry/test/Mesh_Geometry_Tester.cu
//---------------------------------------------------------------------------//

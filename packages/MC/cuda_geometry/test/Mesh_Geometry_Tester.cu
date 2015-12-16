//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/cuda_geometry/test/Mesh_Geometry_Tester.cu
 * \author Steven Hamilton
 * \date   Mon Dec 14 13:28:26 2015
 * \brief  Mesh_Geometry_Tester class definitions.
 * \note   Copyright (c) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "utils/View_Field.hh"
#include "cuda_utils/Definitions.hh"
#include "cuda_utils/Device_Vector_device.t.hh"

#include "../Mesh_Geometry.hh"
#include "Mesh_Geometry_Tester.hh"

namespace cuda_profugus
{

typedef profugus::geometry::cell_type  cell_type;
typedef profugus::geometry::matid_type matid_type;
typedef cuda::Space_Vector             Point;

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
                                      matid_type     *matids)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if( tid < num_points )
    {
        // Create and initialize state on each thread
        // We're only testing matids so direction doesn't matter
        Mesh_State state;
        Point dir = {1.0, 0.0, 0.0};
        mesh.initialize(points[tid],dir,state);

        // Get matid
        matids[tid] = mesh.matid(state);
    }
}

//---------------------------------------------------------------------------//
// Constructor
//---------------------------------------------------------------------------//
Mesh_Geometry_Tester::Mesh_Geometry_Tester( const Vec_Dbl &x_edges,
                                            const Vec_Dbl &y_edges,
                                            const Vec_Dbl &z_edges )
{
    d_mesh = std::make_shared<Mesh_Geometry>(x_edges,y_edges,z_edges);
}

//---------------------------------------------------------------------------//
// Compute volumes of specified cells
//---------------------------------------------------------------------------//
void Mesh_Geometry_Tester::compute_volumes(
        const Vec_Cell_Type &host_cells,
              Vec_Dbl       &host_volumes) const
{
    REQUIRE( host_cells.size() == host_volumes.size() );

    int num_points = host_cells.size();

    // Create memroy on device
    typedef cuda::arch::Device Arch;
    cuda::Device_Vector<Arch,cell_type> device_cells(
        profugus::make_view(host_cells));
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
// Compute matids of specified cells
//---------------------------------------------------------------------------//
void Mesh_Geometry_Tester::compute_matids(
        const Vec_Matid_Type &matids,
        const Vec_Point      &host_points,
              Vec_Matid_Type &host_cell_matids) const
{
    REQUIRE( host_points.size() == host_cell_matids.size() );
    REQUIRE( matids.size() == d_mesh->num_cells() );

    // Set matids with mesh (host call)
    d_mesh->set_matids(matids);

    int num_points = host_points.size();

    // Create memroy on device
    typedef cuda::arch::Device Arch;
    cuda::Device_Vector<Arch,Point> device_points(
        profugus::make_view(host_points));
    cuda::Device_Vector<Arch,matid_type> device_cell_matids(num_points);

    // Execute kernel
    compute_matids_kernel<<<1,num_points>>>(
            *d_mesh,
             num_points,
             device_points.data(),
             device_cell_matids.data());

    REQUIRE( cudaGetLastError() == cudaSuccess );

    // Copy matids back to host
    device_cell_matids.to_host(profugus::make_view(host_cell_matids));
}


//---------------------------------------------------------------------------//
} // end namespace cuda_profugus

//---------------------------------------------------------------------------//
// end of MC/cuda_geometry/test/Mesh_Geometry_Tester.cu
//---------------------------------------------------------------------------//

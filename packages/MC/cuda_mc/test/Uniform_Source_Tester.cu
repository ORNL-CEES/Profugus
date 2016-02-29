//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Uniform_Source_Tester.cc
 * \author Steven Hamilton
 * \date   Wed Jan 20 14:57:16 2016
 * \brief  Uniform_Source_Tester member definitions.
 * \note   Copyright (C) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#include "Uniform_Source_Tester.hh"
#include "../Uniform_Source.cuh"

#include "utils/View_Field.hh"
#include "cuda_utils/Device_Vector.hh"
#include "cuda_geometry/Mesh_Geometry.hh"

namespace cuda_mc
{

typedef cuda_profugus::Mesh_Geometry        Geometry;
typedef cuda_mc::Uniform_Source<Geometry>   Uniform_Src;

__global__ void compute_source_kernel( Uniform_Src           source,
                                       cuda::Space_Vector   *pts,
                                       cuda::Space_Vector   *dirs,
                                       int                   num_vals )
{
     int tid = threadIdx.x + blockIdx.x * blockDim.x;
     if( tid < num_vals )
     {
         curandState_t rng;
         curand_init( 345, tid, 0, &rng );

         auto p = source.get_particle(tid,rng);

         auto geo_state = p.geo_state();
         auto r = geo_state.d_r;
         pts[tid].x = r.x;
         pts[tid].y = r.y;
         pts[tid].z = r.z;

         auto omega = geo_state.d_dir;
         dirs[tid].x = omega.x;
         dirs[tid].y = omega.y;
         dirs[tid].z = omega.z;
     }
}

void Uniform_Source_Tester::test_source( const Vec_Dbl &geom_bounds,
                                         const Vec_Dbl &src_bounds,
                                               Vec_Dbl &x_locs,
                                               Vec_Dbl &y_locs,
                                               Vec_Dbl &z_locs,
                                               Vec_Dbl &mu,
                                               Vec_Dbl &eta,
                                               Vec_Dbl &xi)
{
    // Copy values to device
    int Np = x_locs.size();
    REQUIRE( Np == y_locs.size() );
    REQUIRE( Np == z_locs.size() );
    REQUIRE( Np == mu.size() );
    REQUIRE( Np == eta.size() );
    REQUIRE( Np == xi.size() );
    typedef cuda::arch::Device Arch;
    cuda::Device_Vector<Arch,cuda::Space_Vector> pts_device(Np);
    cuda::Device_Vector<Arch,cuda::Space_Vector> dirs_device(Np);

    REQUIRE( geom_bounds.size() == 6 );

    // Build geometry edges, mesh doesn't really matter
    int Nx = 10;
    double xlo = geom_bounds[0];
    double xhi = geom_bounds[1];
    double dx = (xhi - xlo) / static_cast<double>(Nx);
    std::vector<double> x_edges(Nx+1,0.0);
    for( int ix = 0; ix < Nx+1; ++ix )
        x_edges[ix] = xlo + static_cast<double>(ix) * dx;

    int Ny = 10;
    double ylo = geom_bounds[2];
    double yhi = geom_bounds[3];
    double dy = (yhi - ylo) / static_cast<double>(Ny);
    std::vector<double> y_edges(Ny+1,0.0);
    for( int iy = 0; iy < Ny+1; ++iy )
        y_edges[iy] = ylo + static_cast<double>(iy) * dy;

    int Nz = 10;
    double zlo = geom_bounds[4];
    double zhi = geom_bounds[5];
    double dz = (zhi - zlo) / static_cast<double>(Nz);
    std::vector<double> z_edges(Nz+1,0.0);
    for( int iz = 0; iz < Nz+1; ++iz )
        z_edges[iz] = zlo + static_cast<double>(iz) * dz;

    // Build Geometry
    auto geom = cuda::shared_device_ptr<Geometry>(x_edges,y_edges,z_edges);

    auto db = Teuchos::rcp( new Teuchos::ParameterList() );
    db->set("Np",Np);
    db->set("num_groups",2);

    // Build box shape for source
    REQUIRE( src_bounds.size() == 6 );
    auto src_shape = cuda::shared_device_ptr<cuda_mc::Box_Shape>(
            src_bounds[0], src_bounds[1],
            src_bounds[2], src_bounds[3],
            src_bounds[4], src_bounds[5]);

    // Build source
    Uniform_Src source(db,geom);
    source.build_source(src_shape);

    // Launch kernel
    std::cout << "Launching kernel with " << Np << " particles" << std::endl;
    compute_source_kernel<<<1,Np>>>(source, pts_device.data(),
                                    dirs_device.data(), Np);
    REQUIRE( cudaGetLastError() == cudaSuccess );

    // Copy back to host
    std::vector<cuda::Space_Vector> pts_host(Np);
    std::vector<cuda::Space_Vector> dirs_host(Np);
    pts_device.to_host(profugus::make_view(pts_host));
    dirs_device.to_host(profugus::make_view(dirs_host));

    for( int ip = 0; ip < Np; ++ip )
    {
        x_locs[ip] = pts_host[ip].x;
        y_locs[ip] = pts_host[ip].y;
        z_locs[ip] = pts_host[ip].z;
        mu[ip]     = dirs_host[ip].x;
        eta[ip]    = dirs_host[ip].y;
        xi[ip]     = dirs_host[ip].z;
    }
}

} // end namespace cuda_mc

//---------------------------------------------------------------------------//
//                 end of Uniform_Source_Tester.cc
//---------------------------------------------------------------------------//

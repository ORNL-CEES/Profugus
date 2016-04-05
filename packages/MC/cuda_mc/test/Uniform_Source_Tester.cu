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
#include "../Definitions.hh"

#include "utils/View_Field.hh"
#include "cuda_utils/Device_Vector.hh"
#include "cuda_geometry/Mesh_Geometry.hh"

namespace cuda_mc
{

typedef cuda_profugus::Mesh_Geometry        Geometry;
typedef cuda_mc::Uniform_Source<Geometry>   Uniform_Src;
typedef cuda_mc::Particle<Geometry>         Particle_t;
typedef cuda_mc::RNG_State_t                RNG_State;

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

         auto p = source.get_particle(tid,&rng);

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

// Initialize RNG
__global__ void init_rng( cuda_mc::RNG_State_t *state, int N )
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if( tid < N )
        curand_init( 345, tid, 0, &state[tid] );
}

// Extract data from particles
__global__ void extract_data( Particle_t         *parts,
                              double             *x,
                              double             *y,
                              double             *z,
                              double             *mu,
                              double             *eta,
                              double             *xi,
                              int                 N )
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if( tid < N )
    {
        auto p = parts[tid];

        auto geo_state = p.geo_state();
        auto r = geo_state.d_r;
        x[tid] = r.x;
        y[tid] = r.y;
        z[tid] = r.z;

        auto omega = geo_state.d_dir;
        mu[tid]  = omega.x;
        eta[tid] = omega.y;
        xi[tid]  = omega.z;
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

void Uniform_Source_Tester::test_host_api( const Vec_Dbl &geom_bounds,
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
    auto geom_host = std::make_shared<Geometry>(x_edges,y_edges,z_edges);
    std::vector<int> matids(Nx*Ny*Nz,0);
    geom_host->set_matids(matids);
    cuda::Shared_Device_Ptr<Geometry> geom(geom_host);
    //auto geom = cuda::shared_device_ptr<Geometry>(x_edges,y_edges,z_edges);

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
    auto source_host = std::make_shared<Uniform_Src>(db,geom);
    source_host->build_source(src_shape);
    auto sdp_source = cuda::Shared_Device_Ptr<Uniform_Src>(source_host);

    // Initialize RNG
    thrust::device_vector<RNG_State> rngs( Np );
    init_rng<<<1,Np>>>( rngs.data().get(), Np );

    // Get vector of particles from source
    auto particles = cuda_mc::get_particles(sdp_source,rngs);
    REQUIRE( particles.size() == source_host->num_to_transport() );

    // Extract data from particles
    thrust::device_vector<double> x_dev(   Np );
    thrust::device_vector<double> y_dev(   Np );
    thrust::device_vector<double> z_dev(   Np );
    thrust::device_vector<double> mu_dev(  Np );
    thrust::device_vector<double> eta_dev( Np );
    thrust::device_vector<double> xi_dev(  Np );

    extract_data<<<1,Np>>>( particles.data().get(),
                            x_dev.data().get(),
                            y_dev.data().get(),
                            z_dev.data().get(),
                            mu_dev.data().get(),
                            eta_dev.data().get(),
                            xi_dev.data().get(),
                            Np );

    // Copy back to host
    thrust::host_vector<double> x_host(Np);
    thrust::host_vector<double> y_host(Np);
    thrust::host_vector<double> z_host(Np);
    thrust::host_vector<double> mu_host(Np);
    thrust::host_vector<double> eta_host(Np);
    thrust::host_vector<double> xi_host(Np);
    x_host   = x_dev;
    y_host   = y_dev;
    z_host   = z_dev;
    mu_host  = mu_dev;
    eta_host = eta_dev;
    xi_host  = xi_dev;

    for( int ip = 0; ip < Np; ++ip )
    {
        x_locs[ip] = x_host[ip];
        y_locs[ip] = y_host[ip];
        z_locs[ip] = z_host[ip];
        mu[ip]     = mu_host[ip];
        eta[ip]    = eta_host[ip];
        xi[ip]     = xi_host[ip];
    }
}

} // end namespace cuda_mc

//---------------------------------------------------------------------------//
//                 end of Uniform_Source_Tester.cc
//---------------------------------------------------------------------------//

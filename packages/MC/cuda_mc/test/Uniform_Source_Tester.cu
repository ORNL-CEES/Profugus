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
#include "../Source_Provider.cuh"

#include "gtest/Gtest_Functions.hh"
#include "utils/View_Field.hh"
#include "cuda_geometry/Mesh_Geometry.hh"

using namespace cuda_mc;

typedef cuda_utils::Space_Vector            Space_Vector;
typedef cuda_profugus::Mesh_Geometry        Geometry;
typedef cuda_profugus::Mesh_Geometry_DMM    Geometry_DMM;
typedef Uniform_Source<Geometry>            Uniform_Src;
typedef Uniform_Source_DMM<Geometry>        Uniform_Src_DMM;
typedef Particle<Geometry>                  Particle_t;

__global__ void compute_source_kernel( Uniform_Src   source,
                                       Space_Vector *pts,
                                       Space_Vector *dirs,
                                       int           num_vals )
{
     using def::I; using def::J; using def::K;

     int tid = threadIdx.x + blockIdx.x * blockDim.x;
     if( tid < num_vals )
     {
         curandState_t rng;
         curand_init( 345, tid, 0, &rng );

         auto p = source.get_particle(tid,&rng);

         auto geo_state = p.geo_state();
         pts[tid]  = geo_state.d_r;
         dirs[tid] = geo_state.d_dir;
     }
}

// Extract data from particles
__global__ void extract_data( Particle_t         *parts,
                              Space_Vector       *pts,
                              Space_Vector       *dirs,
                              int                 N )
{
    using def::I; using def::J; using def::K;

    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if( tid < N )
    {
        auto p = parts[tid];

        auto geo_state = p.geo_state();
        pts[tid]  = geo_state.d_r;
        dirs[tid] = geo_state.d_dir;
    }
}

std::shared_ptr<Geometry_DMM> get_geometry()
{
    // Build geometry edges, mesh doesn't really matter
    std::vector<double> geom_bounds = {-4.0, 2.0, -1.0, 1.0, 1.0, 12.0};
    double xlo_geom = -4.0;
    double xhi_geom =  2.0;
    double ylo_geom = -1.0;
    double yhi_geom =  1.0;
    double zlo_geom =  1.0;
    double zhi_geom = 12.0;
    int Nx = 10;
    double dx = (xhi_geom - xlo_geom) / static_cast<double>(Nx);
    std::vector<double> x_edges(Nx+1,0.0);
    for( int ix = 0; ix < Nx+1; ++ix )
        x_edges[ix] = xlo_geom + static_cast<double>(ix) * dx;

    int Ny = 10;
    double dy = (yhi_geom - ylo_geom) / static_cast<double>(Ny);
    std::vector<double> y_edges(Ny+1,0.0);
    for( int iy = 0; iy < Ny+1; ++iy )
        y_edges[iy] = ylo_geom + static_cast<double>(iy) * dy;

    int Nz = 10;
    double dz = (zhi_geom - zlo_geom) / static_cast<double>(Nz);
    std::vector<double> z_edges(Nz+1,0.0);
    for( int iz = 0; iz < Nz+1; ++iz )
        z_edges[iz] = zlo_geom + static_cast<double>(iz) * dz;

    // Build Geometry
    auto geom_host = std::make_shared<Geometry_DMM>(x_edges,y_edges,z_edges);
    std::vector<int> matids(Nx*Ny*Nz,0);
    geom_host->set_matids(matids);
    return geom_host;
}

void test_data(const thrust::device_vector<Space_Vector> &pts,
               const thrust::device_vector<Space_Vector> &dirs,
                     std::vector<double>                  bounds)
{
    using def::I; using def::J; using def::K;

    def::size_type Np = pts.size();
    EXPECT_EQ(dirs.size(),Np);

    // Copy data to host
    thrust::host_vector<Space_Vector> host_pts  = pts;
    thrust::host_vector<Space_Vector> host_dirs = dirs;

    // Now compute mean of position in each direction
    double x_expected = 0.5 * (bounds[0] + bounds[1]);
    double y_expected = 0.5 * (bounds[2] + bounds[3]);
    double z_expected = 0.5 * (bounds[4] + bounds[5]);
    double x_mean = 0.0;
    double y_mean = 0.0;
    double z_mean = 0.0;
    for( auto pt : host_pts )
    {
        x_mean += pt[I];
        y_mean += pt[J];
        z_mean += pt[K];
    }
    x_mean /= static_cast<double>(Np);
    y_mean /= static_cast<double>(Np);
    z_mean /= static_cast<double>(Np);

    std::cout << "Expecting mean position of (" << x_expected
        << ", " << y_expected << ", " << z_expected << ")" << std::endl;
    std::cout << "Actual mean x = " << x_mean << std::endl;
    std::cout << "Actual mean y = " << y_mean << std::endl;
    std::cout << "Actual mean z = " << z_mean << std::endl;

    double sqrt_N = std::sqrt(static_cast<double>(Np));
    EXPECT_TRUE(std::abs(x_mean - x_expected) <
                (bounds[1] - bounds[0])/sqrt_N );
    EXPECT_TRUE(std::abs(y_mean - y_expected) <
                (bounds[3] - bounds[2])/sqrt_N );
    EXPECT_TRUE(std::abs(z_mean - z_expected) <
                (bounds[5] - bounds[4])/sqrt_N );

    // Now compute mean of each direction component
    double mu_mean = 0.0;
    double eta_mean = 0.0;
    double xi_mean = 0.0;
    for( auto dir : host_dirs)
    {
        mu_mean  += dir[I];
        eta_mean += dir[J];
        xi_mean  += dir[K];
    }
    mu_mean  /= static_cast<double>(Np);
    eta_mean /= static_cast<double>(Np);
    xi_mean  /= static_cast<double>(Np);

    std::cout << "Expecting mean angle of (0,0,0)" << std::endl;
    std::cout << "Actual mean mu = " << mu_mean << std::endl;
    std::cout << "Actual mean eta = " << eta_mean << std::endl;
    std::cout << "Actual mean xi = " << xi_mean << std::endl;

    // All angle sampling is isotropic so expected mean is zero
    double mean_expected = 0.0;
    EXPECT_TRUE( std::abs(mu_mean  - mean_expected) < 1.0/sqrt_N );
    EXPECT_TRUE( std::abs(eta_mean - mean_expected) < 1.0/sqrt_N );
    EXPECT_TRUE( std::abs(xi_mean  - mean_expected) < 1.0/sqrt_N );
}

void Uniform_Source_Tester::test_source()
{
    def::size_type Np = 1024;

    auto geom_dmm = get_geometry();
    auto geom = 
        cuda::shared_device_ptr<Geometry>(geom_dmm->device_instance());

    auto db = Teuchos::rcp( new Teuchos::ParameterList() );
    db->set("Np",Np);
    db->set("num_groups",2);

    // Build box shape for source
    double xlo =  0.0;
    double xhi =  1.0;
    double ylo = -1.0;
    double yhi =  1.0;
    double zlo =  3.0;
    double zhi =  5.0;
    auto src_shape = cuda::shared_device_ptr<cuda_mc::Box_Shape>(
            xlo, xhi, ylo, yhi, zlo, zhi);

    // Build source
    Uniform_Src_DMM source(db,geom);
    source.build_source(src_shape);

    // Allocate device vectors
    thrust::device_vector<Space_Vector> device_pts(Np);
    thrust::device_vector<Space_Vector> device_dirs(Np);

    // Launch kernel
    std::cout << "Launching kernel with " << Np << " particles" << std::endl;
    compute_source_kernel<<<1,Np>>>(source.device_instance(),
                                    device_pts.data().get(),
                                    device_dirs.data().get(), Np);
    REQUIRE( cudaGetLastError() == cudaSuccess );

    test_data(device_pts,device_dirs,{xlo, xhi, ylo, yhi, zlo, zhi});
}

void Uniform_Source_Tester::test_host_api()
{

    // Copy values to device
    def::size_type Np = 1024;

    auto geom_dmm = get_geometry();
    auto geom = 
        cuda::shared_device_ptr<Geometry>(geom_dmm->device_instance());

    auto db = Teuchos::rcp( new Teuchos::ParameterList() );
    db->set("Np",Np);
    db->set("num_groups",2);

    // Build box shape for source
    double xlo =  0.0;
    double xhi =  1.0;
    double ylo = -1.0;
    double yhi =  1.0;
    double zlo =  3.0;
    double zhi =  5.0;
    auto src_shape = cuda::shared_device_ptr<cuda_mc::Box_Shape>(
            xlo, xhi, ylo, yhi, zlo, zhi);

    // Build source
    auto source_host = std::make_shared<Uniform_Src_DMM>(db,geom);
    source_host->build_source(src_shape);

    // Initialize RNG
    auto rng_control = std::make_shared<cuda_mc::RNG_Control>(1234);

    // Get vector of particles from source
    thrust::device_vector<Particle<Geometry>> particles;

    cuda_mc::Source_Provider<Geometry> provider;
    thrust::device_vector<int> indices(Np);
    thrust::counting_iterator<int> cnt(0);
    thrust::copy(cnt,cnt+indices.size(),indices.begin());
    provider.get_particles(source_host,rng_control,particles,indices);

    EXPECT_EQ( particles.size(), Np );
    EXPECT_EQ( source_host->num_to_transport(), Np );

    // Extract data from particles
    thrust::device_vector<Space_Vector> pts(Np);
    thrust::device_vector<Space_Vector> dirs(Np);

    extract_data<<<1,Np>>>( particles.data().get(),
                            pts.data().get(),
                            dirs.data().get(),
                            Np );


    test_data(pts,dirs,{xlo, xhi, ylo, yhi, zlo, zhi});

}

//---------------------------------------------------------------------------//
//                 end of Uniform_Source_Tester.cc
//---------------------------------------------------------------------------//

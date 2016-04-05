//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cuda_mc/Uniform_Source.t.cuh
 * \author Steven Hamilton
 * \date   Tue May 06 16:43:26 2014
 * \brief  Uniform_Source member definitions.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef cuda_mc_Uniform_Source_t_cuh
#define cuda_mc_Uniform_Source_t_cuh

#include <numeric>

#include "Teuchos_Array.hpp"

#include "utils/View_Field.hh"
#include "cuda_utils/CudaDBC.hh"
#include "cuda_utils/Utility_Functions.hh"
#include "cuda_utils/Launch_Args.t.cuh"
#include "comm/Timing.hh"
#include "comm/global.hh"
#include "Sampler.cuh"
#include "Uniform_Source.cuh"

namespace cuda_mc
{
//---------------------------------------------------------------------------//
// Functor to populate vector with particles from source
//---------------------------------------------------------------------------//
template <class Geometry>
class Compute_Source
{
  public:

    typedef Uniform_Source<Geometry> UniSource;
    typedef cuda_mc::RNG_State_t     RNG_State;
    typedef Particle<Geometry>       Particle_t;

    Compute_Source( const UniSource *source,
                    RNG_State       *rngs,
                    Particle_t      *particles)
        : d_source(source)
        , d_rngs(rngs)
        , d_particles(particles)
    {
    }

    __device__ void operator()( std::size_t tid ) const
    {
        d_particles[tid] = d_source->get_particle(tid,&d_rngs[tid]);
        ENSURE( d_particles[tid].alive() );
    }

  private:

    const Uniform_Source<Geometry> *d_source;
    RNG_State *d_rngs;
    Particle_t *d_particles;
};

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 *
 * \param db
 * \param geometry
 * \param physics
 */
template <class Geometry>
Uniform_Source<Geometry>::Uniform_Source(RCP_Std_DB     db,
                                         SDP_Geometry   geometry)
    : Base(geometry)
    , d_erg_cdf(db->get<int>("num_groups"))
    , d_np_total(0)
{
    REQUIRE(!db.is_null());

    REQUIRE( profugus::nodes() == 1 );

    // store the total number of requested particles
    d_np_total = static_cast<size_type>(db->get("Np", 1000));
    INSIST(d_np_total > 0., "Number of source particles must be positive");
    d_np_domain = d_np_total / profugus::nodes();

    // get the spectral shape
    d_num_groups = db->get<int>("num_groups");
    const auto &shape = db->get(
        "spectral_shape", Teuchos::Array<double>(d_num_groups, 1.0));
    CHECK(shape.size() == d_num_groups);

    // calculate the normalization
    double norm = std::accumulate(shape.begin(), shape.end(), 0.0);
    CHECK(norm > 0.0);

    // assign to the shape cdf
    double sum = 0.0;
    norm  = 1.0 / norm;
    int n = 0;
    Vec_Dbl erg_cdf(d_num_groups,0.0);
    for (double &c : erg_cdf)
    {
        double val = shape[n] * norm;
        sum += val;
        c = sum;
        ++n;
    }
    ENSURE(cuda::utility::soft_equiv(sum, 1.0));

    // Copy to device
    d_erg_cdf.assign(profugus::make_view(erg_cdf));

    // initialize timers in this class, which may be necessary because domains
    // with no source will not make this timer otherwise
#if UTILS_TIMING > 0
    profugus::Timing_Diagnostics::update_timer(
        "profugus::Uniform_Source.get_particle", 0.0);
#endif
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Build the initial source.
 * \param geometric_shape
 */
template <class Geometry>
void Uniform_Source<Geometry>::build_source(SDP_Shape geometric_shape)
{
    REQUIRE(geometric_shape.get_host_ptr());
    REQUIRE(geometric_shape.get_device_ptr());

    SCOPED_TIMER("profugus::Uniform_Source.build_source");

    // store the spatial shape
    d_geo_shape = geometric_shape.get_device_ptr();

    profugus::global_barrier();
}

//---------------------------------------------------------------------------//
//! \brief Build vector of particles
//---------------------------------------------------------------------------//

template <class Geometry>
thrust::device_vector<Particle<Geometry> > get_particles(
    cuda::Shared_Device_Ptr<Uniform_Source<Geometry>> &source,
    thrust::device_vector<cuda_mc::RNG_State_t> &rngs )
{
    int Np = source.get_host_ptr()->num_to_transport();
    REQUIRE( rngs.size() > 0 );
    REQUIRE( rngs.size() >= Np );

    // Build particle vector
    thrust::device_vector<Particle<Geometry>> particles(Np);

    // Functor to populate vector
    Compute_Source<Geometry> f( source.get_device_ptr(),
            rngs.data().get(), particles.data().get() );

    // Launch kernel
    cuda::Launch_Args<cuda::arch::Device> launch_args;
    launch_args.set_num_elements( Np );

    cuda::parallel_launch( f, launch_args );

    return particles;
}

} // end namespace cuda_mc

#endif // cuda_mc_Uniform_Source_t_cuh

//---------------------------------------------------------------------------//
//                 end of Uniform_Source.t.cuh
//---------------------------------------------------------------------------//

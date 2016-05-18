//---------------------------------*-C++-*-----------------------------------//
/*!
 * \file   MC/cuda_mc/Source_Provider.t.hh
 * \author Steven Hamilton
 * \date   Wed Apr 13 08:38:06 2016
 * \brief  Source_Provider template method definitions.
 * \note   Copyright (c) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef MC_cuda_mc_Source_Provider_t_hh
#define MC_cuda_mc_Source_Provider_t_hh

#include "cuda_utils/Launch_Args.t.cuh"

#include "Source_Provider.cuh"
#include "Uniform_Source.cuh"
#include "Fission_Source.cuh"

namespace cuda_mc
{

//---------------------------------------------------------------------------//
// Functor to populate vector with particles from source
//---------------------------------------------------------------------------//
template <class Geometry, class Src_Type>
class Compute_Source
{
  public:

    typedef cuda_mc::RNG_State_t     RNG_State;
    typedef Particle<Geometry>       Particle_t;

    Compute_Source( const Src_Type  *source,
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
        DEVICE_ENSURE( d_particles[tid].alive() );
    }

  private:

    const Src_Type  *d_source;
    RNG_State       *d_rngs;
    Particle_t      *d_particles;
};
    

//---------------------------------------------------------------------------//
// Get starting particles from source
//---------------------------------------------------------------------------//
template <class Geometry>
void Source_Provider<Geometry>::get_particles(
    SP_Source source, SP_RNG_Control rng_control, 
    Particle_Vector &particles ) const
{
    // Determine type of source
    typedef Uniform_Source<Geometry> Uni_Source;
    typedef Fission_Source<Geometry> Fisn_Source;
    std::shared_ptr<Uni_Source> uni_source;
    std::shared_ptr<Fisn_Source> fisn_source;
    uni_source  = std::dynamic_pointer_cast<Uni_Source>(source);
    fisn_source = std::dynamic_pointer_cast<Fisn_Source>(source);

    // Call to implementation for appropriate source type
    if( uni_source )
    {
        REQUIRE( !fisn_source );

        get_particles_impl(uni_source,rng_control,particles);
    }
    else if( fisn_source )
    {
        REQUIRE( !uni_source );

        get_particles_impl(fisn_source,rng_control,particles);
    }
}

template <class Geometry>
template <class Src_Type>
void Source_Provider<Geometry>::get_particles_impl(
        std::shared_ptr<Src_Type>   source,
        SP_RNG_Control              rng_control,
        Particle_Vector             &particles ) const
{
    int Np = source->num_to_transport();

    rng_control->initialize(Np);

    particles.resize(Np);

    cuda::Shared_Device_Ptr<Src_Type> sdp_source(source);

    // Build functor to populate vector
    auto &rngs = rng_control->get_states();
    Compute_Source<Geometry,Src_Type> f( sdp_source.get_device_ptr(),
                                         rngs.data().get(),
                                         particles.data().get() );

    // Launch kernel
    cuda::Launch_Args<cuda::arch::Device> launch_args;
    launch_args.set_num_elements( Np );

    cuda::parallel_launch( f, launch_args );

    cudaDeviceSynchronize();
}


//---------------------------------------------------------------------------//
} // end namespace cuda_mc

#endif // MC_cuda_mc_Source_Provider_t_hh

//---------------------------------------------------------------------------//
// end of MC/cuda_mc/Source_Provider.t.hh
//---------------------------------------------------------------------------//

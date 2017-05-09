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

    typedef Particle_Vector<Geometry>   Particle_Vector_t;
    typedef def::size_type              size_type;

    Compute_Source( const Src_Type    source,
                    Particle_Vector_t particles,
                    const int        *indices,
                    size_type         num_particles)
        : d_source(source)
        , d_particles(particles)
        , d_indices(indices)
        , d_num_particles(num_particles)
    {
    }

    __device__ void operator()( std::size_t tid )
    {
        int pid = d_indices[tid];
        DEVICE_REQUIRE(pid < d_num_particles);
        d_source.build_particle(pid,&d_particles);
        DEVICE_ENSURE( d_particles.alive(pid) );
    }

  private:

    const Src_Type     d_source;
    Particle_Vector_t  d_particles;
    const int         *d_indices;
    size_type          d_num_particles;
};
    

//---------------------------------------------------------------------------//
// Get starting particles from source
//---------------------------------------------------------------------------//
template <class Geometry>
void Source_Provider<Geometry>::get_particles(
    SP_Source source,
    SP_Particle_Vector_DMM particles, const Index_Vector &indices) const
{
    // Get particle count and set up RNG, particles, source batching
    size_type Np = std::min(source->num_left(),indices.size());
    particles->initialize(Np);
    source->begin_batch(Np);

    // Determine type of source
    typedef Uniform_Source_DMM<Geometry> Uni_Source_DMM;
    typedef Fission_Source_DMM<Geometry> Fisn_Source_DMM;
    std::shared_ptr<Uni_Source_DMM> uni_source_dmm;
    std::shared_ptr<Fisn_Source_DMM> fisn_source_dmm;
    uni_source_dmm  = std::dynamic_pointer_cast<Uni_Source_DMM>(source);
    fisn_source_dmm = std::dynamic_pointer_cast<Fisn_Source_DMM>(source);

    // Call to implementation for appropriate source type
    if( uni_source_dmm )
    {
        REQUIRE( !fisn_source_dmm );

        auto uni_source = uni_source_dmm->device_instance();

        get_particles_impl(uni_source,particles,indices,Np);
    }
    else if( fisn_source_dmm )
    {
        REQUIRE( !uni_source_dmm );

        auto fisn_source = fisn_source_dmm->device_instance();

        get_particles_impl(fisn_source,particles,indices,Np);
    }
    else
    {
        VALIDATE(false,"Unknown source type.");
    }

    // Update remaining particles in source
    source->end_batch(Np);
}

//---------------------------------------------------------------------------//
// Geometry-templated implementation
//---------------------------------------------------------------------------//
template <class Geometry>
template <class Src_Type>
void Source_Provider<Geometry>::get_particles_impl(
        Src_Type                   &source,
        SP_Particle_Vector_DMM      particles,
        const Index_Vector          &indices,
        size_type                    Np) const
{
    // Build functor to populate vector
    Compute_Source<Geometry,Src_Type> f( source,
                                         particles->device_instance(),
                                         indices.data().get(),
                                         Np);

    // Launch kernel
    cuda::Launch_Args<cuda::arch::Device> launch_args;
    launch_args.set_num_elements(Np);

    cuda::parallel_launch( f, launch_args );

    cudaDeviceSynchronize();
}

//---------------------------------------------------------------------------//
} // end namespace cuda_mc

#endif // MC_cuda_mc_Source_Provider_t_hh

//---------------------------------------------------------------------------//
// end of MC/cuda_mc/Source_Provider.t.hh
//---------------------------------------------------------------------------//
